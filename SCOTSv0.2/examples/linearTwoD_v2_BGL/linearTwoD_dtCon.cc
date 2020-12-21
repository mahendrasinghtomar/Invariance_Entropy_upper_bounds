

#include <iostream>
#include <array>
#include <cmath>

/* SCOTS header */
#include "scots.hh"
/* ode solver */
#include "RungeKutta4.hh"
/* time profiling */
#include "TicToc.hh"

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
std::string OSname = "Windows environment";
std::string MATLABDIR = "C:/Users/mahen/GoogleDrive/Currently_working_on_these/lyf_MST_documents/lyf_pc_mac/cpp_packages/SCOTSv0.2/examples/";
#define NOTVISUALSTUDIO 0
#elif __APPLE__
std::string OSname = "Apple Mac";
// Mac
std::string MATLABDIR = "/Users/mst/GoogleDrive/Currently_working_on_these/lyf_MST_documents/lyf_pc_mac/cpp_packages/SCOTSv0.2/examples/";
#define NOTVISUALSTUDIO 1
#elif __linux__
std::string OSname = "Linux environment";
std::string MATLABDIR = "/mnt/c/Users/mahen/GoogleDrive/Currently_working_on_these/lyf_MST_documents/lyf_pc_mac/cpp_packages/SCOTSv0.2/examples/";
#define NOTVISUALSTUDIO 1
#else
    // unknown compiler
#endif

#if NOTVISUALSTUDIO
    /* memory profiling */
    #include <sys/time.h>
    #include <sys/resource.h>
    struct rusage usage;
#endif

#include <stdio.h>
#include <math.h>
#include <algorithm>

#include "Helper_MST_BGL.h"

const std::string Example = "linearTwoD_v2_BGL"; // Name of folder of example
/* state space dim */
const int state_dim=2;
/* input space dim */
const int input_dim=1;

/*
 * data types for the elements of the state space
 * and input space used by the ODE solver
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;
#include "Helper_MST.h"

#define ReachSet 0 // 0 for system_post, radius_post. 1 for reach_set
double eta_s = 0.01; //std::pow(10,-1);
double eta_i = 0.5;

// double eta_s = std::pow(10,-2);
// double eta_i = 0.5;

#if RUN_INDEX == 2
    #include "dtControlClassify.h"
    //#include "Helper_MST.h"
#endif

auto reach_set = [](state_type &r, state_type &x, const input_type &u ) noexcept {
};

/* we integrate the dcdc ode by 0.5 sec (the result is stored in x)  */
auto system_post = [](state_type &x, const input_type &u) noexcept {
   /* the ode describing the dcdc converter */
   state_type z=x;
   x[0]=2*z[0]+u[0];
   x[1]=0.5*z[1]+u[0];
};

auto radius_post = [](state_type &r, const state_type&, const input_type &u) noexcept {
    /* the ode for the growth bound */
    state_type zz=r;
    r[0]=2*zz[0];
    r[1]=0.5*zz[1];
};

int main() {
	MATLABDIR = MATLABDIR + Example + "/";
  /* to measure time */
  TicToc tt, nettime;
  nettime.tic();

  /* setup the workspace of the synthesis problem and the uniform grid */
   /* grid node distance diameter */
    state_type eta={{eta_s,eta_s}};
 /* lower bounds of the hyper-rectangle */
  state_type lb={{-1,-2}};
  /* upper bounds of the hyper-rectangle */
  state_type ub={{1,2}};
  scots::UniformGrid ss(state_dim,lb,ub,eta);
  #if (RUN_INDEX==1)
	std::cout << "Uniform grid details:\n";
	ss.print_info();
  #endif
  
  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
  input_type i_lb={{-1}};
  /* upper bounds of the hyper rectangle */
  input_type i_ub={{1}};
   /* grid node distance diameter */
  input_type i_eta={{eta_i}};
  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  #if (RUN_INDEX==1) 
	is.print_info();
  #endif

  /* compute transition function of symbolic model */
  /* transition function of symbolic model */
  scots::TransitionFunction tf;
  scots::Abstraction<state_type,input_type> abs(ss,is);
  //abs.verbose_off();

  tt.tic();
#if (RUN_INDEX == 1)
	std::cout << "Computing the transition function:\n";
    abs.compute_gb(tf,system_post, radius_post);
    std::cout << "Number of transitions: " << tf.get_no_transitions() <<"\n";
    #if NOTVISUALSTUDIO
      if(!getrusage(RUSAGE_SELF, &usage))
        std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << "\n" << "Memory = " << usage.ru_maxrss << "\n";
    #endif
#else	// RUN_INDEX == 2 :
    scots::WinningDomain win_r;
    if (!read_from_file(win_r, Example+"winDomain")){
        std::cout << "\nFailed to read the file containing winning domain.\n";
        return 1;   // When couldn't read file, 'make' will interrupt and show error
    }
    abs_type Bsize = win_r.get_size();
    std::vector<abs_type> Bpartition(Bsize);
    win_r.get_winning_domain_MST(Bpartition);
    std::vector<int> B_labels(Bsize);   // the labels of Bpartition as per dtc partitioning
    std::vector<std::vector<abs_type>> Bi(Bsize);
    abs.compute_BiBjBlabels(tf,system_post, radius_post, Bpartition, Bsize, Bi, B_labels, classify);
#endif
  tt.toc();

#if (RUN_INDEX == 1)
  /* continue with synthesis */
  /* define function to check if the cell is in the safe set  */
  auto safeset = [&lb, &ub, &ss, &eta](const scots::abs_type& idx) noexcept {
    state_type x;
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (lb[0] <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0)<= ub[0] && lb[1] <= (x[1]-eta[1]/2.0) && (x[1]+eta[1]/2.0)<= ub[1])
      return true;
    return false;
  };
  /* compute winning domain (contains also valid inputs) */
  std::cout << "\nSynthesis: \n";
  tt.tic();
  scots::WinningDomain win = scots::solve_invariance_game(tf,safeset);
  tt.toc();
  abs_type Bsize = win.get_size();
  std::cout << "Winning domain size: " << Bsize << "\n";
  std::vector<abs_type> Bpartition(Bsize);
  win.get_winning_domain_MST(Bpartition);   // winning states copied to Bpartition
  
  std::vector<std::vector<abs_type>> Xj(Bsize);
  std::vector<abs_type> DControllerMinPost(Bsize); // deterministic controller with first such input selected for which no. of posts is minimum
  std::vector<int> B_NPosts(Bsize);   // the number of post states for states in Bpartition
  
  int int_ReachSet = (ReachSet == 1) ? 1 : 0;
  DeterminiseController_Xj_NPosts(int_ReachSet, system_post, radius_post, reach_set, win, ss, is, abs, Bpartition, Bsize, Xj, DControllerMinPost, B_NPosts);

  Graph g;
  /*number of edges in the graph*/
  abs_type num_edges = 0;
  for (const auto& e : B_NPosts)
      num_edges += e;
  /*compute edge weights*/
  std::unique_ptr<double[]> edge_weights = std::make_unique<double[]>(num_edges);
  abs_type count = 0;
  for (abs_type i = 0; i < Xj.size(); i++) {
      for (abs_type j = 0; j < Xj[i].size(); j++) {
          edge_weights[count] = std::log2(B_NPosts[i]);
          count++;
      }
  }

  /*computes max edge weight*/
  // /*double max_edge_weight = *std::max_element(edge_weights.get(), edge_weights.get() + num_edges); */
  // double max_edge_weight = 0;
  // for(auto it=edge_weights.get(); it!= edge_weights.get()+num_edges; it++)
	  // if(*it>max_edge_weight)
		  // max_edge_weight = *it;
  // std::cout << "\nMaximum edge weight, wmax = " << max_edge_weight << "\n";

  /*constructs graph.*/
  constructDirectedGraph(Xj, B_NPosts, Bpartition, g, edge_weights);
  /*prints the maximum cycle mean*/
  print_max_cycle_mean(g);

  /*creates graph, and computes and prints the maximum cyle mean.*/
  //createDirectedGraph(Xj, B_NPosts, Bpartition, g);
  /*computes max edge weight over SCCs with at least one edge.*/
  //compute_wmax....
  
  // createDirectedGraph(Xj, B_NPosts, Bpartition, g);
  // wfile_DController(Bpartition, DControllerMinPost, Example, 1, 0);

  // write_to_file(win, Example+"winDomain");
  // std::cout << "\nWrite controller to " << Example << ".scs \n";
  // if(write_to_file(scots::StaticController(ss,is,std::move(win)),Example))
        // std::cout << "Done. \n";
  
  // TicToc writingB;
    // writingB.tic();
    // std::cout<<"Writing B_NPosts...";
    // wfile(B_NPosts, "B_NPosts", 1, 0, Example);
    // std::cout<<"Done.\nWriting Xi Xj...";
	// wfile(Xj, Bpartition, 1, 0, Example);
    // std::cout<<"Done.\n";
  // writingB.toc();
  // writeGraph_dot(Bpartition,g);

#endif

    #if RUN_INDEX == 2
        TicToc writingB;
        writingB.tic();
        std::cout<<"Writing B_labels...";
        wfile(B_labels, "B_labels", 1, 0, Example);
        std::cout<<"Done.\nWriting Bi...";
        wfile(Bi, 1, 1, Example);
        // std::cout<<"Done.\nWriting Bj...";
        // wfile(Bj, "Bj", 1, 1, Example);
        std::cout<<"Done.\n";
        writingB.toc();   
		
    #endif
	
	#if NOTVISUALSTUDIO
      if(!getrusage(RUSAGE_SELF, &usage))
        std::cout <<  "Memory = " << usage.ru_maxrss << "\n";
    #endif

	std::cout << "\nTotal time:";
	nettime.toc();
	
   if(Bsize == 0) return 1;
      return 0;
    }



