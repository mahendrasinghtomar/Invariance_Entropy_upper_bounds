/*
*with SCOTS_FOCUS_ONLY_INTERSECTION = 0
*	Bypass computation of SCOTS-abstraction during the computation of
*	invariant controller with intersection-domain as the safe set.
*
*with SCOTS_FOCUS_ONLY_INTERSECTION = 1
*	avoid set defined as the set X\domain_intersection.
*/

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
std::string MATLABDIR = "C:/Users/..../SCOTSv0.2/examples/";
#define NOTVISUALSTUDIO 0
#elif __APPLE__
std::string OSname = "Apple Mac";
// Mac
std::string MATLABDIR = "/Users/..../SCOTSv0.2/examples/";
#define NOTVISUALSTUDIO 1
#elif __linux__
std::string OSname = "Linux environment";
std::string MATLABDIR = "/mnt/c/Users/..../SCOTSv0.2/examples/";
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
#include <boost/numeric/interval.hpp>
#include <algorithm>

#include "Helper_MST.h"
#include "Helper_MST_BGL.h"

//typedef boost::numeric::interval<double> interval;
namespace bnumeric = boost::numeric;
namespace Ilib =  boost::numeric::interval_lib;
typedef boost::numeric::interval<double, Ilib::policies<Ilib::save_state<Ilib::rounded_transc_std<double> >,Ilib::checking_base<double>>> interval;

double sup(interval& I){
    return I.upper();
}

double inf(interval& I){
    return I.lower();
}

double mag(interval& I){
    return boost::numeric::norm(I);
}

double midpoint(interval& I){
    return boost::numeric::median(I);
}

double radius(interval& I){
    return 0.5 * boost::numeric::width(I);
}

const std::string Example = "hyperbolic_v2_BGL"; // Name of folder of example
/* state space dim */
const int state_dim = 2;
/* input space dim */
const int input_dim = 2;

#define ReachSet 0 // 0 for system_post, radius_post. 1 for reach_set
double epsil = 0.08;
double eta_s = 0.009; 	//std::pow(10,-2);
double x_rad = 5.9574/2.0;
double eta_i = 0.01;

/*
 * data types for the elements of the state space
 * and input space used by the ODE solver
 */
using state_type = std::array<double, state_dim>;
using input_type = std::array<double, input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

#if RUN_INDEX == 2
    #include "dtControlClassify.h"
#endif

auto reach_set = [](state_type &r, state_type &x, const input_type &u ) noexcept {
    interval xt0(x[0]-r[0], x[0]+r[0]);
    interval xt1(x[1]-r[1], x[1]+r[1]);
    interval xx0, xx1;
	
    xx0 = 5.0 - 0.3 * xt1 - bnumeric::square(xt0) + u[0];
    xx1 = xt0 + u[1];
	
    x[0] = midpoint(xx0);
    x[1] = midpoint(xx1);
    r[0] = radius(xx0);
    r[1] = radius(xx1);
};

auto reach_set_inv = [](state_type &r, state_type &x, const input_type &u ) noexcept {
    interval xt0(x[0]-r[0], x[0]+r[0]);
    interval xt1(x[1]-r[1], x[1]+r[1]);
    interval xx0, xx1;
	
	/* for the inverse system */
	xx0 = xt1 - u[1];
	xx1 = (5.0 - bnumeric::square(xt1 - u[1]) + u[0] - xt0 )/0.3;
	
    x[0] = midpoint(xx0);
    x[1] = midpoint(xx1);
    r[0] = radius(xx0);
    r[1] = radius(xx1);
};

auto system_post = [](state_type &x, const input_type &u) noexcept {
    state_type z = x;
	
    x[0]= 5.0 - 0.3 * z[1] - z[0] * z[0] + u[0];
    x[1] = z[0] + u[1];
		
   };
   
   auto system_post_inv = [](state_type &x, const input_type &u) noexcept {
    state_type z = x;
 	
	// for the inverse system
	x[0] = z[1] - u[1];
	x[1] = (5 - (z[1] - u[1])*(z[1] - u[1]) + u[0] - z[0])/0.3;
	
   };


auto radius_post = [](state_type &r, const state_type &x, const input_type &u) noexcept {
    state_type z = r;
	
    r[0] = 0.3 * z[1] + z[0] * z[0] + 2 * std::abs(x[0]) * z[0];
    r[1] = z[0];
	
};
 
 auto radius_post_inv = [](state_type &r, const state_type &x, const input_type &u) noexcept {
    state_type z = r;
	
	// for the inverse system
	r[0] = z[1];
	//r[1] = (z[1]*z[1] + 2*std::abs(x[1])*z[1] + 2*std::abs(u[1])*z[1] + z[0] )/0.3; //commented on 15July2020 in the ..._minPost folder
    r[1] = (z[1] * z[1] + 2 * std::abs(x[1] - u[1]) * z[1] + z[0]) / 0.3;
 };
 

int main() {
	MATLABDIR = MATLABDIR + Example + "/";
	std::vector<abs_type> domain_fwd, domain_inv;
	
  /* to measure time */
  TicToc tt, nettime;
  nettime.tic();

  /* setup the workspace of the synthesis problem and the uniform grid */
   /* grid node distance diameter */
  state_type eta={{eta_s, eta_s}};
    
 /* lower bounds of the hyper-rectangle */
  state_type lb={{-x_rad, -x_rad}};
 /* upper bounds of the hyper-rectangle */
  state_type ub={{x_rad, x_rad}};
    
  scots::UniformGrid ss(state_dim,lb,ub,eta);
  #if (RUN_INDEX==1)
  std::cout << "Uniform grid details:\n";
  ss.print_info();
  #endif
  
  /* construct grid for the input space */
    input_type i_lb, i_ub, i_eta;
    for(int i=0;i<input_dim;i++){
        i_lb[i] = -epsil;
        i_ub[i] = epsil;
        i_eta[i] = eta_i;
    }
    
  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  #if (RUN_INDEX==1)
  is.print_info();
  #endif
   
  /* compute transition function of symbolic model */
  /* transition function of symbolic model */
  scots::TransitionFunction tf_fwd;
  
  /* for fwd system */
  {
  scots::Abstraction<state_type,input_type> abs(ss,is);
  //abs.verbose_off();

#if (RUN_INDEX == 1)
	tt.tic();
    std::cout << "Computing the transition function:\n";
    #if (ReachSet == 0)
        abs.compute_gb(tf_fwd,system_post, radius_post);
    #else
        abs.compute_gb2(tf_fwd,system_post, radius_post, reach_set);
    #endif
      
      
//      state_type xtemp;
//      ss.itox(7578331, xtemp);
//      input_type utemp;
//      is.itox(43, utemp);
//      /* as per the transition function */
//      std::vector<state_type> post_list = abs.get_post(reach_set, xtemp, utemp, "MST");
//      std::cout << "Posts for 7578331 under 43 as per reach_set: ";
//      for(int ii =0; ii<post_list.size();ii++)
//          std::cout << ss.xtoi(post_list[ii]) << ", ";
////      /*std::cout <<"Posts for 7578331 under 43 as per transition function: ";
////      abs.print_post(tf_fwd, xtemp, utemp);*/
      
      
	std::cout << "Number of transitions: " << tf_fwd.get_no_transitions() <<"\n";
    #if NOTVISUALSTUDIO
      if(!getrusage(RUSAGE_SELF, &usage))
        std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf_fwd.get_no_transitions() << "\n" << "Memory = " << usage.ru_maxrss << "\n";
    #endif
	tt.toc();
#endif
 
  
#if (RUN_INDEX == 1)
  /* continue with synthesis */
  /* define function to check if the cell is in the safe set  */
  auto safeset = [&lb, &ub, &ss, &eta](const scots::abs_type& idx) noexcept {
    state_type x;
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if ((lb[0] <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0) <= ub[0]) && (lb[1] <= (x[1]-eta[1]/2.0) && (x[1]+eta[1]/2.0) <= ub[1]))
      return true;
    return false;
  };
  /* compute winning domain (contains also valid inputs) */
  std::cout << "\nSynthesis: \n";
  tt.tic();
  scots::WinningDomain win = scots::solve_invariance_game(tf_fwd,safeset);
  tt.toc();
  abs_type Bsize = win.get_size();
  std::cout << "Winning domain size: " << Bsize << "\n";
  domain_fwd = win.get_winning_domain();
   // std::cout << "\nWrite controller to " << Example << "_fwd.scs \n";
  // if(write_to_file(scots::StaticController(ss,is,std::move(win)),Example+"_fwd"))
    // std::cout << "Done. \n";
  
  if(Bsize == 0) return 1;
#endif

	}
	
	
	
	
	{
  /* for inverse system */
  /* compute transition function of symbolic model */
  /* transition function of symbolic model */
  scots::TransitionFunction tf;
  scots::Abstraction<state_type,input_type> abs(ss,is);
  //abs.verbose_off();

#if (RUN_INDEX == 1)
	tt.tic();
    std::cout << "\nComputing the transition function for inverse system:\n";
    #if (ReachSet == 0)
        abs.compute_gb(tf,system_post_inv, radius_post_inv);
    #else
        abs.compute_gb2(tf,system_post_inv, radius_post_inv, reach_set_inv);
    #endif
	std::cout << "Number of transitions: " << tf.get_no_transitions() <<"\n";
    #if NOTVISUALSTUDIO
      if(!getrusage(RUSAGE_SELF, &usage))
        std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << "\n"<< "Memory = " << usage.ru_maxrss << "\n";
    #endif
  tt.toc();
#endif
  
#if (RUN_INDEX == 1)
  /* continue with synthesis */
  /* define function to check if the cell is in the safe set  */
  auto safeset = [&lb, &ub, &ss, &eta](const scots::abs_type& idx) noexcept {
    state_type x;
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if ((lb[0] <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0) <= ub[0]) && (lb[1] <= (x[1]-eta[1]/2.0) && (x[1]+eta[1]/2.0) <= ub[1]))
      return true;
    return false;
  };
  /* compute winning domain (contains also valid inputs) */
  std::cout << "\nSynthesis for inverse sys: \n";
  tt.tic();
  scots::WinningDomain win = scots::solve_invariance_game(tf,safeset);
  tt.toc();
  abs_type Bsize = win.get_size();
  std::cout << "Winning domain size: " << Bsize << "\n";

    // std::vector<abs_type> winningDomain(Bsize);
    // win.get_winning_domain_MST(winningDomain);
    // wfile(winningDomain, "WinningDomain.txt");
  domain_inv = win.get_winning_domain();
  // write_to_file(win, Example+"winDomain");
  // std::cout << "\nWrite controller to " << Example << "_inv.scs \n";
  // if(write_to_file(scots::StaticController(ss,is,std::move(win)),Example+"_inv"))
    // std::cout << "Done. \n";
    
#endif
}

#if (RUN_INDEX==1)
std::vector<abs_type> domain_intersection;
unsigned long fwd_size, inv_size;
fwd_size = domain_fwd.size();
inv_size = domain_inv.size();
for(unsigned long i=0; i<fwd_size; i++)
{
	for(unsigned long j=0; j<inv_size; j++)
	{
		if(domain_fwd[i] == domain_inv[j])
		{
			domain_intersection.push_back(domain_fwd[i]);
			break;
		}
	}
}
// wfile(domain_intersection, "domain_intersection", 1, 0, Example);	// redundant; just for the sake of cross-checking
std::cout << "\nDomain intersection size: "<< domain_intersection.size() << "\n";
#endif	
	
abs_type Bsize;
{
  /* With intersection of fwd and inv domains as the safe set*/
  scots::Abstraction<state_type,input_type> abs(ss,is);
  //abs.verbose_off();
  
#if (RUN_INDEX == 1)
  /* continue with synthesis */
  /* define function to check if the cell is in the safe set  */
#if 0
  auto safeset = [&domain_intersection](const scots::abs_type& idx) noexcept {
      for (const auto& e : domain_intersection)
          if (idx == e)
              return true;
      return false;
  };
#else
  auto safeset = [&domain_intersection](const scots::abs_type& idx) noexcept {
      if (std::binary_search(domain_intersection.begin(), domain_intersection.end(), idx))
          return true;
      else
          return false;
  };
#endif  // end safeset selection
  /* compute winning domain (contains also valid inputs) */
  std::cout << "\nSynthesis for intersection: \n";
  tt.tic();
  scots::WinningDomain win = scots::solve_invariance_game(tf_fwd,safeset);
  tt.toc();
  Bsize = win.get_size();
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

  // write_to_file(win, Example+"with_intersection_winDomain");
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
		// tt.tic();
	// std::cout<<"Computing Bi, Bj, Blabels\n";
	// scots::WinningDomain win_r;
    // if (!read_from_file(win_r, Example+"with_intersection_winDomain")){
         // std::cout << "\nFailed to read the file containing winning domain.\n";
         // return 1;   // When couldn't read file, 'make' will interrupt and show error
     // }
	// Bsize = win_r.get_size();
	// std::vector<abs_type> Bpartition(Bsize);
    // win_r.get_winning_domain_MST(Bpartition);
	// std::vector<int> B_labels(Bsize);   // the labels of Bpartition as per dtc partitioning
	// std::vector<std::vector<abs_type>> Bi(Bsize);

    // #if (ReachSet == 0)
        
        // if(!abs.compute_BiBjBlabels(tf_fwd,system_post, radius_post, Bpartition, Bsize, Bi, B_labels, classify))
			// return 1;
    // #else
        // if(!abs.compute_BiBjBlabels2(tf_fwd,system_post, radius_post, reach_set, Bpartition, Bsize, Bi, B_labels, classify))
            // return 1;
    // #endif
    // TicToc writingB;
    // writingB.tic();
    // std::cout<<"Writing B_labels...";
    // wfile(B_labels, "B_labels", 1, 0, Example);
    // std::cout<<"Done.\nWriting Bi...";
    // wfile(Bi, 1, 1, Example);
    // std::cout<<"Done.\n";
    // writingB.toc();
    
#endif	// end RUN_INDEX == 2
	}
	#if NOTVISUALSTUDIO
      if(!getrusage(RUSAGE_SELF, &usage))
        std::cout <<  "Memory = " << usage.ru_maxrss << "\n";
    #endif

	std::cout << "\nTotal time:";
	nettime.toc();
	if(Bsize == 0) return 1;
  return 0;
}
