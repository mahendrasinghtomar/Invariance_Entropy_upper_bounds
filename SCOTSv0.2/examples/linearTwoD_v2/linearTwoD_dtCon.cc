/*
 *   author: MST
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

std::string MATLABDIR = "";

    /* memory profiling */
    #include <sys/time.h>
    #include <sys/resource.h>
    struct rusage usage;

#include <stdio.h>
#include <math.h>
#include <algorithm>

#include "Helper_MST.h"

const std::string Example = "linearTwoD_dtCon";
/* state space dim */
const int state_dim=2;
/* input space dim */
const int input_dim_a=1;	//actual input dimension
const int input_length=1;
const int input_dim=input_dim_a * input_length; // extended input dimension for longer length inputs

/*
 * data types for the elements of the state space
 * and input space used by the ODE solver
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

#define ReachSet 0 
double eta_s = 0.01; 
double eta_i = 0.5;

#if RUN_INDEX == 2
    #include "dtControlClassify.h"   
#endif

auto reach_set = [](state_type& r, state_type& x, const input_type& u) noexcept {
};

auto system_post = [](state_type &x, const input_type &u) noexcept {
   /* the ode describing the dcdc converter */
   state_type z=x;
   for(int i=0;i<input_length;i++){
	   x[0]=2*z[0]+u[i];
	   x[1]=0.5*z[1]+u[i];
	   z=x;
   }
};

auto radius_post = [](state_type &r, const state_type&, const input_type &u) noexcept {
    /* the ode for the growth bound */
    state_type zz=r;
    for(int i=0;i<input_length;i++){
		r[0]=2*zz[0];
        r[1]=0.5*zz[1];
		zz=r;
	}
};

int main() {
  /* to measure time */
  TicToc tt;

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
  input_type i_lb, i_ub, i_eta;
  for(int i=0;i<input_length;i++){
	  /* lower bounds of the hyper rectangle */
	  i_lb[i]=-1;
      /* upper bounds of the hyper rectangle */
      i_ub[i]=1;
      /* grid node distance diameter */
      i_eta[i]=eta_i;
  }

  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  #if (RUN_INDEX==1)
  is.print_info();
  #endif

  /* compute transition function of symbolic model */
  std::cout << "Computing the transition function:\n";
  /* transition function of symbolic model */
  scots::TransitionFunction tf;
  scots::Abstraction<state_type,input_type> abs(ss,is);
  abs.verbose_off();


#if (RUN_INDEX == 1)
    tt.tic();
    abs.compute_gb(tf,system_post, radius_post);
    std::cout << "Number of transitions: " << tf.get_no_transitions() <<"\n";
      if(!getrusage(RUSAGE_SELF, &usage))
        std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << "\n" << "Memory = " << usage.ru_maxrss << "\n";
  tt.toc();

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
#if MINPOST == 0
    write_to_file(win, Example+"winDomain");
      std::cout << "\nWrite controller to " << Example << ".scs \n";
      if(write_to_file(scots::StaticController(ss,is,std::move(win)),Example))
        std::cout << "Done. \n";
#endif  // end MINPOST

#endif  // end RUN_INDEX == 1

    #if RUN_INDEX == 2
      scots::WinningDomain win_r;
      if (!read_from_file(win_r, Example + "winDomain")) {
          std::cout << "\nFailed to read the file containing winning domain.\n";
          return 1;   // When couldn't read file, 'make' will interrupt and show error
      }
      abs_type Bsize = win_r.get_size();
      std::vector<abs_type> Bpartition(Bsize);
      win_r.get_winning_domain_MST(Bpartition);
      std::vector<int> B_labels(Bsize);   // the labels of Bpartition as per dtc partitioning
      std::vector<std::vector<abs_type>> Bi(Bsize);
      abs.compute_BiBjBlabels(tf, system_post, radius_post, Bpartition, Bsize, Bi, B_labels, classify);

        TicToc writingB;
        writingB.tic();
        std::cout<<"Writing B_labels...";
        wfile(B_labels, "B_labels", 1, 0, Example);
        std::cout<<"Done.\nWriting Bi...";
        wfile(Bi, 1, 1, Example);
        std::cout<<"Done.\n";
        writingB.toc();

    
    
    #endif  // end RUN_INDEX == 2
	
      if(!getrusage(RUSAGE_SELF, &usage))
        std::cout <<  "Memory = " << usage.ru_maxrss << "\n";
	
   if(Bsize == 0) return 1;
      return 0;
    }



