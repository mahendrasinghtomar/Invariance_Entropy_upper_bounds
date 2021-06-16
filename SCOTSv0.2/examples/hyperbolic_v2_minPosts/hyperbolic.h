// #define DebugPrint 1

#include <iostream>
#include <array>
#include <cmath>

/* SCOTS header */
#include "scots.hh"
/* ode solver */
#include "RungeKutta4.hh"


/* time profiling */
#include "TicToc.hh"

#define NOTVISUALSTUDIO 0   //0: visual studio. 1: else than visual studio
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

//typedef boost::numeric::interval<double> interval;
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


const std::string Example = "hyperbolic";
/* state space dim */
const int state_dim=2;
/* input space dim */
const int input_dim=2;
/* sampling time */

/*
 * data types for the elements of the state space
 * and input space used by the ODE solver
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;


#define ReachSet 0 // 0 for system_post, radius_post. 1 for reach_set
double epsil = 0.08;
double eta_s = std::pow(10,-2);
double x_rad = 5.9574/2.0;
double eta_i = 0.01;

#if RUN_INDEX == 2
    #include "dtControlClassify.h"
	#include "Helper_MST.h"
#endif

auto reach_set = [](state_type &r, state_type &x, const input_type &u ) noexcept {
    interval xt0(x[0]-r[0], x[0]+r[0]);
    interval xt1(x[1]-r[1], x[1]+r[1]);
    interval xx0, xx1;
    xx0 = 5.0 - 0.3 * xt1 - xt0*xt0 + u[0];
    xx1 = xt0 + u[1];
	
	// for the inverse system
	// xx0 = xt1 - u[1];
	// xx1 = (5.0 - (xt1-u[1])*(xt1-u[1]) + u[0] - xt0 )/0.3;
	
    x[0] = midpoint(xx0);
    x[1] = midpoint(xx1);
    r[0] = radius(xx0);
    r[1] = radius(xx1);
};

auto system_post = [](state_type &x, const input_type &u) noexcept {
    state_type z = x;
    x[0]= 5.0 - 0.3 * z[1] - z[0] * z[0] + u[0];
    x[1] = z[0] + u[1];
	
	// for the inverse system
	// x[0] = z[1] - u[1];
	// x[1] = (5 - (z[1] - u[1])*(z[1] - u[1]) + u[0] - z[0])/0.3;
	
   };

auto radius_post = [](state_type &r, const state_type &x, const input_type &u) noexcept {
    state_type z = r;
    r[0] = 0.3 * z[1] + z[0] * z[0] + 2 * std::abs(x[0]) * z[0];
    r[1] = z[0];
	
	// for the inverse system
	// r[0] = z[1];
	// r[1] = (z[1]*z[1] + 2*std::abs(x[1])*z[1] + 2*std::abs(u[1])*z[1] + z[0] )/0.3;
 };

int main() {
  /* to measure time */
  TicToc tt;

  /* setup the workspace of the synthesis problem and the uniform grid */
   /* grid node distance diameter */
  state_type eta={{eta_s, eta_s}};
    
 /* lower bounds of the hyper-rectangle */
  state_type lb={{-x_rad, -x_rad}};
 /* upper bounds of the hyper-rectangle */
  state_type ub={{x_rad, x_rad}};
    
//      state_type lb= {{-1.1071}};
//      state_type ub= {{1.1071}};
    
  scots::UniformGrid ss(state_dim,lb,ub,eta);
  std::cout << "Uniform grid details:\n";
  ss.print_info();
  
  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
//  /* upper bounds of the hyper rectangle */
//  /* grid node distance diameter */
    input_type i_lb, i_ub, i_eta;
    for(int i=0;i<input_dim;i++){
        i_lb[i] = -epsil;
        i_ub[i] = epsil;
        i_eta[i] = eta_i;
    }
    
  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  is.print_info();

  /* construct grid for the input alphabet */
  /* hyper-rectangle [1,2] with grid node distance 1 */
  /*scots::UniformGrid is(input_dim,input_type{{1}},input_type{{2}},input_type{{1}});
  is.print_info();*/

  /* compute transition function of symbolic model */
  //std::cout << "Computing the transition function:\n";
  /* transition function of symbolic model */
  scots::TransitionFunction tf;
  scots::Abstraction<state_type,input_type> abs(ss,is);
  abs.verbose_off();

  tt.tic();
#if (RUN_INDEX == 1)
    std::cout << "Computing the transition function:\n";
    #if (ReachSet == 0)
        abs.compute_gb(tf,system_post, radius_post);
    #else
        abs.compute_gb2(tf,system_post, radius_post, reach_set);
    #endif
	std::cout << "Number of transitions: " << tf.get_no_transitions() <<"\n";
    #if NOTVISUALSTUDIO
      if(!getrusage(RUSAGE_SELF, &usage))
        std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << "\n";
    #endif
#else   // RUN_INDEX == 2 :
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

    #if (ReachSet == 0)
        
        abs.compute_BiBjBlabels(tf,system_post, radius_post, Bpartition, Bsize, Bi, B_labels, classify);
    #else
        abs.compute_BiBjBlabels2(tf,system_post, radius_post, reach_set, Bpartition, Bsize, Bi, B_labels, classify);
    #endif
#endif
  tt.toc();
  
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
  scots::WinningDomain win = scots::solve_invariance_game(tf,safeset);
  tt.toc();
  abs_type Bsize = win.get_size();
  std::cout << "Winning domain size: " << Bsize << "\n";

    // std::vector<abs_type> winningDomain(Bsize);
    // win.get_winning_domain_MST(winningDomain);
    // wfile(winningDomain, "WinningDomain.txt");
   write_to_file(win, Example+"winDomain");
  std::cout << "\nWrite controller to " << Example << ".scs \n";
  if(write_to_file(scots::StaticController(ss,is,std::move(win)),Example))
    std::cout << "Done. \n";
    
#endif

#if RUN_INDEX == 2
//    std::vector<abs_type> Bpartition(Bsize);
//    win.get_winning_domain_MST(Bpartition);
//
//    std::vector<scots::abs_type> Bi, Bj;    // sparse(Bi(k),Bj(k)) = 1
//    std::vector<int> B_labels(Bsize);   // the labels of Bpartition as per dtc partitioning
//    createBiBjBlabels(ss, is, tf, Bsize, Bpartition, Bi, Bj, B_labels);
    TicToc writingB;
    writingB.tic();
//    wfile(tau, "T_row_b", 1, 0, Example);
//    wfile(row, "T_row_b", 2, 0, Example);
//    wfile(b, "T_row_b", 2, 0, Example);
    std::cout<<"Writing B_labels...";
    wfile(B_labels, "B_labels", 1, 0, Example);
    std::cout<<"Done.\nWriting Bi...";
    wfile(Bi, 1, 1, Example);
    // std::cout<<"Done.\nWriting Bj...";
    // wfile(Bj, "Bj", 1, 1, Example);
    std::cout<<"Done.\n";
    writingB.toc();
    
#endif
  if(Bsize == 0) return 1;
  return 0;
}



