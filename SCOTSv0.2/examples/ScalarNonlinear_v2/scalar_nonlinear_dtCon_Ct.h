/*  mst: continuous time sampled system
*
* With decimal in dtcontrol, use 'double' in state_type and input_type,
* while with float in dtcontrol, use 'float'.
*
*/

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

// #define NOTVISUALSTUDIO 1   //0: visual studio. 1: else than visual studio

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
std::string OSname = "Windows environment";
std::string MATLABDIR = "C:/Users/..../SCOTSv0.2/examples/ScalarNonlinear_v2/";
#define NOTVISUALSTUDIO 0
#elif __APPLE__
std::string OSname = "Apple Mac";
// Mac
std::string MATLABDIR = "/Users/..../SCOTSv0.2/examples/ScalarNonlinear_v2/";
#define NOTVISUALSTUDIO 1
#elif __linux__
std::string OSname = "Linux environment";
std::string MATLABDIR = "/mnt/c/Users/..../SCOTSv0.2/examples/ScalarNonlinear_v2/";
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

#include "Helper_MST.h"

const std::string Example = "scalar_nonlinear_dtCon_Ct";
/* state space dim */
const int state_dim=1;
/* input space dim */
const int input_dim_a=1;	//actual input dimension
const int input_length=4;
const int input_dim=input_dim_a * input_length; // extended input dimension for longer length inputs

/*
 * data types for the elements of the state space
 * and input space used by the ODE solver
 */
using state_type = std::array<double,state_dim>;
using input_type = std::array<double,input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;


#define ReachSet 0 // 0 for system_post, radius_post. 1 for reach_set
const double tau = 0.01;
const double row = 1;  // where row < b*b+1
const double b = 1;
double eta_s = std::pow(10,-6);
double eta_i = 2*row*0.1;

#if RUN_INDEX == 2
    #include "dtControlClassify.h"
#endif

auto reach_set = [](state_type &r, state_type &x, const input_type &u ) noexcept {
//    interval xt(x[0]-r[0], x[0]+r[0]);
//    interval xt2;
////    xt2 = xt + T * sqrt(1+b*b) * sin(atan2(1, b) - 2.0 * xt) + T * u[0] * cos(xt) * cos(xt);
//    xt2 = xt + T * (-2*b*sin(xt)*cos(xt)-sin(xt)*sin(xt)+cos(xt)*cos(xt)+u[0]*cos(xt)*cos(xt));
//    x[0] = midpoint(xt2);
//    r[0] = radius(xt2);
};

auto system_post = [](state_type &x, const input_type &u) noexcept {
   auto rhs =[](state_type& xx,  const state_type &x, const double &u) noexcept {
         xx[0]=-2*b*sin(x[0])*cos(x[0])-sin(x[0])*sin(x[0])+cos(x[0])*cos(x[0]) + u * cos(x[0])*cos(x[0]);
       };
    for(int i=0;i<input_length;i++)
        scots::runge_kutta_fixed4(rhs,x,u[i],state_dim,tau,5);
   };

auto radius_post = [](state_type &r, const state_type&, const input_type &u) noexcept {
	for(int i=0;i<input_length;i++)
		r[0] = exp((2 * sqrt(1 + b * b) + abs(u[i])) * tau) * r[0];

//    /* the ode for the growth bound */
//    auto rhs =[](state_type& rr,  const state_type &r, const double &u) noexcept {
//        rr[0]=(2*sqrt(1+b*b) + abs(u))*r[0];
////        rr[0] = std::abs(-2*b-4-u[0])*r[0];
//      };
//    for(int i=0;i<input_dim;i++)
//        scots::runge_kutta_fixed4(rhs,r,u[i],state_dim,tau,5);
};

int main() {
  /* to measure time */
  TicToc tt;

  /* setup the workspace of the synthesis problem and the uniform grid */
   /* grid node distance diameter */
  state_type eta={{eta_s}};
    
 /* lower bounds of the hyper-rectangle */
  state_type lb={{atan(-b-sqrt(b*b+1+row))}};
 /* upper bounds of the hyper-rectangle */
  state_type ub={{atan(-b-sqrt(b*b+1-row))}};
    
//      state_type lb= {{-1.1071}};
//      state_type ub= {{1.1071}};

    
//  state_type lb={{-1+atan(-b-sqrt(b*b+1+row))}}; // Q enlarged, margin of 1
//  state_type ub={{1+atan(-b-sqrt(b*b+1-row))}};
    
  scots::UniformGrid ss(state_dim,lb,ub,eta);
  #if (RUN_INDEX==1)
  std::cout << "Uniform grid details:\n";
  ss.print_info();
  #endif
  
  /* construct grid for the input space */
  /* lower bounds of the hyper rectangle */
//  input_type i_lb={{-row}};
//  //input_type i_lb={{-1}};
//  /* upper bounds of the hyper rectangle */
//  input_type i_ub={{row}};
//  //input_type i_ub={{1}};
//  /* grid node distance diameter */
//  input_type i_eta={{eta_i}};
  //input_type i_eta={{0.1}};
    input_type i_lb, i_ub, i_eta;
    for(int i=0;i<input_length;i++){
        i_lb[i] = -row;
        i_ub[i] = row;
        i_eta[i] = eta_i;
    }
 
  /* construct grid for the input alphabet */ 
  scots::UniformGrid is(input_dim,i_lb,i_ub,i_eta);
  #if (RUN_INDEX==1)
  is.print_info();
  #endif

  /* compute transition function of symbolic model */
  // std::cout << "Computing the transition function:\n";
  /* transition function of symbolic model */
  scots::TransitionFunction tf;
  scots::Abstraction<state_type,input_type> abs(ss,is);
  // abs.verbose_off();

  
#if (RUN_INDEX == 1)
   tt.tic();
    std::cout << "Computing the transition function:\n";
    #if (ReachSet == 0)
        abs.compute_gb(tf,system_post, radius_post);
    #else
        abs.compute_gb2(tf,system_post, radius_post, reach_set);
    #endif
	std::cout << "Number of transitions: " << tf.get_no_transitions() <<"\n";
    #if NOTVISUALSTUDIO
      if(!getrusage(RUSAGE_SELF, &usage))
        std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << "\n" << "Memory = " << usage.ru_maxrss << "\n";
    #endif
   tt.toc();

  /* continue with synthesis */
  /* define function to check if the cell is in the safe set  */
  auto safeset = [&lb, &ub, &ss, &eta](const scots::abs_type& idx) noexcept {
    state_type x;
    ss.itox(idx,x);
    /* function returns 1 if cell associated with x is in target set  */
    if (lb[0] <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0)<= ub[0]) //&& lb[1] <= (x[1]-eta[1]/2.0) && (x[1]+eta[1]/2.0)<= ub[1])
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
  write_to_file(win, Example + "winDomain");
  std::cout << "\nWrite controller to " << Example << ".scs \n";
  if (write_to_file(scots::StaticController(ss, is, std::move(win)), Example))
      std::cout << "Done. \n";
#else   // MINPOST == 1
  /*Controller updated to contain only those inputs which have the smallest number of posts.*/
  std::vector<scots::abs_type> Bpartition(Bsize);
  win.get_winning_domain_MST(Bpartition);
  std::vector<std::vector<scots::abs_type>> controllerMinPost(Bsize);
  std::vector<int> B_NPosts(Bsize);
  int int_ReachSet = (ReachSet == 1) ? 1 : 0;
  minPost_controller(int_ReachSet, system_post, radius_post, reach_set, win, ss, is, abs, Bpartition, Bsize, controllerMinPost, B_NPosts);
  wfile(B_NPosts, "B_NPosts", 1, 0, Example);
  wfile_winDomain_b(ss, is, Bsize, Bpartition, controllerMinPost, Example + "winDomain.scs", 1);
  wfile_controller_scs_b(ss, is, Bsize, Bpartition, controllerMinPost, Example + ".scs");
#endif  // end MINPOST
#endif  //end RUN_INDEX == 1

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

#if (ReachSet == 0)

  if (!abs.compute_BiBjBlabels(tf, system_post, radius_post, Bpartition, Bsize, Bi, B_labels, classify))
      return 1; // post went outside the winning domain
#else
  abs.compute_gb2(tf, system_post, radius_post, reach_set);
#endif

    // std::vector<abs_type> Bpartition(Bsize);
    // win.get_winning_domain_MST(Bpartition);
    
    // std::vector<scots::abs_type> Bi, Bj;    // sparse(Bi(k),Bj(k)) = 1
    // std::vector<int> B_labels(Bsize);   // the labels of Bpartition as per dtc partitioning
    // createBiBjBlabels(ss, is, tf, Bsize, Bpartition, Bi, Bj, B_labels);
    TicToc writingB;
    writingB.tic();
    wfile(tau, "T_row_b", 1, 0, Example);
    wfile(row, "T_row_b", 2, 0, Example);
    wfile(b, "T_row_b", 2, 0, Example);
    std::cout<<"Writing B_labels...";
    wfile(B_labels, "B_labels", 1, 0, Example);
    std::cout<<"Done.\nWriting Bi...";
	wfile(Bi, 1, 1, Example);
    // wfile(Bi, "Bi", 1, 1, Example);
    // std::cout<<"Done.\nWriting Bj...";
    // wfile(Bj, "Bj", 1, 1, Example);
    std::cout<<"Done.\n";
    writingB.toc();
    
#endif

	#if NOTVISUALSTUDIO
      if(!getrusage(RUSAGE_SELF, &usage))
        std::cout <<  "Memory = " << usage.ru_maxrss << "\n";
    #endif

  if(Bsize == 0) return 1;
  return 0;
}



