/*
* edited 10 March 2021: Iteration for all-time invariant set.
*
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

#if (DTCONTROL == 0 || RUN_INDEX == 2)
#include <lemon/time_measure.h>
#include <lemon/smart_graph.h>
#include <lemon/static_graph.h>
#include <lemon/hartmann_orlin_mmc.h>
#include <lemon/karp_mmc.h>
#endif

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
#include <boost/numeric/interval.hpp>
#include <algorithm>

#include "Helper_MST.h"
#include "Helper_MST_BGL.h"

//typedef boost::numeric::interval<double> interval;
namespace bnumeric = boost::numeric;
namespace Ilib = boost::numeric::interval_lib;
typedef boost::numeric::interval<double, Ilib::policies<Ilib::save_state<Ilib::rounded_transc_std<double> >, Ilib::checking_base<double>>> interval;

double sup(interval& I) {
	return I.upper();
}

double inf(interval& I) {
	return I.lower();
}

double mag(interval& I) {
	return boost::numeric::norm(I);
}

double midpoint(interval& I) {
	return boost::numeric::median(I);
}

double radius(interval& I) {
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
double x_rad = 5.9574 / 2.0;
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

auto reach_set = [](state_type& r, state_type& x, const input_type& u) noexcept {
	interval xt0(x[0] - r[0], x[0] + r[0]);
	interval xt1(x[1] - r[1], x[1] + r[1]);
	interval xx0, xx1;

	xx0 = 5.0 - 0.3 * xt1 - bnumeric::square(xt0) + u[0];
	xx1 = xt0 + u[1];

	x[0] = midpoint(xx0);
	x[1] = midpoint(xx1);
	r[0] = radius(xx0);
	r[1] = radius(xx1);
};

auto reach_set_inv = [](state_type& r, state_type& x, const input_type& u) noexcept {
	interval xt0(x[0] - r[0], x[0] + r[0]);
	interval xt1(x[1] - r[1], x[1] + r[1]);
	interval xx0, xx1;

	/* for the inverse system */
	xx0 = xt1 - u[1];
	xx1 = (5.0 - bnumeric::square(xt1 - u[1]) + u[0] - xt0) / 0.3;

	x[0] = midpoint(xx0);
	x[1] = midpoint(xx1);
	r[0] = radius(xx0);
	r[1] = radius(xx1);
};

auto system_post = [](state_type& x, const input_type& u) noexcept {
	state_type z = x;

	x[0] = 5.0 - 0.3 * z[1] - z[0] * z[0] + u[0];
	x[1] = z[0] + u[1];

};

auto system_post_inv = [](state_type& x, const input_type& u) noexcept {
	state_type z = x;

	// for the inverse system
	x[0] = z[1] - u[1];
	x[1] = (5 - (z[1] - u[1]) * (z[1] - u[1]) + u[0] - z[0]) / 0.3;

};


auto radius_post = [](state_type& r, const state_type& x, const input_type& u) noexcept {
	state_type z = r;

	r[0] = 0.3 * z[1] + z[0] * z[0] + 2 * std::abs(x[0]) * z[0];
	r[1] = z[0];

};

auto radius_post_inv = [](state_type& r, const state_type& x, const input_type& u) noexcept {
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
	TicToc tt, h_step1, h_step2, timePostAllTimeInvariantSet;
	

	/* setup the workspace of the synthesis problem and the uniform grid */
	 /* grid node distance diameter */
	state_type eta = { {eta_s, eta_s} };

	/* lower bounds of the hyper-rectangle */
	state_type lb = { {-x_rad, -x_rad} };
	/* upper bounds of the hyper-rectangle */
	state_type ub = { {x_rad, x_rad} };

	scots::UniformGrid ss(state_dim, lb, ub, eta);
#if (HYPERBOLIC_STEP==1)
	std::cout << "Uniform grid details:\n";
	ss.print_info();
#endif

	/* construct grid for the input space */
	input_type i_lb, i_ub, i_eta;
	for (int i = 0; i < input_dim; i++) {
		i_lb[i] = -epsil;
		i_ub[i] = epsil;
		i_eta[i] = eta_i;
	}

	scots::UniformGrid is(input_dim, i_lb, i_ub, i_eta);
#if (HYPERBOLIC_STEP==1)
	is.print_info();
#endif


#if (HYPERBOLIC_STEP==1)
	h_step1.tic();
	{


		/* compute transition function of symbolic model */
		/* transition function of symbolic model */
		scots::TransitionFunction tf_fwd;

		/* for fwd system */
		{
			scots::Abstraction<state_type, input_type> abs(ss, is);
			//abs.verbose_off();


			tt.tic();
			std::cout << "Computing the transition function:\n";
#if (ReachSet == 0)
			abs.compute_gb(tf_fwd, system_post, radius_post);
#else
			abs.compute_gb2(tf_fwd, system_post, radius_post, reach_set);
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


			std::cout << "Number of transitions: " << tf_fwd.get_no_transitions() << "\n";
#if NOTVISUALSTUDIO
			if (!getrusage(RUSAGE_SELF, &usage))
				std::cout << "Memory per transition: " << usage.ru_maxrss / (double)tf_fwd.get_no_transitions() << "\n" << "Memory = " << usage.ru_maxrss << "\n";
#endif	// NOTVISUALSTUDIO
			tt.toc();




			/* continue with synthesis */
			/* define function to check if the cell is in the safe set  */
			auto safeset = [&lb, &ub, &ss, &eta](const scots::abs_type& idx) noexcept {
				state_type x;
				ss.itox(idx, x);
				/* function returns 1 if cell associated with x is in target set  */
				if ((lb[0] <= (x[0] - eta[0] / 2.0) && (x[0] + eta[0] / 2.0) <= ub[0]) && (lb[1] <= (x[1] - eta[1] / 2.0) && (x[1] + eta[1] / 2.0) <= ub[1]))
					return true;
				return false;
			};
			/* compute winning domain (contains also valid inputs) */
			std::cout << "\nSynthesis: \n";
			tt.tic();
			scots::WinningDomain win = scots::solve_invariance_game(tf_fwd, safeset);
			tt.toc();
			abs_type Bsize = win.get_size();
			std::cout << "Winning domain size: " << Bsize << "\n";
			domain_fwd = win.get_winning_domain();
			// std::cout << "\nWrite controller to " << Example << "_fwd.scs \n";
		   // if(write_to_file(scots::StaticController(ss,is,std::move(win)),Example+"_fwd"))
			 // std::cout << "Done. \n";

			if (Bsize == 0) return 1;
		}

		/*Iteration with previous winning domain as the safe set;
			alternates between forward-time and backward-time system;
			Uses synthesis_invariance_mst() instead of SCOTS::solve_invariance_game(). */
		bool forwar = false;
		int iteration_FwdBwd_count = 1;
		bool needFurther_FwdBwd_iteration = true;
		while (needFurther_FwdBwd_iteration) {

			if (!forwar) {  //time inverse system
				scots::Abstraction<state_type, input_type> abs(ss, is);
				std::vector<abs_type> domain_intersection = domain_fwd;
				//rfile(MATLABDIR + "DataMatlab/S_dtCon_domain_intersection.txt", domain_intersection);
				abs_type Dsize = domain_intersection.size();
				abs_type Bsize;

				std::vector<bool> still_safe(Dsize, true);  // initialized to true
				//bool* input_list_status = new bool[Dsize * is.size()];
				std::unique_ptr<bool[]> input_list_status = std::make_unique<bool[]>(Dsize * is.size());
				/* continue with synthesis */
				tt.tic();
				synthesis_invariance_mst(domain_intersection, still_safe, input_list_status, Bsize, is, ss, abs, system_post_inv, radius_post_inv, reach_set_inv);
				tt.toc();
				std::cout << "Inverse-time winning domain size: " << Bsize << "\n";
				//wfile_winDomain(still_safe, domain_intersection, input_list_status, ss, is, Example + "_timeInv_winDomain.scs", 1);
				std::cout << "\nWrite controller to " << Example << "_timeInv.scs \n";
				wfile_controller_scs(still_safe, domain_intersection, input_list_status, ss, is, Example + "_timeInv.scs");
				std::cout << "Done. \n";
				if (Bsize == Dsize)
					needFurther_FwdBwd_iteration = false;
			}
			else {  // time forward system
				scots::WinningDomain win_temp;
				if (!scots::read_from_file(win_temp, Example + "_timeInv")) {
					std::cout << "\nFailed to read the controller\n";
					return 1;
				}
				domain_inv = win_temp.get_winning_domain();

				scots::Abstraction<state_type, input_type> abs(ss, is);
				std::vector<abs_type> domain_intersection = domain_inv;
				//rfile(MATLABDIR + "DataMatlab/S_dtCon_domain_intersection.txt", domain_intersection);
				abs_type Dsize = domain_intersection.size();
				abs_type Bsize;

				std::vector<bool> still_safe(Dsize, true);  // initialized to true
				//bool* input_list_status = new bool[Dsize * is.size()];
				std::unique_ptr<bool[]> input_list_status = std::make_unique<bool[]>(Dsize * is.size());
				/* continue with synthesis */
				tt.tic();
				synthesis_invariance_mst(domain_intersection, still_safe, input_list_status, Bsize, is, ss, abs, system_post, radius_post, reach_set);
				tt.toc();
				std::cout << "Forward-time winning domain size: " << Bsize << "\n";
				//wfile_winDomain(still_safe, domain_intersection, input_list_status, ss, is, Example + "_timeFwd_winDomain.scs", 1);
				std::cout << "\nWrite controller to " << Example << "_timeFwd.scs \n";
				wfile_controller_scs(still_safe, domain_intersection, input_list_status, ss, is, Example + "_timeFwd.scs");
				std::cout << "Done. \n";

				/*Can be replaced with a for loop that works with still_safe and Bsize.*/
				scots::WinningDomain win_temp2;
				if (!scots::read_from_file(win_temp2, Example + "_timeFwd")) {
					std::cout << "\nFailed to read the controller\n";
					return 1;
				}
				domain_fwd = win_temp2.get_winning_domain();

			}
			forwar = !forwar;
			iteration_FwdBwd_count++;
		}
		std::cout << "Iteration forward-backward-domain count = " << iteration_FwdBwd_count << "\n";
		std::remove((Example + "_timeFwd.scs").c_str());
		std::remove((Example + "_timeInv.scs").c_str());


		//	{
		//  /* for inverse system */
		//  /* compute transition function of symbolic model */
		//  /* transition function of symbolic model */
		//  scots::TransitionFunction tf;
		//  scots::Abstraction<state_type,input_type> abs(ss,is);
		//  //abs.verbose_off();
		//
		//
		//	tt.tic();
		//    std::cout << "\nComputing the transition function for inverse system:\n";
		//    #if (ReachSet == 0)
		//        abs.compute_gb(tf,system_post_inv, radius_post_inv);
		//    #else
		//        abs.compute_gb2(tf,system_post_inv, radius_post_inv, reach_set_inv);
		//    #endif
		//	std::cout << "Number of transitions: " << tf.get_no_transitions() <<"\n";
		//    #if NOTVISUALSTUDIO
		//      if(!getrusage(RUSAGE_SELF, &usage))
		//        std::cout << "Memory per transition: " << usage.ru_maxrss/(double)tf.get_no_transitions() << "\n"<< "Memory = " << usage.ru_maxrss << "\n";
		//    #endif
		//  tt.toc();
		//
		//  
		//
		//  /* continue with synthesis */
		//  /* define function to check if the cell is in the safe set  */
		//  auto safeset = [&lb, &ub, &ss, &eta](const scots::abs_type& idx) noexcept {
		//    state_type x;
		//    ss.itox(idx,x);
		//    /* function returns 1 if cell associated with x is in target set  */
		//    if ((lb[0] <= (x[0]-eta[0]/2.0) && (x[0]+eta[0]/2.0) <= ub[0]) && (lb[1] <= (x[1]-eta[1]/2.0) && (x[1]+eta[1]/2.0) <= ub[1]))
		//      return true;
		//    return false;
		//  };
		//  /* compute winning domain (contains also valid inputs) */
		//  std::cout << "\nSynthesis for inverse sys: \n";
		//  tt.tic();
		//  scots::WinningDomain win = scots::solve_invariance_game(tf,safeset);
		//  tt.toc();
		//  abs_type Bsize = win.get_size();
		//  std::cout << "Winning domain size: " << Bsize << "\n";
		//
		//    // std::vector<abs_type> winningDomain(Bsize);
		//    // win.get_winning_domain_MST(winningDomain);
		//    // wfile(winningDomain, "WinningDomain.txt");
		//  domain_inv = win.get_winning_domain();
		//  // write_to_file(win, Example+"winDomain");
		//  // std::cout << "\nWrite controller to " << Example << "_inv.scs \n";
		//  // if(write_to_file(scots::StaticController(ss,is,std::move(win)),Example+"_inv"))
		//    // std::cout << "Done. \n";
		//    
		//
		//}


		std::cout << "\nAll time invariant set size: " << domain_inv.size() << "\n";
		wfile(domain_inv, "allTimeInvariantSet", 1, 0, Example);
		//wfile_pFacesSafeset(ss , eta, domain_intersection);
		if (domain_inv.size() == 0) return 1;
	}
	h_step1.toc();

#else // (HYPERBOLIC_STEP==2)
		
		h_step2.tic();
		abs_type Bsize;
		std::vector<abs_type> domain_intersection;
		rfile(MATLABDIR + "DataMatlab/S_dtCon_allTimeInvariantSet.txt", domain_intersection);
		/* With intersection of fwd and inv domains as the safe set*/
		scots::Abstraction<state_type, input_type> abs(ss, is);
		//abs.verbose_off();

		std::cout << "Computing the invariant controller for the forward system with higher dimensional input space:\n";

	  /* continue with synthesis */
		std::vector<abs_type> Bpartition;

#if (RUN_INDEX == 1)

#if 0	// 1: SCOTS invariance game, 0: my invariance game

		/*SCOTS invariance game*/
		scots::TransitionFunction tf_fwd;
		tt.tic();
		std::cout << "Computing the transition function:\n";
#if (ReachSet == 0)
		abs.compute_gb(tf_fwd, system_post, radius_post);
#else
		abs.compute_gb2(tf_fwd, system_post, radius_post, reach_set);
#endif
		std::cout << "Number of transitions: " << tf_fwd.get_no_transitions() << "\n";
#if NOTVISUALSTUDIO
		if (!getrusage(RUSAGE_SELF, &usage))
			std::cout << "Memory per transition: " << usage.ru_maxrss / (double)tf_fwd.get_no_transitions() << "\n" << "Memory = " << usage.ru_maxrss << "\n";
#endif
		tt.toc();


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
		scots::WinningDomain win = scots::solve_invariance_game(tf_fwd, safeset);
		tt.toc();
		Bsize = win.get_size();
		std::cout << "All-time invariant winning domain size: " << Bsize << "\n";
		Bpartition.resize(Bsize);
		win.get_winning_domain_MST(Bpartition);   // winning states copied to Bpartition
#else	
	  // my invariance game
		abs_type Dsize = domain_intersection.size();

		std::vector<bool> still_safe(Dsize, true);  // initialized to true
		//bool* input_list_status = new bool[Dsize * is.size()];
		std::unique_ptr<bool[]> input_list_status = std::make_unique<bool[]>(Dsize * is.size());
		/* continue with synthesis */
		tt.tic();
		synthesis_invariance_mst(domain_intersection, still_safe, input_list_status, Bsize, is, ss, abs, system_post, radius_post, reach_set);
		tt.toc();
		std::cout << "Winning domain size: " << Bsize << "\n";
		Bpartition.resize(Bsize);

#if 1	
		/*To get winDomain data structure: write then read controller file.*/
		wfile_winDomain(still_safe, domain_intersection, input_list_status, ss, is, Example + "with_intersection_winDomain.scs", 1);
		std::cout << "\nWrite controller to " << Example << ".scs \n";
		wfile_controller_scs(still_safe, domain_intersection, input_list_status, ss, is, Example + ".scs");
		std::cout << "Done. \n";
		scots::WinningDomain win;
		if (!read_from_file(win, Example + "with_intersection_winDomain")) {
			std::cout << "\nFailed to read the file containing winning domain.\n";
			return 1;   // When couldn't read file, 'make' will interrupt and show error
		}
		win.get_winning_domain_MST(Bpartition);
#else
		/*To get winDomain data structure: define a new class winDomain_mst.*/

		class winDomain_mst {
			std::unordered_map<abs_type, std::vector<abs_type>> storedmap;
		public:
			winDomain_mst(const std::vector<abs_type>& domain, const std::vector<std::vector<abs_type>>& input_list) {
				for (abs_type i = 0; i < domain.size(); i++)
					storedmap[domain[i]] = input_list[i];
			}
			std::vector<abs_type> get_inputs(abs_type xind) {
				return storedmap[xind];
			}
		};
		
		std::vector<std::vector<abs_type>> domain_input_list(Bsize);
		{
			auto M = is.size();
			abs_type count = 0;
			for (scots::abs_type i = 0; i < still_safe.size(); i++)
				if (still_safe[i]) {
					Bpartition[count] = domain_intersection[i];
					for (scots::abs_type j = 0; j < M; j++)
					{
						if (input_list_status[M * i + j])
							domain_input_list[count].push_back(j);
					}
					count++;
				}
		}
		winDomain_mst win(Bpartition, domain_input_list);
#endif // end To get winDomain data structure:...

#endif // end SCOTS vs my invariance game

#if (DTCONTROL == 0)
	std::vector<std::vector<abs_type>> Xj(Bsize);
	std::vector<abs_type> DControllerMinPost(Bsize); // deterministic controller with first such input selected for which no. of posts is minimum
	std::vector<int> B_NPosts(Bsize);   // the number of post states for states in Bpartition

	int int_ReachSet = (ReachSet == 1) ? 1 : 0;
	DeterminiseController_Xj_NPosts(int_ReachSet, system_post, radius_post, reach_set, win, ss, is, abs, Bpartition, Bsize, Xj, DControllerMinPost, B_NPosts);

	std::unordered_map<abs_type, abs_type> Bpartition2ind;
	for (abs_type i = 0; i < Bsize; i++)
		Bpartition2ind[Bpartition[i]] = i;
#endif  // end DTCONTROL == 0

#if (DTCONTROL == 1)    // Save controller to .scs file in RunIndex 1.
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
#endif  // end DTCONTROL = 1

#endif  // end RUN_INDEX == 1

#if RUN_INDEX == 2
	std::cout << "\ndtControl number of partition elements = " << dtcontrol_nPElements << "\n";
	scots::WinningDomain win_r;
	if (!read_from_file(win_r, Example + "winDomain")) {
		std::cout << "\nFailed to read the file containing winning domain.\n";
		return 1;   // When couldn't read file, 'make' will interrupt and show error
	}
	Bsize = win_r.get_size();
	//std::vector<abs_type> Bpartition(Bsize);
	Bpartition.resize(Bsize);
	win_r.get_winning_domain_MST(Bpartition);
	std::vector<int> B_labels(Bsize);   // the labels of Bpartition as per dtc partitioning
	std::vector<std::vector<abs_type>> Bi(Bsize);
	scots::TransitionFunction tf_temp;
	abs.compute_BiBjBlabels(tf_temp, system_post, radius_post, Bpartition, Bsize, Bi, B_labels, classify);
	std::vector<std::vector<abs_type>> Ai(dtcontrol_nPElements);
	/*Smallest number in Bi is 0, in Ai is 0, in B_labels is 1. */
	for (abs_type i = 0; i < Bsize; i++) {
		for (abs_type j = 0; j < Bi[i].size(); j++) {
			if (Ai[B_labels[i] - 1].end() == std::find(Ai[B_labels[i] - 1].begin(), Ai[B_labels[i] - 1].end(), B_labels[Bi[i][j]] - 1))
				Ai[B_labels[i] - 1].push_back(B_labels[Bi[i][j]] - 1);
		}
	}
	std::vector<int> A_NPosts(dtcontrol_nPElements);
	for (abs_type i = 0; i < dtcontrol_nPElements; i++)
		A_NPosts[i] = Ai[i].size();
	std::vector<abs_type> Apartition(dtcontrol_nPElements);
	abs_type count_temp = 0;
	for (abs_type& e : Apartition) {
		e = count_temp;
		count_temp++;
	}
	auto B_NPosts = A_NPosts;
	auto Xj = Ai; // This Xj (unlike the one from RUN_INDEX = 1, which has elements with indexes same as the SCOTS grid) has elements in range 0 to Apartition.size().
	Bpartition = Apartition;

	/* Here Bpartition2ind is redundant. Introduced to mimic Bpartition2ind introduced because of Xj (when DTCONTROL = 0).
	   Defined to be identity, since Bi stores values in range 0 to Bsize-1, and not same as SCOTS grid.*/
	std::unordered_map<abs_type, abs_type> Bpartition2ind;
	for (abs_type i = 0; i < Bsize; i++)
		Bpartition2ind[i] = i;

#endif  // end RUN_INDEX == 2

#if (DTCONTROL == 0 || RUN_INDEX == 2)
#if 0 // BGL
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

	///*computes max edge weight*/
	// //double max_edge_weight = *std::max_element(edge_weights.get(), edge_weights.get() + num_edges);  // Slower than the below for loop.
	// double max_edge_weight = 0;
	// for(auto it=edge_weights.get(); it!= edge_weights.get()+num_edges; it++)
	//	 if(*it>max_edge_weight)
	//		 max_edge_weight = *it;
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
#else   // Lemon
	abs_type num_edges = 0;
	for (const auto& e : B_NPosts)
		num_edges += e;
	lemon::SmartDigraph g;
	auto num_nodes = Bpartition.size();
	g.reserveNode(num_nodes);
	g.reserveArc(num_edges);
	std::vector<lemon::SmartDigraph::Node> nodes(num_nodes);
	for (auto& e : nodes)
		e = g.addNode();
	lemon::SmartDigraph::ArcMap<double> edge_weights(g);
	std::vector<lemon::SmartDigraph::Arc> edges(num_edges);

	abs_type count = 0;
	for (scots::abs_type i = 0; i < Xj.size(); i++) {
		for (scots::abs_type j = 0; j < Xj[i].size(); j++) {
			/*edges[count] = g.addArc(Bpartition[i], Xj[i][j]);*/
			edges[count] = g.addArc(nodes[i], nodes[Bpartition2ind[Xj[i][j]]]);   // Bpartition2ind needed with Xj but not with Bi.
			edge_weights[edges[count]] = -std::log2(B_NPosts[i]); // With negative sign in edge_weights.
			count++;
		}
	}

	//#################################################
#if (MAXEDGEWT_MCM == 1)	// 1=max edge weight, 0=MCM
/*Max edge weight*/
	double maxedgewt = 0;
	for (auto& edge : edges) {
		if (edge_weights[edge] < maxedgewt)
			maxedgewt = edge_weights[edge];	// stores the most negative
	}
	std::cout << "\nMaximum edge weight = " << -maxedgewt << "\n";

	// auto max_edge_it =  std::min_element(edges.begin(), edges.end(),
	// 	[&edge_weights](lemon::SmartDigraph::Arc e1, lemon::SmartDigraph::Arc e2) {return (edge_weights[e1] < edge_weights[e2]); });
	// std::printf("\nMaximum edge weight = %f\n", -edge_weights[*max_edge_it]);

#else
//#################################################
/*Max cycle mean*/
/*{
	lemon::TimeReport trH("Running time of HartmannOrlinMmc: ");
	lemon::HartmannOrlinMmc<lemon::SmartDigraph, lemon::SmartDigraph::ArcMap<double>> Smmc(g, edge_weights);
	Smmc.findCycleMean();
	std::cout << "HartmannOrlin Maximum cycle mean = " << -Smmc.cycleMean() << std::endl;
}*/
	{
		//lemon::TimeReport trH("Running time of KarpMmc: ");
		lemon::KarpMmc<lemon::SmartDigraph, lemon::SmartDigraph::ArcMap<double>> Kmmc(g, edge_weights);
		Kmmc.findCycleMean();
		std::cout << "\nKarp Maximum cycle mean = " << -Kmmc.cycleMean() << std::endl;
	}
#endif	// MAXEDGEWT_MCM
	//#################################################


#endif // end BGL/Lemon

	std::cout << "Hyperbolic STEP2 time ";
	h_step2.toc();

#endif // end (DTCONTROL == 0 || RUN_INDEX == 2)



#if NOTVISUALSTUDIO
	if (!getrusage(RUSAGE_SELF, &usage))
		std::cout << "Memory = " << usage.ru_maxrss << "\n";
#endif
	
	if (Bsize == 0) return 1;
	return 0;
#endif // HYPERBOLIC_STEP==2
}
