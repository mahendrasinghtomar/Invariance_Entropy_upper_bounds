/*
* 
*
*with DTCONTROL = 0:
*	Controller determinised with minPosts control selection.
*       (first such input selected for which no. of posts is minimum)
*	All the cells with same control value considered as a single
*	partition element.
*
*with SCOTS_FOCUS_ONLY_INTERSECTION = 0
*	Bypass computation of SCOTS-abstraction during the computation of
*	invariant controller with intersection-domain as the safe set.
*
*with SCOTS_FOCUS_ONLY_INTERSECTION = 1
*	avoid set defined as the set X\domain_intersection.
*
*
*/

#if (RUN_INDEX == 2)
#define HYPERBOLIC_STEP 2
#endif

#include <iostream>
#include <array>
#include <cmath>
#include <unordered_map>

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

typedef double douORfloa;	// double or float to troubleshoot numerical error with dtControl.

//typedef boost::numeric::interval<douORfloa> interval;
namespace bnumeric = boost::numeric;
namespace Ilib = boost::numeric::interval_lib;
typedef boost::numeric::interval<double, Ilib::policies<Ilib::save_state<Ilib::rounded_transc_std<double> >, Ilib::checking_base<double>>> interval;

douORfloa sup(interval& I) {
	return I.upper();
}

douORfloa inf(interval& I) {
	return I.lower();
}

douORfloa mag(interval& I) {
	return boost::numeric::norm(I);
}

douORfloa midpoint(interval& I) {
	return boost::numeric::median(I);
}

douORfloa radius(interval& I) {
	return 0.5 * boost::numeric::width(I);
}

#define SCOTS_FOCUS_ONLY_INTERSECTION 0 //if 0 then SCOTS-abstraction bypassed during longer length control sequence

const std::string ExampleFolder = "hyperbolic_v2_minPosts";
const std::string Example = "hyperbolic";
/* state space dim */
const int state_dim = 2;
/* input space dim */
const int input_dim_a = 2;    //actual input dimension
const int input_length = 1;

#define ReachSet 0 // 0 for system_post, radius_post. 1 for reach_set
douORfloa epsil = 0.08;
douORfloa eta_s = 0.009;     //std::pow(10,-2);
douORfloa x_rad = 5.9574 / 2.0;  // 2.9787
douORfloa eta_i = 0.01;

/*
 * data types for the elements of the state space
 * and input space used by the ODE solver
 */
using state_type = std::array<douORfloa, state_dim>;
//using input_type = std::array<douORfloa,input_dim>;
#if (HYPERBOLIC_STEP==1)
const int input_dim = input_dim_a;
using input_type = std::array<douORfloa, input_dim>;
const int uLength = 1;
#else   // (HYPERBOLIC_STEP==2)
const int input_dim = input_dim_a * input_length;
using input_type = std::array<douORfloa, input_dim>;
const int uLength = input_length;
#endif  // end HYPERBOLIC_STEP

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;


#if (RUN_INDEX == 2) && (DTCONTROL == 1)
#include "dtControlClassify.h"
#endif

auto reach_set = [](state_type& r, state_type& x, const input_type& u) noexcept {
	//    interval xt0(x[0]-r[0], x[0]+r[0]);
	//    interval xt1(x[1]-r[1], x[1]+r[1]);
	//    interval xx0, xx1;
	//
	//	for(abs_type i=0;i< uLength;i++)
	//	{
	//        abs_type u_begin = 2*i;	// the start index of the control vector for the current iteration		
	//		xx0 = 5.0 - 0.3 * xt1 - bnumeric::square(xt0) + u[u_begin];
	//		xx1 = xt0 + u[u_begin+1];
	//		xt0 = xx0;
	//		xt1 = xx1;
	//	}
	//	
	//    x[0] = midpoint(xx0);
	//    x[1] = midpoint(xx1);
	//    r[0] = radius(xx0);
	//    r[1] = radius(xx1);
};

auto reach_set_inv = [](state_type& r, state_type& x, const input_type& u) noexcept {
	//    interval xt0(x[0]-r[0], x[0]+r[0]);
	//    interval xt1(x[1]-r[1], x[1]+r[1]);
	//    interval xx0, xx1;
	//
	//	xx0 = xt1 - u[1];
	//	xx1 = (5.0 - bnumeric::square(xt1-u[1]) + u[0] - xt0 )/0.3;
	//	
	//    x[0] = midpoint(xx0);
	//    x[1] = midpoint(xx1);
	//    r[0] = radius(xx0);
	//    r[1] = radius(xx1);
};

//template<class Fs, class Fi>
auto system_post = [](state_type& x, const input_type& u) noexcept {
	state_type z = x;
	//x[0]= 5.0 - 0.3 * z[1] - z[0] * z[0] + u[0];
	//x[1] = z[0] + u[1];

	//int uLength = u.size()/input_dim_a;
	for (abs_type i = 0; i < uLength; i++)
	{
		abs_type u_begin = 2 * i;	// the start index of the control vector for the current iteration		
		x[0] = 5.0 - 0.3 * z[1] - z[0] * z[0] + u[u_begin];
		x[1] = z[0] + u[u_begin + 1];
		z = x;
	}
};

//template<class Fs, class Fi>
auto system_post_inv = [](state_type& x, const input_type& u) noexcept {
	state_type z = x;
	// x[0]= 5.0 - 0.3 * z[1] - z[0] * z[0] + u[0];
	// x[1] = z[0] + u[1];

	// for the inverse system
	x[0] = z[1] - u[1];
	x[1] = (5 - (z[1] - u[1]) * (z[1] - u[1]) + u[0] - z[0]) / 0.3;

};


//template<class Fs, class Fi>
auto radius_post = [](state_type& r, const state_type& x, const input_type& u) noexcept {
	state_type z = r;
	//r[0] = 0.3 * z[1] + z[0] * z[0] + 2 * std::abs(x[0]) * z[0];
	//r[1] = z[0];

	state_type zx = x;
	state_type xx;
	//abs_type uLength = u.size()/input_dim_a;
	for (abs_type i = 0; i < uLength; i++)
	{
		abs_type u_begin = 2 * i;	// the start index of the control vector for the current iteration		
		r[0] = 0.3 * z[1] + z[0] * z[0] + 2 * std::abs(zx[0]) * z[0]; // + u[u_begin];
		r[1] = z[0]; // + u[u_begin+1];
		z = r;

		xx[0] = 5.0 - 0.3 * zx[1] - zx[0] * zx[0] + u[u_begin];
		xx[1] = zx[0] + u[u_begin + 1];
		zx = xx;

	}
};

//template<class Fs, class Fi>
auto radius_post_inv = [](state_type& r, const state_type& x, const input_type& u) noexcept {
	state_type z = r;
	// r[0] = 0.3 * z[1] + z[0] * z[0] + 2 * std::abs(x[0]) * z[0];
	// r[1] = z[0];

	// for the inverse system
	r[0] = z[1];
	//r[1] = (z[1]*z[1] + 2*std::abs(x[1])*z[1] + 2*std::abs(u[1])*z[1] + z[0] )/0.3; commented on 15July2020
	r[1] = (z[1] * z[1] + 2 * std::abs(x[1] - u[1]) * z[1] + z[0]) / 0.3;
};


int main() {
	MATLABDIR = MATLABDIR + ExampleFolder + "/";
	/* to measure time */
	TicToc tt;
	scots::TransitionFunction* tf_fwd_pt;    // pointer to forward system transition function
	std::vector<abs_type>* domain_intersection_pt;  // pointer to the domain_intersection
	scots::UniformGrid* ss_pt; // pointer to state space grid
	state_type* lb_pt, * ub_pt, * eta_pt;

	/* state space grid */
	/* grid node distance diameter */
	state_type eta = { {eta_s, eta_s} };
	/* lower bounds of the hyper-rectangle */
	state_type lb = { {-x_rad, -x_rad} };
	/* upper bounds of the hyper-rectangle */
	state_type ub = { {x_rad, x_rad} };
	lb_pt = &lb;
	ub_pt = &ub;
	eta_pt = &eta;
	scots::UniformGrid ss(state_dim, lb, ub, eta);
	ss_pt = &ss;

#if (HYPERBOLIC_STEP==1) 
	{
#if (RUN_INDEX==1)
		std::vector<abs_type> domain_fwd, domain_inv;

		std::cout << "Uniform grid details (state space):\n";
		ss.print_info();

		/* construct grid for the input space */
		input_type i_lb, i_ub, i_eta;
		for (int i = 0; i < input_dim; i++) {
			i_lb[i] = -epsil;
			i_ub[i] = epsil;
			i_eta[i] = eta_i;
		}

		scots::UniformGrid is(input_dim, i_lb, i_ub, i_eta);
		is.print_info();


		/* compute transition function of symbolic model */
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
#endif  // end ReachSet 
			std::cout << "Number of transitions: " << tf_fwd.get_no_transitions() << "\n";
			tf_fwd_pt = &tf_fwd;
#if NOTVISUALSTUDIO
			if (!getrusage(RUSAGE_SELF, &usage))
				std::cout << "Memory per transition: " << usage.ru_maxrss / (douORfloa)tf_fwd.get_no_transitions() << "\n" << "Memory = " << usage.ru_maxrss << "\n";
#endif  // end NOTVISUALSTUDIO
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
			std::cout << "Forward-time winning domain size: " << Bsize << "\n";
			domain_fwd = win.get_winning_domain();
			// std::cout << "\nWrite controller to " << Example << "_fwd.scs \n";
			//if(write_to_file(scots::StaticController(ss,is,std::move(win)),Example+"_fwd"))
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
				domain_intersection_pt = &domain_intersection;
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
				if (!scots::read_from_file(win_temp, "hyperbolic_timeInv")) {
					std::cout << "\nFailed to read the controller\n";
					return 1;
				}
				domain_inv = win_temp.get_winning_domain();

				scots::Abstraction<state_type, input_type> abs(ss, is);
				std::vector<abs_type> domain_intersection = domain_inv;
				//rfile(MATLABDIR + "DataMatlab/S_dtCon_domain_intersection.txt", domain_intersection);
				domain_intersection_pt = &domain_intersection;
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
				if (!scots::read_from_file(win_temp2, "hyperbolic_timeFwd")) {
					std::cout << "\nFailed to read the controller\n";
					return 1;
				}
				domain_fwd = win_temp2.get_winning_domain();

			}
			forwar = !forwar;
			iteration_FwdBwd_count++;
		}
		std::cout << "Iteration forward-backward-domain count = " << iteration_FwdBwd_count << "\n";
		std::remove("hyperbolic_timeFwd.scs");
		std::remove("hyperbolic_timeInv.scs");

		//            {
		//                /* for inverse system */
		//                /* compute transition function of symbolic model */
		//                /* transition function of symbolic model */
		//                scots::TransitionFunction tf;
		//                scots::Abstraction<state_type, input_type> abs(ss, is);
		//                //abs.verbose_off();
		//
		//                tt.tic();
		//                std::cout << "\nComputing the transition function for inverse system:\n";
		//#if (ReachSet == 0)
		//                abs.compute_gb(tf, system_post_inv, radius_post_inv);
		//#else
		//                abs.compute_gb2(tf, system_post_inv, radius_post_inv, reach_set_inv);
		//#endif
		//                std::cout << "Number of transitions: " << tf.get_no_transitions() << "\n";
		//#if NOTVISUALSTUDIO
		//                if (!getrusage(RUSAGE_SELF, &usage))
		//                    std::cout << "Memory per transition: " << usage.ru_maxrss / (douORfloa)tf.get_no_transitions() << "\n" << "Memory = " << usage.ru_maxrss << "\n";
		//#endif
		//                tt.toc();
		//
		//                /* continue with synthesis */
		//                /* define function to check if the cell is in the safe set  */
		//                auto safeset = [&lb, &ub, &ss, &eta](const scots::abs_type& idx) noexcept {
		//                    state_type x;
		//                    ss.itox(idx, x);
		//                    /* function returns 1 if cell associated with x is in target set  */
		//                    if ((lb[0] <= (x[0] - eta[0] / 2.0) && (x[0] + eta[0] / 2.0) <= ub[0]) && (lb[1] <= (x[1] - eta[1] / 2.0) && (x[1] + eta[1] / 2.0) <= ub[1]))
		//                        return true;
		//                    return false;
		//                };
		//                /* compute winning domain (contains also valid inputs) */
		//                std::cout << "\nSynthesis for inverse sys: \n";
		//                tt.tic();
		//                scots::WinningDomain win = scots::solve_invariance_game(tf, safeset);
		//                tt.toc();
		//                abs_type Bsize = win.get_size();
		//                std::cout << "Winning domain size: " << Bsize << "\n";
		//
		//                // std::vector<abs_type> winningDomain(Bsize);
		//                // win.get_winning_domain_MST(winningDomain);
		//                // wfile(winningDomain, "WinningDomain.txt");
		//                domain_inv = win.get_winning_domain();
		//                // write_to_file(win, Example+"winDomain");
		//                // std::cout << "\nWrite controller to " << Example << "_inv.scs \n";
		//                // if(write_to_file(scots::StaticController(ss,is,std::move(win)),Example+"_inv"))
		//                  // std::cout << "Done. \n";
		//
		//
		//        }




				/*unsigned long fwd_size, inv_size;
				fwd_size = domain_fwd.size();
				inv_size = domain_inv.size();
				for (unsigned long i = 0; i < fwd_size; i++)
				{
					for (unsigned long j = 0; j < inv_size; j++)
					{
						if (domain_fwd[i] == domain_inv[j])
						{
							domain_intersection.push_back(domain_fwd[i]);
							break;
						}
					}
				}*/

		std::cout << "\nAll time invariant set size: " << domain_inv.size() << "\n";
		wfile(domain_inv, "allTimeInvariantSet", 1, 0, Example);
		//wfile_pFacesSafeset(ss , eta, domain_intersection);
		if (domain_inv.size() == 0) return 1;
		domain_intersection_pt = &domain_inv;
#endif	// end RUN_INDEX==1
	}
#endif  // end HYPERBOLIC_STEP==1



#if (HYPERBOLIC_STEP==2)
	abs_type Bsize;
	{
		/* With intersection of fwd and inv domains as the safe set*/
		/* construct grid for the input space */
		input_type i_lb, i_ub, i_eta;
		for (int i = 0; i < input_dim; i++) {
			i_lb[i] = -epsil;
			i_ub[i] = epsil;
			i_eta[i] = eta_i;
		}
		scots::UniformGrid is(input_dim, i_lb, i_ub, i_eta);
		//#if (RUN_INDEX==1)
		std::cout << "\nInput grid\n";
		is.print_info();
		//#endif

		  /* compute transition function of symbolic model */
		  /* transition function of symbolic model */
		scots::Abstraction<state_type, input_type> abs(*ss_pt, is);

#if (RUN_INDEX == 1)
		// /*If input length!=1, only then tf_fwd recomputed.*/
		// if (input_length != 1) 


		std::cout << "Computing the invariant controller for the forward system with higher dimensional input space:\n";
#if (SCOTS_FOCUS_ONLY_INTERSECTION == 1)
		{
			tt.tic();

			std::vector<abs_type> domain_intersection;
			rfile(MATLABDIR + "DataMatlab/S_dtCon_allTimeInvariantSet.txt", domain_intersection);
			domain_intersection_pt = &domain_intersection;
			abs_type Dsize = domain_intersection.size();
			abs_type xGridSize = ss.size();
			std::unique_ptr<bool[]> cells_outside_intersection(new bool[xGridSize]()); //avoid set for compute_gb
			for (abs_type i = 0; i < xGridSize; i++)
				cells_outside_intersection[i] = true;
			for (abs_type i = 0; i < Dsize; i++)
				cells_outside_intersection[domain_intersection[i]] = false;
			auto avoid_lambda = [&cells_outside_intersection](const abs_type& xid) noexcept {
				return cells_outside_intersection[xid]; };

			scots::TransitionFunction tf_fwd;
#if (ReachSet == 0)
			abs.compute_gb(tf_fwd, system_post, radius_post, avoid_lambda);
#else
			abs.compute_gb2(tf_fwd, system_post, radius_post, reach_set, avoid_lambda);
#endif  // end ReachSet
			std::cout << "Number of transitions: " << (tf_fwd).get_no_transitions() << "\n";
			// tf_fwd_pt = &tf_fwd;
#if NOTVISUALSTUDIO
			if (!getrusage(RUSAGE_SELF, &usage))
				std::cout << "Memory = " << usage.ru_maxrss << "\n";
#endif  // end NOTVISUALSTUDIO

			std::vector<bool> still_safe(Dsize, true);  // initialized to true
	  /* continue with synthesis */
	  /* define function to check if the cell is in the safe set  */
			std::unordered_map<abs_type, abs_type> domain_intersection_xind_to_ind; // maps state-index to the index in the domain_intersection vector.
			for (abs_type i = 0; i < Dsize; i++)
				domain_intersection_xind_to_ind[domain_intersection[i]] = i;
			auto safeset = [&domain_intersection_pt, &still_safe, &domain_intersection_xind_to_ind](const scots::abs_type& idx) noexcept {
				/*unsigned long d_int_size = (*domain_intersection_pt).size();
				for(unsigned long i=0; i<d_int_size; i++)
				{
					if (idx == (*domain_intersection_pt)[i])
					{
						if (still_safe[i])
							return true;
						else
							return false;
					}
				}
				return false;*/
				if (domain_intersection_xind_to_ind.count(idx) > 0)
				{
					if (still_safe[domain_intersection_xind_to_ind[idx]])
						return true;
					else
						return false;
				}
				else
					return false;
			};
			/* compute winning domain (contains also valid inputs) */
			std::cout << "\nSynthesis for intersection: \n";
			tt.tic();
			scots::WinningDomain win = scots::solve_invariance_game(tf_fwd, safeset);
			Bsize = win.get_size();
			std::cout << "Winning domain size: " << Bsize << "\n";
			write_to_file(win, Example + "with_intersection_winDomain");
			std::cout << "\nWrite controller to " << Example << ".scs \n";
			if (write_to_file(scots::StaticController(ss, is, std::move(win)), Example))
				std::cout << "Done. \n";

			tt.toc();
		}

#else // SCOTS_FOCUS_ONLY_INTERSECTION==0

		/*my invariance game*/
		lb_pt = NULL;
		ub_pt = NULL;
		eta_pt = NULL;
		std::vector<abs_type> domain_intersection;
		rfile(MATLABDIR + "DataMatlab/S_dtCon_allTimeInvariantSet.txt", domain_intersection);
		domain_intersection_pt = &domain_intersection;
		abs_type Dsize = domain_intersection.size();

		std::vector<bool> still_safe(Dsize, true);  // initialized to true
		//bool* input_list_status = new bool[Dsize * is.size()];
		std::unique_ptr<bool[]> input_list_status = std::make_unique<bool[]>(Dsize * is.size());
		/* continue with synthesis */
		tt.tic();
		synthesis_invariance_mst(domain_intersection, still_safe, input_list_status, Bsize, is, ss, abs, system_post, radius_post, reach_set);
		tt.toc();
		std::cout << "Winning domain size: " << Bsize << "\n";
		wfile_winDomain(still_safe, domain_intersection, input_list_status, ss, is, Example + "with_intersection_winDomain.scs", 1);
		std::cout << "\nWrite controller to " << Example << ".scs \n";
		wfile_controller_scs(still_safe, domain_intersection, input_list_status, ss, is, Example + ".scs");
		std::cout << "Done. \n";

#endif // SCOTS_FOCUS_ONLY_INTERSECTION




#if 0
		//// begin 09 March 2021
			/*Checks for controlled invariance, for the backward-time system, of the domain from hyperbolic.scs.
			An alternative: compute invariant controller with domain from hyperbolic.scs as the state-space,
			then compare the size of the domain with the size of the winning domain of the computed controller.*/
		{
			input_type i_lb_temp, i_ub_temp, i_eta_temp;
			for (int i = 0; i < input_dim; i++) {
				i_lb_temp[i] = -epsil;
				i_ub_temp[i] = epsil;
				i_eta_temp[i] = eta_i;
			}
			scots::UniformGrid is_temp(input_dim, i_lb_temp, i_ub_temp, i_eta_temp);
			//#if (RUN_INDEX==1)
			std::cout << "\nInput grid\n";
			is_temp.print_info();
			scots::Abstraction<state_type, input_type> abs_temp(ss, is_temp);
			scots::WinningDomain win_temp;
			if (!scots::read_from_file(win_temp, "hyperbolic")) {
				std::cout << "\nFailed to read the controller\n";
				return 1;
			}
			std::vector<abs_type> domain_temp = win_temp.get_winning_domain();
			state_type x;
			input_type u;
			for (const abs_type& xind : domain_temp) {
				ss.itox(xind, x);
				abs_type M = is_temp.size();
				for (abs_type j = 0; j < M; j++) {
					bool invariantInput = true;
					/* current input */
					is_temp.itox(j, u);
					std::vector<state_type> posts = abs_temp.get_post(system_post_inv, radius_post_inv, x, u);
					//std::vector<state_type> posts = abs_temp.get_post(system_post, radius_post, x, u);
					for (const state_type& apost : posts) {
						abs_type apostind = ss.xtoi(apost);
						if (!win_temp.is_winning(apostind)) {
							/*The post is outside the domain_temp, go to next input.*/
							invariantInput = false;
							if (j == M - 1) {
								/*It was the last element in our considered control input set.*/
								std::cout << "\nThe state-cell index =" << xind << " has no input in our considered input set for invariance of  the backward-time system.\n";
								return 1;   // Show error.
							}
							break;
							//std::cout << "\nPost of state = " << xind << " (" << x[0] << ", " << x[1] << ") under the dtcontrol assigned input = " << is.xtoi(u) << " (" << u[0] << ", " << u[1] << ") went outside the controller domain\n";
						}
					}
					if (invariantInput == true) {
						/*The control input is admissible for the state-cell. Go to next state-cell.*/
						break;
					}
				}

			}
			std::cout << "\nThe domain is controlled invariant for the time-backward system.\n";
		}
		//// end 09 March 2021
#endif




#if MINPOST == 0
 /* write_to_file(win, Example + "winDomain");
  std::cout << "\nWrite controller to " << Example << ".scs \n";
  if (write_to_file(scots::StaticController(ss, is, std::move(win)), Example))
	  std::cout << "Done. \n";*/
#else   // MINPOST == 1
  /*Controller updated to contain only those inputs which have the smallest number of posts.*/

#if 0   // Switch read of winDomain between SCOTS-read and mine.
		scots::WinningDomain win;
		if (!read_from_file(win, Example + "with_intersection_winDomain")) {
			std::cout << "\nFailed to read the file containing winning domain.\n";
			return 1;   // When couldn't read file, 'make' will interrupt and show error
		}
		Bsize = win.get_size();
		std::vector<abs_type> Bpartition(Bsize);
		win.get_winning_domain_MST(Bpartition);
#else
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

		std::vector<abs_type> Bpartition;
		std::vector<std::vector<abs_type>> domain_input_list;
		if (!read_winDomain_mst(Example + "with_intersection_winDomain.scs", domain_input_list, Bpartition)) {
			std::cout << "\nFailed to read the file: " << Example + "with_intersection_winDomain.scs\n";
			return 1; // Error signal for makefile.
		}
		Bsize = Bpartition.size();
		winDomain_mst win(Bpartition, domain_input_list);

#endif  // end 0


		//std::vector<scots::abs_type> Bpartition(Bsize);
		//win.get_winning_domain_MST(Bpartition);
		std::vector<std::vector<scots::abs_type>> controllerMinPost(Bsize);
		std::vector<int> B_NPosts(Bsize);
		int int_ReachSet = (ReachSet == 1) ? 1 : 0;
		minPost_controller(int_ReachSet, system_post, radius_post, reach_set, win, ss, is, abs, Bpartition, Bsize, controllerMinPost, B_NPosts);
		wfile(B_NPosts, "B_NPosts", 1, 0, Example);
		wfile_winDomain_b(ss, is, Bsize, Bpartition, controllerMinPost, Example + "winDomain.scs", 1);
		wfile_controller_scs_b(ss, is, Bsize, Bpartition, controllerMinPost, Example + ".scs");
#endif  // end MINPOST

#endif   // end RUN_INDEX == 1

#if (RUN_INDEX == 2)   // RUN_INDEX == 2 :
		tt.tic();
		std::cout << "Computing Bi, Bj, Blabels\n";

#if 0
		scots::WinningDomain win_r;
		if (!read_from_file(win_r, Example + "with_intersection_winDomain")) {
			std::cout << "\nFailed to read the file containing winning domain.\n";
			return 1;   // When couldn't read file, 'make' will interrupt and show error
		}
		Bsize = win_r.get_size();
		std::vector<abs_type> Bpartition(Bsize);
		win_r.get_winning_domain_MST(Bpartition);
#else
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

		std::vector<abs_type> Bpartition;
		std::vector<std::vector<abs_type>> domain_input_list;
		if (!read_winDomain_mst(Example + "with_intersection_winDomain.scs", domain_input_list, Bpartition)) {
			std::cout << "\nFailed to read the file: " << Example + "with_intersection_winDomain.scs\n";
			return 1; // Error signal for makefile.
		}
		Bsize = Bpartition.size();
		winDomain_mst win_r(Bpartition, domain_input_list);

#endif  // end 0



#if (DTCONTROL == 0)
		std::vector<std::vector<abs_type>> Xj(Bsize);
		std::vector<abs_type> DControllerMinPost(Bsize); // deterministic controller with first such input selected for which no. of posts is minimum
		std::vector<int> B_NPosts(Bsize);   // the number of post states for states in Bpartition

		int int_ReachSet = (ReachSet == 1) ? 1 : 0;
		DeterminiseController_Xj_NPosts(int_ReachSet, system_post, radius_post, reach_set, win_r, ss, is, abs, Bpartition, Bsize, Xj, DControllerMinPost, B_NPosts);
		wfile_DController(Bpartition, DControllerMinPost, Example, 1, 0);
		std::vector<int> Bpartition_coarsePno(Bsize, 0);    // initialized to zero. Coarse partition index taken to start from 1. Same purpose as B_labels. 
		std::vector<int> CoarsePartition;
		compute_coarsePartition(DControllerMinPost, CoarsePartition, Bpartition_coarsePno);

		//write_to_file(win, Example + "winDomain");
		//std::cout << "\nWrite controller to " << Example << ".scs \n";
		//if (write_to_file(scots::StaticController(ss, is, std::move(win_r)), Example))
		//    std::cout << "Done. \n";

		TicToc writingB;
		writingB.tic();
		std::cout << "Writing B_labels...";
		wfile(Bpartition_coarsePno, "B_labels", 1, 0, Example);
		std::cout << "Done.\nWriting Bi...";
		wfile(Xj, Bpartition);
		std::cout << "Done.\n";
		writingB.toc();
#endif // DTCONTROL == 0

#if (DTCONTROL == 1)
		std::vector<int> B_labels(Bsize);   // the labels of Bpartition as per dtc partitioning
		std::vector<std::vector<abs_type>> Bi(Bsize);

		scots::TransitionFunction tf_fwd_d; // dummy
#if (ReachSet == 0)

		if (!abs.compute_BiBjBlabels(tf_fwd_d, system_post, radius_post, Bpartition, Bsize, Bi, B_labels, classify))
			return 1;
#else
		if (!abs.compute_BiBjBlabels2(tf_fwd_d, system_post, radius_post, reach_set, Bpartition, Bsize, Bi, B_labels, classify))
			return 1;
#endif  // end ReachSet
		tt.toc();

		TicToc writingB;
		writingB.tic();
		std::cout << "Writing B_labels...";
		wfile(B_labels, "B_labels", 1, 0, Example);
		std::cout << "Done.\nWriting Bi...";
		wfile(Bi, 1, 1, Example);
		std::cout << "Done.\n";
		writingB.toc();
#endif // end DTCONTROL == 1

#endif  // end RUN_INDEX == 2
	}
#if NOTVISUALSTUDIO
	if (!getrusage(RUSAGE_SELF, &usage))
		std::cout << "Memory = " << usage.ru_maxrss << "\n";
#endif  // END NOTVISUALSTUDIO


	if (Bsize == 0) return 1;
	return 0;
#endif  // end HYPERBOLIC_STEP==2
}
