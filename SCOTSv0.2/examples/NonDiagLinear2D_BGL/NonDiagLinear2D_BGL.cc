

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

std::string MATLABDIR = "";

/* memory profiling */
#include <sys/time.h>
#include <sys/resource.h>
struct rusage usage;


#include <stdio.h>
#include <math.h>
#include <algorithm>

#include "Helper_MST_BGL.h"
#include "Helper_MST.h"

const std::string Example = "NonDiagLinear2D_BGL"; // Name of folder of example
/* state space dim */
const int state_dim = 2;
/* input space dim */
const int input_dim = 1;

/*
 * data types for the elements of the state space
 * and input space used by the ODE solver
 */
using state_type = std::array<double, state_dim>;
using input_type = std::array<double, input_dim>;

/* abbrev of the type for abstract states and inputs */
using abs_type = scots::abs_type;

#define ReachSet 0 // 0 for system_post, radius_post. 1 for reach_set
double eta_s = 0.2; //std::pow(10,-1);
double eta_i = 0.05;
double w_dist = 0.1; // disturbance [-w,w]

#if RUN_INDEX == 2
#include "dtControlClassify.h"
#endif

auto reach_set = [](state_type& r, state_type& x, const input_type& u) noexcept {
};

auto system_post = [](state_type& x, const input_type& u) noexcept {
	/* system dynamics */
	state_type z = x;
	x[0] = 2 * z[0] + z[1] + u[0];
	x[1] = 0.5 * z[1] - 0.4 * z[0] + u[0];
};

auto radius_post = [](state_type& r, const state_type&, const input_type& u, int CLongerU = 0) noexcept {
	/* the growth bound */
	state_type zz = r;
	r[0] = 2 * zz[0] + zz[1] + w_dist;
	r[1] = 0.5 * zz[1] - 0.4 * zz[0] + w_dist;
};

int main() {
	MATLABDIR = MATLABDIR + Example + "/";
	/* to measure time */
	TicToc tt;

	/* setup the workspace of the synthesis problem and the uniform grid */
	 /* grid node distance diameter */
	state_type eta = { {eta_s,eta_s} };
	/* lower bounds of the hyper-rectangle */
	state_type lb = { {-1,-2} };
	/* upper bounds of the hyper-rectangle */
	state_type ub = { {1,2} };
	scots::UniformGrid ss(state_dim, lb, ub, eta);
#if (RUN_INDEX==1)
	std::cout << "Uniform grid details:\n";
	ss.print_info();
#endif

	/* construct grid for the input space */
	/* lower bounds of the hyper rectangle */
	input_type i_lb = { {-1} };
	/* upper bounds of the hyper rectangle */
	input_type i_ub = { {1} };
	/* grid node distance diameter */
	input_type i_eta = { {eta_i} };
	scots::UniformGrid is(input_dim, i_lb, i_ub, i_eta);
#if (RUN_INDEX==1) 
	is.print_info();
#endif

	/* compute transition function of symbolic model */
	/* transition function of symbolic model */
	scots::TransitionFunction tf;
	scots::Abstraction<state_type, input_type> abs(ss, is);
	//abs.verbose_off();

	tt.tic();
#if (RUN_INDEX == 1)

	
	std::cout << "Computing the transition function:\n";
	abs.compute_gb(tf, system_post, radius_post);
	std::cout << "Number of transitions: " << tf.get_no_transitions() << "\n";

	if (!getrusage(RUSAGE_SELF, &usage))
		std::cout << "Memory per transition: " << usage.ru_maxrss / (double)tf.get_no_transitions() << "\n" << "Memory = " << usage.ru_maxrss << "\n";


	
	/*The win domain with eta_s=0.2 and eta_i=0.05 as the safe set.*/
	double toler = 0.0000000001;
	std::vector<double> llhsy{ 0.1, -0.3, -0.7, -1.1, -1.5, -1.7, -1.9 };
	std::vector<double> llhsx{ -0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3 };
	std::vector<double> urhsy{ -0.1, 0.3, 0.7, 1.1, 1.5, 1.7, 1.9 };
	std::vector<double> urhsx{ 0.9, 0.7, 0.5, 0.3, 0.1, -0.1, -0.3 };
	for (int i = 0; i < 7; i++) {
		llhsy[i] = llhsy[i] - toler;
		llhsx[i] = llhsx[i] - toler;
		urhsy[i] = urhsy[i] + toler;
		urhsx[i] = urhsx[i] + toler;
	}
	auto safeset = [&ss, &eta, &llhsy, &llhsx, &urhsy, &urhsx, &toler](const scots::abs_type& idx) noexcept {

		state_type x, cell_lb, cell_ub;
		ss.itox(idx, x);
		for (int i = 0; i < state_dim; i++) {
			cell_lb[i] = (x[i] - eta[i] / 2.0);
			cell_ub[i] = (x[i] + eta[i] / 2.0);
		}
		/* function returns 1 if cell associated with x is in target set  */
		if (-0.9 - toler <= cell_lb[0] &&
			cell_ub[0] <= 0.9 + toler &&
			-1.9 - toler <= cell_lb[1] &&
			cell_ub[1] <= 1.9 + toler &&
			cell_lb[1] >= -2 * cell_lb[0] - 1.7 - toler &&
			cell_ub[1] <= -2 * cell_ub[0] + 1.7 + toler) {
			int indt;
			for (int i = 0; i < 7; i++) {
				if (cell_lb[1] >= llhsy[i]) {
					indt = i;
					break;
				}
			}
			if (cell_lb[0] >= llhsx[indt]) {
				for (int i = 0; i < 7; i++) {
					if (cell_ub[1] <= urhsy[i]) {
						indt = i;
						break;
					}
				}
				if (cell_ub[0] <= urhsx[indt])
					return true;
			}
		}
		return false;
	};


/* compute winning domain (contains also valid inputs) */
	std::cout << "\nSynthesis: \n";
	scots::WinningDomain win = scots::solve_invariance_game(tf, safeset);
	abs_type Bsize = win.get_size();
	std::cout << "Winning domain size: " << Bsize << "\n";

#if (DTCONTROL == 0)
	std::vector<abs_type> Bpartition(Bsize);
	win.get_winning_domain_MST(Bpartition);   // winning states copied to Bpartition

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
#if (MINPOST == 0)
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



tt.toc();
#endif  // end RUN_INDEX == 1

#if (RUN_INDEX == 2)
	std::cout<< "\ndtControl number of partition elements = " << dtcontrol_nPElements << "\n";
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

	/*computes max edge weight*/
	double max_edge_weight = *std::max_element(edge_weights.get(), edge_weights.get() + num_edges);
	std::cout << "\nMaximum edge weight, wmax = " << max_edge_weight << "\n";

	/*constructs graph.*/
	constructDirectedGraph_b(Xj, B_NPosts, Bpartition, g, edge_weights, num_edges);
	/*prints the maximum cycle mean*/
	print_max_cycle_mean(g);

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

tt.toc();
#endif // end (DTCONTROL == 0 || RUN_INDEX == 2)


	if (!getrusage(RUSAGE_SELF, &usage))
		std::cout << "Memory = " << usage.ru_maxrss << "\n";


	if (Bsize == 0) return 1;
	return 0;
	}



