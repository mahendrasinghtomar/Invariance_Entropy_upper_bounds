//
//  Helper_MST_BGL.h
//  SCOTSv02MST
//
//  Created by MST on 06.07.20.
//  Copyright © 2020 MST. All rights reserved.
//

#ifndef Helper_MST_BGL_h
#define Helper_MST_BGL_h

#include <boost/config.hpp>
#include <iostream> // for std::cout
#include <utility> // for std::pair
#include <algorithm> // for std::for_each
#include <boost/utility.hpp> // for boost::tie
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/howard_cycle_ratio.hpp>
#include <boost/graph/strong_components.hpp>

// create a typedef for the Graph type
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property,
    boost::property< boost::edge_weight_t, double, boost::property<boost::edge_index_t, int > > >
    Graph;

template<class abs_type>
bool createDirectedGraph(const std::vector<std::vector<abs_type>>& Xj, const std::vector<int>& B_NPosts, const std::vector<abs_type>& Bpartition, Graph& g) {
    /*Creates directed graph, and computes and prints the maximum cycle mean. */
    int num_vertices = Bpartition.size();
    int num_edges = 0;
    /*for (int i = 0; i < B_NPosts.size(); i++)
        num_edges = num_edges + B_NPosts[i];*/
    for (const auto& e : B_NPosts)
        num_edges += e;
    double* edge_weights = new double[num_edges];

    typedef std::pair<abs_type, abs_type> Edge;
    Edge* edge_array = new Edge[num_edges];
    abs_type count = 0;
    for (abs_type i = 0; i < Xj.size(); i++) {
        for (abs_type j = 0; j < Xj[i].size(); j++) {
            edge_array[count] = Edge(Bpartition[i], Xj[i][j]);
            edge_weights[count] = std::log2(B_NPosts[i]);
            count++;
        }
    }
    g = Graph(edge_array, edge_array + num_edges, edge_weights, num_vertices);
    boost::property_map< Graph, boost::edge_weight_t >::type ewm
        = boost::get(boost::edge_weight, g);
    boost::property_map< Graph, boost::vertex_index_t >::type vim
        = boost::get(boost::vertex_index, g);
    boost::property_map< Graph, boost::edge_index_t >::type eim
        = boost::get(boost::edge_index, g);
    double max_cr = boost::maximum_cycle_mean(g, vim, ewm, eim);
    std::cout << "\nMaximum cycle mean = " << max_cr << "\n";

    delete[] edge_weights;
    delete[] edge_array;
    return true;
}

bool constructDirectedGraph(const std::vector<std::vector<scots::abs_type>>& Xj, const std::vector<int>& B_NPosts, const std::vector<scots::abs_type>& Bpartition, Graph& g, std::unique_ptr<double[]>& edge_weights) {
    /*Constructs directed graph.*/
    int num_vertices = Bpartition.size();
    int num_edges = 0;
    for (const auto& e : B_NPosts)
        num_edges += e;
    //double* edge_weights = new double[num_edges];

    typedef std::pair<scots::abs_type, scots::abs_type> Edge;
    Edge* edge_array = new Edge[num_edges];
    scots::abs_type count = 0;
    for (scots::abs_type i = 0; i < Xj.size(); i++) {
        for (scots::abs_type j = 0; j < Xj[i].size(); j++) {
            edge_array[count] = Edge(Bpartition[i], Xj[i][j]);
            //edge_weights[count] = std::log2(B_NPosts[i]);
            count++;
        }
    }
    g = Graph(edge_array, edge_array + num_edges, edge_weights.get(), num_vertices);
    
    //delete[] edge_weights;
    delete[] edge_array;
    return true;
}

bool constructDirectedGraph_b(const std::vector<std::vector<scots::abs_type>>& Xj, const std::vector<int>& B_NPosts, const std::vector<scots::abs_type>& Bpartition, Graph& g, std::unique_ptr<double[]>& edge_weights, const scots::abs_type& num_edges) {
    /*Constructs directed graph.*/
    int num_vertices = Bpartition.size();
    /*int num_edges = 0;
    for (const auto& e : B_NPosts)
        num_edges += e;*/
    //double* edge_weights = new double[num_edges];

    typedef std::pair<scots::abs_type, scots::abs_type> Edge;
    Edge* edge_array = new Edge[num_edges];
    scots::abs_type count = 0;
    for (scots::abs_type i = 0; i < Xj.size(); i++) {
        for (scots::abs_type j = 0; j < Xj[i].size(); j++) {
            edge_array[count] = Edge(Bpartition[i], Xj[i][j]);
            //edge_weights[count] = std::log2(B_NPosts[i]);
            count++;
        }
    }
    g = Graph(edge_array, edge_array + num_edges, edge_weights.get(), num_vertices);

    //delete[] edge_weights;
    delete[] edge_array;
    return true;
}

void print_max_cycle_mean(Graph& g) {
    boost::property_map< Graph, boost::edge_weight_t >::type ewm
        = boost::get(boost::edge_weight, g);
    boost::property_map< Graph, boost::vertex_index_t >::type vim
        = boost::get(boost::vertex_index, g);
    boost::property_map< Graph, boost::edge_index_t >::type eim
        = boost::get(boost::edge_index, g);
    double max_cr = boost::maximum_cycle_mean(g, vim, ewm, eim);
    std::cout << "\nMaximum cycle mean = " << max_cr << "\n";
}

//void compute_wmax_overSCCs(Graph& g) {
//    /*Computes the maximum edge weight in any scc of g. Assumption: for each vertex all the outgoing edges have the same weight.*/
//    boost::property_map< Graph, boost::edge_weight_t >::type edgeWeights
//        = boost::get(boost::edge_weight, g);
//    std::map<boost::graph_traits<Graph>::vertex_descriptor, scots::abs_type> compMap;
//    boost::associative_property_map< std::map<boost::graph_traits<Graph>::vertex_descriptor, scots::abs_type>> componentMap(compMap);
//    scots::abs_type Nscc = boost::strong_components(g, componentMap);
//    std::vector<double> scc_wmax(Nscc,0);
//    boost::graph_traits<Graph>::out_edge_iterator ei, ei_end;
//    for (const auto& e : compMap) {
//        for (scots::abs_type i = 0; i < Nscc; i++) {
//            if (e.second == i) {
//                //if()
//                boost::tie(ei, ei_end) = out_edges(e.first, g);
//                //if(scc_wmax)
//
//                //std::cout<< "Out edge weight = " << boost::out_edges
//            }
//        }
//    }
//
//    std::cout << "\n";
//    for (const auto& e : compMap)
//        std::cout << "e.first = " << e.first << ", e.second = " << e.second << "\n";
//    std::cout << "\n";
//
//}

template<class abs_type>
void writeGraph_dot(const std::vector<abs_type>& Bpartition, Graph& g){
	int num_vertices = Bpartition.size();
	std::string* name = new std::string[num_vertices];
    for (int i = 0; i < num_vertices; i++)
        name[i] = std::to_string(Bpartition[i]);
	
	/* It's the second time ewm is computed. First during maxCycleMean.*/
	boost::property_map< Graph, boost::edge_weight_t >::type ewm
        = boost::get(boost::edge_weight, g);
		
    std::ofstream gvz_file;
    gvz_file.open("gvz_file.txt");
    boost::write_graphviz(gvz_file, g, boost::make_label_writer(name),
        boost::make_label_writer(ewm));
    gvz_file.close();
	delete[] name;
}


#endif /* Helper_MST_BGL_h */
