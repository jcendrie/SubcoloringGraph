#include <iostream> // for std::cout
#include <deque>
#include <queue>
#include <stack>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>
#include "tools.hpp"
#include "upgrade_score.hpp"

#ifndef COLORING
#define COLORING

using namespace boost;

typedef adjacency_list< vecS, vecS, undirectedS, property< vertex_color_t, int > > Graph;
typedef typename graph_traits< Graph >::vertex_descriptor vertex_desc;
typedef typename graph_traits< Graph >::vertex_iterator vertex_iter;
typedef graph_traits< Graph >::adjacency_iterator adj_iter;
typedef graph_traits< Graph >::edge_descriptor edge_desc;
typedef typename std::deque< vertex_desc > Clique;


// void max_uncolor_degree(const Graph& g, vertex_desc& max_vertex, int& max_degree);

// void max_degree(const Graph& g, int& max_degree);

// void sub_graph(Graph& neighbor_graph, vertex_desc center, Graph g);

// void color_and_erase(Graph& g, Clique k, int c, int& uncolored_node);

// void clean(Graph& g, int& uncolored_node);

// void betweenness(const Graph& g, std::vector<float>& cb, vertex_desc& max_vertex, float& max_betweenness);

// void insertion_tri(const Graph& g, vertex_desc vertex, std::vector<float> cb, std::vector<vertex_desc>& ltri);

// void insertion_tri(const Graph& g, vertex_desc vertex, std::vector<int> aud, std::vector<vertex_desc>& ltri);

// void find_clique_max_betweenness(const Graph& g, std::vector<float> cb, vertex_desc center, Clique& k);

// void max_uncolor_degree(const Graph& g, std::vector<int>& actual_uncolor_degree, vertex_desc& max_vertex, int& max_degree);

// void find_clique_max_degree(const Graph& g, std::vector<int>& actual_uncolor_degree, vertex_desc center, Clique& k);

// void maj_uncolor_degree() {

int subcoloring(Graph& g, std::string mode, bool normal = false);

int channel_by_arrival(Graph& g, int k_max);

int coloring_by_Markov(Graph& g, float threshold);

int naive_then_subcoloring(Graph& g, float threshold, int k_naive, int k_max);

int extract_low_degree_vertex(Graph& g, Graph& g_ori, std::stack<vertex_desc>& stack_low_degree, int k_max); // Coloration en k_max+1 des sommets de degrée non colorée < à k_max de manière récursive.

int back_low_degree_vertex(Graph& g, Graph& g_ori, std::stack<vertex_desc>& stack_low_degree, int k_max);

#endif

