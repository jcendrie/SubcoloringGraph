#include <iostream> // for std::cout
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphml.hpp>
#include <vector>

#ifndef TOOLS
#define TOOLS

using namespace boost;

typedef adjacency_list< vecS, vecS, undirectedS, property< vertex_color_t, int > > Graph;
typedef typename graph_traits< Graph >::vertex_descriptor vertex_desc;
typedef typename graph_traits< Graph >::vertex_iterator vertex_iter;
typedef graph_traits< Graph >::adjacency_iterator adj_iter;
typedef graph_traits< Graph >::edge_descriptor edge_desc;

Graph initGraph();

void initGraph(Graph& g, const char* name_in);

void writeGraph(Graph& g, const char* name_out);

void write_Results(char* real_name, int n, float threshold_sub, int k_max_naive, int k_max, std::vector<std::vector<float>> stat);

void write_result(char* real_name, int k, std::vector<float> stat, int n, int m,  bool normal);

void displayGraph(Graph g, bool visible_edge = false);

void blankgraph(Graph& g);

void statColor(Graph g, int k);

int computeLibertyDegree(Graph g, int k, std::vector<std::vector<bool>>& dlColor);

void insertion_tri(std::vector<float>& list, float elem);

#endif

