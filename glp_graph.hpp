#include <glpk.h>
#include <stack>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>

#ifndef GLP_GRAPH
#define GLP_GRAPH

using namespace boost;

typedef adjacency_list< vecS, vecS, undirectedS, property< vertex_color_t, int > > Graph;

typedef typename graph_traits< Graph >::vertex_descriptor vertex_desc;
typedef typename graph_traits< Graph >::vertex_iterator vertex_iter;
typedef typename graph_traits< Graph >::adjacency_iterator adj_iter;
typedef typename graph_traits< Graph >::edge_descriptor edge_desc;


void calcul_P3(Graph g, std::stack<std::vector<int>>& l, int& nb);

void calcul_P3_init(Graph g, std::stack<std::vector<int>>& l, int& nb);

bool glpk_solve_graph(glp_prob* mip, int k, int n, int nb, std::stack<std::vector<int>> l);

bool glpk_solve_graph_k_color(Graph& g, int k);

int verif_sol(glp_prob* mip, int k, int n, int nb, std::stack<std::vector<int>> l);

bool glpk_solve_graph_k_color_by_fragment(Graph& g, int k);

bool glpk_solve_minimize_clique_k_color(Graph& g, int k, int min_d);

bool glpk_solve_minimize_conflicts(Graph& g, int k, int min_conflict);

// Test ss-coloriage pour toygraphe avec 2 couleurs !   
void glpk_toy();
void glpk_toy_v2();

#endif
