#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bron_kerbosch_all_cliques.hpp>
#include "tools.hpp"
#include <stdlib.h>
#include <queue>
using namespace boost;

#ifndef UPGRADE_SCORE
#define UPGRADE_SCORE

typedef adjacency_list< vecS, vecS, undirectedS, property< vertex_color_t, int > > Graph;

typedef typename graph_traits< Graph >::vertex_descriptor vertex_desc;
typedef typename graph_traits< Graph >::vertex_iterator vertex_iter;
typedef graph_traits< Graph >::adjacency_iterator adj_iter;
typedef graph_traits< Graph >::edge_descriptor edge_desc;
typedef typename std::deque< vertex_desc > Clique;

// float end_score(Graph g, float score_min, std::vector<float>& score, int k, const char* name_out);

std::pair<float, float> score_sub_colo(Graph g, std::vector<float>& score, int k, const char* name_out);

std::pair<float, float> score_MIS(Graph g, std::vector<float>& score, int k, const char* name_out);

std::pair<float, float> score_Markov(Graph g, std::vector<float>& score, int& k, const char* name_out, bool max_mis = false);

int taille_clique_max_Bron(Graph g);

int taille_clique_max_Markov(Graph g);

// bool compute_one_color_inv_score(Graph g, int k, std::vector<std::vector<int>>& dl_color_inv_score, vertex_desc v, int i);

// int compute_dl_color_inv_score(Graph g, int k, std::vector<std::vector<int>>& dl_color_inv_score);


void upgrade_with_dl(Graph& g, int& k);

void upgrade_with_dl_v0(Graph& g, int& k);

bool check_sub_coloring(Graph g, int k);

void add_color_divide_biggest_clique(Graph& g, std::vector<float>& score, int& k);

#endif


