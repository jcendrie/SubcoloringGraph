#include <iostream> // for std::cout
#include <fstream>
#include <math.h>
#include <boost/graph/adjacency_list.hpp>

#ifndef GRAPH_GENERATION
#define GRAPH_GENERATION

using namespace boost;

typedef adjacency_list< vecS, vecS, undirectedS, property< vertex_color_t, int > > Graph;
typedef typename graph_traits< Graph >::vertex_descriptor vertex_desc;
typedef typename graph_traits< Graph >::vertex_iterator vertex_iter;
typedef graph_traits< Graph >::adjacency_iterator adj_iter;
typedef graph_traits< Graph >::edge_descriptor edge_desc;

Graph rand_graph_bino(int n, float p);

Graph rand_graph_bino_density(int n, float d);

Graph rand_UDG(int n, float taille_cube);

Graph rand_UDG_density(int n, float d);

Graph rand_quasi_UDG(int n, float taille_cube, float dist);

Graph rand_quasi_UDG_density(int n, float d, float dist);

Graph rand_3D(int n, float taille_cube, float dist);

Graph rand_ubg_density(int n, float d);

Graph rand_3D_density(int n, float d, float dist);

Graph rand_sbm(int n, int r, float p, float q);

Graph rand_sbm_density(int n, int r, float rpq, float d);

Graph stadium(int n, float l, float larg);

Graph stadiumV2(int n, float taille_cube);

Graph stadium_density(int n, float d, float dist);

Graph read_coord(std::string name);

Graph etoile(int b);

#endif

