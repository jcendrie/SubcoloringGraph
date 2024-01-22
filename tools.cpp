#include "tools.hpp"

using namespace boost;

Graph initGraph() {
    const int num_vertices = 9;

    // writing out the edges in the graph
    typedef std::pair< int, int > Edge;
    Edge edge_array[] = {
        Edge(0, 1),
        Edge(1, 2),
        Edge(2, 3),
        Edge(2, 6),
        Edge(2, 8),
        Edge(3, 4),
        Edge(3, 6),
        Edge(4, 5),
        Edge(4, 7),
        Edge(5, 7),
        Edge(6, 7),
        Edge(6, 9),
        Edge(8, 9)
    };

    const int num_edges = sizeof(edge_array) / sizeof(edge_array[0]);

    std::vector<int> base_color(num_edges, 0);

    Graph g(num_vertices);
    property_map< Graph, vertex_color_t >::type colormap = get(vertex_color, g);
    for (std::size_t j = 0; j < num_edges; ++j)
    {
        graph_traits< Graph >::edge_descriptor e;
        bool inserted;
        tie(e, inserted)
            = add_edge(edge_array[j].first, edge_array[j].second, g);
    }
    std::pair< vertex_iter, vertex_iter > vp_temp;
    std::size_t j = 0;
    for (vp_temp = vertices(g); vp_temp.first != vp_temp.second; ++vp_temp.first) {
        colormap[*vp_temp.first] = base_color[j];
        j++;      
    }

    return g;
}

void initGraph(Graph& g, const char* name_in) {
    std::ifstream ifs;
    ifs.open(name_in);
    dynamic_properties dp;
    read_graphml(ifs, g, dp);
    ifs.close();
}

void writeGraph(Graph& g, const char* name_out) {
    std::ofstream ofs;
    ofs.open(name_out);    
    dynamic_properties dp;
    write_graphml(ofs, g, dp);
    ofs.close();
}

void write_Results(char* real_name, int n, float threshold_sub, int k_max_naive, int k_max, std::vector<std::vector<float>> stat) {
    char name_out[1000];
    strcpy(name_out, "graph_avg_perf_nb_naive/");
    strcat(name_out, real_name);
    strcat(name_out, ".txt");
    std::ofstream ofs;
    ofs.open(name_out);
    ofs << n << std::endl;
    ofs << threshold_sub << std::endl;
    ofs << k_max_naive << std::endl;
    ofs << k_max << std::endl;
    for (int i = 0; i <= k_max_naive; i++) {
        ofs << stat[0][i] << " " << stat[1][i] << " " << stat[2][i] << " " << stat[3][i] << " " << stat[4][i] << " " << stat[5][i] << " " << stat[6][i] << " " << std::endl;
    }
    ofs.close();
}

void write_result(char* real_name, int k, std::vector<float> stat, int n, int m, bool normal) {
    char name_out[1000];
    strcpy(name_out, "graph_stat_perf_subcolo/");
    strcat(name_out, real_name);
    strcat(name_out, ".txt");
    std::ofstream ofs;
    ofs.open(name_out, std::ofstream::out | std::ofstream::app);
    ofs << stat[0] << " " << stat[1] << " " << stat[2] << " " << stat[3] << " " << stat[4] << " " << stat[5] << " " << k << " " << n << " " << m << " " << normal << std::endl;
    ofs.close();
}

void displayGraph(Graph g, bool visible_edge) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);

    std::cout << "vertices, degree and color of g = ";
    std::pair< vertex_iter, vertex_iter > vp;
    long long unsigned int max_degree = 0;
    float avg_degree = 0.;

    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        std::cout << "(" << get(vertex_id, *vp.first) << ", " << degree(*vp.first, g) << ", " << get(color, *vp.first) << ") ";
        if (degree(*vp.first, g) > max_degree) {
            max_degree = degree(*vp.first, g);
        }
        avg_degree += degree(*vp.first, g);
    }
    std::cout << std::endl;
    std::cout << std::endl;

    avg_degree /= num_vertices(g);

    std::cout << "G(" << num_vertices(g) << ", " << num_edges(g) << "). " << std::endl << "Le degre max du graphe est : " << max_degree << std::endl << "Le degre moyen du graphe est : " << avg_degree << std::endl;


    if (visible_edge) {
        std::cout << "edges(g) = ";
        graph_traits< Graph >::edge_iterator ei, ei_end;
        for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
            std::cout << "(" << get(vertex_id, source(*ei, g)) << ", "
                      << get(vertex_id, target(*ei, g)) << ") ";
        std::cout << std::endl;
    }
}

void blankgraph(Graph& g) {
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);

    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        if (get(color, *vp.first) != 0) {
            put(vertex_color, g, *vp.first, 0);
            // std::cout << "Mauvaise initialisation auto de la couleur... ?" << std::endl;
        }
    }
}

void statColor(Graph g, int k) {
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);

    std::pair< vertex_iter, vertex_iter > vp;    
    std::vector<int> save_color(k+1, 0);
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        if (degree(*vp.first, g) != 0) {
            save_color[get(color, *vp.first)]++;
        } else {
            save_color[0]++;            
        }
    }
    if (save_color[0] != 0) {
        std::cout << "Number of vertex of color undetermined yet (but well colored with score 1 by low_degree_vertex) : " << save_color[0] << std::endl;
    }
    for (int i = 1; i <= k; i++) {
        std::cout << "Number of vertex of color " << i << " : " << save_color[i] << std::endl;
    }
}

int computeLibertyDegree(Graph g, int k, std::vector<std::vector<bool>>& dlColor) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);
    int dl = 0;

    std::pair< vertex_iter, vertex_iter > vp;    
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        for (int i = 1; i <= k; i++) {
            if (get(color, *vp.first) != i) {
                std::pair<adj_iter, adj_iter> neighbor = adjacent_vertices(*vp.first, g);
                bool possibleColor = true;
                bool foundCandidatClique = false;
                while (!foundCandidatClique and neighbor.first != neighbor.second) {
                    if (get(color, *neighbor.first) == i) {
                        foundCandidatClique = true;
                        std::pair<adj_iter, adj_iter> temp = neighbor;
                        temp.first++;
                        while (possibleColor and temp.first != temp.second) {
                            if (get(color, *temp.first) == i) {
                                std::pair<edge_desc, bool> e = edge(*neighbor.first, *temp.first, g);
                                if (!e.second) {
                                    possibleColor = false;
                                }
                            }
                            temp.first++;
                        }
                        std::pair<adj_iter, adj_iter> neighborCand = adjacent_vertices(*neighbor.first, g);
                        while (possibleColor and neighborCand.first != neighborCand.second) {
                            if (*neighborCand.first != *vp.first and get(color, *neighborCand.first) == i) {
                                std::pair<edge_desc, bool> e = edge(*vp.first, *neighborCand.first, g);
                                if (!e.second) {
                                    possibleColor = false;
                                }                            
                            }
                            neighborCand.first++;
                        }
                    }
                    neighbor.first++;
                }
                if (possibleColor) {
                    dl++;
                    dlColor[get(vertex_id, *vp.first)][i] = true;
                    // std::cout << "Le sommet " << get(vertex_id, *vp.first) << " peut être de la couleur " << i << std::endl;
                }
            }
        }
    }
    return dl;
}

void insertion_tri(std::vector<float>& list, float elem) {
    if (list.size() == 0 or elem <= list.back()) {
        list.push_back(elem);
    } else {
        std::vector<float>::iterator temp = list.begin();
        while (temp != list.end() and elem <= *temp) {
            temp++;
        }
        list.insert(temp, elem);
    }
}


