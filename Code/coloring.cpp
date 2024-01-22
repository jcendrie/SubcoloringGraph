#include "coloring.hpp"

using namespace boost;

void max_degree(const Graph& g, int& max_degree) {
    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        if(int(degree(*vp.first, g)) > max_degree) {
            max_degree = degree(*vp.first, g);
        }
    }
}

void sub_graph(Graph& neighbor_graph, vertex_desc center, Graph g) { // create the subgraph of g with neighbor of center uncolorfull
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);

    std::pair<adj_iter, adj_iter> neighbor = adjacent_vertices(center, g);
    for (std::pair<adj_iter, adj_iter> temp = neighbor; temp.first != temp.second; ++temp.first) {
        if (get(color, *temp.first) == 0) {
            std::pair<adj_iter, adj_iter> neighbor_temp = adjacent_vertices(*temp.first, g);
            for (std::pair<adj_iter, adj_iter> temp2 = neighbor_temp; temp2.first != temp2.second; ++temp2.first) {
                if (get(color, *temp2.first) == 0) {
                    std::pair<edge_desc, bool> e = edge(center, *temp2.first, g);
                    if (e.second and get(vertex_id, *temp.first) < get(vertex_id, *temp2.first)) {
                        graph_traits< Graph >::edge_descriptor e;
                        bool inserted;
                        tie(e, inserted) = add_edge(*temp.first, *temp2.first, neighbor_graph);
                    }
                }
            }
        graph_traits< Graph >::edge_descriptor e;
        bool inserted;
        tie(e, inserted) = add_edge(center, *temp.first, neighbor_graph);
        }
    }
}

void color_and_erase(Graph& g, Clique cl, int c, int& uncolored_node) {
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);

    for (Clique::iterator it = cl.begin(); it != cl.end(); ++it) {
        std::pair<adj_iter, adj_iter> neighbor = adjacent_vertices(*it, g);
        if (get(color, *it) == 0) {
            uncolored_node--;
        }
        put(vertex_color, g, *it, c);
        for (std::pair<adj_iter, adj_iter> temp = neighbor; temp.first != temp.second; ++temp.first) {
            if (get(color, *temp.first) == 0) {
                put(vertex_color, g, *temp.first, -1);
                uncolored_node--;
            }
        }
    }
}

void clean(Graph& g, int& uncolored_node) {
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);
    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        if (get(color, *vp.first) == -1) {
            put(vertex_color, g, *vp.first, 0);
            uncolored_node++;
        }
    }        
}

void betweenness(Graph g, std::vector<float>& cb, vertex_desc& max_vertex, float& max_betweenness) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);

    int n = num_vertices(g);

    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        if (get(color, *vp.first) == 0) {
            int vs = get(vertex_id, *vp.first);
            std::stack<vertex_desc> s; // returns vertices in order of non-increasing distance from s
            std::vector<std::vector<int>> p(n);
            std::vector<float> sigma(n, 0.);
            sigma[vs] = 1;
            std::vector<int> d(n, -1);
            d[vs] = 0;
            std::queue<vertex_desc> q;
            q.push(*vp.first);
            while(!q.empty()) {
                vertex_desc vvertex = q.front();
                q.pop();
                int v = get(vertex_id, vvertex);
                s.push(vvertex);

                std::pair<adj_iter, adj_iter> neighbor = adjacent_vertices(vvertex, g);
                for (std::pair<adj_iter, adj_iter> temp = neighbor; temp.first != temp.second; ++temp.first) {
                    if (get(color, *temp.first) == 0) {
                        int w = get(vertex_id, *temp.first);
                        if (d[w] < 0) {
                            q.push(*temp.first);
                            d[w] = d[v] + 1;
                        }
                        if (d[w] == d[v] + 1) {
                            sigma[w] = sigma[w] + sigma[v];
                            p[w].push_back(v);
                        }
                    
                    }
                }
            }

            std::vector<float> delta(n, 0.);
            while (!s.empty()) {
                vertex_desc wvertex = s.top();
                int w = get(vertex_id, wvertex);
                s.pop();
                for (std::vector<int>::size_type i = 0; i < p[w].size(); i++) {
                    int v = p[w][i];
                    delta[v] = delta[v] + (sigma[v]/sigma[w])*(1+delta[w]);
                }
                if (w != vs) {
                    cb[w] = cb[w] + delta[w];
                    if (cb[w] > max_betweenness) {
                        max_betweenness = cb[w];
                        max_vertex = wvertex;
                    }
                }
            }
        }
    }
    if (max_betweenness == 0) {
        vp = vertices(g);
        bool recherche = true;
        while (recherche and vp.first != vp.second) {
            if (get(color, *vp.first) == 0) {
                max_vertex = *vp.first;
                recherche = false;
            } else {
                vp.first++;
            }
        }
    }
}

void insertion_tri(const Graph& g, vertex_desc vertex, std::vector<float> cb, std::vector<vertex_desc>& ltri) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    if (ltri.size() == 0 or cb[get(vertex_id, vertex)] <= cb[get(vertex_id, ltri.back())]) {
        ltri.push_back(vertex);
    } else {
        std::vector<vertex_desc>::iterator temp = ltri.begin();
        while (temp != ltri.end() and cb[get(vertex_id, vertex)] <= cb[get(vertex_id, *temp)]) {
            temp++;
        }
        ltri.insert(temp, vertex);
    }
}

void insertion_tri(const Graph& g, vertex_desc vertex, std::vector<int> aud, std::vector<vertex_desc>& ltri) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    if (ltri.size() == 0 or aud[get(vertex_id, vertex)] <= aud[get(vertex_id, ltri.back())]) {
        ltri.push_back(vertex);
    } else {
        std::vector<vertex_desc>::iterator temp = ltri.begin();
        while (temp != ltri.end() and aud[get(vertex_id, vertex)] <= aud[get(vertex_id, *temp)]) {
            temp++;
        }
        ltri.insert(temp, vertex);
    }
}

void find_clique_max_betweenness(const Graph& g, std::vector<float> cb, vertex_desc center, Clique& cl) {
    cl.push_back(center);
    std::vector<vertex_desc> ltri;

    for (std::pair<adj_iter, adj_iter> temp = adjacent_vertices(center, g); temp.first != temp.second; ++temp.first) {
        insertion_tri(g, *temp.first, cb, ltri);
    }

    for (std::vector<vertex_desc>::iterator temp = ltri.begin(); temp != ltri.end(); ++temp) {
        vertex_desc candidat = *temp;
        bool candidat_valide = true;
        std::deque<vertex_desc>::iterator temp2 = cl.begin();
        while (candidat_valide and temp2 != cl.end()) {
            std::pair<edge_desc, bool> e = edge(candidat, *temp2, g);
            if (!e.second) {
                candidat_valide = false;
            }
            temp2++;
        }
        if (candidat_valide) {
            cl.push_back(candidat);
        }
    }
}

void max_uncolor_degree(Graph g, std::vector<int>& actual_uncolor_degree, vertex_desc& max_vertex, int& max_degree) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);
    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        if (get(color, *vp.first) == 0) {
            if(actual_uncolor_degree[get(vertex_id, *vp.first)] > max_degree) {
                max_degree = out_degree(*vp.first, g);
                max_vertex = *vp.first;
            }            
        }
    }
}

void find_clique_max_degree(const Graph& g, std::vector<int> actual_uncolor_degree, vertex_desc center, Clique& cl) {
    cl.push_back(center);
    std::vector<vertex_desc> ltri;

    for (std::pair<adj_iter, adj_iter> temp = adjacent_vertices(center, g); temp.first != temp.second; ++temp.first) {
        insertion_tri(g, *temp.first, actual_uncolor_degree, ltri);
    }

    for (std::vector<vertex_desc>::iterator temp = ltri.begin(); temp != ltri.end(); ++temp) {
        vertex_desc candidat = *temp;
        bool candidat_valide = true;
        std::deque<vertex_desc>::iterator temp2 = cl.begin();
        while (candidat_valide and temp2 != cl.end()) {
            std::pair<edge_desc, bool> e = edge(candidat, *temp2, g);
            if (!e.second) {
                candidat_valide = false;
            }
            temp2++;
        }
        if (candidat_valide) {
            cl.push_back(candidat);
        }
    }
}

void maj_uncolor_degree(const Graph& g, std::vector<int>& actual_uncolor_degree, Clique& cl) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    for (std::deque<vertex_desc>::iterator elem = cl.begin(); elem != cl.end(); elem++) {
        for (std::pair<adj_iter, adj_iter> neighbor = adjacent_vertices(*elem, g); neighbor.first != neighbor.second; neighbor.first++) {
            actual_uncolor_degree[get(vertex_id, *neighbor.first)]--;
        }
    }
}

int subcoloring(Graph& g, std::string mode, bool normal) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);

    blankgraph(g);
    int k = 0;
    int n = num_vertices(g);
    int uncolored_node = n;
    int step = 10*n;

    std::vector<int> actual_uncolor_degree(n, 0);
    
    if (mode == "degree") {
        std::pair< vertex_iter, vertex_iter > vp;
        for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
            actual_uncolor_degree[get(vertex_id, *vp.first)] = degree(*vp.first, g);
        }
    }

    int first_max_degree = 0;
    max_degree(g, first_max_degree);

    
    while (uncolored_node != 0 && step > 0) {
        k++;
        while (uncolored_node != 0 && step > 0) {
            vertex_desc max_vertex;
            std::vector<float> cb(n, 0.); // betweenness

            if (mode == "degree") {
                int max_degree = -1;
                max_uncolor_degree(g, actual_uncolor_degree, max_vertex, max_degree);
            } else if (mode == "between") { // max degree
                float max_betweenness = 0;
                betweenness(g, cb, max_vertex, max_betweenness);
            }

            if (normal) {
                Clique cl;
                cl.push_back(max_vertex);
                color_and_erase(g, cl, k, uncolored_node);
            } else {
                Graph neighbor_graph;
                sub_graph(neighbor_graph, max_vertex, g);

                if (num_edges(neighbor_graph) == 0) {
                    put(vertex_color, g, max_vertex, k);
                    uncolored_node--;
                } else {
                    Clique cl;
                    if (mode == "degree") {
                        find_clique_max_degree(neighbor_graph, actual_uncolor_degree, max_vertex, cl);
                        maj_uncolor_degree(g, actual_uncolor_degree, cl);
                    } else if (mode == "between") {
                        find_clique_max_betweenness(neighbor_graph, cb, max_vertex, cl);
                    }
                    color_and_erase(g, cl, k, uncolored_node);
                }                
            }                
            // displayGraph(g);
            // std::cout << uncolored_node << std::endl;   s
            step--;
        }
        clean(g, uncolored_node);
    }
    if (uncolored_node > 0) {
        std::cout << "Warning, uncomplete coloration, increase the number of step !" << std::endl;
    } else {
        std::cout << "Complete coloration of the graph (" << num_vertices(g) << ", " << num_edges(g) << ") of max degree " << first_max_degree << " with " << k << " color(s) !" << std::endl;
    }
    return k;
}

int channel_by_arrival(Graph& g, int k_max) {
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);

    blankgraph(g);

    int k = 0;
    std::vector<vertex_desc> list_vertex;

    for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; vp.first++) {
        list_vertex.push_back(*vp.first);
    }

    std::random_shuffle(list_vertex.begin(), list_vertex.end());

    for (long long unsigned int i = 0; i < list_vertex.size(); i++) {
        std::vector<int> color_neighbor(k_max+1, 0);
        std::pair<adj_iter, adj_iter> neighbor;
        for (neighbor = adjacent_vertices(list_vertex[i], g); neighbor.first != neighbor.second; neighbor.first++) {
            color_neighbor[get(color, *neighbor.first)]++;
        }
        int min_color = 1;
        int min_nb_color = color_neighbor[1];
        for (int j = 2; j < k_max+1; j++) {
            if (color_neighbor[j] < min_nb_color) {
                min_nb_color = color_neighbor[j];
                min_color = j;
            }
        }
        put(vertex_color, g, list_vertex[i], min_color);
        if (min_color > k) {
            k = min_color;
        }
    }
    return k;
}

int coloring_by_Markov(Graph& g, float threshold) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);

    blankgraph(g);

    int k = 0;
    int n = num_vertices(g);
    int uncolored = n;
    int previous_uncolored = 0;
    int step = 10000;
    std::vector<float> score(n, 0);
    while (uncolored != 0 && step > 0) {
        k++;
        previous_uncolored = uncolored;
        for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; vp.first++) {
            if (get(color, *vp.first) == 0) {
                put(vertex_color, g, *vp.first, k);
                uncolored--;
            }
        }

        score_Markov(g, score, k, "");

        for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; vp.first++) {
            if (score[get(vertex_id, *vp.first)] < threshold) {
                put(vertex_color, g, *vp.first, 0);
                uncolored++;
                if (uncolored == previous_uncolored) {
                    put(vertex_color, g, *vp.first, k);
                    uncolored--;
                }
            }
        }
    }
    if (uncolored > 0) {
        std::cout << "Warning, uncomplete coloration, increase the number of step !" << std::endl;
    } else {
        std::cout << "Complete coloration of the graph (" << num_vertices(g) << ", " << num_edges(g) << ") with " << k << " color(s) !" << std::endl;
    }
    return k;
}

int naive_then_subcoloring(Graph& g, float threshold, int k_naive, int k_max) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);

    blankgraph(g);

    int n = num_vertices(g);
    int k = 0;
    int uncolored = n;
    int previous_uncolored = 0;
    std::vector<float> score(n, 0);


    if (k_naive != 0) {
        k = channel_by_arrival(g, k_naive);

        while (uncolored != previous_uncolored) {
            previous_uncolored = uncolored;
            uncolored = 0;
            score_Markov(g, score, k, "");

            for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; vp.first++) {
                if (score[get(vertex_id, *vp.first)] <= threshold) {
                    put(vertex_color, g, *vp.first, 0);
                    uncolored++;
                }
            }
        }
        // std::cout << uncolored << std::endl;
    }
    if (uncolored != 0) {
        Graph subcoloring_graph(n);

        for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; vp.first++) {
            if (get(color, *vp.first) == 0) {
                std::pair<vertex_iter, vertex_iter> vp2 = vp;
                for (vp2.first++; vp2.first != vp2.second; vp2.first++) {
                    if (get(color, *vp2.first) == 0) {
                        std::pair<edge_desc, bool> e = edge(*vp.first, *vp2.first, g);
                        if (e.second) {
                            graph_traits< Graph >::edge_descriptor e;
                            bool inserted;
                            tie(e, inserted) = add_edge(*vp.first, *vp2.first, subcoloring_graph);
                        }
                    }
                }
            }
        }

        int k_sub = subcoloring(subcoloring_graph, "degree");

        score_sub_colo(subcoloring_graph, score, k + k_sub, "");

        while (k_sub < k_max-k) {
            add_color_divide_biggest_clique(subcoloring_graph, score, k_sub);
        }
        upgrade_with_dl(subcoloring_graph, k_sub);

        std::vector<int> new_color(n, 0);
        for (std::pair<vertex_iter, vertex_iter> vp = vertices(subcoloring_graph); vp.first != vp.second; vp.first++) {
            new_color[get(vertex_index, subcoloring_graph, *vp.first)] = get(vertex_color, subcoloring_graph, *vp.first);
        }

        // score_sub_colo(subcoloring_graph, score, k, "Test 3");


        for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; vp.first++) {
            if (get(color, *vp.first) == 0) {
                put(vertex_color, g, *vp.first, new_color[get(vertex_id, *vp.first)] + k);
            }
        }
        k += k_sub;
    }
    return k;
}

int extract_low_degree_vertex(Graph& g, Graph& g_ori, std::stack<vertex_desc>& stack_low_degree, int k_max) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);

    int n = num_vertices(g);
    Graph gprime(n);
    int count = 0;

    std::vector<bool> low_degree_vertex(n, false);
    std::vector<int> actual_degree(n, 0);

    for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; ++vp.first) {
        actual_degree[get(vertex_id, *vp.first)] = degree(*vp.first, g);
    }
    
    bool extract = true;
    while (extract) {
        extract = false;
        for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; vp.first++) {
            int ind = get(vertex_id, *vp.first);
            if (!low_degree_vertex[ind] and actual_degree[ind] < k_max) {
                low_degree_vertex[ind] = true;
                stack_low_degree.push(*vp.first);
                count++;
                extract = true;
                for (std::pair<adj_iter, adj_iter> neighbor = adjacent_vertices(*vp.first, g); neighbor.first != neighbor.second; neighbor.first++) {
                    if (!low_degree_vertex[get(vertex_id, *neighbor.first)]) {
                        actual_degree[get(vertex_id, *neighbor.first)]--;
                    }
                }
            }
        }
    }
    for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; vp.first++) {
        if (!low_degree_vertex[get(vertex_id, *vp.first)]) {
            std::pair<vertex_iter, vertex_iter> vp2 = vp;
            for (vp2.first++; vp2.first != vp2.second; vp2.first++) {
                if (!low_degree_vertex[get(vertex_id, *vp2.first)]) {
                    std::pair<edge_desc, bool> e = edge(*vp.first, *vp2.first, g);
                    if (e.second) {
                        graph_traits< Graph >::edge_descriptor e;
                        bool inserted;
                        tie(e, inserted) = add_edge(*vp.first, *vp2.first, gprime);
                    }
                }
            }
        }
    }

    g_ori = g;
    g = gprime;
    std::cout << "Nombre de sommets retires par analyse du degre des sommets : " << count << std::endl;
    return count;
}

int back_low_degree_vertex(Graph& g, Graph& g_ori, std::stack<vertex_desc>& stack_low_degree, int k_max) {
    int n = num_vertices(g_ori);
    int i_max = 0;
    std::vector<bool> in_g_ori(n, false);
    for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; vp.first++) {
        if (degree(*vp.first, g) != 0) {
            put(vertex_color, g_ori, *vp.first, get(vertex_color, g, *vp.first));
        }
        if (degree(*vp.first, g) != 0) {
            in_g_ori[get(vertex_index, g, *vp.first)] = true;
        }
    }

    while (!stack_low_degree.empty()) {
        vertex_desc cand = stack_low_degree.top();
        std::vector<bool> color_see(k_max+1, false);
        for (std::pair<adj_iter, adj_iter> neighbor = adjacent_vertices(cand, g_ori); neighbor.first != neighbor.second; neighbor.first++) {
            if (in_g_ori[get(vertex_index, g, *neighbor.first)]) {
                color_see[get(vertex_color, g, *neighbor.first)] = true;
            }
        }
        int i = 1;
        while (color_see[i]) {
            i++;
        }
        if (i > i_max) {
            i_max = i;
        }
        put(vertex_color, g_ori, cand, i);
        in_g_ori[get(vertex_index, g, cand)] = true;
        stack_low_degree.pop();
        g = g_ori;
    }
    return i_max;
}


