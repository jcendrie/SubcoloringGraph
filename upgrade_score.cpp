#include "upgrade_score.hpp"




float end_score(Graph g, float score_min, std::vector<float>& score, int k, const char* name_out) {
    std::cout << "Min_score : " << score_min << std::endl;
    int n = num_vertices(g);
    float tot = 0;
    char name[1000];
    strcpy(name, "histo/");
    strcat(name, name_out);
    strcat(name, ".txt");
    std::ofstream histo(name);
    if (histo) {
        histo << k << std::endl;
        for (int i = 0; i < n; i++) {
            // histo << score[i] << " " << get(vertex_color, g, i) << std::endl;
            histo << score[i] << std::endl;
            tot += score[i];
        }
    } else {
        std::cout << "Error with the file" << std::endl;
    }
    histo.close();
    std::cout << "Avg_score : " << tot/n << std::endl;
    return tot/n;
}


std::pair<float, float> score_sub_colo(Graph g, std::vector<float>& score, int k, const char* name_out) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);

    std::pair< vertex_iter, vertex_iter > vp;
    float score_min = 1;
    float avg_score = 0;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        float nb_neighbor_same_color = 0.;
        std::pair<adj_iter, adj_iter> neighbor;
        for (neighbor = adjacent_vertices(*vp.first, g); neighbor.first != neighbor.second; neighbor.first++) {
            if (get(color, *vp.first) == get(color, *neighbor.first)) {
                nb_neighbor_same_color += 1.;
            }
        }
        score[get(vertex_id, *vp.first)] = 1/(nb_neighbor_same_color + 1);
        if (score[get(vertex_id, *vp.first)] < score_min) {
            score_min = score[get(vertex_id, *vp.first)];
        }
    }
    if (strcmp(name_out, "")) {
        avg_score = end_score(g, score_min, score, k, name_out);
    }
    return std::make_pair(score_min, avg_score);
}

struct clique_visitor_score {
    clique_visitor_score(long long unsigned int& max, std::vector<Clique>& mis) : k_max(max), k_mis(mis) {}
    long long unsigned int& k_max;
    std::vector<Clique>& k_mis;
    
    template <typename Clique, typename Graph>
    inline void clique(const Clique& p, const Graph& g) {
        if (p.size() > k_max) {
            k_mis.clear();
            k_max = p.size();
            k_mis.push_back(p);
        } else if (p.size() == k_max) {
            k_mis.push_back(p);
        }
    }
};

inline clique_visitor_score find_clique(long long unsigned int& max, std::vector<Clique>& mis) {
    return clique_visitor_score(max, mis);
}

std::pair<float, float> score_MIS(Graph g, std::vector<float>& score, int k, const char* name_out) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);
    float score_min = 1;
    float avg_score = 0;
    int n = num_vertices(g);
    

    for (int i = 1; i <= k; i++) {
        Graph gprime;

        std::vector<int> good_color;
        std::pair< vertex_iter, vertex_iter > vp;
        for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
            if (get(color, *vp.first) == i) {
                good_color.push_back(get(vertex_id, *vp.first));
                std::pair< vertex_iter, vertex_iter > u = vp;    
                for (u.first++; u.first != u.second; ++u.first) {
                    if (get(color, *u.first) == i) {
                        std::pair<edge_desc, bool> e = edge(*vp.first, *u.first, g);
                        if (!e.second) {
                            add_edge(*vp.first, *u.first, gprime);
                        }
                    }
                }
            }
        }

        std::vector<Clique> mis;
        long long unsigned int max = 0;

        bron_kerbosch_all_cliques(gprime, find_clique(max, mis));

        std::vector<int> in_mis(n, 0);
        for (long long unsigned int j = 0; j < mis.size(); j++) {
            for (long long unsigned int l = 0; l < mis[j].size(); l++) {
                in_mis[mis[j][l]]++;
            }
        }

        for (long long unsigned int j = 0; j < good_color.size(); j++) {
            score[good_color[j]] = float(in_mis[good_color[j]]) / float(mis.size());
            if (score[good_color[j]] < score_min) {
                score_min = score[good_color[j]];
            }
        }
    }
    if (strcmp(name_out, "")) {
        avg_score = end_score(g, score_min, score, k, name_out);
    }
    return std::make_pair(score_min, avg_score);
}

std::pair<float, float> score_Markov(Graph g, std::vector<float>& score, int& k, const char* name_out, bool max_MIS) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);
    float score_min = 1;
    float avg_score = 0;
    bool new_seed = true;
    int n = num_vertices(g);
    int taille_MIS = 0;
    
    for (int i = 0; i < n; i++) {
        score[i] = 0.;
    }

    for (int i = 1; i <= k; i++) {
        std::vector<vertex_desc> tab_sommets;
        std::pair< vertex_iter, vertex_iter > vp;
        int nprime = 0;
        for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
            if (get(color, *vp.first) == i) {
                tab_sommets.push_back(*vp.first);
                nprime++;
            }
        }

        std::vector<int> dis_time(n, 0);
        std::vector<int> fin_time(n, 0);
        // discover_time -> permet de vérifier depuis quelle taille de MIS le sommet est dans le MIS 
        // finish_time -> permet de savoir à quel moment il est entrée dans le MIS
        // donc si, au moment où on retire le sommet, le discover_time = len_max_MIS, on prend en compte nb_MIS - finish_time pour le score, sinon on prend seulement nb_MIS.

        float lambda = nprime; // A observer, implique une forte dépendance dans l'obtention des MIS de grandes tailles mais aussi le changement de MIS...
        if (max_MIS) {
            lambda = nprime/10;
        }
        int len_max_MIS = 0;
        int len_actual_MIS = 0;
        int nb_MIS = 0;
        std::vector<int> in_MIS(n, 0);

        float conv = nprime;
        int step = 10000000;
        int intervalle = nprime*1000;
        int next_stop = step - intervalle;
        while (conv > 0.001*nprime and step > 0) { 
            step--;
            vertex_desc candidat = tab_sommets[rand()%nprime];
            if (dis_time[get(vertex_id, candidat)] == 0) {
                bool possibleCand = true;
                std::pair<adj_iter, adj_iter> neighbor = adjacent_vertices(candidat, g);
                while (possibleCand and neighbor.first != neighbor.second) {
                    if (get(color, *neighbor.first) == i) {
                        if (dis_time[get(vertex_id, *neighbor.first)] != 0) {
                            possibleCand = false;
                        }                        
                    }
                    neighbor.first++;
                }
                float test = ((double) rand() / (RAND_MAX));
                if (possibleCand and test < lambda/(1+lambda)) {
                    len_actual_MIS++;
                    if (len_actual_MIS > len_max_MIS) {
                        len_max_MIS = len_actual_MIS;
                        next_stop = step - intervalle;
                        nb_MIS = 0;
                        for (int alpha = 0; alpha < n; alpha++) {
                            in_MIS[alpha] = 0;
                        }
                    }
                    dis_time[get(vertex_id, candidat)] = len_max_MIS;
                    fin_time[get(vertex_id, candidat)] = nb_MIS;
                }
            } else {
                float test = ((double) rand() / (RAND_MAX));
                if (test < 1/(1+lambda)) {
                    if (dis_time[get(vertex_id, candidat)] == len_max_MIS) {
                        in_MIS[get(vertex_id, candidat)] += nb_MIS - fin_time[get(vertex_id, candidat)];
                    } else {
                        in_MIS[get(vertex_id, candidat)] += nb_MIS;
                    }
                    dis_time[get(vertex_id, candidat)] = 0;
                    fin_time[get(vertex_id, candidat)] = 0;
                    len_actual_MIS--;
                }
            }
            if (len_actual_MIS == len_max_MIS) {
                nb_MIS++;
            }
            if (step == next_stop) {
                // std::cout << i << " " << len_max_MIS << " " << len_actual_MIS << " " << step << " " << nb_MIS << std::endl;
                next_stop = step - intervalle;
                conv = 0;
                for (int j = 0; j < nprime; j++) {
                    if (dis_time[get(vertex_id, tab_sommets[j])] != 0) {
                        if (dis_time[get(vertex_id, tab_sommets[j])] == len_max_MIS) {
                            in_MIS[get(vertex_id, tab_sommets[j])] += nb_MIS - fin_time[get(vertex_id, tab_sommets[j])];
                        } else {
                            in_MIS[get(vertex_id, tab_sommets[j])] += nb_MIS;
                        }
                        if (new_seed) {
                            dis_time[get(vertex_id, tab_sommets[j])] = 0;
                            fin_time[get(vertex_id, tab_sommets[j])] = 0;
                            len_actual_MIS--;
                        } else {
                            dis_time[get(vertex_id, tab_sommets[j])] = len_max_MIS;
                            fin_time[get(vertex_id, tab_sommets[j])] = nb_MIS;                            
                        }
                    }
                    float temp = float(in_MIS[get(vertex_id, tab_sommets[j])]) / float(nb_MIS);
                    conv += (score[tab_sommets[j]]-temp)*(score[tab_sommets[j]]-temp);
                    score[tab_sommets[j]] = temp;
                    // std::cout << temp << std::endl;
                }
                conv = sqrt(conv);
                // std::cout << "conv :" << conv << std::endl;
            }
        }
        // std::cout << "couleurs " << i << " taille max MIS " << len_max_MIS << std::endl;
        if (step == 0) {
            std::cout << "Attention, arrêt par manque d'étapes et non par convergence..." << std::endl;
        }

        for (int j = 0; j < nprime; j++) {
            if (dis_time[get(vertex_id, tab_sommets[j])] != 0) {
                if (dis_time[get(vertex_id, tab_sommets[j])] == len_max_MIS) {
                    in_MIS[get(vertex_id, tab_sommets[j])] += nb_MIS - fin_time[get(vertex_id, tab_sommets[j])];
                } else {
                    in_MIS[get(vertex_id, tab_sommets[j])] += nb_MIS;
                }
            }
            score[tab_sommets[j]] = float(in_MIS[get(vertex_id, tab_sommets[j])]) / float(nb_MIS);
            if (score[tab_sommets[j]] < score_min) {
                score_min = score[tab_sommets[j]];
            }
        }
        if (max_MIS and len_max_MIS > taille_MIS) {
            taille_MIS = len_max_MIS;
        }
    }
    if (strcmp(name_out, "")) {
        avg_score = end_score(g, score_min, score, k, name_out);
    }
    if (max_MIS) {
        k = taille_MIS;
    }
    return std::make_pair(score_min, avg_score);
}

int taille_clique_max_Bron(Graph g) {
    return (int) bron_kerbosch_clique_number(g);
}

int taille_clique_max_Markov(Graph g) {
    int n = num_vertices(g);
    int taille_clique_max = 1;
    int replication = 5;
    
    Graph gprime(n);

    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        std::pair< vertex_iter, vertex_iter > u = vp;    
        put(vertex_color, gprime, *vp.first, 1);
        for (u.first++; u.first != u.second; ++u.first) {
            std::pair<edge_desc, bool> e = edge(*vp.first, *u.first, g);
            if (!e.second) {
                add_edge(*vp.first, *u.first, gprime);
            }
        }
    }


    for (int i = 0; i < replication; i++) {
        int taille_clique = 1;
        std::vector<float> score(n, 0);

        score_Markov(gprime, score, taille_clique, "", true);
        if (taille_clique > taille_clique_max) {
            taille_clique_max = taille_clique;
        }
    }
    // for (int i = 0; i < n; i++) {
    //     std::cout << score[i] << std::endl;
    // }
    return taille_clique_max;
}

int compute_dl_color_inv_score(Graph g, int k, std::vector<std::vector<int>>& dl_color_inv_score) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);
    int max_score = 0;

    for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; vp.first++) {
        for (int i = 1; i <= k; i++) {
            int nb_neighbor_same_color = 0;
            bool possibleColor = true;
            bool foundCandidatClique = false;
            std::pair<adj_iter, adj_iter> neighbor = adjacent_vertices(*vp.first, g);
            while (!foundCandidatClique and neighbor.first != neighbor.second) {
                if (get(color, *neighbor.first) == i) {
                    foundCandidatClique = true;
                    std::pair<adj_iter, adj_iter> temp = neighbor;
                    temp.first++;
                    nb_neighbor_same_color++;
                    while (possibleColor and temp.first != temp.second) {
                        if (get(color, *temp.first) == i) {
                            std::pair<edge_desc, bool> e = edge(*neighbor.first, *temp.first, g);
                            nb_neighbor_same_color++;
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
                dl_color_inv_score[get(vertex_id, *vp.first)][i] = nb_neighbor_same_color;
                if (i == get(color, *vp.first) and nb_neighbor_same_color > max_score) {
                    max_score = nb_neighbor_same_color;
                }        
            } else {
                dl_color_inv_score[get(vertex_id, *vp.first)][i] = -1;
            }
        }
    }
    return max_score;
}

void upgrade_with_dl(Graph& g, int& k) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);
    int n = num_vertices(g);

    std::vector<std::vector<int>> dl_color_inv_score;
    for (int i = 0; i < n; i++) {
        std::vector<int> temp(k+1, -1);
        dl_color_inv_score.push_back(temp);
    }
    int max_score = compute_dl_color_inv_score(g, k, dl_color_inv_score);

    std::vector<std::queue<vertex_desc>> list_file;
    for (int i = 0; i <= max_score; i++) {
        std::queue<vertex_desc> temp;
        list_file.push_back(temp);
    }

    for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; vp.first++) {
        list_file[dl_color_inv_score[get(vertex_id, *vp.first)][get(color, *vp.first)]].push(*vp.first);
    }

    std::vector<bool> vu(n, false);

    // std::cout << max_score << std::endl;

    int step = n*n;

    while (max_score != 0 and step > 0) {
        vertex_desc cand = list_file[max_score].front();
        vu[get(vertex_id, cand)] = true;
        list_file[max_score].pop();

        // for (int i = 16; i < 17; i++) {
        //     for (int j = 1; j <= k; j++) {
        //         std::cout << i << " " << j << " " << dl_color_inv_score[i][j] << " " << vu[i] << std::endl;
        //     }
        // } 
        
        if (max_score == dl_color_inv_score[get(vertex_id, cand)][get(color, cand)]) {
            int min_score = max_score;
            int min_color = get(color, cand);
            for (int i = 1; i <= k; i++) {
                if (dl_color_inv_score[get(vertex_id, cand)][i] != -1 and dl_color_inv_score[get(vertex_id, cand)][i] < min_score) {
                    min_score = dl_color_inv_score[get(vertex_id, cand)][i];
                    min_color = i;
                }
            }
            if (min_color != get(color, cand)) {
                int old_color = get(color, cand);
                // std::cout << "test sommets " << get(vertex_id, cand) << " old " << old_color << " new " << min_color << std::endl;
                put(vertex_color, g, cand, min_color);
                for (std::pair<adj_iter, adj_iter> neighbor = adjacent_vertices(cand, g); neighbor.first != neighbor.second; neighbor.first++) {
                    bool change = false;
                    if (dl_color_inv_score[get(vertex_id, *neighbor.first)][old_color] != -1) {
                        dl_color_inv_score[get(vertex_id, *neighbor.first)][old_color]--;
                        change = true;
                    } else {
                        int nb_neighbor_same_color = 0;
                        bool possibleColor = true;
                        bool foundCandidatClique = false;
                        std::pair<adj_iter, adj_iter> neighneighbor = adjacent_vertices(*neighbor.first, g);
                        while (!foundCandidatClique and neighneighbor.first != neighneighbor.second) {
                            if (get(color, *neighneighbor.first) == old_color) {
                                foundCandidatClique = true;
                                std::pair<adj_iter, adj_iter> neighneighbor2 = neighneighbor;
                                neighneighbor2.first++;
                                nb_neighbor_same_color++;
                                while (possibleColor and neighneighbor2.first != neighneighbor2.second) {
                                    if (get(color, *neighneighbor2.first) == old_color) {
                                        std::pair<edge_desc, bool> e = edge(*neighneighbor.first, *neighneighbor2.first, g);
                                        nb_neighbor_same_color++;
                                        if (!e.second) {
                                            possibleColor = false;
                                        }
                                    }
                                    neighneighbor2.first++;
                                }
                                std::pair<adj_iter, adj_iter> neighneighneighbor = adjacent_vertices(*neighneighbor.first, g);
                                while (possibleColor and neighneighneighbor.first != neighneighneighbor.second) {
                                    if (*neighneighneighbor.first != *neighbor.first and get(color, *neighneighneighbor.first) == old_color) {
                                        std::pair<edge_desc, bool> e = edge(*neighbor.first, *neighneighneighbor.first, g);
                                        if (!e.second) {
                                            possibleColor = false;
                                        }                            
                                    }
                                    neighneighneighbor.first++;
                                }
                            }
                            neighneighbor.first++;
                        }
                        if (possibleColor) {
                            dl_color_inv_score[get(vertex_id, *neighbor.first)][old_color] = nb_neighbor_same_color;
                            change = true;
                        }
                    }
                    if (dl_color_inv_score[get(vertex_id, *neighbor.first)][min_color] != -1) {
                        bool stayPossible = true;
                        std::pair<adj_iter, adj_iter> neighneighbor = adjacent_vertices(*neighbor.first, g);
                        while (neighneighbor.first != neighneighbor.second) {
                            if (*neighneighbor.first != cand and get(color, *neighneighbor.first) == min_color) {
                                std::pair<edge_desc, bool> e = edge(cand, *neighneighbor.first, g);
                                if (!e.second) {
                                    stayPossible = false;
                                }
                            }
                            neighneighbor.first++;
                        }
                        std::pair<adj_iter, adj_iter> neighbor2 = adjacent_vertices(cand, g);
                        while (stayPossible and neighbor2.first != neighbor2.second) {
                            if (*neighbor2.first != *neighbor.first and get(color, *neighbor2.first) == min_color) {
                                std::pair<edge_desc, bool> e = edge(*neighbor.first, *neighbor2.first, g);
                                if (!e.second) {
                                    stayPossible = false;
                                }
                            }
                            neighbor2.first++;
                        }
                        if (!stayPossible) {
                            // if (get(color, *neighbor.first) == min_color) {
                            //     std::cout << "FIRST WARNING" << std::endl;
                            //     exit(-1);
                            // }
                            dl_color_inv_score[get(vertex_id, *neighbor.first)][min_color] = -1;
                        } else {
                            dl_color_inv_score[get(vertex_id, *neighbor.first)][min_color]++;
                            change = true;
                        }
                    }
                    if (get(color, *neighbor.first) == old_color or get(color, *neighbor.first) == min_color or 
                        (change and vu[get(vertex_id, *neighbor.first)] and dl_color_inv_score[get(vertex_id, *neighbor.first)][old_color] != -1 and dl_color_inv_score[get(vertex_id, *neighbor.first)][old_color] < dl_color_inv_score[get(vertex_id, *neighbor.first)][get(color, *neighbor.first)])) {
                        vu[get(vertex_id, *neighbor.first)] = false;
                        // std::cout << dl_color_inv_score[get(vertex_id, *neighbor.first)][get(color, *neighbor.first)] << std::endl;
                        list_file[dl_color_inv_score[get(vertex_id, *neighbor.first)][get(color, *neighbor.first)]].push(*neighbor.first);
                        if (dl_color_inv_score[get(vertex_id, *neighbor.first)][get(color, *neighbor.first)] > max_score) {
                            max_score = dl_color_inv_score[get(vertex_id, *neighbor.first)][get(color, *neighbor.first)];
                        }
                    }


                    if (get(color, *neighbor.first) == old_color) {
                        for (std::pair<adj_iter, adj_iter> neighneighbor = adjacent_vertices(*neighbor.first, g); neighneighbor.first != neighneighbor.second; neighneighbor.first++) {
                            if (*neighneighbor.first != cand) {
                                std::pair<edge_desc, bool> e = edge(cand, *neighneighbor.first, g);
                                if (!e.second) {
                                    bool possibleColor = true;
                                    int nb_neighbor_same_color = 1;
                                    std::pair<adj_iter, adj_iter> neighneighneighbor = adjacent_vertices(*neighneighbor.first, g);
                                    while (possibleColor and neighneighneighbor.first != neighneighneighbor.second) {
                                        if (*neighneighneighbor.first != *neighbor.first and get(color, *neighneighneighbor.first) == old_color) {
                                            std::pair<edge_desc, bool> e = edge(*neighbor.first, *neighneighneighbor.first, g);
                                            nb_neighbor_same_color++;
                                            if (!e.second) {
                                                possibleColor = false;
                                            }
                                        }
                                        neighneighneighbor.first++;
                                    }
                                    std::pair<adj_iter, adj_iter> neighneighbor2 = adjacent_vertices(*neighbor.first, g);
                                    while (possibleColor and neighneighbor2.first != neighneighbor2.second) {
                                        if (*neighneighbor2.first != *neighneighbor.first and get(color, *neighneighbor2.first) == old_color) {
                                            std::pair<edge_desc, bool> e = edge(*neighneighbor.first, *neighneighbor2.first, g);
                                            if (!e.second) {
                                                possibleColor = false;
                                            }                            
                                        }
                                        neighneighbor2.first++;
                                    }
                                    if (possibleColor) {
                                        dl_color_inv_score[get(vertex_id, *neighneighbor.first)][old_color] = nb_neighbor_same_color;
                                        if (vu[get(vertex_id, *neighneighbor.first)] and nb_neighbor_same_color < dl_color_inv_score[get(vertex_id, *neighneighbor.first)][get(color, *neighneighbor.first)]) {
                                            vu[get(vertex_id, *neighneighbor.first)] = false;
                                            list_file[dl_color_inv_score[get(vertex_id, *neighneighbor.first)][get(color, *neighneighbor.first)]].push(*neighneighbor.first);
                                            if (dl_color_inv_score[get(vertex_id, *neighneighbor.first)][get(color, *neighneighbor.first)] > max_score) {
                                                max_score = dl_color_inv_score[get(vertex_id, *neighneighbor.first)][get(color, *neighneighbor.first)];
                                            }                                            
                                        }
                                    }
                                }
                            }
                        }
                    } else if (get(color, *neighbor.first) == min_color) {
                        for (std::pair<adj_iter, adj_iter> neighneighbor = adjacent_vertices(*neighbor.first, g); neighneighbor.first != neighneighbor.second; neighneighbor.first++) {
                            if (*neighneighbor.first != cand) {
                                std::pair<edge_desc, bool> e = edge(cand, *neighneighbor.first, g);
                                if (!e.second) {
                                    // if (get(color, *neighneighbor.first) == min_color) {
                                    //     std::cout << "WARNING" << std::endl;
                                    //     exit(-1);
                                    // }
                                    dl_color_inv_score[get(vertex_id, *neighneighbor.first)][min_color] = -1;
                                }                            
                            }
                        }                        
                    }
                }
            }
        }
        while (list_file[max_score].empty() and max_score != 0) {
            max_score--;
        }
        step--;
    }

    if (step == 0) {
        std::cout << "Warning, uncomplete upgrade, please make more step" << std::endl;
    }

    std::vector<int> save_color(k+1, 0);
    for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; ++vp.first) {
        save_color[get(color, *vp.first)]++;
    }
    int sub_k = 0;
    for (int i = 1; i <= k; i++) {
        if (save_color[i] == 0) {
            sub_k++;
        }
    }
    if (sub_k != 0) {
        k -= sub_k;
        std::cout << "On utilise " << sub_k << " couleurs en moins que le maximum autorisé." << std::endl;
    }
}

// void upgrade_with_dl_v0(Graph& g, int& k) {
//     property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
//     property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);

//     std::vector<std::vector<bool>> dlColor;
//     int n = num_vertices(g);
//     for (int i = 0; i < n; i++) {
//         std::vector<bool> temp(k+1, false);
//         dlColor.push_back(temp);
//     }
//     computeLibertyDegree(g, k, dlColor);

//     bool upgraded = true;
//     int step = 1000;

//     while (upgraded and step > 0) {
//         upgraded = false;
//         std::pair< vertex_iter, vertex_iter > vp;    
//         std::cout << "step : " << step << std::endl;
//         for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
//             std::vector<int> neighborColor(k+1, 0);
//             for (int i = 1; i <= k; i++) {
//                 if (dlColor[get(vertex_id, *vp.first)][i]) {
//                     std::pair<adj_iter, adj_iter> neighbor;
//                     for (neighbor = adjacent_vertices(*vp.first, g); neighbor.first != neighbor.second; neighbor.first++) {
//                         neighborColor[get(color, *neighbor.first)]++;
//                     }
//                 }
//             }
//             int min_nb_neighbor_same_color = neighborColor[get(color, *vp.first)];
//             int min_color = get(color, *vp.first);
//             for (int i = 1; i <= k; i++) {
//                 if (dlColor[get(vertex_id, *vp.first)][i] and neighborColor[i] < min_nb_neighbor_same_color) {
//                     min_nb_neighbor_same_color = neighborColor[i];
//                     min_color = i;
//                 }                
//             }
//             if (min_color != get(color, *vp.first)) {
//                 int old_color = get(color, *vp.first);
//                 upgraded = true;
//                 // std::cout << "test sommets " << get(vertex_id, *vp.first) << " old " << old_color << " new " << min_color << std::endl;
//                 put(vertex_color, g, *vp.first, min_color);
//                 dlColor[get(vertex_id, *vp.first)][min_color] = false;
//                 dlColor[get(vertex_id, *vp.first)][old_color] = true;
//                 std::pair<adj_iter, adj_iter> neighbor;
//                 for (neighbor = adjacent_vertices(*vp.first, g); neighbor.first != neighbor.second; neighbor.first++) {
//                     if (dlColor[*neighbor.first][min_color]) {
//                         bool stayPossible = true;
//                         std::pair<adj_iter, adj_iter> neighneighbor = adjacent_vertices(*neighbor.first, g);
//                         while (neighneighbor.first != neighneighbor.second) {
//                             if (*neighneighbor.first != *vp.first and get(color, *neighneighbor.first) == min_color) {
//                                 std::pair<edge_desc, bool> e = edge(*vp.first, *neighneighbor.first, g);
//                                 if (!e.second) {
//                                     stayPossible = false;
//                                 }
//                             }
//                             neighneighbor.first++;
//                         }
//                         neighneighbor = adjacent_vertices(*vp.first, g);
//                         while (stayPossible and neighneighbor.first != neighneighbor.second) {
//                             if (*neighneighbor.first != *neighbor.first and get(color, *neighneighbor.first) == min_color) {
//                                 std::pair<edge_desc, bool> e = edge(*neighbor.first, *neighneighbor.first, g);
//                                 if (!e.second) {
//                                     stayPossible = false;
//                                 }
//                             }
//                             neighneighbor.first++;
//                         }
//                         if (!stayPossible) {
//                             dlColor[*neighbor.first][min_color] = false;
//                         }
//                     }
//                     if (!dlColor[*neighbor.first][old_color] and get(color, *neighbor.first) != old_color) {
//                         bool becomesPossible = true;
//                         bool foundCandidatClique = false;
//                         std::pair<adj_iter, adj_iter> neighneighbor = adjacent_vertices(*neighbor.first, g);
//                         while (!foundCandidatClique and neighneighbor.first != neighneighbor.second) {
//                             if (get(color, *neighneighbor.first) == old_color) {
//                                 foundCandidatClique = true;
//                                 std::pair<adj_iter, adj_iter> temp = neighneighbor;
//                                 temp.first++;
//                                 while (becomesPossible and temp.first != temp.second) {
//                                     if (get(color, *temp.first) == old_color) {
//                                         std::pair<edge_desc, bool> e = edge(*neighneighbor.first, *temp.first, g);
//                                         if (!e.second) {
//                                             becomesPossible = false;
//                                         }
//                                     }
//                                     temp.first++;
//                                 }
//                                 std::pair<adj_iter, adj_iter> neighborCand;
//                                 for (neighborCand = adjacent_vertices(*neighneighbor.first, g); neighborCand.first != neighborCand.second; neighborCand.first++) {
//                                     if (*neighborCand.first != *neighbor.first and get(color, *neighborCand.first) == old_color) {
//                                         std::pair<edge_desc, bool> e = edge(*neighbor.first, *neighborCand.first, g);
//                                         if (!e.second) {
//                                             becomesPossible = false;
//                                         }                            
//                                     }
//                                 }
//                             }
//                             neighneighbor.first++;
//                         }
//                         if (becomesPossible) {
//                             dlColor[*neighbor.first][old_color] = true;
//                         }
//                     }
//                     if (get(color, *neighbor.first) == min_color) {
//                         std::pair<adj_iter, adj_iter> neighneighbor;
//                         for (neighneighbor = adjacent_vertices(*neighbor.first, g); neighneighbor.first != neighneighbor.second; neighneighbor.first++) {
//                             if (*neighneighbor.first != *vp.first and dlColor[*neighneighbor.first][min_color]) {
//                                 bool stayPossibleCand = true;
//                                 std::pair<adj_iter, adj_iter> neighneighneighbor = adjacent_vertices(*neighneighbor.first, g);
//                                 while (stayPossibleCand and neighneighneighbor.first != neighneighneighbor.second) {
//                                     if (*neighneighneighbor.first != *neighbor.first and get(color, *neighneighneighbor.first) == min_color) {
//                                         std::pair<edge_desc, bool> e = edge(*neighbor.first, *neighneighneighbor.first, g);
//                                         if (!e.second) {
//                                             stayPossibleCand = false;
//                                         }
//                                     }
//                                     neighneighneighbor.first++;
//                                 }
//                                 neighneighneighbor = adjacent_vertices(*neighbor.first, g);
//                                 while (stayPossibleCand and neighneighneighbor.first != neighneighneighbor.second) {
//                                     if (*neighneighneighbor.first != *neighneighbor.first and get(color, *neighneighneighbor.first) == min_color) {
//                                         std::pair<edge_desc, bool> e = edge(*neighneighbor.first, *neighneighneighbor.first, g);
//                                         if (!e.second) {
//                                             stayPossibleCand = false;
//                                         }
//                                     }
//                                     neighneighneighbor.first++;
//                                 }
//                                 if (!stayPossibleCand) {
//                                     dlColor[*neighneighbor.first][min_color] = false;
//                                 }
//                             }                            
//                         }                        
//                     }
//                 }
//             }
//         }
//         step--;
//     }
//     if (step == 0) {
//         std::cout << "Warning, uncomplete upgrade, please make more step" << std::endl;
//     }
//     std::vector<int> save_color(k, 0);
//     for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; ++vp.first) {
//         save_color[get(color, *vp.first)-1]++;
//     }
//     int sub_k = 0;
//     for (int i = 1; i <= k; i++) {
//         if (save_color[i-1] == 0) {
//             sub_k++;
//         }
//     }
//     if (sub_k != 0) {
//         k -= sub_k;
//         std::cout << "On utilise " << sub_k << " couleurs en moins que le maximum autorisé." << std::endl;
//     }
// }

bool check_sub_coloring(Graph g, int k) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);
    bool correct_sub_coloring = true;

    std::pair< vertex_iter, vertex_iter > vp = vertices(g);    
    // while (correct_sub_coloring and vp.first != vp.second) {
    while (vp.first != vp.second) {
        std::pair<adj_iter, adj_iter> neighbor = adjacent_vertices(*vp.first, g);
        bool foundCandidatClique = false;
        while (!foundCandidatClique and neighbor.first != neighbor.second) {
            if (get(color, *neighbor.first) == get(color, *vp.first)) {
                foundCandidatClique = true;
                std::pair<adj_iter, adj_iter> temp = neighbor;
                temp.first++;
                // while (correct_sub_coloring and temp.first != temp.second) {
                while (temp.first != temp.second) {
                    if (get(color, *temp.first) == get(color, *vp.first)) {
                        std::pair<edge_desc, bool> e = edge(*neighbor.first, *temp.first, g);
                        if (!e.second) {
                            correct_sub_coloring = false;
                            std::cout << "Here is the error ! step 1 / v = " << get(vertex_id, *vp.first) << std::endl;
                            std::cout << "Il manque l'arête " << get(vertex_id, *neighbor.first) << " " << get(vertex_id, *temp.first) << std::endl;
                        }
                    }
                    temp.first++;
                }
                std::pair<adj_iter, adj_iter> neighborCand = adjacent_vertices(*neighbor.first, g);
                // while (correct_sub_coloring and neighborCand.first != neighborCand.second) {
                while (neighborCand.first != neighborCand.second) {
                    if (*neighborCand.first != *vp.first and get(color, *neighborCand.first) == get(color, *vp.first)) {
                        std::pair<edge_desc, bool> e = edge(*vp.first, *neighborCand.first, g);
                        if (!e.second) {
                            correct_sub_coloring = false;
                            std::cout << "Here is the error ! step 2 / v = " << get(vertex_id, *vp.first) << std::endl;
                            std::cout << "Il manque l'arête " << get(vertex_id, *vp.first) << " " << get(vertex_id, *neighborCand.first) << std::endl;
                        }                            
                    }
                    neighborCand.first++;
                }
            }
            neighbor.first++;
        }
        vp.first++;
    }
    return correct_sub_coloring;
}

void add_color_divide_biggest_clique(Graph& g, std::vector<float>& score, int& k) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);
    property_map< Graph, vertex_color_t >::type color = get(vertex_color, g);
    k++;
    float score_min = 1;
    vertex_desc vertex_min = *vertices(g).first;
    std::pair< vertex_iter, vertex_iter > vp;    
    for (vp = vertices(g); vp.first != vp.second; vp.first++) {
        if (score[get(vertex_id, *vp.first)] < score_min) {
            score_min = score[get(vertex_id, *vp.first)];
            vertex_min = *vp.first;
        }
    }
    if (score_min != 1) {
        std::pair<adj_iter, adj_iter> neighbor;
        int color_min = get(color, vertex_min);
        int nb_min = 1;
        for (neighbor = adjacent_vertices(vertex_min, g); neighbor.first != neighbor.second; neighbor.first++) {
            if (get(color, *neighbor.first) == color_min) {
                nb_min++;
            }
        }
        
        float score_new_color, score_old_color;
        int nb_new_color, nb_old_color;
        if (nb_min % 2 == 0) {
            nb_new_color = nb_min/2;
            nb_old_color = nb_min/2;
            score_new_color = 1/float(nb_new_color);
            score_old_color = 1/float(nb_old_color);
        } else {
            nb_new_color = nb_min/2 + 1;
            nb_old_color = nb_min/2;
            score_new_color = 1/float(nb_new_color);
            score_old_color = 1/float(nb_old_color);
        }
        for (neighbor = adjacent_vertices(vertex_min, g); neighbor.first != neighbor.second; neighbor.first++) {
            if (get(color, *neighbor.first) == color_min) {
                if (nb_new_color != 0) {
                    if (nb_old_color != 0) {
                        int old_or_new = rand()%2;
                        if (old_or_new) {
                            nb_new_color--;
                            put(vertex_color, g, *neighbor.first, k);
                            score[get(vertex_id, *neighbor.first)] = score_new_color;
                        } else {
                            nb_old_color--;
                            score[get(vertex_id, *neighbor.first)] = score_old_color;
                        }
                    } else {
                        nb_new_color--;
                        put(vertex_color, g, *neighbor.first, k);
                        score[get(vertex_id, *neighbor.first)] = score_new_color;
                    }
                } else {
                    nb_old_color--;
                    score[get(vertex_id, *neighbor.first)] = score_old_color;
                }
            }
        }
        if (nb_new_color != 0) {
            nb_new_color--;
            put(vertex_color, g, vertex_min, k);
            score[get(vertex_id, vertex_min)] = score_new_color;
        } else {
            nb_old_color--;
            score[get(vertex_id, vertex_min)] = score_old_color;
        }
    }
}


