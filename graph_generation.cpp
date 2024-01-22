#include "graph_generation.hpp"

using namespace boost;


Graph rand_graph_bino(int n, float p) {
    Graph g(n);

    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; vp.first++) {
        std::pair< vertex_iter, vertex_iter > vp2 = vp; 
        for (vp2.first++; vp2.first != vp2.second; vp2.first++) {
            if (rand()%100 < p*100) {
                graph_traits< Graph >::edge_descriptor e;
                bool inserted;
                tie(e, inserted) = add_edge(*vp.first, *vp2.first, g);
            }
        }
    }
    std::cout << "Construction of the graph (" << num_vertices(g) << ", " << num_edges(g) << ") complete" << std::endl;
    return g;
}

Graph rand_UDG(int n, float taille_cube) {
    Graph g(n);
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);

    std::vector<std::pair<float, float>> nodes;
    for (int i = 0; i < n; i++) {
        float x = float((rand()%((int)(taille_cube*1000)))) / 1000;
        float y = float((rand()%((int)(taille_cube*1000)))) / 1000;
        nodes.push_back(std::make_pair(x, y));
    }


    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; vp.first++) {
        std::pair< vertex_iter, vertex_iter > vp2 = vp; 
        for (vp2.first++; vp2.first != vp2.second; vp2.first++) {
            std::pair<float, float> xy1 = nodes[get(vertex_id, *vp.first)];
            std::pair<float, float> xy2 = nodes[get(vertex_id, *vp2.first)];
            if ((xy1.first - xy2.first)*(xy1.first - xy2.first) + (xy1.second - xy2.second)*(xy1.second - xy2.second) < 1) {
                graph_traits< Graph >::edge_descriptor e;
                bool inserted;
                tie(e, inserted) = add_edge(*vp.first, *vp2.first, g);
            }
        }
    }
    std::cout << "Construction of the graph (" << num_vertices(g) << ", " << num_edges(g) << ") complete" << std::endl;
    return g;
}

Graph rand_quasi_UDG(int n, float taille_cube, float dist) {
    Graph g(n);
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);

    std::vector<std::pair<float, float>> nodes;
    for (int i = 0; i < n; i++) {
        float x = float((rand()%((int)(taille_cube*1000)))) / 1000;
        float y = float((rand()%((int)(taille_cube*1000)))) / 1000;
        nodes.push_back(std::make_pair(x, y));
    }
    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; vp.first++) {
        std::pair< vertex_iter, vertex_iter > vp2 = vp; 
        for (vp2.first++; vp2.first != vp2.second; vp2.first++) {
            std::pair<float, float> xy1 = nodes[get(vertex_id, *vp.first)];
            std::pair<float, float> xy2 = nodes[get(vertex_id, *vp2.first)];
            float proba = ((double) rand() / (RAND_MAX));
            float temp_dist = sqrt((xy1.first - xy2.first)*(xy1.first - xy2.first) + (xy1.second - xy2.second)*(xy1.second - xy2.second));
            if (temp_dist < 1 - proba + dist * proba) {
                graph_traits< Graph >::edge_descriptor e;
                bool inserted;
                tie(e, inserted) = add_edge(*vp.first, *vp2.first, g);
            }
        }
    }
    std::cout << "Construction of the graph (" << num_vertices(g) << ", " << num_edges(g) << ") complete" << std::endl;
    return g;
}

Graph rand_3D(int n, float taille_cube, float dist) {
    Graph g(n);
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);

    std::vector<std::vector<float>> nodes;
    for (int i = 0; i < n; i++) {
        float x = float((rand()%((int)(taille_cube*1000)))) / 1000;
        float y = float((rand()%((int)(taille_cube*1000)))) / 1000;
        float z = float((rand()%((int)(taille_cube*1000)))) / 1000;
        nodes.push_back({x, y, z});
    }
    for (std::pair<vertex_iter, vertex_iter> vp = vertices(g); vp.first != vp.second; vp.first++) {
        std::pair<vertex_iter, vertex_iter> vp2 = vp; 
        for (vp2.first++; vp2.first != vp2.second; vp2.first++) {
            std::vector<float> xyz1 = nodes[get(vertex_id, *vp.first)];
            std::vector<float> xyz2 = nodes[get(vertex_id, *vp2.first)];
            float proba = ((double) rand() / (RAND_MAX));
            float temp_dist = (xyz1[0] - xyz2[0])*(xyz1[0] - xyz2[0]) + (xyz1[1] - xyz2[1])*(xyz1[1] - xyz2[1]) + (xyz1[2] - xyz2[2])*(xyz1[2] - xyz2[2]);
            if (temp_dist < 1 - proba + dist * proba) {
                graph_traits< Graph >::edge_descriptor e;
                bool inserted;
                tie(e, inserted) = add_edge(*vp.first, *vp2.first, g);
            }
        }
    }
    std::cout << "Construction of the graph (" << num_vertices(g) << ", " << num_edges(g) << ") complete" << std::endl;
    return g;
}

Graph read_coord(std::string name) {
    std::ifstream coordFile;
    coordFile.open(name);
    int n, range;

    std::string nline;
    std::getline(coordFile, nline);
    n = std::stoi(nline);
    std::string rangeline;
    std::getline(coordFile, rangeline);
    range = std::stoi(rangeline);
    Graph g(n);

    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);

    std::vector<std::pair<float, float>> coord; 
    for (int i = 0; i < n; i++) {
        std::string line;
        std::getline(coordFile, line);
        std::istringstream iss(line);
        float x, y;
        if (!(iss >> x >> y)) { break; } // error
        coord.push_back({x, y});
    }
    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; vp.first++) {
        std::pair< vertex_iter, vertex_iter > vp2 = vp; 
        for (vp2.first++; vp2.first != vp2.second; vp2.first++) {
            std::pair<float, float> xy1 = coord[get(vertex_id, *vp.first)];
            std::pair<float, float> xy2 = coord[get(vertex_id, *vp2.first)];
            if ((xy1.first - xy2.first)*(xy1.first - xy2.first) + (xy1.second - xy2.second)*(xy1.second - xy2.second) <= range*range) {
                graph_traits< Graph >::edge_descriptor e;
                bool inserted;
                tie(e, inserted) = add_edge(*vp.first, *vp2.first, g);
            }
        }
    }
    std::cout << "Construction of the graph (" << num_vertices(g) << ", " << num_edges(g) << ") complete" << std::endl;
    coordFile.close();
    return g;
}


Graph rand_sbm(int n, int r, float p, float q) {
    Graph g(n);
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);

    if (r > n/2) {
        std::cout << "Choix de r trop grand par rapport à n..." << std::endl;
        exit(-1);
    }
    float** mat_p = new float*[n];
    for (int i = 0; i < n; i++) {
        mat_p[i] = new float[n];
    }
    float sum = 0.;
    float partition[r];

    for (int i = 0; i < r; i++) {
        partition[i] = ((double) rand() / (RAND_MAX));
        sum += partition[i];
    }
    for (int i = 0; i < r; i++) {
        partition[i] = partition[i]*(n/(2*sum)) + n/(2*(float)r);
        if (i != 0) {
            partition[i] += partition[i-1];
        }
        std::cout << partition[i] << std::endl;
    }
    int compte = 0;
    int ind = 0;
    for (int i = 0; i < n; i++) {
        compte++;
        while (compte > int(partition[ind]) and ind < r - 1) {
            ind++;
        }
        for (int j = i+1; j < n; j++) {
            if (j < int(partition[ind])) {
                mat_p[i][j] = p;
            } else {
                mat_p[i][j] = q;
            }
        }
    }

    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; vp.first++) {
        std::pair< vertex_iter, vertex_iter > vp2 = vp; 
        for (vp2.first++; vp2.first != vp2.second; vp2.first++) {
            if (((double) rand() / (RAND_MAX)) < mat_p[get(vertex_id, *vp.first)][get(vertex_id, *vp2.first)]) {
                graph_traits< Graph >::edge_descriptor e;
                bool inserted;
                tie(e, inserted) = add_edge(*vp.first, *vp2.first, g);
            }
        }
    }
    std::cout << "Construction of the graph (" << num_vertices(g) << ", " << num_edges(g) << ") complete" << std::endl;
    for (int i = 0; i < n; i++) {
        delete mat_p[i];
    }
    return g;
}

Graph stadium(int n, float l, float larg) {
    int x = std::sqrt((l*n)/larg) + 1;
    int y = std::sqrt((larg*n)/l) + 1;
    n = x*y;

    Graph g(n);
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);

    std::vector<std::pair<float, float>> nodes;

    float pasx = l / (float) (x + 1);
    float pasy = larg / (float) (y + 1);

    for (int i = 1; i <= x; i++) {
        for (int j = 1; j <= y; j++) {
            nodes.push_back(std::make_pair(pasx*i, pasy*j));
        }
    }

    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; vp.first++) {
        std::pair< vertex_iter, vertex_iter > vp2 = vp; 
        for (vp2.first++; vp2.first != vp2.second; vp2.first++) {
            std::pair<float, float> xy1 = nodes[get(vertex_id, *vp.first)];
            std::pair<float, float> xy2 = nodes[get(vertex_id, *vp2.first)];
            if ((xy1.first - xy2.first)*(xy1.first - xy2.first) + (xy1.second - xy2.second)*(xy1.second - xy2.second) < 1) {
                graph_traits< Graph >::edge_descriptor e;
                bool inserted;
                tie(e, inserted) = add_edge(*vp.first, *vp2.first, g);
            }
        }
    }
    std::cout << "Construction of the graph (" << num_vertices(g) << ", " << num_edges(g) << ") complete" << std::endl;
    return g;
}

Graph stadiumV2(int n, float taille_cube) {
    int sub_n = n / 8;
    float sub_t = taille_cube / 3;
    int nb = std::sqrt(sub_n) + 1;
    n = 8 * nb * nb;

    Graph g(n);
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);

    std::vector<std::pair<float, float>> nodes;

    float pas = sub_t / (float) (nb);

    for (int x = 0; x < 3; x++) {
        for (int y = 0; y < 3; y++) {
            if (!(x == 1 and y == 1)) {
                for (int i = 0; i < nb; i++) {
                    for (int j = 0; j < nb; j++) {
                        nodes.push_back(std::make_pair(sub_t*x + pas*(i+0.5), sub_t*y + pas*(j+0.5)));
                    }
                }
            }
        }
    }
    
    // for (int i = 0; i < n; i++) {
    //     std::cout << nodes[i].first << " " << nodes[i].second << std::endl;
    // }

    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; vp.first++) {
        std::pair< vertex_iter, vertex_iter > vp2 = vp; 
        for (vp2.first++; vp2.first != vp2.second; vp2.first++) {
            std::pair<float, float> xy1 = nodes[get(vertex_id, *vp.first)];
            std::pair<float, float> xy2 = nodes[get(vertex_id, *vp2.first)];
            if ((xy1.first - xy2.first)*(xy1.first - xy2.first) + (xy1.second - xy2.second)*(xy1.second - xy2.second) < 1) {
                graph_traits< Graph >::edge_descriptor e;
                bool inserted;
                tie(e, inserted) = add_edge(*vp.first, *vp2.first, g);
            }
        }
    }
    std::cout << "Construction of the graph (" << num_vertices(g) << ", " << num_edges(g) << ") complete" << std::endl;
    return g;
}

Graph etoile(int b) {
    Graph g(2*b+1);
    std::pair<vertex_iter, vertex_iter> vp = vertices(g);

    vertex_desc center = *vp.first;
    for (vp.first++; vp.first != vp.second; vp.first += 2) {
        graph_traits< Graph >::edge_descriptor e1;
        bool inserted1;
        tie(e1, inserted1) = add_edge(center, *vp.first, g);
        graph_traits< Graph >::edge_descriptor e2;
        bool inserted2;
        tie(e2, inserted2) = add_edge(*vp.first, *vp.first + 1, g);
    }
    std::cout << "Construction of the graph (" << num_vertices(g) << ", " << num_edges(g) << ") complete" << std::endl;
    return g;
}



