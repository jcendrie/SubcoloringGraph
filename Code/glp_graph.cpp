#include "glp_graph.hpp"

using namespace boost;

void calcul_P3(Graph g, std::stack<std::vector<int>>& l, int& nb) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);

    std::pair< vertex_iter, vertex_iter > s;
    for (s = vertices(g); s.first != s.second; s.first++) {
        std::pair<adj_iter, adj_iter> v;
        for (v = adjacent_vertices(*s.first, g); v.first != v.second; v.first++) {
            std::pair<adj_iter, adj_iter> u = v;
            for (u.first++; u.first != u.second; u.first++) {
                std::pair<edge_desc, bool> e = edge(*v.first, *u.first, g);
                if (!e.second) {
                    std::vector<int> temp = {int(get(vertex_id, *v.first)), int(get(vertex_id, *s.first)), int(get(vertex_id, *u.first))};
                    l.push(temp);
                    nb++;
                }
            }
        }
    }
}

void calcul_P3_init(Graph g, std::stack<std::vector<int>>& l, int& nb) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);

    std::pair< vertex_iter, vertex_iter > s;
    for (s = vertices(g); s.first != s.second; s.first++) {
        int ten_P3 = 0;
        std::pair<adj_iter, adj_iter> v;
        for (v = adjacent_vertices(*s.first, g); v.first != v.second; v.first++) {
            std::pair<adj_iter, adj_iter> u = v;
            for (u.first++; u.first != u.second; u.first++) {
                int r = 0;
                if (ten_P3 > 10) {
                    r = rand()%100;
                }
                std::pair<edge_desc, bool> e = edge(*v.first, *u.first, g);
                if (!e.second and r < 10) {
                    std::vector<int> temp = {int(get(vertex_id, *v.first)), int(get(vertex_id, *s.first)), int(get(vertex_id, *u.first))};
                    l.push(temp);
                    ten_P3++;
                    nb++;
                }                    
            }
        }
    }
}

bool glpk_solve_graph(glp_prob* mip, int k, int n, int nb, std::stack<std::vector<int>> l) {
    bool optimal_found = false;
    int* ia = new int[1+k*(n+nb*3)];
    int* ja = new int[1+k*(n+nb*3)];
    double* ar = new double[1+k*(n+nb*3)];

    glp_add_rows(mip, n + k*nb);
    for (int i = 1; i <= n; i++) {
        glp_set_row_bnds(mip, i, GLP_FX, 1.0, 1.0);
    }
    for (int i = n+1; i <= n + k*nb; i++) {
        glp_set_row_bnds(mip, i, GLP_UP, 0.0, 2.0);
    }

    glp_add_cols(mip, k*n);
    for (int i = 1; i <= k*n; i++) {
        glp_set_col_kind(mip, i, GLP_BV);
        glp_set_obj_coef(mip, i, 1.0);
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            ia[i*k+j+1] = i+1, ja[i*k+j+1] = i*k+j+1, ar[i*k+j+1] = 1.0;
        }
    }

    for (int i = 0; i < nb; i++) {
        if (l.empty()) {
            std::cout << "Error in the count of the P3..." << std::endl;
        } else {
            std::vector<int> elem = l.top();        
            l.pop();
            for (int j = 0; j < k; j++) {
                ia[n*k+1 + (i*k+j)*3] = n+1 + i*k+j, ja[n*k+1 + (i*k+j)*3] = elem[0]*k+j+1, ar[n*k+1 + (i*k+j)*3] = 1.0;            
                ia[n*k+1 + (i*k+j)*3+1] = n+1 + i*k+j, ja[n*k+1 + (i*k+j)*3+1] = elem[1]*k+j+1, ar[n*k+1 + (i*k+j)*3+1] = 1.0;            
                ia[n*k+1 + (i*k+j)*3+2] = n+1 + i*k+j, ja[n*k+1 + (i*k+j)*3+2] = elem[2]*k+j+1, ar[n*k+1 + (i*k+j)*3+2] = 1.0;            
            }
        }
    }

    glp_load_matrix(mip, k*(n+nb*3), ia, ja, ar);
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.presolve = GLP_ON;
    glp_intopt(mip, &parm);
    int error = glp_mip_status(mip);
    if (error == 5) {
        optimal_found = true;
    }
    delete[] ia;
    delete[] ja;
    delete[] ar;
    return optimal_found;
}

bool glpk_solve_graph_k_color(Graph& g, int k) {
    std::stack<std::vector<int>> l;
    bool optimal_found = false;
    int n = num_vertices(g);
    int nb = 0;
    calcul_P3(g, l, nb);

    std::cout << "nb P3 : " << nb << std::endl;

    glp_prob* mip;
    mip = glp_create_prob();
    glp_set_prob_name(mip, "Sub_Coloring");
    glp_set_obj_dir(mip, GLP_MAX);

    optimal_found = glpk_solve_graph(mip, k, n, nb, l);

    double z = glp_mip_obj_val(mip);
    std::cout << std::endl << "z = " << z << std::endl;
    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        for (int j = 1; j <= k; j++) {
            int x = glp_mip_col_val(mip, get(vertex_index, g, *vp.first)*k + j);
            if (x == 1) {
                put(vertex_color, g, *vp.first, j);
            }
        }
    }
    glp_delete_prob(mip);
    return optimal_found;
}

int verif_sol(glp_prob* mip, int k, int n, int nb, std::stack<std::vector<int>> l) {
    int more_frag = 0;
    int num_row = glp_get_num_rows(mip);
    for (int i = 0; i < nb; i++) {
        std::vector<int> elem = l.top();
        l.pop();
        bool valid_P3 = true;
        int j = 1;
        while (j <= k and valid_P3) {
            if (glp_mip_col_val(mip, elem[0]*k+j) + glp_mip_col_val(mip, elem[1]*k+j) + glp_mip_col_val(mip, elem[2]*k+j) > 2) {
                valid_P3 = false;
                more_frag++;
            }
            j++;
        }
        if (!valid_P3) {
            glp_add_rows(mip, k);
            for (j = 1; j <= k; j++) {
                num_row++;
                glp_set_row_bnds(mip, num_row, GLP_UP, 0.0, 2.0);
                int ind[4] = {0, elem[0]*k + j, elem[1]*k + j, elem[2]*k + j};
                double val[4] = {1., 1., 1., 1.};
                glp_set_mat_row(mip, num_row, 3, ind, val);
            }
        }
    }
    return more_frag;
}



bool glpk_solve_graph_k_color_by_fragment(Graph& g, int k) {
    std::stack<std::vector<int>> l;
    std::stack<std::vector<int>> l_init;
    bool optimal_found = false;
    int n = num_vertices(g);
    int nb = 0;
    int nb_init = 0;
    calcul_P3(g, l, nb);
    calcul_P3_init(g, l_init, nb_init);

    std::cout << "nb P3 : " << nb << std::endl;
    std::cout << "nb P3 init : " << nb_init << std::endl;

    glp_prob* mip;
    mip = glp_create_prob();
    glp_set_prob_name(mip, "Sub_Coloring");
    glp_set_obj_dir(mip, GLP_MAX);
    optimal_found = glpk_solve_graph(mip, k, n, nb_init, l_init);
    int more_frag = verif_sol(mip, k, n, nb, l);

    while (optimal_found and more_frag != 0) {
        glp_iocp parm;
        glp_init_iocp(&parm);
        parm.presolve = GLP_ON;
        glp_intopt(mip, &parm);
        int error = glp_mip_status(mip);
        if (error == 5) {
            optimal_found = true;
        }
        more_frag = verif_sol(mip, k, n, nb, l);
        nb_init += more_frag;
    }

    std::cout << "nb P3 end : " << nb_init << std::endl;

    double z = glp_mip_obj_val(mip);
    std::cout << std::endl << "z = " << z << std::endl;
    std::pair< vertex_iter, vertex_iter > vp;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        for (int j = 1; j <= k; j++) {
            int x = glp_mip_col_val(mip, get(vertex_index, g, *vp.first)*k + j);
            if (x == 1) {
                put(vertex_color, g, *vp.first, j);
            }
        }
    }
    glp_delete_prob(mip);
    return optimal_found;
}

bool glpk_solve_minimize_clique_k_color(Graph& g, int k, int min_d) {
    property_map< Graph, vertex_index_t >::type vertex_id = get(vertex_index, g);

    std::stack<std::vector<int>> l;
    bool optimal_found = false;
    int n = num_vertices(g);
    int nb = 0;
    calcul_P3(g, l, nb);

    std::cout << "nb P3 : " << nb << std::endl;

    glp_prob* mip;
    mip = glp_create_prob();
    glp_set_prob_name(mip, "Sub_Coloring_minimize_clique");
    glp_set_obj_dir(mip, GLP_MIN);

    int* ia = new int[1+k*(n+nb*3)+n*n*k];
    int* ja = new int[1+k*(n+nb*3)+n*n*k];
    double* ar = new double[1+k*(n+nb*3)+n*n*k];

    glp_add_cols(mip, k*n + 1);
    for (int i = 1; i <= k*n; i++) {
        glp_set_col_kind(mip, i, GLP_BV);
        glp_set_obj_coef(mip, i, 0.0);
    }
    glp_set_col_name(mip, k*n + 1, "d");
    glp_set_col_bnds(mip, k*n + 1, GLP_LO, min_d+n, 0.0);
    glp_set_obj_coef(mip, k*n + 1, 1.0);

    glp_add_rows(mip, n + k*nb + k*n);
    for (int i = 1; i <= n; i++) {
        glp_set_row_bnds(mip, i, GLP_FX, 1.0, 1.0);
    }
    for (int i = n+1; i <= n + k*nb; i++) {
        glp_set_row_bnds(mip, i, GLP_UP, 0.0, 2.0);
    }
    for (int i = n + k*nb + 1; i <= n + k*nb + k*n; i++) {
        glp_set_row_bnds(mip, i, GLP_UP, 0.0, 0.0);
    }


    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            ia[i*k+j+1] = i+1, ja[i*k+j+1] = i*k+j+1, ar[i*k+j+1] = 1.0;
        }
    }

    for (int i = 0; i < nb; i++) {
        if (l.empty()) {
            std::cout << "Error in the count of the P3..." << std::endl;
        } else {
            std::vector<int> elem = l.top();        
            l.pop();
            for (int j = 0; j < k; j++) {
                ia[n*k+1 + (i*k+j)*3] = n+1 + i*k+j, ja[n*k+1 + (i*k+j)*3] = elem[0]*k+j+1, ar[n*k+1 + (i*k+j)*3] = 1.0;            
                ia[n*k+1 + (i*k+j)*3+1] = n+1 + i*k+j, ja[n*k+1 + (i*k+j)*3+1] = elem[1]*k+j+1, ar[n*k+1 + (i*k+j)*3+1] = 1.0;            
                ia[n*k+1 + (i*k+j)*3+2] = n+1 + i*k+j, ja[n*k+1 + (i*k+j)*3+2] = elem[2]*k+j+1, ar[n*k+1 + (i*k+j)*3+2] = 1.0;            
            }
        }
    }

    int row = n + k*nb;
    int ind = k*(n+nb*3);

    std::pair< vertex_iter, vertex_iter > vp;    
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        for (int i = 0; i < k; i++) {
            row++;
            ind++;
            ia[ind] = row, ja[ind] = get(vertex_id, *vp.first)*k+i+1, ar[ind] = float(n) + 1;
            std::pair<adj_iter, adj_iter> neighbor;
            for (neighbor = adjacent_vertices(*vp.first, g); neighbor.first != neighbor.second; neighbor.first++) {
                ind++;
                ia[ind] = row, ja[ind] = get(vertex_id, *neighbor.first)*k+i+1, ar[ind] = 1.0;
            }
            ind++;            
            ia[ind] = row, ja[ind] = k*n + 1, ar[ind] = -1.0;            
        }
    }

    glp_load_matrix(mip, ind, ia, ja, ar);
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.presolve = GLP_ON;
    glp_intopt(mip, &parm);
    int error = glp_mip_status(mip);
    if (error == 5) {
        optimal_found = true;
    }
    delete[] ia;
    delete[] ja;
    delete[] ar;

    double d = glp_mip_obj_val(mip);
    std::cout << std::endl << "d = " << d - n << std::endl;
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        for (int j = 1; j <= k; j++) {
            int x = glp_mip_col_val(mip, get(vertex_index, g, *vp.first)*k + j);
            if (x == 1) {
                put(vertex_color, g, *vp.first, j);
            }
        }
    }
    glp_delete_prob(mip);
    return optimal_found;
}


bool glpk_solve_minimize_conflicts(Graph& g, int k, int min_conflict) {
    bool optimal_found = false;
    int n = num_vertices(g);
    int m = num_edges(g);

    glp_prob* mip;
    mip = glp_create_prob();
    glp_set_prob_name(mip, "minimize_conflict");
    glp_set_obj_dir(mip, GLP_MIN);

    int* ia = new int[1+k*(n+m*4)+1];
    int* ja = new int[1+k*(n+m*4)+1];
    double* ar = new double[1+k*(n+m*4)+1];

    glp_add_cols(mip, k*(n+m)+1);
    for (int i = 1; i <= k*n; i++) {
        glp_set_col_kind(mip, i, GLP_BV);
        glp_set_obj_coef(mip, i, 0.0);
    }

    for (int i = k*n+1; i <= k*(n+m); i++) {
        glp_set_col_kind(mip, i, GLP_BV);
        glp_set_obj_coef(mip, i, 0.0);
    }
    glp_set_col_bnds(mip, k*(n+m)+1, GLP_LO, min_conflict, 0.0);
    glp_set_obj_coef(mip, k*(n+m)+1, 1.0);


    glp_add_rows(mip, n + k*m + 1);
    for (int i = 1; i <= n; i++) {
        glp_set_row_bnds(mip, i, GLP_FX, 1.0, 1.0);
    }
    for (int i = n+1; i <= n + k*m; i++) {
        glp_set_row_bnds(mip, i, GLP_UP, 0.0, 1.0);
    }
    glp_set_row_bnds(mip, n + k*m + 1, GLP_FX, 0.0, 0.0);


    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            ia[i*k+j+1] = i+1, ja[i*k+j+1] = i*k+j+1, ar[i*k+j+1] = 1.0;
        }
    }

    int i = 0;
    typedef property_map<Graph, vertex_index_t>::type IndexMap;
    IndexMap index = get(vertex_index, g);
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ei++) {
        int u = index[source(*ei, g)];
        int v = index[target(*ei, g)];
        if (i > m) {
            std::cout << "Error in the number of edges..." << std::endl;
        } else {
            for (int j = 0; j < k; j++) {
                ia[n*k+1 + (i*k+j)*4] = n+1 + i*k+j, ja[n*k+1 + (i*k+j)*4] = u*k+j+1, ar[n*k+1 + (i*k+j)*4] = 1.0;
                ia[n*k+1 + (i*k+j)*4+1] = n+1 + i*k+j, ja[n*k+1 + (i*k+j)*4+1] = v*k+j+1, ar[n*k+1 + (i*k+j)*4+1] = 1.0;
                ia[n*k+1 + (i*k+j)*4+2] = n+1 + i*k+j, ja[n*k+1 + (i*k+j)*4+2] = n*k+i*k+j+1, ar[n*k+1 + (i*k+j)*4+2] = -1.0;
                ia[n*k+1 + (i*k+j)*4+3] = n + k*m + 1, ja[n*k+1 + (i*k+j)*4+3] = n*k+i*k+j+1, ar[n*k+1 + (i*k+j)*4+3] = 1.0;
            }
            i++;
        }
    }
    ia[1+k*(n+m*4)] = n + k*m + 1, ja[1+k*(n+m*4)] = k*(n+m)+1, ar[1+k*(n+m*4)] = -1.0;
    glp_load_matrix(mip, 1+k*(n+m*4), ia, ja, ar);
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.presolve = GLP_ON;
    glp_intopt(mip, &parm);
    int error = glp_mip_status(mip);
    if (error == 5) {
        optimal_found = true;
    }
    delete[] ia;
    delete[] ja;
    delete[] ar;

    double nbc = glp_mip_obj_val(mip);
    std::cout << std::endl << "nbc = " << nbc << std::endl;
    std::pair< vertex_iter, vertex_iter > vp;    
    for (vp = vertices(g); vp.first != vp.second; ++vp.first) {
        for (int j = 1; j <= k; j++) {
            int x = glp_mip_col_val(mip, get(vertex_index, g, *vp.first)*k + j);
            if (x == 1) {
                put(vertex_color, g, *vp.first, j);
            }
        }
    }
    glp_delete_prob(mip);
    return optimal_found;
}






// Test ss-coloriage pour toygraphe avec 2 couleurs !   
void glpk_toy() {
    glp_prob *lp;
    int ia[1+1000], ja[1+1000];
    double ar[1+1000], z, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12;
    lp = glp_create_prob();
    glp_set_prob_name(lp, "Sub_Coloring");
    glp_set_obj_dir(lp, GLP_MAX);
    glp_add_rows(lp, 20);
    glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 2, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 3, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 4, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 5, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 6, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 7, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 8, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 9, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 10, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 11, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 12, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 13, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 14, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 15, GLP_FX, 1.0, 1.0);
    glp_set_row_bnds(lp, 16, GLP_FX, 1.0, 1.0);
    glp_set_row_bnds(lp, 17, GLP_FX, 1.0, 1.0);
    glp_set_row_bnds(lp, 18, GLP_FX, 1.0, 1.0);
    glp_set_row_bnds(lp, 19, GLP_FX, 1.0, 1.0);
    glp_set_row_bnds(lp, 20, GLP_FX, 1.0, 1.0);

    // glp_set_row_bnds(lp, 21, GLP_LO, 0.0, 0.0);
    glp_add_cols(lp, 12);
    glp_set_col_name(lp, 1, "x1");
    glp_set_col_kind(lp, 1, GLP_BV);
    glp_set_obj_coef(lp, 1, 1.0);
    glp_set_col_name(lp, 2, "x2");
    glp_set_col_kind(lp, 2, GLP_BV);
    glp_set_obj_coef(lp, 2, 1.0);
    glp_set_col_name(lp, 3, "x3");
    glp_set_col_kind(lp, 3, GLP_BV);
    glp_set_obj_coef(lp, 3, 1.0);
    glp_set_col_name(lp, 4, "x4");
    glp_set_col_kind(lp, 4, GLP_BV);
    glp_set_obj_coef(lp, 4, 1.0);
    glp_set_col_name(lp, 5, "x5");
    glp_set_col_kind(lp, 5, GLP_BV);
    glp_set_obj_coef(lp, 5, 1.0);
    glp_set_col_name(lp, 6, "x6");
    glp_set_col_kind(lp, 6, GLP_BV);
    glp_set_obj_coef(lp, 6, 1.0);
    glp_set_col_name(lp, 7, "x7");
    glp_set_col_kind(lp, 7, GLP_BV);
    glp_set_obj_coef(lp, 7, 1.0);
    glp_set_col_name(lp, 8, "x8");
    glp_set_col_kind(lp, 8, GLP_BV);
    glp_set_obj_coef(lp, 8, 1.0);
    glp_set_col_name(lp, 9, "x9");
    glp_set_col_kind(lp, 9, GLP_BV);
    glp_set_obj_coef(lp, 9, 1.0);
    glp_set_col_name(lp, 10, "x10");
    glp_set_col_kind(lp, 10, GLP_BV);
    glp_set_obj_coef(lp, 10, 1.0);
    glp_set_col_name(lp, 11, "x11");
    glp_set_col_kind(lp, 11, GLP_BV);
    glp_set_obj_coef(lp, 11, 1.0);
    glp_set_col_name(lp, 12, "x12");
    glp_set_col_kind(lp, 12, GLP_BV);
    glp_set_obj_coef(lp, 12, 1.0);

    // glp_set_col_name(lp, 13, "d");
    // glp_set_col_kind(lp, 13, GLP_IV);
    // glp_set_obj_coef(lp, 13, 1.0);
    ia[1] = 1, ja[1] = 1, ar[1] = 1.0;
    ia[2] = 1, ja[2] = 2, ar[2] = 1.0;
    ia[3] = 1, ja[3] = 3, ar[3] = 1.0;
    ia[4] = 2, ja[4] = 1, ar[4] = 1.0;
    ia[5] = 2, ja[5] = 5, ar[5] = 1.0;
    ia[6] = 2, ja[6] = 6, ar[6] = 1.0;
    ia[7] = 3, ja[7] = 2, ar[7] = 1.0;
    ia[8] = 3, ja[8] = 5, ar[8] = 1.0;
    ia[9] = 3, ja[9] = 6, ar[9] = 1.0;
    ia[10] = 4, ja[10] = 2, ar[10] = 1.0;
    ia[11] = 4, ja[11] = 5, ar[11] = 1.0;
    ia[12] = 4, ja[12] = 3, ar[12] = 1.0;
    ia[13] = 5, ja[13] = 6, ar[13] = 1.0;
    ia[14] = 5, ja[14] = 5, ar[14] = 1.0;
    ia[15] = 5, ja[15] = 3, ar[15] = 1.0;
    ia[16] = 6, ja[16] = 3, ar[16] = 1.0;
    ia[17] = 6, ja[17] = 5, ar[17] = 1.0;
    ia[18] = 6, ja[18] = 4, ar[18] = 1.0;
    ia[19] = 7, ja[19] = 4, ar[19] = 1.0;
    ia[20] = 7, ja[20] = 5, ar[20] = 1.0;
    ia[21] = 7, ja[21] = 6, ar[21] = 1.0;
    ia[22] = 8, ja[22] = 7, ar[22] = 1.0;
    ia[23] = 8, ja[23] = 8, ar[23] = 1.0;
    ia[24] = 8, ja[24] = 9, ar[24] = 1.0;
    ia[25] = 9, ja[25] = 7, ar[25] = 1.0;
    ia[26] = 9, ja[26] = 11, ar[26] = 1.0;
    ia[27] = 9, ja[27] = 12, ar[27] = 1.0;
    ia[28] = 10, ja[28] = 8, ar[28] = 1.0;
    ia[29] = 10, ja[29] = 11, ar[29] = 1.0;
    ia[30] = 10, ja[30] = 12, ar[30] = 1.0;
    ia[31] = 11, ja[31] = 8, ar[31] = 1.0;
    ia[32] = 11, ja[32] = 11, ar[32] = 1.0;
    ia[33] = 11, ja[33] = 9, ar[33] = 1.0;
    ia[34] = 12, ja[34] = 12, ar[34] = 1.0;
    ia[35] = 12, ja[35] = 11, ar[35] = 1.0;
    ia[36] = 12, ja[36] = 9, ar[36] = 1.0;
    ia[37] = 13, ja[37] = 9, ar[37] = 1.0;
    ia[38] = 13, ja[38] = 11, ar[38] = 1.0;
    ia[39] = 13, ja[39] = 10, ar[39] = 1.0;
    ia[40] = 14, ja[40] = 10, ar[40] = 1.0;
    ia[41] = 14, ja[41] = 11, ar[41] = 1.0;
    ia[42] = 14, ja[42] = 12, ar[42] = 1.0;

    ia[43] = 15, ja[43] = 1, ar[43] = 1.0;
    ia[44] = 15, ja[44] = 7, ar[44] = 1.0;
    ia[45] = 16, ja[45] = 2, ar[45] = 1.0;
    ia[46] = 16, ja[46] = 8, ar[46] = 1.0;
    ia[47] = 17, ja[47] = 3, ar[47] = 1.0;
    ia[48] = 17, ja[48] = 9, ar[48] = 1.0;
    ia[49] = 18, ja[49] = 4, ar[49] = 1.0;
    ia[50] = 18, ja[50] = 10, ar[50] = 1.0;
    ia[51] = 19, ja[51] = 5, ar[51] = 1.0;
    ia[52] = 19, ja[52] = 11, ar[52] = 1.0;
    ia[53] = 20, ja[53] = 6, ar[53] = 1.0;
    ia[54] = 20, ja[54] = 12, ar[54] = 1.0;


    // ia[55] = 21, ja[55] = 1, ar[55] = 1.0;
    // ia[56] = 21, ja[56] = 2, ar[56] = 1.0;
    // ia[57] = 21, ja[57] = 6, ar[57] = 1.0;
    // ia[58] = 22, ja[58] = 1, ar[58] = 1.0;
    // ia[59] = 22, ja[59] = 2, ar[59] = 1.0;
    // ia[60] = 22, ja[60] = 3, ar[60] = 1.0;
    // ia[61] = 22, ja[61] = 6, ar[61] = 1.0;
    // ia[62] = 21, ja[62] = 11, ar[62] = 1.0;
    // ia[63] = 21, ja[63] = 6, ar[63] = 1.0;
    // ia[64] = 21, ja[64] = 12, ar[64] = 1.0;
    // ia[65] = 21, ja[65] = 6, ar[65] = 1.0;
    // ia[66] = 21, ja[66] = 12, ar[66] = 1.0;
    glp_load_matrix(lp, 54, ia, ja, ar);
    glp_simplex(lp, NULL);
    z = glp_get_obj_val(lp);
    x1 = glp_get_col_prim(lp, 1);
    x2 = glp_get_col_prim(lp, 2);
    x3 = glp_get_col_prim(lp, 3);
    x4 = glp_get_col_prim(lp, 4);
    x5 = glp_get_col_prim(lp, 5);
    x6 = glp_get_col_prim(lp, 6);
    x7 = glp_get_col_prim(lp, 7);
    x8 = glp_get_col_prim(lp, 8);
    x9 = glp_get_col_prim(lp, 9);
    x10 = glp_get_col_prim(lp, 10);
    x11 = glp_get_col_prim(lp, 11);
    x12 = glp_get_col_prim(lp, 12);
    printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g; x4 = %g; x5 = %g; x6 = %g; x7 = %g; x8 = %g; x9 = %g; x10 = %g; x11 = %g; x12 = %g\n", z, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12);
    glp_delete_prob(lp);
}


void glpk_toy_v2() {
    glp_prob *lp;
    int ia[1+1000], ja[1+1000];
    double ar[1+1000], z, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12;
    lp = glp_create_prob();
    glp_set_prob_name(lp, "Sub_Coloring");
    glp_set_obj_dir(lp, GLP_MIN);
    glp_add_cols(lp, 13);
    glp_set_col_name(lp, 1, "x1");
    glp_set_col_kind(lp, 1, GLP_BV);
    glp_set_obj_coef(lp, 1, 0.0);
    glp_set_col_name(lp, 2, "x2");
    glp_set_col_kind(lp, 2, GLP_BV);
    glp_set_obj_coef(lp, 2, 0.0);
    glp_set_col_name(lp, 3, "x3");
    glp_set_col_kind(lp, 3, GLP_BV);
    glp_set_obj_coef(lp, 3, 0.0);
    glp_set_col_name(lp, 4, "x4");
    glp_set_col_kind(lp, 4, GLP_BV);
    glp_set_obj_coef(lp, 4, 0.0);
    glp_set_col_name(lp, 5, "x5");
    glp_set_col_kind(lp, 5, GLP_BV);
    glp_set_obj_coef(lp, 5, 0.0);
    glp_set_col_name(lp, 6, "x6");
    glp_set_col_kind(lp, 6, GLP_BV);
    glp_set_obj_coef(lp, 6, 0.0);
    glp_set_col_name(lp, 7, "x7");
    glp_set_col_kind(lp, 7, GLP_BV);
    glp_set_obj_coef(lp, 7, 0.0);
    glp_set_col_name(lp, 8, "x8");
    glp_set_col_kind(lp, 8, GLP_BV);
    glp_set_obj_coef(lp, 8, 0.0);
    glp_set_col_name(lp, 9, "x9");
    glp_set_col_kind(lp, 9, GLP_BV);
    glp_set_obj_coef(lp, 9, 0.0);
    glp_set_col_name(lp, 10, "x10");
    glp_set_col_kind(lp, 10, GLP_BV);
    glp_set_obj_coef(lp, 10, 0.0);
    glp_set_col_name(lp, 11, "x11");
    glp_set_col_kind(lp, 11, GLP_BV);
    glp_set_obj_coef(lp, 11, 0.0);
    glp_set_col_name(lp, 12, "x12");
    glp_set_col_kind(lp, 12, GLP_BV);
    glp_set_obj_coef(lp, 12, 0.0);


    glp_set_col_name(lp, 13, "d");
    glp_set_col_bnds(lp, 13, GLP_LO, 1.0, 10.0);
    glp_set_obj_coef(lp, 13, 1.0);



    glp_add_rows(lp,32);
    glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 2, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 3, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 4, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 5, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 6, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 7, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 8, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 9, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 10, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 11, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 12, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 13, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 14, GLP_UP, 0.0, 2.0);
    glp_set_row_bnds(lp, 15, GLP_FX, 1.0, 1.0);
    glp_set_row_bnds(lp, 16, GLP_FX, 1.0, 1.0);
    glp_set_row_bnds(lp, 17, GLP_FX, 1.0, 1.0);
    glp_set_row_bnds(lp, 18, GLP_FX, 1.0, 1.0);
    glp_set_row_bnds(lp, 19, GLP_FX, 1.0, 1.0);
    glp_set_row_bnds(lp, 20, GLP_FX, 1.0, 1.0);

    glp_set_row_bnds(lp, 21, GLP_UP, 0.0, 0.0);
    glp_set_row_bnds(lp, 22, GLP_UP, 0.0, 0.0);
    glp_set_row_bnds(lp, 23, GLP_UP, 0.0, 0.0);
    glp_set_row_bnds(lp, 24, GLP_UP, 0.0, 0.0);
    glp_set_row_bnds(lp, 25, GLP_UP, 0.0, 0.0);
    glp_set_row_bnds(lp, 26, GLP_UP, 0.0, 0.0);
    glp_set_row_bnds(lp, 27, GLP_UP, 0.0, 0.0);
    glp_set_row_bnds(lp, 28, GLP_UP, 0.0, 0.0);
    glp_set_row_bnds(lp, 29, GLP_UP, 0.0, 0.0);
    glp_set_row_bnds(lp, 30, GLP_UP, 0.0, 0.0);
    glp_set_row_bnds(lp, 31, GLP_UP, 0.0, 0.0);
    glp_set_row_bnds(lp, 32, GLP_UP, 0.0, 0.0);

    ia[1] = 1, ja[1] = 1, ar[1] = 1.0;
    ia[2] = 1, ja[2] = 2, ar[2] = 1.0;
    ia[3] = 1, ja[3] = 3, ar[3] = 1.0;
    ia[4] = 2, ja[4] = 1, ar[4] = 1.0;
    ia[5] = 2, ja[5] = 5, ar[5] = 1.0;
    ia[6] = 2, ja[6] = 6, ar[6] = 1.0;
    ia[7] = 3, ja[7] = 2, ar[7] = 1.0;
    ia[8] = 3, ja[8] = 5, ar[8] = 1.0;
    ia[9] = 3, ja[9] = 6, ar[9] = 1.0;
    ia[10] = 4, ja[10] = 2, ar[10] = 1.0;
    ia[11] = 4, ja[11] = 5, ar[11] = 1.0;
    ia[12] = 4, ja[12] = 3, ar[12] = 1.0;
    ia[13] = 5, ja[13] = 6, ar[13] = 1.0;
    ia[14] = 5, ja[14] = 5, ar[14] = 1.0;
    ia[15] = 5, ja[15] = 3, ar[15] = 1.0;
    ia[16] = 6, ja[16] = 3, ar[16] = 1.0;
    ia[17] = 6, ja[17] = 5, ar[17] = 1.0;
    ia[18] = 6, ja[18] = 4, ar[18] = 1.0;
    ia[19] = 7, ja[19] = 4, ar[19] = 1.0;
    ia[20] = 7, ja[20] = 5, ar[20] = 1.0;
    ia[21] = 7, ja[21] = 6, ar[21] = 1.0;
    ia[22] = 8, ja[22] = 7, ar[22] = 1.0;
    ia[23] = 8, ja[23] = 8, ar[23] = 1.0;
    ia[24] = 8, ja[24] = 9, ar[24] = 1.0;
    ia[25] = 9, ja[25] = 7, ar[25] = 1.0;
    ia[26] = 9, ja[26] = 11, ar[26] = 1.0;
    ia[27] = 9, ja[27] = 12, ar[27] = 1.0;
    ia[28] = 10, ja[28] = 8, ar[28] = 1.0;
    ia[29] = 10, ja[29] = 11, ar[29] = 1.0;
    ia[30] = 10, ja[30] = 12, ar[30] = 1.0;
    ia[31] = 11, ja[31] = 8, ar[31] = 1.0;
    ia[32] = 11, ja[32] = 11, ar[32] = 1.0;
    ia[33] = 11, ja[33] = 9, ar[33] = 1.0;
    ia[34] = 12, ja[34] = 12, ar[34] = 1.0;
    ia[35] = 12, ja[35] = 11, ar[35] = 1.0;
    ia[36] = 12, ja[36] = 9, ar[36] = 1.0;
    ia[37] = 13, ja[37] = 9, ar[37] = 1.0;
    ia[38] = 13, ja[38] = 11, ar[38] = 1.0;
    ia[39] = 13, ja[39] = 10, ar[39] = 1.0;
    ia[40] = 14, ja[40] = 10, ar[40] = 1.0;
    ia[41] = 14, ja[41] = 11, ar[41] = 1.0;
    ia[42] = 14, ja[42] = 12, ar[42] = 1.0;

    ia[43] = 15, ja[43] = 1, ar[43] = 1.0;
    ia[44] = 15, ja[44] = 7, ar[44] = 1.0;
    ia[45] = 16, ja[45] = 2, ar[45] = 1.0;
    ia[46] = 16, ja[46] = 8, ar[46] = 1.0;
    ia[47] = 17, ja[47] = 3, ar[47] = 1.0;
    ia[48] = 17, ja[48] = 9, ar[48] = 1.0;
    ia[49] = 18, ja[49] = 4, ar[49] = 1.0;
    ia[50] = 18, ja[50] = 10, ar[50] = 1.0;
    ia[51] = 19, ja[51] = 5, ar[51] = 1.0;
    ia[52] = 19, ja[52] = 11, ar[52] = 1.0;
    ia[53] = 20, ja[53] = 6, ar[53] = 1.0;
    ia[54] = 20, ja[54] = 12, ar[54] = 1.0;


    ia[55] = 21, ja[55] = 1, ar[55] = 1.0;
    ia[56] = 21, ja[56] = 2, ar[56] = 1.0;
    ia[57] = 21, ja[57] = 6, ar[57] = 1.0;
    ia[58] = 21, ja[58] = 13, ar[58] = -1.0;
    ia[59] = 22, ja[59] = 1, ar[59] = 1.0;
    ia[60] = 22, ja[60] = 2, ar[60] = 1.0;
    ia[61] = 22, ja[61] = 3, ar[61] = 1.0;
    ia[62] = 22, ja[62] = 6, ar[62] = 1.0;
    ia[63] = 22, ja[63] = 13, ar[63] = -1.0;
    ia[64] = 23, ja[64] = 2, ar[64] = 1.0;
    ia[65] = 23, ja[65] = 3, ar[65] = 1.0;
    ia[66] = 23, ja[66] = 5, ar[66] = 1.0;
    ia[67] = 23, ja[67] = 13, ar[67] = -1.0;
    ia[68] = 24, ja[68] = 4, ar[68] = 1.0;
    ia[69] = 24, ja[69] = 5, ar[69] = 1.0;
    ia[70] = 24, ja[70] = 13, ar[70] = -1.0;
    ia[71] = 25, ja[71] = 3, ar[71] = 1.0;
    ia[72] = 25, ja[72] = 4, ar[72] = 1.0;
    ia[73] = 25, ja[73] = 5, ar[73] = 1.0;
    ia[74] = 25, ja[74] = 6, ar[74] = 1.0;
    ia[75] = 25, ja[75] = 13, ar[75] = -1.0;
    ia[76] = 26, ja[76] = 1, ar[76] = 1.0;
    ia[77] = 26, ja[77] = 2, ar[77] = 1.0;
    ia[78] = 26, ja[78] = 5, ar[78] = 1.0;
    ia[79] = 26, ja[79] = 6, ar[79] = 1.0;
    ia[80] = 26, ja[80] = 13, ar[80] = -1.0;

    ia[81] = 27, ja[81] = 7, ar[81] = 1.0;
    ia[82] = 27, ja[82] = 8, ar[82] = 1.0;
    ia[83] = 27, ja[83] = 12, ar[83] = 1.0;
    ia[84] = 27, ja[84] = 13, ar[84] = -1.0;
    ia[85] = 28, ja[85] = 7, ar[85] = 1.0;
    ia[86] = 28, ja[86] = 8, ar[86] = 1.0;
    ia[87] = 28, ja[87] = 9, ar[87] = 1.0;
    ia[88] = 28, ja[88] = 12, ar[88] = 1.0;
    ia[89] = 28, ja[89] = 13, ar[89] = -1.0;
    ia[90] = 29, ja[90] = 8, ar[90] = 1.0;
    ia[91] = 29, ja[91] = 9, ar[91] = 1.0;
    ia[92] = 29, ja[92] = 11, ar[92] = 1.0;
    ia[93] = 29, ja[93] = 13, ar[93] = -1.0;
    ia[94] = 30, ja[94] = 10, ar[94] = 1.0;
    ia[95] = 30, ja[95] = 11, ar[95] = 1.0;
    ia[96] = 30, ja[96] = 13, ar[96] = -1.0;
    ia[97] = 31, ja[97] = 9, ar[97] = 1.0;
    ia[98] = 31, ja[98] = 10, ar[98] = 1.0;
    ia[99] = 31, ja[99] = 11, ar[99] = 1.0;
    ia[100] = 31, ja[100] = 12, ar[100] = 1.0;
    ia[101] = 31, ja[101] = 13, ar[101] = -1.0;
    ia[102] = 32, ja[102] = 7, ar[102] = 1.0;
    ia[103] = 32, ja[103] = 8, ar[103] = 1.0;
    ia[104] = 32, ja[104] = 11, ar[104] = 1.0;
    ia[105] = 32, ja[105] = 12, ar[105] = 1.0;
    ia[106] = 32, ja[106] = 13, ar[106] = -1.0;
    
    glp_load_matrix(lp, 106, ia, ja, ar);
    glp_iocp parm;
    glp_init_iocp(&parm);
    parm.presolve = GLP_ON;
    glp_intopt(lp, &parm);
    z = glp_mip_obj_val(lp);
    x1 = glp_mip_col_val(lp, 1);
    x2 = glp_mip_col_val(lp, 2);
    x3 = glp_mip_col_val(lp, 3);
    x4 = glp_mip_col_val(lp, 4);
    x5 = glp_mip_col_val(lp, 5);
    x6 = glp_mip_col_val(lp, 6);
    x7 = glp_mip_col_val(lp, 7);
    x8 = glp_mip_col_val(lp, 8);
    x9 = glp_mip_col_val(lp, 9);
    x10 = glp_mip_col_val(lp, 10);
    x11 = glp_mip_col_val(lp, 11);
    x12 = glp_mip_col_val(lp, 12);
    float d = glp_mip_col_val(lp, 13);
    printf("\nz = %g; x1 = %g; x2 = %g; x3 = %g; x4 = %g; x5 = %g; x6 = %g; x7 = %g; x8 = %g; x9 = %g; x10 = %g; x11 = %g; x12 = %g; d = %g\n", z, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, d);
    glp_delete_prob(lp);
}
