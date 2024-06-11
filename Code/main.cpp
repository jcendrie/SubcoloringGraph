#include <iostream> // for std::cout
#include <stdlib.h>
#include <time.h>
#include <boost/graph/adjacency_list.hpp>
#include <chrono>
#include <cstring>
#include <string>
#include "coloring.hpp"
// #include "glp_graph.hpp"
#include "upgrade_score.hpp"
#include "graph_generation.hpp"
#include "tools.hpp"


using namespace boost;

typedef adjacency_list< vecS, vecS, undirectedS, property< vertex_color_t, int > > Graph;

void usage() {
    std::cout << "Aide :" << std::endl;
    std::cout << "  -h / --help : affiche cette aide !" << std::endl;
    std::cout << "  -i / --input nom_fichier : utilise le graphe en .graphml a l'adresse nom_fichier." << std::endl;
    std::cout << "  -c / --coord nom_fichier : utilise le graphe en .txt a l'adresse nom_fichier avec des coordonnés (x y) des ap et une portée." << std::endl;
    std::cout << "  -w / --write nom_fichier : le graphe precedement construit est ecrit dans un fichier nom_fichier.graphml dans le dossier graph." << std::endl;
    std::cout << std::endl;
    std::cout << "  -u / --udg n dist : utilise un graphe avec n sommets construit en positionnant aleatoirement n points dans un espace de cote taille et en créant le unit disk graph correspondant." << std::endl;
    std::cout << "  -qu / --quasi_udg n d dist : utilise un graphe avec n sommets construit en positionnant aleatoirement n points dans un carre et en créant le quasi unit disk graph correspondant de densite d, avec des arêtes à proba 1 si la distance entre deux sommets est inférieure à dist et des arêtes avec proba proba qui décroit de manière linéaire" << std::endl;
    std::cout << "  -r / --rand n p : utilise un graph d'Erdos Renyi avec n sommets et une probabilité de p (entre 0 et 1) pour les aretes." << std::endl;
    std::cout << "  -sbm / --stochastic_block_model n r rpq d : utilise un graphe avec n sommets construit par le Stochastic block model avec r le nombre de communautés, p la proba pour une arête entre membre d'une même communauté et q sinon, trouvée à partir de rpq le rapport de p par q et avec la densité d du graphe." << std::endl;
    std::cout << "  -st / --stadium n long larg : utilise un graphe avec n sommets construits en positionnant de manière régulière des points dans un espace de cote long et larg." << std::endl;
    std::cout << "  -st2 / --stadiumV2 n d dist : utilise un graphe avec n sommets construits en positionnant de manière régulière des points dans un carre avec un trou au centre, de densite d." << std::endl;
    std::cout << "  -et / --etoile b : graphe sous forme d'étoile avec b branches, les branches sont de longueur 2 et sont toutes reliés à un même centre. (n = 2*b+1)" << std::endl;
    std::cout << "  -ubg / --unit_ball_graph n d : utilise un graphe avec n sommets construit en positionnant aleatoirement n points dans un cube et en créant le quasi unit circle graph correspondant de densite d, avec des arêtes à proba 1 si la distance entre deux sommets est inférieure à 1" << std::endl;
    std::cout << "  -3d / --3d_quasi_ucg n d dist : utilise un graphe avec n sommets construit en positionnant aleatoirement n points dans un cube et en créant le quasi unit circle graph correspondant de densite d, avec des arêtes à proba 1 si la distance entre deux sommets est inférieure à dist et des arêtes avec proba proba qui décroit de manière linéaire" << std::endl;
    std::cout << std::endl;
    std::cout << "  -g / --display_graph : affiche le graphe lors des étapes clés." << std::endl;
    std::cout << "  -e / --display_edges : affiche les arêtes lors de l'affichage des graphes." << std::endl;
    std::cout << "  -nv / --non_verbose : empêche l'affichage de multiples informations." << std::endl;
    std::cout << "  -dl / --degreeofliberty : affiche les degre de liberté des sommets dans le graphe produit." << std::endl;
    std::cout << "  -clmc / --compute_len_max_clique : affiche une approximation de la taille de la clique max." << std::endl;
    std::cout << std::endl;
    std::cout << "  -n / --normal : l'algo construit une coloration simple (et pas une sous-coloration)." << std::endl;
    std::cout << "  -m / --mode nom_mode : l'algo suit l'heuristique indiquee, degree ou between, pour favoriser le choix des sommets avec un haut degre ou avec une grande valeur de betweenness." << std::endl;
    std::cout << "  -k / --color_max k_max : donne le droit d'utiliser jusqu'à k_max couleurs pour la partie upgrade du résultat." << std::endl;
    std::cout << std::endl;
    std::cout << "  -cm / --coloration_by_markov threshold_mar : donne une coloration compléte en 1-ext en utilisant uniquement le score donnée par les chaines de Markov avec un threshold_mar sur les sommets que l'on garde avant de mettre tous les autres dans une deuxieme étape" << std::endl;
    std::cout << "  -ns / --naive_then_subcol threshold_sub k_naive k_max : donne une coloration compléte en 1-ext en utilisant la technique naive pour k_naive couleurs puis subcoloration sur les sommets de scores inf à threshold_sub avec une upgrade sur k_max couleurs (peut se retrouver à en utiliser plus car subcolo glouton)." << std::endl;
    std::cout << "  -nr / --naive_random k_max : donne le droit d'utiliser jusqu'à k_max couleurs pour un résultat d'allocation naif avec de la minimisation du nombre de conflit." << std::endl;
    std::cout << "  -exact / --exact : utilise la recherche de tous les MIS pour obtenir la perf des sommets." << std::endl;
    std::cout << "  -markov / --markov : utilise les chaine de Markov pour obtenir la perf des sommets." << std::endl;
    std::cout << std::endl;
    // std::cout << "  -s / --solve k : résolution par programmation linéaire (algo exponentielle) avec k couleurs max." << std::endl;
    // std::cout << "  -f / --fragment : résolution par programmation linéaire avec frament de contraintes." << std::endl;
    // std::cout << "  -min / --minimalize min_d : utilise la résolution par programmation linéaire en minimisant le nombre de sommets des grandes cliques (borne avec min_d)." << std::endl;
    // std::cout << "  -mc / --min_conflict min_c : utilise la résolution par programmation linéaire en minimisant le nombre de conflits (borne avec min_d)." << std::endl;
    // std::cout << std::endl;
    std::cout << "  -t / --test replica threshold_sub k_max_naive k_max : construit la courbe des performance moyenne des sommets du graphe lors du ns avec les parametres threshold_sub i allant de 0 à k_max_naive et k_max, le tout avec replica replications pour construire chaque points." << std::endl;
    std::cout << "  -save / --save_result : enregistre les résultats à l'endroit prévu pour cela" << std::endl;
    exit(-1);
}

//sbn graph pour cluster séparés par un pont 

int main(int argc, char** argv) {
    int seed = time(NULL);
    srand(seed);
    auto start = std::chrono::high_resolution_clock::now();
    Graph g;
    std::string mode = "degree";
    int k = 1;
    int k_max = 1;
    float threshold_mar = 0.;
    float threshold_sub = 0.;
    int k_naive = 1;
    int k_max_naive = 0;
    int replica = 1;
    // int min_d = 1;
    // int min_c = 1;
    char name_in[1000] = "";
    char name_out[1000];
    std::string type = "";

    bool init_graph = false;
    bool write_graph = false;
    bool display_graph = false;
    bool display_edges = false;
    bool verbose = true;
    bool dl = false;
    bool clmc = false;
    bool normal = false;
    bool approx = false;
    bool upgrade = false;
    bool color_by_markov = false;
    bool naive = false;
    bool naive_then_subcol = false;
    bool exact = false;
    bool markov = false;
    bool solve_graph = false;
    // bool fragment = false;
    // bool minimalize = false;
    // bool min_conflict = false;
    bool test = false;
    bool save = false;

    bool optimal_found = true;

    for (int i=1; i<argc; i++) {
        if ((strcmp(argv[i],"-h")==0) || (strcmp(argv[i],"--help")==0) ) {
            usage();
        } else if ((strcmp(argv[i],"-i")==0) || (strcmp(argv[i],"--input")==0) ) {
            if (i==argc-1) {
                usage();
            }
            sscanf(argv[++i], "%s", name_in);
            if (!init_graph) {
                init_graph = true;
                initGraph(g, name_in);
            }
        } else if ((strcmp(argv[i],"-c")==0) || (strcmp(argv[i],"--coord")==0) ) {
            if (i==argc-1) {
                usage();
            }
            i++;
            if (!init_graph) {            
                init_graph = true;
                g = read_coord(argv[i]);
            }
        } else if ((strcmp(argv[i],"-w")==0) || (strcmp(argv[i],"--write")==0) ) {
            if (i==argc-1) {
                usage();
            }
            write_graph = true;
            sscanf(argv[++i], "%s", name_out);
        } else if ((strcmp(argv[i],"-u")==0) || (strcmp(argv[i],"--udg")==0) ) {
            if (i==argc-2) {
                usage();
            }
            int n;
            float d;
            sscanf(argv[++i], "%d", &n);
            sscanf(argv[++i], "%f", &d);
            if (!init_graph) {            
                init_graph = true;
                g = rand_UDG_density(n, d);
                type = "UDG";
            }
        } else if ((strcmp(argv[i],"-qu")==0) || (strcmp(argv[i],"--quasi_udg")==0) ) {
            if (i==argc-3) {
                usage();
            }
            int n;
            float d;
            float dist;
            sscanf(argv[++i], "%d", &n);
            sscanf(argv[++i], "%f", &d);
            sscanf(argv[++i], "%f", &dist);
            if (!init_graph) {
                init_graph = true;
                g = rand_quasi_UDG_density(n, d, dist);
                type = "QUDG";
            }
        } else if ((strcmp(argv[i],"-r")==0) || (strcmp(argv[i],"--rand")==0) ) {
            if (i==argc-2) {
                usage();
            }
            int n;
            float d;
            sscanf(argv[++i], "%d", &n);
            sscanf(argv[++i], "%f", &d);
            if (!init_graph) {            
                init_graph = true;
                g = rand_graph_bino_density(n, d);
                type = "ER";
            }
        } else if ((strcmp(argv[i],"-sbm")==0) || (strcmp(argv[i],"--stochastic_block_model")==0) ) {
            if (i==argc-4) {
                usage();
            }
            int n;
            int r;
            float rpq;
            float d;
            sscanf(argv[++i], "%d", &n);
            sscanf(argv[++i], "%d", &r);
            sscanf(argv[++i], "%f", &rpq);
            sscanf(argv[++i], "%f", &d);
            if (!init_graph) {            
                init_graph = true;
                g = rand_sbm_density(n, r, rpq, d);
                type = "SBM";
            }
        } else if ((strcmp(argv[i],"-st")==0) || (strcmp(argv[i],"--stadium")==0) ) {
            if (i==argc-3) {
                usage();
            }
            int n;
            float l;
            float larg;
            sscanf(argv[++i], "%d", &n);
            sscanf(argv[++i], "%f", &l);
            sscanf(argv[++i], "%f", &larg);
            if (!init_graph) {            
                init_graph = true;
                g = stadium(n, l, larg);
            }
        } else if ((strcmp(argv[i],"-st2")==0) || (strcmp(argv[i],"--stadiumV2")==0) ) {
            if (i==argc-3) {
                usage();
            }
            int n;
            float d;
            float dist;
            sscanf(argv[++i], "%d", &n);
            sscanf(argv[++i], "%f", &d);
            sscanf(argv[++i], "%f", &dist);
            if (!init_graph) {            
                init_graph = true;
                g = stadium_density(n, d, dist);
                type = "ST";
            }
        } else if ((strcmp(argv[i],"-et")==0) || (strcmp(argv[i],"--etoile")==0) ) {
            if (i==argc-1) {
                usage();
            }
            int b;
            sscanf(argv[++i], "%d", &b);
            if (!init_graph) {            
                init_graph = true;
                g = etoile(b);
            }
        } else if ((strcmp(argv[i],"-ubg")==0) || (strcmp(argv[i],"--unit_ball_graph")==0) ) {
            if (i==argc-2) {
                usage();
            }
            int n;
            float d;
            sscanf(argv[++i], "%d", &n);
            sscanf(argv[++i], "%f", &d);
            if (!init_graph) {
                init_graph = true;
                g = rand_ubg_density(n, d);
                type = "UBG";
            }
        } else if ((strcmp(argv[i],"-3d")==0) || (strcmp(argv[i],"--3d_quasi_ucg")==0) ) {
            if (i==argc-3) {
                usage();
            }
            int n;
            float d;
            float dist;
            sscanf(argv[++i], "%d", &n);
            sscanf(argv[++i], "%f", &d);
            sscanf(argv[++i], "%f", &dist);
            if (!init_graph) {
                init_graph = true;
                g = rand_3D_density(n, d, dist);
                type = "QUBG";
            }
        } else if ((strcmp(argv[i],"-g")==0) || (strcmp(argv[i],"--display_graph")==0) ) {
            display_graph = true;
        } else if ((strcmp(argv[i],"-e")==0) || (strcmp(argv[i],"--display_edges")==0) ) {
            display_edges = true;
        } else if ((strcmp(argv[i],"-nv")==0) || (strcmp(argv[i],"--non_verbose")==0) ) {
            verbose = false;
        } else if ((strcmp(argv[i],"-dl")==0) || (strcmp(argv[i],"--degreeofliberty")==0) ) {
            dl = true;
        } else if ((strcmp(argv[i],"-clmc")==0) || (strcmp(argv[i],"--compute_len_max_clique")==0) ) {
            clmc = true;
        } else if ((strcmp(argv[i],"-n")==0) || (strcmp(argv[i],"--normal")==0) ) {
            std::cout << "Mode with only a normal coloring" << std::endl;
            normal = true;
        } else if ((strcmp(argv[i],"-m")==0) || (strcmp(argv[i],"--mode")==0) ) {
            if (i==argc-1) {
                usage();
            }
            i++;
            approx = true;
            if (strcmp(argv[i], "degree") != 0 and strcmp(argv[i], "between") != 0) {
                std::cout << "Bad mode in input. Use of the degree mode" << std::endl;
            } else {
               mode = argv[i];
            }
        } else if ((strcmp(argv[i],"-k")==0) || (strcmp(argv[i],"--color_max")==0) ) {
            if (i==argc-1) {
                usage();
            }
            sscanf(argv[++i], "%d", &k_max);
            upgrade = true;
        } else if ((strcmp(argv[i],"-cm")==0) || (strcmp(argv[i],"--coloration_by_markov")==0) ) {
            if (i==argc-1) {
                usage();
            }
            sscanf(argv[++i], "%f", &threshold_mar);
            color_by_markov = true;
        } else if ((strcmp(argv[i],"-nr")==0) || (strcmp(argv[i],"--naive_random")==0) ) {
            if (i==argc-1) {
                usage();
            }
            sscanf(argv[++i], "%d", &k_max);
            naive = true;
        } else if ((strcmp(argv[i],"-ns")==0) || (strcmp(argv[i],"--naive_then_subcol")==0) ) {
            if (i==argc-3) {
                usage();
            }
            sscanf(argv[++i], "%f", &threshold_sub);
            sscanf(argv[++i], "%d", &k_naive);
            sscanf(argv[++i], "%d", &k_max);
            naive_then_subcol = true;
        } else if ((strcmp(argv[i],"-exact")==0) || (strcmp(argv[i],"--exact")==0) ) {
            exact = true;
        } else if ((strcmp(argv[i],"-markov")==0) || (strcmp(argv[i],"--markov")==0) ) {
            markov = true;
        // } else if ((strcmp(argv[i],"-s")==0) || (strcmp(argv[i],"--solve")==0) ) {
        //     if (i==argc-1) {
        //         usage();
        //     }
        //     solve_graph = true;
        //     sscanf(argv[++i], "%d", &k);
        // } else if ((strcmp(argv[i],"-f")==0) || (strcmp(argv[i],"--fragment")==0) ) {
        //     fragment = true;
        // } else if ((strcmp(argv[i],"-min")==0) || (strcmp(argv[i],"--minimalize")==0) ) {
        //     if (i==argc-1) {
        //         usage();
        //     }
        //     sscanf(argv[++i], "%d", &min_d);
        //     minimalize = true;
        // } else if ((strcmp(argv[i],"-mc")==0) || (strcmp(argv[i],"--min_conflict")==0) ) {
        //     if (i==argc-1) {
        //         usage();
        //     }
        //     sscanf(argv[++i], "%d", &min_c);
        //     min_conflict = true;
        } else if ((strcmp(argv[i],"-t")==0) || (strcmp(argv[i],"--test")==0) ) {
            if (i==argc-4) {
                usage();
            }
            sscanf(argv[++i], "%d", &replica);
            sscanf(argv[++i], "%f", &threshold_sub);
            sscanf(argv[++i], "%d", &k_max_naive);
            sscanf(argv[++i], "%d", &k_max);
            test = true;
        } else if ((strcmp(argv[i],"-save")==0) || (strcmp(argv[i],"--save_result")==0) ) {
            save = true;
        }
    }

    if (!init_graph) {
        if (verbose) {
            std::cout << "No graph in input. Use of the default toygraph" << std::endl;
        }
        g = initGraph();
        init_graph = true;
    }

    if (display_graph) {
        displayGraph(g, display_edges);
    }

    if (clmc) {
        // std::cout << taille_clique_max_Bron(g) << std::endl;
        std::cout << "Approximation de la borne inf de la taille de la clique max du graphe : " << taille_clique_max_Markov(g) << std::endl;
    }

    int n = num_vertices(g);
    int m = num_edges(g);
    std::vector<float> score(n, 0);
    Graph g_ori;
    std::stack<vertex_desc> stack_low_degree;
    // int num_change = extract_low_degree_vertex(g, g_ori, stack_low_degree, k_max);
    int num_change = 0;
    if (num_change != n) {
        if (naive) {
            k = channel_by_arrival(g, k_max);
            if (verbose) {
                statColor(g, k);
            }
            if (exact) {
                score_MIS(g, score, k, "via_min_conflict(exact)");        
            }
            if (markov) {
                score_Markov(g, score, k, "via_min_conflict(markov)");        
            }
            if (display_graph) {
                displayGraph(g, display_edges);
            }
            result_article(seed, type, m/(float)(n*(n-1)/2), "GCA", 0, score);
        }

        if (color_by_markov) {
            k = coloring_by_Markov(g, threshold_mar);
            if (verbose) {
                statColor(g, k);
            }
            score_Markov(g, score, k, "via_markov_coloring");
            if (display_graph) {
                displayGraph(g, display_edges);
            }
        }
     

        if (solve_graph or normal or approx) {
            // if (solve_graph) {
            //     if (minimalize) {
            //         optimal_found = glpk_solve_minimize_clique_k_color(g, k, min_d);
            //     } else if (fragment) {
            //         optimal_found = glpk_solve_graph_k_color_by_fragment(g, k);
            //     } else if (min_conflict) {
            //         optimal_found = glpk_solve_minimize_conflicts(g, k, min_c);
            //     } else {
            //         optimal_found = glpk_solve_graph_k_color(g, k);
            //     }
            //     if (verbose) {
            //         statColor(g, k);
            //     }
            //     if (save) {
            //         score_sub_colo(g, score, k, "");
            //     } else {
            //         score_sub_colo(g, score, k, "avant_upgrade");
            //         // score_Markov(g, score, k, "avant_upgrade(markov)");
            //     }
            // }
            
            if (approx) {
                k = subcoloring(g, mode);
            }

            if (normal) {
                k = subcoloring(g, mode, true);        
            }

            if (verbose) {
                statColor(g, k);
            }

            if (save) {
                score_sub_colo(g, score, k, "");
            } else {
                score_sub_colo(g, score, k, "avant_upgrade");
                // score_Markov(g, score, k, "avant_upgrade(markov)");
            }

            if (display_graph) {
                displayGraph(g, display_edges);
            }

            if (optimal_found and dl) {
                std::vector<std::vector<bool>> dlColor;
                for (int i = 0; i < n; i++) {
                    std::vector<bool> temp(k, false);
                    dlColor.push_back(temp);
                }
                int ddl = computeLibertyDegree(g, k, dlColor);
                std::cout << "Total degree of liberty : " << ddl << std::endl;
            }
        }

        if (upgrade) {
            int k_avant_upgrade = k;

            while (k < k_max) {
                add_color_divide_biggest_clique(g, score, k);
            }

            upgrade_with_dl(g, k);

            if (verbose) {
                statColor(g, k);
            }

            if (save) {
                score_sub_colo(g, score, k, "");
            } else {            
                score_sub_colo(g, score, k, "apres_upgrade");
                  // score_Markov(g, score, k, "apres_upgrade(markov)");
            }

            if (display_graph) {
                displayGraph(g, display_edges);
            }

            if (save) {
                char real_name[1000] = "";
                strcpy(real_name, name_in);
                std::vector<float> score_tri;
                float sum = 0;
                for (int i = 0; i < n; i++) {
                    insertion_tri(score_tri, score[i]);
                    sum += score[i];
                }
                std::vector<float> stat = {sum/float(n), score_tri[0], score_tri[n/4 - 1], score_tri[n/2 - 1], score_tri[3*n/4 - 1], score_tri[n - 1]};
                write_result(real_name, k_avant_upgrade, stat, n, m, normal);
            }
            result_article(seed, type, m/(float)(n*(n-1)/2), "SCA", k_avant_upgrade, score);
        }


        // if (check_sub_coloring(g, k)) {
        //     std::cout << "Le graphe a un sous-coloriage valide." << std::endl;
        // } else {
        //     std::cout << "WARNING : Le graphe n'a pas un sous-coloriage valide !!!" << std::endl;
        // }



        if (naive_then_subcol) {
            k = naive_then_subcoloring(g, threshold_sub, k_naive, k_max);
            if (verbose) {
                statColor(g, k);
            }
            if (exact) {
                score_MIS(g, score, k, "via_ns(exact)");        
            }
            if (markov) {
                score_Markov(g, score, k, "via_ns(markov)");        
            }
            if (display_graph) {
                displayGraph(g, display_edges);
            }
        }



        if (test) {
            char real_name[1000] = "";
            char out_name[1000] = "";
            // for (long long unsigned int i = 12; i < strlen(name_in)-8; i++) {
            //     real_name[i-12] = name_in[i];
            // }
            strcpy(real_name, name_in);
            std::vector<float> tab_max_avg_score(k_max_naive + 1);
            std::vector<float> tab_max_score(k_max_naive + 1);
            std::vector<float> tab_3rd_quartile_score(k_max_naive + 1);
            std::vector<float> tab_med_score(k_max_naive + 1);
            std::vector<float> tab_1st_quartile_score(k_max_naive + 1);
            std::vector<float> tab_min_score(k_max_naive + 1);
            std::vector<float> tab_rep(k_max_naive + 1);
            strcpy(out_name, real_name);
            // strcat(real_name, "_naive_then_sub_colo_");
            for (int k_naive = 0; k_naive <= k_max_naive; k_naive++) {
                std::pair<float, float> max_stat_score(0,  0);
                int rep = 0;
                int k_save = 0;
                char score_name[1000];
                strcpy(score_name, "histo/");
                strcat(score_name, real_name);
                for (int i = 0; i < replica; i++) {
                    char temp_name[1000];
                    strcpy(temp_name, real_name);
                    std::vector<float> score(n, 0);
                    k = naive_then_subcoloring(g, threshold_sub, k_naive, k_max);
                    std::string k_str = std::to_string(k);
                    strcat(temp_name, k_str.c_str());
                    strcat(temp_name, "_");
                    std::string k_naive_str = std::to_string(k_naive);
                    strcat(temp_name, k_naive_str.c_str());
                    strcat(temp_name, "_");
                    std::string tresh_str = std::to_string(threshold_sub);
                    strcat(temp_name, tresh_str.c_str());
                    strcat(temp_name, "_");
                    std::string i_str = std::to_string(i);
                    strcat(temp_name, i_str.c_str());
                    std::pair<float, float> stat_score = score_Markov(g, score, k, temp_name);
                    std::cout << temp_name << " " << k_naive << " " << i << " " << stat_score.first << " " << stat_score.second << " " << k << std::endl;
                    if (stat_score.second > max_stat_score.second) {
                        max_stat_score = stat_score;
                        rep = i;
                        k_save = k;
                    }
                }
                std::vector<float> score;

                std::string k_save_str = std::to_string(k_save);
                strcat(score_name, k_save_str.c_str());
                strcat(score_name, "_");
                std::string k_naive_str = std::to_string(k_naive);
                strcat(score_name, k_naive_str.c_str());
                strcat(score_name, "_");
                std::string tresh_str = std::to_string(threshold_sub);
                strcat(score_name, tresh_str.c_str());
                strcat(score_name, "_");
                std::string rep_str = std::to_string(rep);
                strcat(score_name, rep_str.c_str());
                strcat(score_name, ".txt");
                std::string s;
                std::ifstream histo(score_name);
                if (histo.is_open()) {
                    std::getline(histo, s);
                    while (std::getline(histo, s)) {
                        insertion_tri(score, std::stof(s));
                    }
                    histo.close();
                }
                tab_max_avg_score[k_naive] = max_stat_score.second;
                tab_max_score[k_naive] = score[0];
                tab_3rd_quartile_score[k_naive] = score[n/4 - 1];
                tab_med_score[k_naive] = score[n/2 - 1];
                tab_1st_quartile_score[k_naive] = score[3*n/4 - 1];
                tab_min_score[k_naive] = score[n - 1];
                tab_rep[k_naive] = rep;
            }
            std::vector<std::vector<float>> stat;
            stat.push_back(tab_max_avg_score);
            stat.push_back(tab_max_score);
            stat.push_back(tab_3rd_quartile_score);
            stat.push_back(tab_med_score);
            stat.push_back(tab_1st_quartile_score);
            stat.push_back(tab_min_score);
            stat.push_back(tab_rep);
            write_Results(real_name, n, threshold_sub, k_max_naive, k_max, stat);
        }
    }

    if (num_change != 0) {
        int k_use = back_low_degree_vertex(g, g_ori, stack_low_degree, k_max);
        if (k_use > k) {
            k = k_use;
        }

        if (verbose) {
            statColor(g, k);
        }
        if (display_graph) {
            displayGraph(g, display_edges);
        }
    }

    // char real_name[1000] = "";

    // if (strcmp(name_in, "")) {
    //     for (long long unsigned int i = 12; i < strlen(name_in)-8; i++) {
    //         real_name[i-12] = name_in[i];
    //     }
    //     std::string k_str = std::to_string(k);
    //     strcat(real_name, "_");
    //     strcat(real_name, k_str.c_str());
    //     if (solve_graph or normal or approx) {
    //         strcat(real_name, "_subcoloring");
    //         score_sub_colo(g, score, k, real_name);
    //     } else if (naive_then_subcol) {
    //         strcat(real_name, "_naive_then_subcoloring");
    //         score_Markov(g, score, k, real_name);        
    //     } else if (naive) {
    //         strcat(real_name, "_naive");
    //         score_Markov(g, score, k, real_name);        
    //     } else if (color_by_markov) {
    //         strcat(real_name, "_treshold");
    //         score_Markov(g, score, k, real_name);        
    //     }
    // }






    if (write_graph) {
        char name[1000];
        strcpy(name, "fixed_graph/");
        strcat(name, name_out);
        strcat(name, ".graphml");
        writeGraph(g, name);
    }


    auto stop = std::chrono::high_resolution_clock::now();
 
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);


 
    std::cout << "Time taken by program : " << (float) duration.count()/1000000 << " seconds" << std::endl;

    return 0;
}
