# SubcoloringGraph


Aide temporaire pour l'utilisation du programme. Il y a besoin de la librairie BGL, dont la localisation est à indiquer au début du Makefile.

Commande pour le programme: ./subcol.out ...

Aide :
-h / --help : affiche cette aide !
-i / --input nom_fichier : utilise le graphe en .graphml a l'adresse nom_fichier.
-c / --coord nom_fichier : utilise le graphe en .txt a l'adresse nom_fichier avec des coordonnés (x y) des ap et une portée.
-w / --write nom_fichier : le graphe precedement construit est ecrit dans un fichier nom_fichier.graphml dans le dossier graph.

-u / --udg n taille : utilise un graphe avec n sommets construit en positionnant aleatoirement n points dans un espace de cote taille et en créant le unit disk graph correspondant.
-qu / --quasi_udg n taille dist : utilise un graphe avec n sommets construit en positionnant aleatoirement n points dans un espace de cote taille et en créant le quasi unit disk graph correspondant, avec des arêtes à proba 1 si la distance entre deux sommets est inférieure à dist et des arêtes avec proba proba qui décroit de manière linéaire
-r / --rand n p : utilise un graph d'Erdos Renyi avec n sommets et une probabilité de p (entre 0 et 1) pour les aretes.
-sbm / --stochastic_block_model n r p q : utilise un graphe avec n sommets construit par le Stochastic block model avec r le nombre de communautés, p la proba pour une arête entre membre d'une même communauté et q sinon.
-st / --stadium n long larg : utilise un graphe avec n sommets construits en positionnant de manière régulière des points dans un espace de cote long et larg.
-st2 / --stadiumV2 n taille_cube : utilise un graphe avec n sommets construits en positionnant de manière régulière des points dans un espace de cote taille_cube avec un trou au centre.

-et / --etoile b : graphe sous forme d'étoile avec b branches, les branches sont de longueur 2 et sont toutes reliés à un même centre.
-3d / --3d_quasi_ucg n taille dist : utilise un graphe avec n sommets construit en positionnant aleatoirement n points dans un cube de cote taille et en créant le quasi unit circle graph correspondant, avec des arêtes à proba 1 si la distance entre deux sommets est inférieure à dist et des arêtes avec proba proba qui décroit de manière linéaire

-g / --display_graph : affiche le graphe lors des étapes clés.
-e / --display_edges : affiche les arêtes lors de l'affichage des graphes.
-nv / --non_verbose : empêche l'affichage de multiples informations.
-dl / --degreeofliberty : affiche les degre de liberté des sommets dans le graphe produit.
-clmc / --compute_len_max_clique : affiche une approximation de la taille de la clique max.

-n / --normal : l'algo construit une coloration simple (et pas une sous-coloration).
-m / --mode nom_mode : l'algo suit l'heuristique indiquee, degree ou between, pour favoriser le choix des sommets avec un haut degre ou avec une grande valeur de betweenness.
-k / --color_max k_max : donne le droit d'utiliser jusqu'à k_max couleurs pour la partie upgrade du résultat.

-cm / --coloration_by_markov threshold_mar : donne une coloration compléte en 1-ext en utilisant uniquement le score donnée par les chaines de Markov avec un threshold_mar sur les sommets que l'on garde avant de mettre tous les autres dans une deuxieme étape
-ns / --naive_then_subcol threshold_sub k_naive k_max : donne une coloration compléte en 1-ext en utilisant la technique naive pour k_naive couleurs puis subcoloration sur les sommets de scores inf à threshold_sub avec une upgrade sur k_max couleurs (peut se retrouver à en utiliser plus car subcolo glouton).
-nr / --naive_random k_max : donne le droit d'utiliser jusqu'à k_max couleurs pour un résultat d'allocation naif avec de la minimisation du nombre de conflit.
-exact / --exact : utilise la recherche de tous les MIS pour obtenir la perf des sommets.
-markov / --markov : utilise les chaine de Markov pour obtenir la perf des sommets.

-t / --test replica threshold_sub k_max_naive k_max : construit la courbe des performance moyenne des sommets du graphe lors du ns avec les parametres threshold_sub i allant de 0 à k_max_naive et k_max, le tout avec replica replications pour construire chaque points.
-save / --save_result : enregistre les résultats à l'endroit prévu pour cela
