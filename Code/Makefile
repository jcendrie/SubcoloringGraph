warn = -pedantic -Wall -Wextra -Wold-style-cast -Woverloaded-virtual -Wfloat-equal -Wwrite-strings -Wpointer-arith -Wcast-qual -Wcast-align -Wconversion -Wshadow -Weffc++ -Wredundant-decls -Wdouble-promotion -Winit-self -Wswitch-default -Wswitch-enum -Wundef -Wlogical-op -Winline
localisation_boost = src/
# localisation_glpk = usr/local/lib/

# ifeq (test,$(firstword $(MAKECMDGOALS)))
#   TEST_ARGS := $(wordlist 2,$(words $(MAKECMDGOALS)),$(MAKECMDGOALS))
#   $(eval $(TEST_ARGS):;@:)
# endif

all: subcol.out

subcol.out: main.o coloring.o upgrade_score.o graph_generation.o tools.o readml.o 
	g++ -Wall main.o coloring.o upgrade_score.o graph_generation.o tools.o readml.o -o subcol.out

main.o: main.cpp coloring.hpp upgrade_score.hpp graph_generation.hpp tools.hpp
	g++ -Wall -I "$(localisation_boost)boost_1_81_0/" -c main.cpp -o main.o

coloring.o: coloring.cpp coloring.hpp upgrade_score.o tools.o
	g++ -Wall -I "$(localisation_boost)boost_1_81_0/" -c coloring.cpp -o coloring.o

upgrade_score.o: upgrade_score.cpp upgrade_score.hpp tools.o
	g++ -Wall -I "$(localisation_boost)boost_1_81_0/" -c upgrade_score.cpp -o upgrade_score.o

graph_generation.o: graph_generation.cpp graph_generation.hpp
	g++ -Wall -I "$(localisation_boost)boost_1_81_0/" -c graph_generation.cpp -o graph_generation.o

tools.o: tools.cpp tools.hpp
	g++ -Wall -I "$(localisation_boost)boost_1_81_0/" -c tools.cpp -o tools.o

readml.o: 
	g++ -Wall -I "$(localisation_boost)boost_1_81_0/" -c "$(localisation_boost)boost_1_81_0/libs/graph/src/graphml.cpp" -o readml.o


target_qudg: subcol.out
	for densite1 in 0 ; do \
		for densite2 in 0 1 2 3 4 5 6 7 8 9 ; do \
			./subcol.out -nv -m degree -k 23 -nr 23 -markov -qu 1000 0.1 0.7 ; \
			for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
				./subcol.out -nv -m degree -k 23 -nr 23 -markov -qu 1000 0.$$densite1$$densite2 0.7 ; \
		   		echo $$densite1$$densite2$$number1; \
		   		sleep 1; \
		   	done \
		done \
	done 

target_qubg: subcol.out
	for densite1 in 0 ; do \
		for densite2 in 0 1 2 3 4 5 6 7 8 9 ; do \
			./subcol.out -nv -m degree -k 23 -nr 23 -markov -3d 1000 0.1 0.7 ; \
			for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
				./subcol.out -nv -m degree -k 23 -nr 23 -markov -3d 1000 0.$$densite1$$densite2 0.7 ; \
		   		echo $$densite1$$densite2$$number1; \
		   		sleep 1; \
		   	done \
		done \
	done 

target_er: subcol.out
	for densite1 in 0 ; do \
		for densite2 in 0 1 2 3 4 5 6 7 8 9 ; do \
			./subcol.out -nv -m degree -k 23 -nr 23 -markov -r 1000 0.1 ; \
			for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
				./subcol.out -nv -m degree -k 23 -nr 23 -markov -r 1000 0.$$densite1$$densite2 ; \
		   		echo $$densite1$$densite2$$number1; \
		   		sleep 1; \
		   	done \
		done \
	done 

target_sbm: subcol.out
	for densite1 in 0 ; do \
		for densite2 in 0 1 2 3 4 5 6 7 8 9 ; do \
			./subcol.out -nv -m degree -k 23 -nr 23 -markov -sbm 1000 50 10 0.1 ; \
			for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
				./subcol.out -nv -m degree -k 23 -nr 23 -markov -sbm 1000 50 10 0.$$densite1$$densite2 ; \
		   		echo $$densite1$$densite2$$number1; \
		   		sleep 1; \
		   	done \
		done \
	done 

target_stadium: subcol.out
	for densite1 in 0 ; do \
		for densite2 in 0 1 2 3 4 5 6 7 8 9 ; do \
			./subcol.out -nv -m degree -k 23 -nr 23 -markov -st2 1000 0.1 0.7 ; \
			for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
				./subcol.out -nv -m degree -k 23 -nr 23 -markov -st2 1000 0.$$densite1$$densite2 0.7 ; \
		   		echo $$densite1$$densite2$$number1; \
		   		sleep 1; \
		   	done \
		done \
	done 






# all: subcol.out

# subcol.out: main.o coloring.o upgrade_score.o graph_generation.o tools.o glp_graph.o readml.o 
# 	g++ -Wall -I "$(localisation_glpk)" main.o coloring.o upgrade_score.o graph_generation.o tools.o glp_graph.o readml.o -o subcol.out

# main.o: main.cpp coloring.hpp upgrade_score.hpp graph_generation.hpp tools.hpp glp_graph.hpp
# 	g++ -Wall -I "$(localisation_boost)boost_1_81_0/" -I "$(localisation_glpk)" -c main.cpp -o main.o

# coloring.o: coloring.cpp coloring.hpp upgrade_score.o tools.o
# 	g++ -Wall -I "$(localisation_boost)boost_1_81_0/" -c coloring.cpp -o coloring.o

# upgrade_score.o: upgrade_score.cpp upgrade_score.hpp tools.o
# 	g++ -Wall -I "$(localisation_boost)boost_1_81_0/" -c upgrade_score.cpp -o upgrade_score.o

# graph_generation.o: graph_generation.cpp graph_generation.hpp
# 	g++ -Wall -I "$(localisation_boost)boost_1_81_0/" -c graph_generation.cpp -o graph_generation.o

# tools.o: tools.cpp tools.hpp
# 	g++ -Wall -I "$(localisation_boost)boost_1_81_0/" -c tools.cpp -o tools.o
	
# glp_graph.o: glp_graph.cpp glp_graph.hpp
# 	g++ -Wall -I "$(localisation_boost)boost_1_81_0/" -I "$(localisation_glpk)" -c glp_graph.cpp -o glp_graph.o

# readml.o: 
# 	g++ -Wall -I "$(localisation_boost)boost_1_81_0/" -c "$(localisation_boost)boost_1_81_0/libs/graph/src/graphml.cpp" -o readml.o

# .PHONY: test
# test: all
# 	./subcol.out -i $(TEST_ARGS) -t 4 0.2 10 23

# target_ran: all
# 	for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 		for number2 in 0 1 ; do \
# 	   		./subcol.out -r 1000 0.08 -i "r_1000_0.08_$$number1$$number2" -t 1 0.2 20 23 ; \
# 		done \
# 	done

# target_sbm: all
# 	for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 		for number2 in 0 1 ; do \
# 	   		./subcol.out -sbm 1000 50 0.95 0.06 -i "sbm_1000_50_0.95_0.06_$$number1$$number2" -t 1 0.2 20 23 ; \
# 		done \
# 	done

# target_udg: all
# 	for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 		for number2 in 0 1 ; do \
# 			./subcol.out -u 1000 5.9 -i "udg_1000_5.9_$$number1$$number2" -t 1 0.2 20 23 ; \
# 		done \
# 	done

# target_qudg: all
# 	for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 		for number2 in 0 1 ; do \
# 	   		./subcol.out -qu 1000 5.3 0.7 -i "qudg_1000_5.3_0.7_$$number1$$number2" -t 1 0.2 20 23 ; \
# 		done \
# 	done

# target_sub_ran: all
# 	for densite in 0.06 0.07 0.08 0.09 0.1 0.15 0.2 0.3 0.5 0.9 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			for number2 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 				./subcol.out -r 1000 $$densite -i "r_1000_$$densite" -m degree -k 23 -save -nv ; \
# 				./subcol.out -r 1000 $$densite -i "r_1000_$$densite" -n -k 0 -save -nv ; \
# 			done \
# 		done \
# 	done

# target_sub_sbm: all
# 	for densite in 0.04 0.05 0.06 0.07 0.8 0.1 0.15 0.2 0.25 0.29 0.49; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			for number2 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 				./subcol.out -sbm 1000 50 0.95 $$densite -i "sbm_1000_50_0.95_$$densite" -m degree -k 23 -save -nv ; \
# 				./subcol.out -sbm 1000 50 0.95 $$densite -i "sbm_1000_50_0.95_$$densite" -n -k 0 -save -nv ; \
# 			done \
# 		done \
# 	done

# target_sub_udg: all
# 	for densite in 8.5 7.5 6.8 6.2 5.9 5.5 5.1 4.8 4.5 4 3.5 2.8 2 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			for number2 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 				./subcol.out -u 1000 $$densite -i "udg_1000_$$densite" -m degree -k 23 -save -nv ; \
# 				./subcol.out -u 1000 $$densite -i "udg_1000_$$densite" -n -k 0 -save -nv ; \
# 			done \
# 		done \
# 	done

# target_sub_qudg: all
# 	for densite in 7.8 7 6.25 5.8 5.3 5 4.7 4.3 3.8 3.2 2.5 1.8 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			for number2 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 				./subcol.out -qu 1000 $$densite 0.7 -i "qudg_1000_0.7_$$densite" -m degree -k 23 -save -nv ; \
# 				./subcol.out -qu 1000 $$densite 0.7 -i "qudg_1000_0.7_$$densite" -n -k 0 -save -nv ; \
# 			done \
# 		done \
# 	done

# target_sub_stadium: all
# 	for densite in 9 8.1 7.2 7 6.3 5.95 5.65 5 4.4 3.7 3.2 2.8 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			for number2 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 				./subcol.out -st2 1000 $$densite -i "stadium_1000_$$densite" -m degree -k 23 -save -nv ; \
# 				./subcol.out -st2 1000 $$densite -i "stadium_1000_$$densite" -n -k 0 -save -nv ; \
# 			done \
# 		done \
# 	done

# target_sub_3D: all
# 	for densite in 3.6 3.4 3.2 3.05 2.9 2.8 2.4 2.1 1.9 1.4 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			for number2 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 				./subcol.out -3d 1000 $$densite 0.7 -i "3d_1000_0.7_$$densite" -m degree -k 23 -save -nv ; \
# 				./subcol.out -3d 1000 $$densite 0.7 -i "3d_1000_0.7_$$densite" -n -k 0 -save -nv ; \
# 			done \
# 		done \
# 	done

# target_sub_ran_100: all
# 	for densite in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			./subcol.out -r 100 $$densite -i "r_10_100_$$densite" -m degree -k 10 -save -nv ; \
# 			./subcol.out -r 100 $$densite -i "r_23_100_$$densite" -m degree -k 23 -save -nv ; \
# 		done \
# 	done \

# target_sub_ran_200: all
# 	for densite in 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			./subcol.out -r 200 $$densite -i "r_10_200_$$densite" -m degree -k 10 -save -nv ; \
# 			./subcol.out -r 200 $$densite -i "r_23_200_$$densite" -m degree -k 23 -save -nv ; \
# 		done \
# 	done \

# target_sub_ran_500: all
# 	for densite in 0.02 0.04 0.06 0.08 0.1 0.12 0.14 0.16 0.18 0.2 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			./subcol.out -r 500 $$densite -i "r_10_500_$$densite" -m degree -k 10 -save -nv ; \
# 			./subcol.out -r 500 $$densite -i "r_23_500_$$densite" -m degree -k 23 -save -nv ; \
# 		done \
# 	done \

# target_sub_ran_1000: all
# 	for densite in 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			./subcol.out -r 1000 $$densite -i "r_10_1000_$$densite" -m degree -k 10 -save -nv ; \
# 			./subcol.out -r 1000 $$densite -i "r_23_1000_$$densite" -m degree -k 23 -save -nv ; \
# 		done \
# 	done \

# target_sub_ran_2000: all
# 	for densite in 0.005 0.01 0.015 0.02 0.025 0.03 0.035 0.04 0.045 0.05 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			./subcol.out -r 2000 $$densite -i "r_10_2000_$$densite" -m degree -k 10 -save -nv ; \
# 			./subcol.out -r 2000 $$densite -i "r_23_2000_$$densite" -m degree -k 23 -save -nv ; \
# 		done \
# 	done \

# target_sub_ran_5000: all
# 	for densite in 0.002 0.004 0.006 0.008 0.01 0.012 0.014 0.016 0.018 0.02 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			./subcol.out -r 5000 $$densite -i "r_10_5000_$$densite" -m degree -k 10 -save -nv ; \
# 			./subcol.out -r 5000 $$densite -i "r_23_5000_$$densite" -m degree -k 23 -save -nv ; \
# 		done \
# 	done \

# target_sub_qudg_100: all
# 	for densite in 4.3 3 2.3 1.9 1.67 1.5 1.32 1.1 0.9 0.5 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			./subcol.out -qu 100 $$densite 0.7 -i "qudg_10_100_0.7_$$densite" -m degree -k 10 -save -nv ; \
# 			./subcol.out -qu 100 $$densite 0.7 -i "qudg_23_100_0.7_$$densite" -m degree -k 23 -save -nv ; \
# 		done \
# 	done \

# target_sub_qudg_200: all
# 	for densite in 6.3 4.4 3.5 3 2.6 2.3 2.1 1.9 1.8 1.65 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			./subcol.out -qu 200 $$densite 0.7 -i "qudg_10_200_0.7_$$densite" -m degree -k 10 -save -nv ; \
# 			./subcol.out -qu 200 $$densite 0.7 -i "qudg_23_200_0.7_$$densite" -m degree -k 23 -save -nv ; \
# 		done \
# 	done \

# target_sub_qudg_500: all
# 	for densite in 10 7.2 5.8 5 4.3 3.9 3.6 3.4 3.15 2.9 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			./subcol.out -qu 500 $$densite 0.7 -i "qudg_10_500_0.7_$$densite" -m degree -k 10 -save -nv ; \
# 			./subcol.out -qu 500 $$densite 0.7 -i "qudg_23_500_0.7_$$densite" -m degree -k 23 -save -nv ; \
# 		done \
# 	done \

# target_sub_qudg_1000: all
# 	for densite in 15 10.3 8.3 7.3 6.4 5.8 5.3 5 4.7 4.4 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			./subcol.out -qu 1000 $$densite 0.7 -i "qudg_10_1000_0.7_$$densite" -m degree -k 10 -save -nv ; \
# 			./subcol.out -qu 1000 $$densite 0.7 -i "qudg_23_1000_0.7_$$densite" -m degree -k 23 -save -nv ; \
# 		done \
# 	done \

# target_sub_qudg_2000: all
# 	for densite in 23 15 12 10 9.4 8.5 7.7 7.25 6.8 6.35 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			./subcol.out -qu 2000 $$densite 0.7 -i "qudg_10_2000_0.7_$$densite" -m degree -k 10 -save -nv ; \
# 			./subcol.out -qu 2000 $$densite 0.7 -i "qudg_23_2000_0.7_$$densite" -m degree -k 23 -save -nv ; \
# 		done \
# 	done \

# target_sub_qudg_5000: all
# 	for densite in 40 26 21 17 15.3 14 12.8 12 11 10.3 ; do \
# 		for number1 in 0 1 2 3 4 5 6 7 8 9 ; do \
# 			./subcol.out -qu 5000 $$densite 0.7 -i "qudg_10_5000_0.7_$$densite" -m degree -k 10 -save -nv ; \
# 			./subcol.out -qu 5000 $$densite 0.7 -i "qudg_23_5000_0.7_$$densite" -m degree -k 23 -save -nv ; \
# 		done \
# 	done \




test: all
# 	./subcol.out -i "fixed_graph/hexagraph.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/r_1000_0.03.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/r_1000_0.04.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/r_1000_0.05.graphml" -t 3 0.2 15 23 
# 	./subcol.out -i "fixed_graph/r_1000_0.06.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/r_1000_0.07.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/r_1000_0.08.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/r_1000_0.09.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/udg_1000_9.8.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/udg_1000_8.5.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/udg_1000_7.5.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/udg_1000_6.8.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/udg_1000_6.2.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/udg_1000_5.9.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/udg_1000_5.5.graphml" -t 3 0.2 15 23

# 	./subcol.out -i "fixed_graph/sbm_1000_50_0.95_0.01.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/sbm_1000_50_0.95_0.02.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/sbm_1000_50_0.95_0.03.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/sbm_1000_50_0.95_0.04.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/sbm_1000_50_0.95_0.05.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/sbm_1000_50_0.95_0.06.graphml" -t 3 0.2 15 23
# 	./subcol.out -i "fixed_graph/sbm_1000_50_0.95_0.07.graphml" -t 3 0.2 15 23

clean:
	rm -f *.o

veryclean: clean
	rm -f *.out
