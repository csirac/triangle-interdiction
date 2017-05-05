CPPFLAGS=-std=c++11 -Wall -O3
CPLEX_INCLUDE=-I/opt/ibm/ILOG/CPLEX_Studio1263/cplex/include -I/opt/ibm/ILOG/CPLEX_Studio1263/concert/include
CPLEX_LIB=-L/opt/ibm/ILOG/CPLEX_Studio1263/cplex/lib/x86-64_linux/static_pic -L/opt/ibm/ILOG/CPLEX_Studio1263/concert/lib/x86-64_linux/static_pic -lilocplex -lcplex -lconcert -lm -lpthread
CPLEX_FLAGS=${CPLEX_INCLUDE} ${CPLEX_LIB} -DNDEBUG -DILOSTRICTPOD -DIL_STD

all: main_glpk.cpp
	g++ main_glpk.cpp -o grec ${CPPFLAGS} -lglpk -lpthread
test: main_glpk.cpp
	g++ main_glpk.cpp -o grec_test -std=c++11 -O0 -g -lglpk -lpthread
cplex: main_cplex.cpp
	g++ main_cplex.cpp -o grec_cplex ${CPPFLAGS} -lglpk ${CPLEX_FLAGS}
er: er.cpp
	g++ -std=c++11 er.cpp -o er `pkg-config --cflags --libs igraph`
simplify: simplify.cpp
	g++ -std=c++11 simplify.cpp -o simplify `pkg-config --cflags --libs igraph` 
el2bin: el2bin.cpp
	g++ -std=c++11 el2bin.cpp -o el2bin
bin2adj: bin2adj.cpp
	g++ -std=c++11 bin2adj.cpp -o bin2adj -l pthread
small: smallGraph.cpp
	g++ smallGraph.cpp -o smg ${CPPFLAGS} -lglpk -lpthread
