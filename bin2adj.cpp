#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <sstream>
#include "mygraph.cpp"

using namespace std;
using namespace mygraph;

int main(int argc, char ** argv) {
  if (argc == 1) {
     cerr << "Usage: " << argv[0] << " <input graph in binary edge list>" << endl;
     return 1;
  }
   
  string fname( argv[1] );
  
  string outfile = fname.substr(0, fname.find_last_of( '.' )) + ".adj";

  tinyGraph g;
  g.read_edge_list_bin( fname );

  ofstream ofile( outfile.c_str(), ios::out | ios::binary );
  ofile.write( (char*) &(g.n), sizeof( g.n ) );
  unsigned ss;
  for (unsigned i = 0; i < g.n; ++i) {
     ss = g.adjList[i].neis.size();
     ofile.write( (char*) (&ss), sizeof( ss ) );
     for (unsigned j = 0; j < ss; ++j) {
	ofile.write( (char*) (&(g.adjList[i].neis[j] )), sizeof( uint32_t ) );
     }
  }

  ofile.close();
  
  return 0;
}
