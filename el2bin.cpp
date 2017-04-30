#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <sstream>
using namespace std;

int main(int argc, char ** argv) {
  if (argc == 1) {
     cerr << "Usage: " << argv[0] << " <input graph in plain text edge list>" << endl;
     return 1;
  }
   
  string fname( argv[1] );
  
  string outfile = fname.substr(0, fname.find_last_of( '.' )) + ".bin";

  ifstream ifile( fname.c_str() );
  ofstream ofile( outfile.c_str(), ios::out | ios::binary );
  unsigned from,to;
  vector< unsigned > v_edges;
  unsigned m = 0;
  unsigned n = 0;
  string line;
  istringstream iss;
  while ( getline( ifile, line ) ) {
     if ( line[0]  != '#' && line[0] != ' ') {
	if (line.size() == 0)
	   break;
	  
	iss.clear();
	iss.str( line );
	//need to add an edge 
	iss >> from;
	iss >> to;

	if (from >= n)
	   n = from + 1;
	if (to >= n)
	   n = to + 1;
	
	v_edges.push_back( from );
	v_edges.push_back( to );

	++m;
     }
  }

  //Write number of edges
  ofile.write( (char*) &n, sizeof( n ) );
  ofile.write( (char*) &m, sizeof( m ) );
  for (size_t i = 0; i < m; ++i) {
     //Write each edge
     ofile.write( (char*) &(v_edges[2*i]), sizeof( unsigned ) );
     ofile.write( (char*) &(v_edges[2*i + 1]), sizeof( unsigned ) );
  }

  ifile.close();
  ofile.close();
  
  return 0;
}
