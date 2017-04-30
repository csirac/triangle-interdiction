#include <igraph.h>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;

int main( int argc, char** argv ) {
   unsigned n = stoi( argv[1] );
   unsigned m = stoi( argv[2] );

   igraph_t g;
   igraph_erdos_renyi_game( &g, IGRAPH_ERDOS_RENYI_GNM, n, m, false, false );

   string name = "er_" + to_string( n ) + "_" + to_string( m ) + ".el";
   
   FILE *fout = fopen(name.c_str(), "w");
      
   igraph_write_graph_edgelist( &g, fout );

   igraph_destroy( &g );
   return 0;
}
