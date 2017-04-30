#include <string>
#include "mygraph.cpp"

using namespace mygraph;

int main(int argc, char ** argv) {
   smallGraph G;
   string fname( argv[1] );
   G.read_edge_list_bin( fname );
   G.logg(INFO, "Listing triangles..." );
   G.list_triangles();
   G.logg(INFO, "Starting DART..." );
   G.dart_base_free();
   
}
   
