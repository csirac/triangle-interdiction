#include <string>
#include "mygraph.cpp"
#include <iostream>

using namespace mygraph;

int main(int argc, char ** argv) {
     tinyGraph G;
     string fname( argv[1] );
     G.read_adj_list_bin( fname );
     //     G.logg(INFO, "Listing triangles..." );
     //     G.list_triangles();
     G.logg(INFO, "Starting DART..." );
     G.dart_base_free();

     G.logg(INFO, "Starting pruning..." );
     G.free_prune();
     G.logg(INFO, "Free prune: " + to_string( G.countS() ) );
     //     G.dart_base_free_smaller();

     // G.logg(INFO, "Checking feasibility...");
     // Graph H;
     // H.read_edge_list_bin( fname );
     // H.logg(INFO, "Listing triangles..." );
     // H.list_triangles();
     // H.logg(INFO, "Initializing S..." );
     // for (auto it1 = H.E.begin();
     // 	  it1 != H.E.end();
     // 	  ++it1 ) {
     // 	if ( (*it1).from < (*it1).to ) {
     // 	   auto it2 = G.adjList[ (*it1).from ].begin();
     // 	   while ((*it2).getId() != (*it1).to)
     // 	      ++it2;
     // 	   if ((*it2).inS())
     // 	      (*it1).in_S = true;
     // 	} else {
     // 	   auto it2 = G.adjList[ (*it1).to ].begin();
     // 	   while ((*it2).getId() != (*it1).from)
     // 	      ++it2;
     // 	   if ((*it2).inS())
     // 	      (*it1).in_S = true;
     // 	}
     // }

     // // H.clear_edges();
     // // H.dart_base();
     // H.logg(INFO, "Checking feasibility..." );
     // cout << H.ensure_feasibility() << endl;

     // G.logg(INFO, "Checking validity..." );
     // cout << G.check_validity() << endl;
     return 0;
}
   
