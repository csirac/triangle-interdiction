#include <iostream>
#include "mygraph.cpp"
#include "glpk_solver.cpp"
#include <string>
#include <iomanip>
#include <ctime>
#include <unistd.h>
#include <random>
#include "proc_size.cpp"

double randomAddEdges( Graph& G, unsigned mAdd);
double randomAddAndRemoveEdges( Graph& G, unsigned mAdd);

void outputFile( ostream& os,
		 string& fname,
		 string algName,
		 unsigned n,
		 unsigned m,
		 double preprocessTime,
		 unsigned solSize,
		 double runningTime ) {
  os << fname;
  os << ' ';
  os << algName << ' ';
  os << n << ' ';
  os << m << ' ';
  os << preprocessTime << ' ';
  os << solSize << ' ';
  os << runningTime << ' ';
  os << getPeakRSS() / (1024.0 * 1024.0) << endl; //peak resources used in Mb
}
		 

void print_help() {
  cout << "Options: " << endl;
  cout << "-G <graph filename in edge list format>" << endl
       << "-g <graph filename in binary edge list format>" << endl
       << "-D [run DART]" << endl
       << "-E [run DART (second implementation)]" << endl
       << "-T [run TARL]" << endl
       << "-K [run 2-approx. of Kortsarz et al.]" << endl
       << "-O [run optimal solution via GNU GLPK]" << endl
       << "-A <m_add> (run DART, then adaptively add <m_add> random edges to the network)" << endl
       << "-t [time limit in hours, default is 4]" << endl
       << "-x [Max number of threads, default is 18]" << endl;
  
   
}

using namespace std;
using namespace mygraph;

int main(int argc, char ** argv) {
  int c;
  extern char *optarg;
  //  extern int optint, optopt;

  if (argc == 1) {
     print_help();
     return 1;
  }

  string fname;
  bool bKortsarz = false;
  bool bOpt = false;
  bool bDart = false;
  bool bDart2 = false;
  bool bTarl = false;
  bool bBinaryFormat = false;
  
  string s_arg;
  unsigned nThreads = 18;
  double max_hours = 4.0;
  unsigned mAdd;
  bool bAdd = false;
  string outfilename;
  bool bOut = false;
  
  while ((c = getopt( argc, argv, ":G:OKTDEt:x:A:g:o:") ) != -1) {
    switch(c) {
    case 'o':
       s_arg.assign( optarg );
       outfilename = s_arg;
       bOut = true;
       break;
    case 'A':
       bAdd = true;
       s_arg.assign( optarg );
       mAdd = stoi( s_arg );
       break;
    case 'G':
      //graph specification
      fname.assign( optarg );
      break;
    case 'g':
      //graph specification
      fname.assign( optarg );
      bBinaryFormat = true;
      break;
    case 'O':
       bOpt = true;
       break;
    case 'K':
       bKortsarz = true;
       break;
    case 'T':
       bTarl = true;
       break;
    case 'x':
       s_arg.assign( optarg );
       nThreads = stoi( s_arg );
       break;
    case 't':
       s_arg.assign( optarg );
       max_hours = stod( s_arg );
       break;
    case 'D':
       bDart = true;
       break;
    case 'E':
       bDart2 = true;
       break;
    case '?':
      print_help();
      return 1;
      break;
    }
  }
  

  if (fname.size() == 0) {
     print_help();
     cerr << "Input graph is required.\n";
     return 1;
  }

  //  string outfile = "run-" + fname.substr(0, fname.find_last_of( '.' )) + "-" + to_string( time(0) ) + ".txt";
  string outfile = "log" + to_string( time(0) ) + ".txt";
  ofstream of( outfile.c_str() );
  Graph G(DEBUG, of, true);

  ofstream ofile;
  
  if (bOut) {
     ofile.open( outfilename.c_str(), ios::app );
  }

  if (bKortsarz || bOpt || bTarl || bDart) {
     G.logg(INFO, "Reading graph...");
     if (bBinaryFormat) {
	G.read_edge_list_bin( fname );
     } else {
	G.read_edge_list( fname );
     }
     G.logg(INFO, "Basic graph info (n, m): " + to_string( G.V.size() ) + " "  + to_string( G.E.size() ) );
  }

  double t_triangle = 0.0;
  
  if (bKortsarz || bTarl || bOpt ) {
     G.logg(INFO, "Starting triangle-listing (single threaded)..." );
     clock_t t_start = clock();
     G.list_triangles();
     double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
     G.logg(INFO, "Triangle-listing: " + to_string( G.T.size() ) + " "  + to_string( t_elapsed ) );
     //     G.logg(INFO, "Starting triangle-listing (multi-threaded)..." );
     //     clock_t t_start = clock();
     //     list_triangles_multi( G, nThreads );
     //     double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
    
     //     G.logg(INFO, "Triangle-listing: " + to_string( G.T.size() ) + " "  + to_string( t_elapsed ) );
     t_triangle = t_elapsed;
  }
  
  if (bDart) {
     G.logg( INFO, "Starting dart-base (first implementation)..." );
     G.dart_base_free();
     G.countS();
     G.logg( INFO, "dart-base: " + to_string(G.sizeS) + " " + to_string(G.runningTime)
	     + " " + to_string( getPeakRSS() / (1024.0 * 1024.0)));

     if (bOut) {
       outputFile( ofile,
		   fname,
		   "dart-base1",
		   G.V.size(),
		   G.E.size(),
		   G.preprocessTime,
		   G.sizeS,
		   G.runningTime );
     }
     
     G.logg(INFO, "Starting free_prune()...");

     G.free_prune();
     G.countS();
     G.logg( INFO, "After pruning: " + to_string(G.sizeS) + " " + to_string(G.runningTime)
	     + " " + to_string( getPeakRSS() / (1024.0 * 1024.0)));

     if (bOut) {
       outputFile( ofile,
		   fname,
		   "prune1",
		   G.V.size(),
		   G.E.size(),
		   G.preprocessTime,
		   G.sizeS,
		   G.runningTime );
     }
     
     if (bAdd) {
	G.init_dynamic();
	if( !G.verify_graph() ) {
	   G.logg(ERROR, "Graph structure is incorrect." );
	   exit(1);
	} else {
	   G.logg(INFO, "Graph structure is correct." );
	}
	
	G.logg(INFO, "Adding and removing edges...");
	G.runningTime += randomAddAndRemoveEdges( G, mAdd );
	G.countS();
	G.logg( INFO, "After adding and removing " +to_string(mAdd) + " edges, " + to_string(G.sizeS) + " " + to_string(G.runningTime) );

	G.logg(INFO, "Basic graph info (n, m): " + to_string( G.V.size() ) + " "  + to_string( G.E.size() ) );

	G.logg(DEBUG, "Comprehensive feasibility check...");
	G.T.clear();
	G.init_static();
	if( !G.verify_graph() ) {
	   G.logg(ERROR, "Graph structure is incorrect." );
	   exit(1);
	} else {
	   G.logg(INFO, "Graph structure is correct." );
	}

	G.list_triangles();
	G.logg(INFO, "Triangle-listing: " + to_string( G.T.size() ) );
	if (G.ensure_feasibility()) {
	   G.logg(DEBUG, "Dart_add has maintained feasibility...");
	}  else {
	   G.logg(ERROR, "Dart_add has violated feasibility.");
	}
     }


     
     // G.logg( INFO, "Starting dart-base-free..." );
     // t_start = clock();

     // size = G.dart_base_free();
     // t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
     
     // G.logg( INFO, "dart-base-free: " + to_string(size) + " " + to_string(t_elapsed) +
     // 	     " " + to_string(G.ensure_feasibility() ));
     // //     G.dart_base();
     // G.logg(INFO, "Starting better_prune()...");
     // t_start = clock();
     // G.better_prune_S();
     // t_elapsed = t_elapsed + double (clock() - t_start) / CLOCKS_PER_SEC;
     // size = G.countS();
     // G.logg( INFO, "After pruning: " + to_string(size) + " " + to_string(t_elapsed) + " " + to_string(G.count_infeasible() ));

     // auto e1 = G.E.begin();
     // auto e2 = H.E.begin();

     // while (e1 != G.E.end()) {
     // 	if (e1->in_S != e2->in_S) {
     // 	   cerr << "Soultions disagree on edge " << e1-> to << ' ' << e1->from << endl;
     // 	   cerr << "H: " << e2->to << ' ' << e2->from << endl;
     // 	   cerr << e1->in_S << ' ' << e2->in_S << endl;
     // 	}
     // 	++e1; ++e2;
     // }
     
     G.clear_edges();

     // G.logg( INFO, "Starting dart-base-integrated..." );
     // t_start = clock();

     // G.dart_base_integrated();
     // t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
     // size = G.countS();
     // G.logg( INFO, "dart-base-integrated: " + to_string(size) + " " + to_string(t_elapsed) + " " + to_string(G.ensure_feasibility() ));

     // G.logg( INFO, "Starting dart-base-free..." );
     // t_start = clock();

     // size = G.dart_base_free();
     // t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;

     // G.logg( INFO, "dart-base-free: " + to_string(size) + " " + to_string(t_elapsed) +
     // 	     " " + to_string(G.ensure_feasibility() ));
     // //     G.dart_base();
     // G.logg(INFO, "Starting free_prune_old()...");
     // t_start = clock();
     // G.free_prune_old();
     // t_elapsed = t_elapsed + double (clock() - t_start) / CLOCKS_PER_SEC;
     // size = G.countS();
     // G.logg( INFO, "After pruning: " + to_string(size) + " " + to_string(t_elapsed) + " " + to_string(G.ensure_feasibility() ));
     
     G.clear_edges();
  }
  

  if (bDart2) {
     G.logg( INFO, "Reading graph into tinyGraph structure..." );
     tinyGraph g;
     if (bBinaryFormat) {
	g.read_edge_list_bin( fname );
     } else {
	//	G.read_edge_list( fname );
     }
     G.logg( INFO, "Starting dart-base (second implementation)..." );
     clock_t t_start = clock();

     size_t size = g.dart_base_free();
     double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;

     G.logg( INFO, "dart-base2: " + to_string(size) + " " +
	     to_string(t_elapsed) + " " + to_string( getPeakRSS() / (1024.0 * 1024.0)));

     if (bOut) {
       outputFile( ofile,
		   fname,
		   "dart2",
		   g.n ,
		   g.m,
		   g.preprocessTime,
		   size,
		   t_elapsed );
     }
     
     G.logg(INFO, "Starting free_prune()...");
     t_start = clock();
     g.free_prune();
     t_elapsed = t_elapsed + double (clock() - t_start) / CLOCKS_PER_SEC;
     size = g.countS();
     G.logg( INFO, "After pruning: " + to_string(size) + " " + to_string(t_elapsed)
	     + " " + to_string( getPeakRSS() / (1024.0 * 1024.0)));

     if (bOut) {
       outputFile( ofile,
		   fname,
		   "prune2",
		   g.n ,
		   g.m,
		   g.preprocessTime,
		   size,
		   t_elapsed );
     }
     
     if (bAdd) {
	G.init_dynamic();
	if( !G.verify_graph() ) {
	   G.logg(ERROR, "Graph structure is incorrect." );
	   exit(1);
	} else {
	   G.logg(INFO, "Graph structure is correct." );
	}
	
	G.logg(INFO, "Adding and removing edges...");
	t_elapsed += randomAddAndRemoveEdges( G, mAdd );
	size = G.countS();
	G.logg( INFO, "After adding and removing " +to_string(mAdd) + " edges, " + to_string(size) + " " + to_string(t_elapsed) );

	G.logg(INFO, "Basic graph info (n, m): " + to_string( G.V.size() ) + " "  + to_string( G.E.size() ) );

	G.logg(DEBUG, "Comprehensive feasibility check...");
	G.T.clear();
	G.init_static();
	if( !G.verify_graph() ) {
	   G.logg(ERROR, "Graph structure is incorrect." );
	   exit(1);
	} else {
	   G.logg(INFO, "Graph structure is correct." );
	}

	G.list_triangles();
	G.logg(INFO, "Triangle-listing: " + to_string( G.T.size() ) );
	if (G.ensure_feasibility()) {
	   G.logg(DEBUG, "Dart_add has maintained feasibility...");
	}  else {
	   G.logg(ERROR, "Dart_add has violated feasibility.");
	}
     }


     
     // G.logg( INFO, "Starting dart-base-free..." );
     // t_start = clock();

     // size = G.dart_base_free();
     // t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
     
     // G.logg( INFO, "dart-base-free: " + to_string(size) + " " + to_string(t_elapsed) +
     // 	     " " + to_string(G.ensure_feasibility() ));
     // //     G.dart_base();
     // G.logg(INFO, "Starting better_prune()...");
     // t_start = clock();
     // G.better_prune_S();
     // t_elapsed = t_elapsed + double (clock() - t_start) / CLOCKS_PER_SEC;
     // size = G.countS();
     // G.logg( INFO, "After pruning: " + to_string(size) + " " + to_string(t_elapsed) + " " + to_string(G.count_infeasible() ));

     // auto e1 = G.E.begin();
     // auto e2 = H.E.begin();

     // while (e1 != G.E.end()) {
     // 	if (e1->in_S != e2->in_S) {
     // 	   cerr << "Soultions disagree on edge " << e1-> to << ' ' << e1->from << endl;
     // 	   cerr << "H: " << e2->to << ' ' << e2->from << endl;
     // 	   cerr << e1->in_S << ' ' << e2->in_S << endl;
     // 	}
     // 	++e1; ++e2;
     // }
     
     G.clear_edges();

     // G.logg( INFO, "Starting dart-base-integrated..." );
     // t_start = clock();

     // G.dart_base_integrated();
     // t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
     // size = G.countS();
     // G.logg( INFO, "dart-base-integrated: " + to_string(size) + " " + to_string(t_elapsed) + " " + to_string(G.ensure_feasibility() ));

     // G.logg( INFO, "Starting dart-base-free..." );
     // t_start = clock();

     // size = G.dart_base_free();
     // t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;

     // G.logg( INFO, "dart-base-free: " + to_string(size) + " " + to_string(t_elapsed) +
     // 	     " " + to_string(G.ensure_feasibility() ));
     // //     G.dart_base();
     // G.logg(INFO, "Starting free_prune_old()...");
     // t_start = clock();
     // G.free_prune_old();
     // t_elapsed = t_elapsed + double (clock() - t_start) / CLOCKS_PER_SEC;
     // size = G.countS();
     // G.logg( INFO, "After pruning: " + to_string(size) + " " + to_string(t_elapsed) + " " + to_string(G.ensure_feasibility() ));
     
     G.clear_edges();
     
  }

  if (bTarl) {
     G.logg( INFO, "Starting TARL..." );
     clock_t t_start = clock();
     if (glpk_tarl( G, max_hours )) {
	double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
	unsigned size = G.countS();
	G.logg( INFO, "TARL: " + to_string(size) + " " + to_string(t_elapsed) + " " + to_string(G.ensure_feasibility() ));
	if (G.ensure_feasibility())
	   G.logg(INFO, "TARL solution is feasible.");
	else
	   G.logg(WARN, "TARL solution is infeasible!");

	if (bOut) {
	  outputFile( ofile,
		      fname,
		      "tarl",
		      G.V.size(),
		      G.E.size(),
		      G.preprocessTime,
		      size,
		      t_elapsed + t_triangle );

	}
     } else {
	G.logg(INFO, "TARL (GLPK) exceeded time limit!");
	if (bOut) {
	  outputFile( ofile,
		      fname,
		      "tarl",
		      G.V.size(),
		      G.E.size(),
		      G.preprocessTime,
		      0,
		      0.0 );
	}
     }
     G.clear_edges();
  }

  if (bKortsarz) {
     G.logg( INFO, "Starting Kortsarz..." );
     clock_t t_start = clock();
     if (glpk_kortsarz( G, max_hours ) ) {
	double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
	unsigned size = G.countS();
	G.logg( INFO, "Kortsarz (GLPK): " + to_string(size) + " " + to_string(t_elapsed) + " " + to_string(G.ensure_feasibility() ));

	if (bOut) {
	  outputFile( ofile,
		      fname,
		      "kortsarz",
		      G.V.size(),
		      G.E.size(),
		      G.preprocessTime,
		      size,
		      t_elapsed + t_triangle );

	}
	
	if (G.ensure_feasibility())
	   G.logg(INFO, "Kortsarz solution is feasible.");
	else
	   G.logg(WARN, "Kortsarz solution is infeasible!");
     } else {
	G.logg(INFO, "Kortsarz (GLPK) exceeded time limit!");
	if (bOut) {
	  outputFile( ofile,
		      fname,
		      "kortsarz",
		      G.V.size(),
		      G.E.size(),
		      G.preprocessTime,
		      0,
		      0.0);

	}
     }
     
     G.clear_edges();
  }

  if (bOpt) {
     G.logg(INFO, "Starting GLPK IP solver..." );
     clock_t t_start = clock();
     GLPK_solver GLPK( G, max_hours );
     vector< unsigned > vout;
     if (GLPK.MIP_solve( G )) {
	double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
	unsigned size = G.countS();

	G.logg(INFO, "OPT (GLPK): " + to_string(size) + " " + to_string(t_elapsed) + " " + to_string(G.ensure_feasibility() ));
	if (bOut) {
	  outputFile( ofile,
		      fname,
		      "opt",
		      G.V.size(),
		      G.E.size(),
		      G.preprocessTime,
		      size,
		      t_elapsed + t_triangle );

	}
	
     } else {
	G.logg(INFO, "OPT (GLPK) exceeded time limit!");
	if (bOut) {
	  outputFile( ofile,
		      fname,
		      "opt",
		      G.V.size(),
		      G.E.size(),
		      G.preprocessTime,
		      0,
		      0 );
	}
     }

     G.clear_edges();
  }

  of.close();
  return 0;
}

/*
 * Add random edges to G and update the Dart solution
 * Returns the total time elapsed
 */

double randomAddEdges( Graph& G, unsigned mAdd) {
   random_device rd;
   mt19937 gen( rd() );
   uniform_int_distribution<> vdist(0, G.V.size() - 1);
   double t_elapsed = 0.0;
   
   unsigned eAdded = 0;
   while (eAdded < mAdd) {
      node_id from = vdist( gen );
      node_id to = vdist( gen );
      pedge e;
      if ( G.dynamic_add_edge( from, to, e ) ) {
	 t_elapsed += G.dart_add_edge( from, to, e );
	 ++eAdded;
      }
   }
   return t_elapsed;
}

/*
 * Add random edges to G and update the Dart solution
 * Returns the total time elapsed
 */

double randomAddAndRemoveEdges( Graph& G, unsigned mAdd) {
   random_device rd;
   mt19937 gen( rd() );
   uniform_int_distribution<> vdist(0, G.V.size() - 1);
   double t_elapsed = 0.0;
   
   unsigned eAdded = 0;
   vector< pedge > edges_added;
   while (eAdded < mAdd) {
      node_id from = vdist( gen );
      node_id to = vdist( gen );
      pedge e;
      if ( G.dynamic_add_edge( from, to, e ) ) {
	 edges_added.push_back( e );
	 t_elapsed += G.dart_add_edge( from, to, e );
	 ++eAdded;
      }
   }

   size_t size = G.countS();
   G.logg( INFO, "After adding " +to_string(mAdd) + " edges, " + to_string(size) + " " + to_string(t_elapsed) );
   
   for (auto it = edges_added.begin();
	it != edges_added.end();
	++it ) {
      t_elapsed += G.dart_remove_edge( *it, true );
   }
   
   return t_elapsed;
}
