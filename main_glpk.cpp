#include <iostream>
#include "mygraph.cpp"
#include "glpk_solver.cpp"
#include <string>
#include <iomanip>
#include <ctime>
#include <unistd.h>
#include <random>
#include "proc_size.cpp"

enum Algo {DART1, DART2, OPT};

void randomRemoveEdges( Graph& G_in,
			unsigned mRemove,
			vector< Algo >& A,
			ostream& os, unsigned checkpoint);

void randomRemoveEdges( tinyGraph& g,
			unsigned mRemove,
			double& graphUpdateTime,
			double& dart2RemTime );

void randomRemoveEdges( Graph& G,
			unsigned mRemove,
			double& graphUpdateTime,
			double& dart1RemTime );

void randomAddEdges( Graph& G,
		     unsigned mAdd, vector< Algo >& A, ostream& os, unsigned checkpoint);
void randomAddEdges( Graph& G, unsigned mAdd, double& graphUpdateTime, double& dart1AddTime);
void randomAddEdges( tinyGraph& G, unsigned mAdd, double& , double& );
double randomAddAndRemoveEdges( Graph& G, unsigned mAdd);
double randomAddAndRemoveEdges( tinyGraph& G, unsigned mAdd);


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
       << "-P [run primal-dual algorithm]" << endl
       << "-p [Preprocess graph and count triangles only]" << endl
       << "-A <m_add> (run DART, then adaptively add <m_add> random edges to the network)" << endl
       << "-R <m_remove> (run DART, then adaptively remove <m_remove> random edges from the network)" << endl
       << "-S <m_addremove> (run DART, then adaptively add <m_addremove> random edges to the network, and then remove the same edges)" << endl
       << "-t [time limit in hours, default is 4]" << endl
       << "-C [Perform in-depth comparison of dynamic algorithms]" << endl;
   
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
  bool bPD = false;
  bool bBinaryFormat = false;
  
  string s_arg;
  unsigned nThreads = 18;
  double max_hours = 4.0;
  unsigned mAdd;
  bool bAdd = false;
  string outfilename;
  bool bOut = false;
  bool bRemove = false;
  unsigned mRemove;
  bool bAddRemove = false;
  unsigned Nreps = 1;
  unsigned erN = 0;
  double erP;
  bool bDynCompare = false;
  unsigned checkpoint = 1;
  bool bVisualize = false;
  string VisFname = "";
  bool bPreprocessOnly = false;
  
  while ((c = getopt( argc, argv, ":G:OKTDEt:x:A:R:S:g:o:N:PC:V:p") ) != -1) {
    switch(c) {
    case 'p':
      bPreprocessOnly = true;
      break;
    case 'V':
      bVisualize = true;
      VisFname.assign( optarg );
      break;
    case 'C':
       bDynCompare = true;
       s_arg.assign( optarg );
       checkpoint = stoi( s_arg );
       break;
    case 'P':
       bPD = true;
       break;
    case 'N':
       s_arg.assign( optarg );
       Nreps = stoi( s_arg );
       break;
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
    case 'R':
       bRemove = true;
       s_arg.assign( optarg );
       mRemove = stoi( s_arg );
       break;
    case 'S':
       bAddRemove = true;
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
  } else {
    if (fname.substr(0, 2) == "ER") {
      size_t pos_colon = fname.find_first_of( ':' );
      string s_n = fname.substr(2, (pos_colon - 2));
      string s_p = fname.substr( pos_colon + 1 );
      erP = stod( s_p );
      erN = stoi( s_n );
    }
  }

  //  string outfile = "run-" + fname.substr(0, fname.find_last_of( '.' )) + "-" + to_string( time(0) ) + ".txt";
  //  string outfile = "log" + to_string( time(0) ) + ".txt";
  //  ofstream of( outfile.c_str() );
  Graph G(INFO, cout);

  resultsHandler myResults;

  for (unsigned iter = 0; iter < Nreps; ++iter) {
     G.clear_graph();
     myResults.data.clear();
     
     if (bKortsarz || bOpt || bTarl || bDart || bPD || bPreprocessOnly ) {
	G.logg(INFO, "Reading graph...");
	if (bBinaryFormat) {
	   G.read_edge_list_bin( fname );
	} else {
	  if (erN > 0) {
	    //generate ER graph
	    G.genER( erN, erP );
	    myResults.add( "erP", erP );
	  } else {
	    G.read_edge_list( fname );
	  }
	}
	
	G.logg(INFO, "Basic graph info (n, m): " + to_string( G.V.size() ) + " "  + to_string( G.E.size() ) );
	myResults.add( "GraphNodes", G.V.size() );
	myResults.add( "GraphEdges", G.E.size() );
	myResults.add( "GraphPreprocess", G.preprocessTime );
	myResults.add( "GraphName", fname );

	if (bPreprocessOnly) {
	  size_t nTriangles = G.count_triangles();
	  myResults.add( "NumberTriangles", nTriangles );
	  iter = Nreps;
	  bKortsarz = bOpt = bTarl = bDart = bPD = bDart2 = false;
	}
     }

     if (bDynCompare) {
	if (bAdd) {
	   //Perform comprehensive dynamic addition of edges test
	   vector< Algo > A;
	   if (bOpt)
	      A.push_back( OPT );
	   if (bDart)
	      A.push_back( DART1 );
	   if (bDart2)
	      A.push_back( DART2 );
	   
	   for (unsigned iter = 0; iter < Nreps; ++iter) {
	      ofstream ofile;
	      ofile.open( outfilename.c_str(), ios::app );

	      randomAddEdges( G, 
			      mAdd, A, ofile, checkpoint );
     
	      ofile.close();
	   }
	}

	if (bRemove) {
	   //Perform comprehensive dynamic removal of edges test
	   vector< Algo > A;
	   if (bOpt)
	      A.push_back( OPT );
	   if (bDart)
	      A.push_back( DART1 );
	   if (bDart2)
	      A.push_back( DART2 );
	   
	   for (unsigned iter = 0; iter < Nreps; ++iter) {
	      ofstream ofile;
	      ofile.open( outfilename.c_str(), ios::app );

	      randomRemoveEdges( G, 
				 mRemove, A, ofile, checkpoint );
     
	      ofile.close();
	   }
	}
	
	exit( 0 );
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
	myResults.add( "TriangleSize", G.T.size() );
	myResults.add( "TriangleTime", t_triangle );
     }

     if (bPD) {
       G.clear_edges();
	G.primal_dual();
	myResults.add( "PrimalDualSize", G.sizeS );
	myResults.add( "PrimalDualTime", G.runningTime );
	myResults.add( "PrimalDualMem" , getPeakRSS() / (1024.0 * 1024.0));
     }
     
     if (bDart) {
       	G.clear_edges();

	G.logg( INFO, "Starting dart-base (first implementation)..." );
	G.dart_base();
	G.countS();
	G.logg( INFO, "dart-base: " + to_string(G.sizeS) + " " + to_string(G.runningTime)
		+ " " + to_string( getPeakRSS() / (1024.0 * 1024.0)));
	
	myResults.add( "Dart1Size", G.sizeS );
	myResults.add( "Dart1Time", G.runningTime );
	myResults.add( "Dart1Mem" , getPeakRSS() / (1024.0 * 1024.0));
     
	if (bAdd) {
	   G.init_dynamic(true);
	   G.logg(INFO, "Adding edges...");
	   double GraphUpdateTime;
	   double dart1AddTime;
	   randomAddEdges( G, mAdd, GraphUpdateTime, dart1AddTime );
	   G.countS();
	   G.logg( INFO, "After adding " +to_string(mAdd) + " edges, " + to_string(G.sizeS) + " " + to_string(G.runningTime) );

	   G.logg(INFO, "Basic graph info (n, m): " + to_string( G.V.size() ) + " "  + to_string( G.E.size() ) );

	   myResults.add( "GraphUpdateTime", GraphUpdateTime );
	   myResults.add( "Dart1AddTime", dart1AddTime );
	   myResults.add( "Dart1AddSize", G.sizeS );
	   myResults.add( "NumberAdded", mAdd );
	   
	//    G.logg(DEBUG, "Comprehensive feasibility check...");
	//    G.T.clear();
	//    G.init_static();
	//    if( !G.verify_graph() ) {
	//       G.logg(ERROR, "Graph structure is incorrect." );
	//       exit(1);
	//    } else {
	//       G.logg(INFO, "Graph structure is correct." );
	//    }

	//    G.list_triangles();
	//    G.logg(INFO, "Triangle-listing: " + to_string( G.T.size() ) );
	//    if (G.ensure_feasibility()) {
	//       G.logg(DEBUG, "Dart_add has maintained feasibility...");
	//    }  else {
	//       G.logg(ERROR, "Dart_add has violated feasibility.");
	//    }
	// 
	}

	if (bRemove) {
	   G.init_dynamic(true);
	   G.logg(INFO, "Removing edges...");
	   double GraphUpdateTime;
	   double dart1RemTime;
	   randomRemoveEdges( G, mRemove, GraphUpdateTime, dart1RemTime );
	   G.countS();
	   G.logg( INFO, "After removing " +to_string(mRemove) + " edges, " + to_string(G.sizeS) + " " + to_string( dart1RemTime) );

	   G.logg(INFO, "Basic graph info (n, m): " + to_string( G.V.size() ) + " "  + to_string( G.E.size() ) );

	   myResults.add( "GraphUpdateTime", GraphUpdateTime );
	   myResults.add( "Dart1RemoveTime", dart1RemTime );
	   myResults.add( "Dart1RemoveSize", G.sizeS );
	   myResults.add( "NumberRemoved", mRemove );
	   
	//    G.logg(DEBUG, "Comprehensive feasibility check...");
	//    G.T.clear();
	//    G.init_static();
	//    if( !G.verify_graph() ) {
	//       G.logg(ERROR, "Graph structure is incorrect." );
	//       exit(1);
	//    } else {
	//       G.logg(INFO, "Graph structure is correct." );
	//    }

	//    G.list_triangles();
	//    G.logg(INFO, "Triangle-listing: " + to_string( G.T.size() ) );
	//    if (G.ensure_feasibility()) {
	//       G.logg(DEBUG, "Dart_add has maintained feasibility...");
	//    }  else {
	//       G.logg(ERROR, "Dart_add has violated feasibility.");
	//    }
	// 
	}
	
     }
  

     if (bDart2) {
	G.logg( INFO, "Reading graph into tinyGraph structure..." );
	G.clear_edges();
	tinyGraph g;
	if (bBinaryFormat) {
	   g.read_edge_list_bin( fname );
	} else {
	   //	G.read_edge_list( fname );
	}
	G.logg( INFO, "Starting dart-base (second implementation)..." );
	clock_t t_start = clock();
	g.dart_base();
	double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
	size_t size = g.countS();
	G.logg( INFO, "dart-base2: " + to_string(size) + " " +
		to_string(t_elapsed) + " " + to_string( getPeakRSS() / (1024.0 * 1024.0)));

	myResults.add( "Dart2Size", size );
	myResults.add( "Dart2Time", t_elapsed );
	myResults.add( "Dart2Mem" , getPeakRSS() / (1024.0 * 1024.0));
	myResults.add( "tinyGraphPreprocess", g.preprocessTime );
	myResults.add( "tinyGraphNodes", g.n );
	myResults.add( "tinyGraphEdges", g.m );

	if (bAdd) {
	   G.logg(INFO, "Adding edges...");
	   double tinyGraphUpdateTime;
	   double dart2AddTime;
	   randomAddEdges( g, mAdd, tinyGraphUpdateTime, dart2AddTime );
	   size = g.countS();
	   G.logg( INFO, "After adding " + to_string(mAdd) + " edges, " + to_string(size) + " " + to_string(dart2AddTime) );

	   G.logg(INFO, "Basic graph info (n, m): " + to_string( g.n ) + " "  + to_string( g.m ) );

	   myResults.add( "tinyGraphUpdateTime", tinyGraphUpdateTime );
	   myResults.add( "Dart2AddTime", dart2AddTime );
	   myResults.add( "Dart2AddSize", size );
	   myResults.add( "NumberAdded", mAdd );
	}

	if (bRemove) {
	   G.logg(INFO, "Removing edges...");
	   double tinyGraphUpdateTime;
	   double dart2RemoveTime;
	   randomRemoveEdges( g, mRemove, tinyGraphUpdateTime, dart2RemoveTime );
	   size = g.countS();
	   G.logg( INFO, "After removing " + to_string(mRemove) + " edges, " + to_string(size) + " " + to_string(dart2RemoveTime) );

	   G.logg(INFO, "Basic graph info (n, m): " + to_string( g.n ) + " "  + to_string( g.m ) );

	   myResults.add( "tinyGraphUpdateTime", tinyGraphUpdateTime );
	   myResults.add( "Dart2RemoveTime", dart2RemoveTime );
	   myResults.add( "Dart2RemoveSize", size );
	   myResults.add( "NumberRemoved", mRemove );
	}

	if (bAddRemove) {
	   G.logg(INFO, "Adding and removing edges...");
	   t_elapsed += randomAddAndRemoveEdges( g, mAdd );
	   size = g.countS();
	   G.logg( INFO, "After adding and removing " +to_string(mAdd) + " edges, " + to_string(size) + " " + to_string(t_elapsed) );

	   G.logg(INFO, "Basic graph info (n, m): " + to_string( g.n ) + " "  + to_string( g.m ) );
	}

     }

     if (bTarl) {
	G.logg( INFO, "Starting TARL..." );
	G.clear_edges();

	clock_t t_start = clock();
	if (glpk_tarl2( G, max_hours )) {
	   double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
	   unsigned size = G.countS();
	   G.logg( INFO, "TARL: " + to_string(size) + " " + to_string(t_elapsed) + " " + to_string(G.ensure_feasibility() ));
	   if (G.ensure_feasibility())
	      G.logg(INFO, "TARL solution is feasible.");
	   else
	      G.logg(WARN, "TARL solution is infeasible!");

	   myResults.add( "TarlSize", size );
	   myResults.add( "TarlTime", t_elapsed + t_triangle);

	} else {
	   G.logg(INFO, "TARL (GLPK) exceeded time limit!");
	   bTarl = false;
	   myResults.set( "TarlSize", 0 );
	   myResults.set( "TarlTime", 0 );
	   return 1;
	}

     }

     if (bKortsarz) {
	G.logg( INFO, "Starting Kortsarz..." );
	G.clear_edges();
	clock_t t_start = clock();
	if (glpk_kortsarz( G, max_hours ) ) {
	   double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
	   unsigned size = G.countS();
	   G.logg( INFO, "Kortsarz (GLPK): " + to_string(size) + " " + to_string(t_elapsed) + " " + to_string(G.ensure_feasibility() ));

	   myResults.add( "KortsarzSize", size );
	   myResults.add( "KortsarzTime", t_elapsed + t_triangle);
	
	   if (G.ensure_feasibility())
	      G.logg(INFO, "Kortsarz solution is feasible.");
	   else
	      G.logg(WARN, "Kortsarz solution is infeasible!");
	} else {
	   G.logg(INFO, "Kortsarz (GLPK) exceeded time limit!");
	   bKortsarz = false;
	   myResults.set( "KortsarzSize", 0 );
	   myResults.set( "KortsarzTime", 0 );
	   return 1;
	}
     
     }

     if (bOpt) {
	G.logg(INFO, "Starting GLPK IP solver..." );
	G.clear_edges();
	clock_t t_start = clock();
	GLPK_solver GLPK( G, max_hours );
	vector< unsigned > vout;
	if (GLPK.MIP_solve( G )) {
	   double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
	   unsigned size = G.countS();

	   G.logg(INFO, "OPT (GLPK): " + to_string(size) + " " + to_string(t_elapsed) + " " + to_string(G.ensure_feasibility() ));

	   myResults.add( "OptSize", size );
	   myResults.add( "OptTime", t_elapsed + t_triangle);
	   
	} else {
	   G.logg(INFO, "OPT (GLPK) exceeded time limit!");
	   bOpt = false;

	   myResults.set( "OptSize", 0 );
	   myResults.set( "OptTime", 0 );
	   return 1;
	}


     }

     if (bOut) {
       ofstream ofile;
       ofile.open( outfilename.c_str(), ios::app );
       myResults.print_xml( ofile );
       ofile.close();
     }

     if (bVisualize) {
       ofstream ofile;
       ofile.open( VisFname.c_str() );
       G.write_graphml( ofile );
     }
     //     myResults.print_xml( cout );
  }


  
  
  return 0;
}

void outputResult( Graph& G, string algName, double time, ostream& os, bool xml = true ) {
   resultsHandler mr;
   mr.add( algName + "Size", G.countS() );
   mr.add( algName + "Time", time );
   mr.add( "GraphNodes", G.V.size() );
   mr.add( "GraphEdges", G.E.size() );
   mr.add( "GraphPreprocess", G.preprocessTime );
   if (xml)
      mr.print_xml( os );
   else
      mr.print( os );

}

void outputResult( Graph& G, string algName, unsigned solSize, double time, ostream& os ) {
   resultsHandler mr;
   mr.add( algName + "Size", solSize );
   mr.add( algName + "Time", time );
   mr.add( "GraphNodes", G.V.size() );
   mr.add( "GraphEdges", G.E.size() );
   mr.add( "GraphPreprocess", G.preprocessTime );
   mr.print_xml( os );
}

void outputResult( tinyGraph& G, string algName, double time, ostream& os, bool xml = true ) {
   resultsHandler mr;
   mr.add( algName + "Size", G.countS() );
   mr.add( algName + "Time", time );
   mr.add( "GraphNodes", G.n );
   mr.add( "GraphEdges", G.m );
   mr.add( "GraphPreprocess", G.preprocessTime );
   if (xml)
      mr.print_xml( os );
   else
      mr.print( os );
}

void outputResult( tinyGraph& G, string algName, unsigned solSize, double time, ostream& os ) {
   resultsHandler mr;
   mr.add( algName + "Size", solSize );
   mr.add( algName + "Time", time );
   mr.add( "GraphNodes", G.n );
   mr.add( "GraphEdges", G.m );
   mr.add( "GraphPreprocess", G.preprocessTime );
   mr.print_xml( os );
}

void runStaticAlgs( vector< Algo >& A,
		    tinyGraph& h,
		    Graph& H,
		    ostream& os,
		    Logger& logg) {
   clock_t t_start;
   double t_elapsed;
   for (unsigned i = 0; i < A.size(); ++i) {
      switch( A[i] ) {
      case DART1:
	 logg( INFO, "Starting Dart1..." );
	 t_start = clock();
	 H.dart_base();
	 t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
	 outputResult( H, "Dart1", t_elapsed, os );
	 logg( INFO, "Dart1 finished." );
	 outputResult( H, "Dart1", t_elapsed, logg.of, false );	 
	 break;

      case DART2:
	 logg( INFO, "Starting Dart2..." );
	 t_start = clock();
	 h.dart_base();
	 t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
	 outputResult( h, "Dart2", t_elapsed, os );
	 logg( INFO, "Dart2 finished." );
	 outputResult( h, "Dart2", t_elapsed, logg.of, false );
	 break;

      case OPT:
	 logg( INFO, "Starting OPT(GLPK)..." );
	 t_start = clock();
	 H.clear_edges();
	 H.list_triangles();
	 GLPK_solver GLPK( H, 1 );
	 GLPK.MIP_solve( H );
	 t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
	 outputResult( H, "Opt", t_elapsed, os );
	 logg( INFO, "OPT(GLPK) finished." );	       
	 outputResult( H, "Opt", t_elapsed, logg.of, false );
	 
	 H.clear_edges();

	 break;
      }
   }	       
}
		   

void randomAddEdges( Graph& G_in,
		     unsigned mAdd, vector< Algo >& A, ostream& os, unsigned checkpoint = 1) {
   G_in.logg( INFO, "Starting in-depth comparison, addition of edges..." );
   Graph G( G_in );
   //DartAdd, Dart2Add need an initial static run.
   G_in.logg( INFO, "Starting initial static solutions..." );
   clock_t t_start = clock();
   G.dart_base();
   double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
   unsigned solSize1 = G.countS();
   outputResult( G, "Dart1Add", solSize1, t_elapsed, os );

   G.init_dynamic();
   
   tinyGraph g( G );
   t_start = clock();
   g.dart_base();
   t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
   unsigned solSize2 = g.countS();
   outputResult( g, "Dart2Add", solSize2, t_elapsed, os );
   //Begin edge addition procedure
   uniform_int_distribution<> vdist(0, g.n - 1);
   unsigned eAdded = 0;
   while (eAdded < mAdd) {
      node_id from = vdist( gen );
      node_id to = vdist( gen );
      vector< tinyEdge >::iterator e;
      t_start = clock();
      if ( g.add_edge_half( from, to, e ) ) {
	 g.add_edge_half(to, from, e );
	 g.m += 1;
	 g.preprocessTime = double (clock() - t_start) / CLOCKS_PER_SEC;
	 t_start = clock();
	 pedge f;
	 G.dynamic_add_edge( from, to, f );
         G.preprocessTime = double (clock() - t_start) / CLOCKS_PER_SEC;
	 //Run dynamic algorithms
	 t_start = clock();
	 solSize1 += G.dart_add_edge( from, to, f );
	 t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
	 outputResult( G, "Dart1Add", solSize1, t_elapsed, os );
	 t_start = clock();
	 solSize2 += g.dart_add_edge( to, e );
	 t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
	 outputResult( g, "Dart2Add", solSize2, t_elapsed, os );
	 if ((eAdded + 1) % checkpoint == 0) {
	    G_in.logg( INFO, "Reached checkpoint..." );
	    //Create static graphs
	    Graph H( G );
	    tinyGraph h( g );
	 
	    runStaticAlgs( A, h, H, os, G_in.logg );
	 }
	 
	 ++eAdded;
      }
   }
}

void randomRemoveEdges( Graph& G_in,
			unsigned mRemove,
			vector< Algo >& A,
			ostream& os, unsigned checkpoint = 1) {
   G_in.logg( INFO, "Starting in-depth comparison, removal of edges..." );
   Graph G( G_in );
   //Need an initial static run.
   G_in.logg( INFO, "Starting initial static solutions..." );
   clock_t t_start = clock();
   G.dart_base();
   double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
   int solSize1 = G.countS();
   outputResult( G, "Dart1Static", solSize1, t_elapsed, os );

   G.init_dynamic();
   
   tinyGraph g( G );
   t_start = clock();
   g.dart_base();
   t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
   int solSize2 = g.countS();
   outputResult( g, "Dart2Static", solSize2, t_elapsed, os );

   //Begin edge removal procedure
   uniform_int_distribution<> vdist(0, g.n - 1);
   unsigned numRemoved = 0;
   double dartTime2, dartTime1;
   double updateTime2, updateTime1;
   while (numRemoved < mRemove) {
      node_id from;
      do {
	 from = vdist( gen );
      } while (g.adjList[from].neis.size() == 0);

      uniform_int_distribution<> neiDist(0, g.adjList[from].neis.size() - 1);
      size_t toIndex = neiDist( gen );
      vector< tinyEdge >::iterator e = g.adjList[ from ].neis.begin() + toIndex;

      //remove edge
      dartTime2 = 0.0;
      updateTime2 = 0.0;
      solSize2 += g.dart_remove_edge( from, e, dartTime2, updateTime2 );

      g.preprocessTime = updateTime2;
      outputResult( g, "Dart2Remove", solSize2, dartTime2, os );

      Node& From = G.V[ from ];
      list< pedge >::iterator f = From.neighbors.begin();
      size_t tmpST = 0;
      while (tmpST < toIndex) {
	 ++tmpST;
	 ++f;
      }

      dartTime1 = 0.0;
      updateTime1 = 0.0;
      solSize1 += G.dart_remove_edge( *f, dartTime1, updateTime1 );
      G.preprocessTime = updateTime1;
      outputResult( G, "Dart1Remove", solSize1, dartTime1, os );
      if ((numRemoved + 1) % checkpoint == 0) {
	 G_in.logg( INFO, "Reached checkpoint..." );
	 //Create static graphs
	 Graph H( G );
	 tinyGraph h( g );
	 
	 runStaticAlgs( A, h, H, os, G_in.logg );
      }
	 
      ++numRemoved;
   }
}

void randomRemoveEdges( tinyGraph& g,
			unsigned mRemove,
			double& graphUpdateTime,
			double& dart2RemTime ) {

   //Begin edge removal procedure
   uniform_int_distribution<> vdist(0, g.n - 1);
   unsigned numRemoved = 0;
   graphUpdateTime = 0.0;
   dart2RemTime = 0.0;
   while (numRemoved < mRemove) {
      node_id from;
      do {
	 from = vdist( gen );
      } while (g.adjList[from].neis.size() == 0);

      uniform_int_distribution<> neiDist(0, g.adjList[from].neis.size() - 1);
      size_t toIndex = neiDist( gen );
      vector< tinyEdge >::iterator e = g.adjList[ from ].neis.begin() + toIndex;

      //remove edge
      g.dart_remove_edge( from, e, dart2RemTime, graphUpdateTime );
	 
      ++numRemoved;
   }
}

void randomRemoveEdges( Graph& G,
			unsigned mRemove,
			double& graphUpdateTime,
			double& dart1RemTime ) {

   //Begin edge removal procedure
   uniform_int_distribution<> vdist(0, G.V.size() - 1);
   unsigned numRemoved = 0;
   graphUpdateTime = 0.0;
   dart1RemTime = 0.0;
   while (numRemoved < mRemove) {
      node_id from;
      do {
	 from = vdist( gen );
      } while (G.V[from].neighbors.size() == 0);

      uniform_int_distribution<> neiDist(0, G.V[from].neighbors.size() - 1);
      size_t toIndex = neiDist( gen );

      Node& From = G.V[ from ];
      auto f = From.neighbors.begin();
      size_t tmpST = 0;
      while (tmpST < toIndex) {
	 ++tmpST;
	 ++f;
      }

      //remove edge
      G.dart_remove_edge( *f, dart1RemTime, graphUpdateTime );
	 
      ++numRemoved;
   }
}

/*
 * Add random edges to tinygraph g and update the Dart2 solution
 * Returns the time to update the tinyGraph structure and the time
 * taken by Dart2Add
 */
void randomAddEdges( tinyGraph& g, unsigned mAdd, double& graphUpdateTime, double& dart2AddTime) {
   uniform_int_distribution<> vdist(0, g.n - 1);

   unsigned eAdded = 0;
   clock_t t_start;
   graphUpdateTime = 0.0;
   dart2AddTime = 0.0;
   bool bAdded;
   while (eAdded < mAdd) {
      node_id from = vdist( gen );
      node_id to = vdist( gen );
      vector< tinyEdge >::iterator e;
      t_start = clock();
      bAdded = g.add_edge_half( from, to, e );
      if ( bAdded ) {
	 g.add_edge_half(to, from, e );
	 g.m += 1;
	 graphUpdateTime += double (clock() - t_start) / CLOCKS_PER_SEC;
	 t_start = clock();
	 g.dart_add_edge( to, e );
	 dart2AddTime += double (clock() - t_start) / CLOCKS_PER_SEC;
	 ++eAdded;
      }
   }
}

double randomAddAndRemoveEdges( tinyGraph& G, unsigned mAdd) {
   random_device rd;
   mt19937 gen( rd() );
   uniform_int_distribution<> vdist(0, G.n - 1);
   double t_elapsed = 0.0;
   
   unsigned eAdded = 0;
   vector< node_id > edgesAdded;
   while (eAdded < mAdd) {
      node_id from = vdist( gen );
      node_id to = vdist( gen );
      vector< tinyEdge >::iterator e;
      if ( G.add_edge_half( from, to, e ) ) {
	G.add_edge_half(to, from, e );
	edgesAdded.push_back( to );
	edgesAdded.push_back( from );
	G.m += 1;
	t_elapsed += G.dart_add_edge( to, e );
	++eAdded;
      }
   }

   for (unsigned i = 0; i < edgesAdded.size(); i += 2) {
     vector< tinyEdge >::iterator ee = G.findEdgeInList( edgesAdded[ i ], edgesAdded[ i + 1 ] );
     //     t_elapsed += G.dart_remove_edge( edgesAdded[i], ee );
   }
   
   return t_elapsed;
}

/*
 * Add random edges to G and update the Dart solution
 * Returns the total time elapsed of dart1, and total time to update Graph
 */

void randomAddEdges( Graph& G, unsigned mAdd, double& graphUpdateTime, double& dart1AddTime) {
   random_device rd;
   mt19937 gen( rd() );
   uniform_int_distribution<> vdist(0, G.V.size() - 1);
   clock_t t_start;
   unsigned eAdded = 0;
   graphUpdateTime = 0.0;
   dart1AddTime = 0.0;
   while (eAdded < mAdd) {
      node_id from = vdist( gen );
      node_id to = vdist( gen );
      pedge e;

      t_start = clock();
      if ( G.dynamic_add_edge( from, to, e ) ) {
	 graphUpdateTime += double (clock() - t_start) / CLOCKS_PER_SEC;
	 t_start = clock();
	 G.dart_add_edge( from, to, e );
	 dart1AddTime += double (clock() - t_start) / CLOCKS_PER_SEC;
	 ++eAdded;
      }
   }
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
      //      t_elapsed += G.dart_remove_edge( *it, true );
   }
   
   return t_elapsed;
}

		     
