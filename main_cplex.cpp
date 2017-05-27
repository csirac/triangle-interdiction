#include <iostream>
#include "mygraph.cpp"
#include "lp_solver.cpp"
#include "glpk_solver.cpp"
#include <string>
#include <iomanip>
#include <ctime>

void print_help() {
  cout << "Options: " << endl;
  cout << "-G <graph filename in edge list format>" << endl
       << "-D [run DART]" << endl
       << "-T [run TARL]" << endl
       << "-K [run 2-approx. of Kortsarz et al.]" << endl
       << "-O [run optimal solution via IBM CPLEX]" << endl
       << "-t [time limit in hours, default is 4]" << endl
       << "-x [Max number of threads, default is 18]" << endl;
   
}

using namespace std;
using namespace mygraph;

int main(int argc, char ** argv) {
  int c;
  extern char *optarg;
  extern int optint, optopt;

  if (argc == 1) {
     print_help();
     return 1;
  }

  string fname;
  bool bKortsarz = false;
  bool bOpt = false;
  bool bDart = false;
  bool bTarl = false;
  bool bGLPK = false;
  string s_arg;
  unsigned nThreads = 18;
  double max_hours = 4.0;
  
  while ((c = getopt( argc, argv, ":G:LOKTDt:x:") ) != -1) {
    switch(c) {
    case 'G':
      //graph specification
      fname.assign( optarg );
      break;
    case 'O':
       bOpt = true;
       break;
    case 'L':
       bGLPK = true;
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

  string outfile = "run-" + fname.substr(0, fname.find_last_of( '.' )) + "-" + to_string( time(0) ) + ".txt";
  ofstream of( outfile.c_str() );
  
  // of << "# time Triangle_listing, (size,time of) ";
  // if (bDart)
  //    of << "Dart ";
  // if (bTarl)
  //    of << "Tarl ";
  // if (bKortsarz)
  //    of << "Kortsarz ";
  // if (bOpt)
  //    of << "Opt ";
  //  ofstream of( "log.txt" );

  //  of << endl;
  Graph G(INFO, cout);

  G.logg(INFO, "Reading graph...");
  G.read_edge_list( fname );

  //  if (bOpt || bTarl || bKortsarz || bDart ) {
  of << "# Triangle-listing:" << endl;
  clock_t t_start = clock();
  G.list_triangles();
  
  double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
  //  G.check_triangles();
    
  of << G.T.size() << ' ' << t_elapsed << endl;
     //  }
     //  else {
     //     cout << "No algorithms specified.\n";
     //     of.close();
     //     return 0;
     //  }

  if (bDart) {
     of << "# Dart_base" << endl;
     clock_t t_start = clock();

     unsigned size = G.dart_base();
     double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
     of << size << ' ' << t_elapsed << ' ' << G.ensure_feasibility() << endl;
     G.clear_edges();
     of << "# Dart_base with pruning" << endl;
     t_start = clock();

     size = G.dart_base_heu();
     t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;

     of << size << ' ' << t_elapsed << ' ' << G.ensure_feasibility() << endl;
     //     G.dart_base();

     if (G.ensure_feasibility())
	G.logg(INFO, "Pruned solution remains feasible.");
     else
	G.logg(WARN, "Pruned solution is infeasible!");
     
     G.clear_edges();

  }

  if (bTarl) {
     of << "# TARL" << endl;
     clock_t t_start = clock();
     unsigned size = tarl( G, nThreads, max_hours );
     double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
     of << size << ' ' << t_elapsed << ' ' << G.ensure_feasibility() << endl;
     if (G.ensure_feasibility())
	G.logg(INFO, "TARL solution is feasible.");
     else
	G.logg(WARN, "TARL solution is infeasible!");
     
     G.clear_edges();
  }

  if (bKortsarz) {
     of << "# Kortsarz" << endl;
     G.info_check( of );
     clock_t t_start = clock();
     kortsarz( G, nThreads, max_hours );
     double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
     unsigned size = G.countS();
     of << size << ' ' << t_elapsed << ' ' << G.ensure_feasibility() << endl;

     if (G.ensure_feasibility())
	G.logg(INFO, "Kortsarz solution is feasible.");
     else
	G.logg(WARN, "Kortsarz solution is infeasible!");
	  
     G.clear_edges();

  }

  if (bOpt) {
     of << "# Opt (CPLEX)" << endl;
     G.info_check( of );
     
     vector< double > v_sol;
     bool status;
     clock_t t_start = clock();
     status = solve_ip( G, v_sol, max_hours, nThreads );
     double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
     unsigned size;
     if (status)
	size = G.countS();
     else
	size = 0;
     of << size << ' ' << t_elapsed << ' ' << G.ensure_feasibility() << endl;
     cout << size << ' ' << t_elapsed << ' ' << G.ensure_feasibility() << endl;
     if (G.ensure_feasibility())
	G.logg(INFO, "CPLEX solution is feasible.");
     else
	G.logg(WARN, "CPLEX solution is infeasible!");

     G.clear_edges();
  }

  if (bGLPK) {
     of << "# Opt (GLPK)" << endl;
     clock_t t_start = clock();
     GLPK_solver GLPK( G );
     vector< unsigned > vout;
     bool status;
     status = GLPK.MIP_solve( G );
     double t_elapsed = double (clock() - t_start) / CLOCKS_PER_SEC;
     unsigned size;
     if (status)
	size = G.countS();
     else
	size = 0;
     
     of << size << ' ' << t_elapsed << ' ' << G.ensure_feasibility() << endl;

     G.clear_edges();
  }

  of.close();
  return 0;
}
