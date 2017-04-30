#ifndef GLPK_SOLVER_CPP
#define GLPK_SOLVER_CPP

#include <glpk.h>
#include "mygraph.cpp"
using namespace std;
using namespace mygraph;

class GLPK_solver {
public:
   glp_prob* lp;
   int* ia;
   int* ja;
   double* ar;
   unsigned k;
   unsigned n;
   glp_smcp parmLP; //parameters for LP solver
   glp_iocp parmMIP; //parameters for MIP solver
   
   GLPK_solver( Graph& G, double max_hours = 4.0 ) {
      
      G.logg( DEBUG, "Initializing GLPK..." );
      lp = glp_create_prob();
      glp_set_obj_dir(lp, GLP_MIN);
      k = G.countTvalid();
      
      glp_add_rows( lp, k ); //For each triangle, we have a constraint
      for (size_t i = 1; i <= k; ++i) {
	glp_set_row_bnds( lp, i, GLP_LO, 1.0, 0.0 ); 
      }

      n = G.E.size();
      glp_add_cols( lp, n );
      for (size_t i = 1; i <= n; ++i) {
	glp_set_col_kind( lp, i, GLP_BV ); //whether we choose edge
	glp_set_obj_coef( lp, i, 1.0 ); //unweighted for now
      }
      G.logg(DEBUG, "Allocating matrices...");

      ia = new int[ k * 3 + 1];
      ja = new int[ k * 3 + 1];
      ar = new double[ k * 3 + 1];
      size_t index = 1; //index of current element in the constraint matrix
      size_t tid = 0;
      for (auto t = G.T.begin();
	   t != G.T.end();
	   ++t ) {
	 if (!G.triangle_valid(t))
	    continue;
	 ia[ index ] = tid + 1;
	 ja[ index ] = t->e1->eid + 1;
	 ar[ index ] = 1.0;
	 ++index;
	 ia[ index ] = tid + 1;
	 ja[ index ] = t->e2->eid + 1;
	 ar[ index ] = 1.0;
	 ++index;
	 ia[ index ] = tid + 1;
	 ja[ index ] = t->e3->eid + 1;
	 ar[ index ] = 1.0;
	 ++index;

	 if (t->e1 == t->e2 || t->e2 == t->e3 || t->e3 == t->e1) {
	    G.logg(ERROR, "Triangle " + to_string(tid) + " does not have distinct edges!");

	 }
	 
	 ++tid;
      }

      G.logg(DEBUG, "Creating GLPK matrix...");
      //all path constraints added, create the matrix
      glp_load_matrix( lp, (index - 1), ia, ja, ar );

      //Set LP paramaters
      glp_init_smcp( &parmLP );
      parmLP.msg_lev = GLP_MSG_OFF;
      //      parmLP.meth = GLP_PRIMAL;
      //      parmLP.pricing = GLP_PT_PSE;
      //      parmLP.r_test = GLP_RT_HAR;
      //      parmLP.tol_bnd = 1e-7;
      //      parmLP.tol_dj = 1e-7;
      //      parmLP.tol_priv = 1e-9;
      //      parmLP.obj_ll = -1.0 * DBL_MAX;
      //      parmLP.obj_ul = DBL_MAX;
      //      parmLP.it_lim = INT_MAX;
      parmLP.tm_lim = static_cast< size_t >( 1000 * 3600 * max_hours ); //convert hours to milliseconds
      //      parmLP.out_frq = 500;
      //      parmLP.out_dly = 0;
      parmLP.presolve = GLP_ON;

      //set MIP parameters
      //start with defaults
      glp_init_iocp( &parmMIP );
      parmMIP.presolve = GLP_ON;
      parmMIP.binarize = GLP_ON;
      parmMIP.msg_lev = GLP_MSG_OFF;
      parmMIP.tm_lim = static_cast< size_t >( 1000 * 3600 * max_hours ); //convert hours to milliseconds
      
      G.logg(DEBUG, "GLPK instance constructed." );

   }

   /*
    * Results in v_out, must already have size n
    */
   
   bool LP_solve( vector< double >& v_out ) {
      //      G.logg(DEBUG, "Starting GLPK simplex..." );
      //      v_out.clear();
      //      v_out.reserve( n );
      glp_simplex( lp, &parmLP );
      for (size_t i = 0; i < n; ++i) {
	 v_out[i] = glp_get_col_prim( lp, i + 1 );
	 //	 v_out.push_back( glp_get_col_prim( lp, i + 1 ) );
      }

      return ( glp_get_status( lp ) == GLP_OPT );
   }

   bool MIP_solve(Graph& G) {
      G.logg(INFO, "Starting GLPK simplex..." );
      glp_simplex( lp, &parmLP );
      G.logg(INFO, "Starting GLPK integer optimization...");
      glp_intopt( lp, &parmMIP );

      if (glp_mip_status( lp ) != GLP_OPT ) {
	 G.logg(ERROR, "Optimal solution not found...");
	 return false;
      }

      size_t eid = 0;
      for (pedge e = G.E.begin();
	   e != G.E.end();
	   ++e ) {
	 if (glp_mip_col_val( lp, eid + 1 ) == 1) {
	    e->in_S = true;
	 }
	 ++eid;
      }

      return true;

   }

   ~GLPK_solver() {
      delete [] ia;
      delete [] ja;
      delete [] ar;
   }
};

bool glpk_kortsarz( Graph& G, double max_hours = 4.0 ) {
   G.logg(INFO, "Starting Kortsarz et al....");
   vector< double > lpsol( G.E.size(), 0.0 );

   bool b_solve;
   unsigned eid;
   unsigned sizeS = 0;
   double elapsed_time = 0.0;
   double max_secs = max_hours * 3600;
   
   clock_t t_start = clock();
   
   do {
      if (G.countTvalid() == 0)
	 break;
      elapsed_time = double (clock() - t_start) / CLOCKS_PER_SEC;
      if (elapsed_time > max_secs) {
	 return false;
      }
      
      GLPK_solver KSolver( G, max_hours );
      bool status = KSolver.LP_solve( lpsol );
      if (status == false)
	 return false;
      
      b_solve = false;
      eid = 0;
      G.logg(DEBUG, "Size of E: " + to_string( G.E.size() ));
      for (auto e = G.E.begin();
	   e != G.E.end();
	   ++e) {
	 if (lpsol[eid] > 0.5) {
	    b_solve = true;
	    e->in_S = true;
	    //	    G.remove_edge( e );
	    //	    cerr << eid << " found, breaking...\n";
	    ++sizeS;
	    break;
	 }	    
	 ++eid;
      }
      
   } while (b_solve);

   //Next step: remove all Triangle-redundant edges
   for (auto e = G.E.begin();
	e != G.E.end();
	++e) {
      if ( !e->in_S ) {
	 if ( G.is_triangle_redundant( e ) ) {
	    e->in_W = true;
	 }
      }
   }

   for (auto i = G.V.begin();
	i != G.V.end();
	++i) {
      i->in_A = false;
      i->in_C = false;
   }
   
   //Run Bipartite-complement
   vector< pedge > R;
   G.bipartite_complement( R );
   G.logg(DEBUG, "Size of R:" + to_string(R.size()));

   // for (auto i = G.V.begin();
   // 	i != G.V.end();
   // 	++i) {
   //    i->in_A = false;
   //    i->in_C = false;
   // }
   
   // G.logg(DEBUG, "Starting bipartite...");
   // vector< pedge > B;
   // G.bipartite( B );
   // G.logg(DEBUG, "Size of G':" + to_string(B.size() + R.size()));

   for (auto i = R.begin();
	i != R.end();
	++i) {
      (*i)->in_S = true;
   }
   
   //   G.logg(INFO, "Solution size:" + to_string(R.size() + sizeS));

   return true; //R.size() + sizeS;

}

bool glpk_tarl( Graph& G, double max_hours = 4.0 ) {
   vector< double > lpsol( G.E.size(), 0.0 );

   GLPK_solver KSolver( G, max_hours );
   bool status = KSolver.LP_solve( lpsol );
   if (status == false)
      return false;
   
   unsigned eid = 0;
   vector< pedge > W;
   unsigned nWprime = 0;
   unsigned nedges_added = 0;
   double weight = 0.0;
   double sWeight = 0.0;
   for (pedge e = G.E.begin();
	e != G.E.end();
	++e ) {
      if ( lpsol[eid] < 1.0 / 5 ) {
	 e-> in_W = true;
	 W.push_back( e );
	 //	 lpsol[eid] = 0.0;
	    
	 if ( !(e-> in_S) ) {
	    ++nWprime;
	    //ADD other edges of triangles containing e
	    
	    for ( auto it = e->Delta.begin(); it != e->Delta.end(); ++it ) {
	       if (G.triangle_valid(*it)) {
		  nedges_added = 0;
		  weight = 0.0;
		  if ((*it)->e1 != e) {
		     (*it)->e1->in_S = true;
		     //		  lpsol[ (*it)->e1->eid ] = 1.0;
		     ++nedges_added;
		     weight +=  lpsol[ (*it)->e1->eid ];

		  }
		  if ((*it)->e2 != e) {
		     (*it)->e2->in_S = true;
		     //		  lpsol[ (*it)->e2->eid ] = 1.0;
		     ++nedges_added;
		     weight +=  lpsol[ (*it)->e2->eid ];

		  }
		  if ((*it)->e3 != e) {
		     (*it)->e3->in_S = true;
		     //		  lpsol[ (*it)->e3->eid ] = 1.0;
		     ++nedges_added;
		     weight +=  lpsol[ (*it)->e3->eid ];

		  }
		  
		  if (nedges_added > 2)
		     G.logg(DEBUG, "More than 2 edges added from triangle...");
		  if (weight < 4.0 / 5)
		     G.logg(DEBUG, "Weight added is too low...");

		  sWeight += weight;
	       }
	    }
	 }
      }
     
      ++eid;
   }
   G.logg(DEBUG, "sWeight: " + to_string(sWeight));
   G.logg(DEBUG, "size of W: " + to_string(W.size()));
   G.logg(DEBUG, "size of Wprime: " + to_string(nWprime));
   //   for (auto i = W.begin();
   //   	i != W.end();
   //   	++i) {
   //         (*i)->in_S = false;
   //   }

   // G.logg(DEBUG, "Removing S...");

   unsigned sizeS = 0;
   double weightS = 0.0;
   auto i = G.E.begin();
   while ( i != G.E.end() ) {
      if ( i->in_S ) {
	 //	 i = G.remove_edge( i );
	 ++sizeS;
	 weightS += lpsol[ i->eid ];
      }
      ++i;
   }
   //   else {
   //     ++i;
   //   }
   // }
   G.logg(DEBUG, "Size of S: " + to_string( sizeS ) );
   G.logg(DEBUG, "Weight of S: " + to_string( weightS ) );
   // G.logg(DEBUG, "Removing W...");
   // for (auto i = W.begin(); i != W.end(); ++i) {
   //   G.remove_edge( *i );
   // }

   G.logg(DEBUG, "Starting bipartite (complement)...");
   vector< pedge > R;
   G.bipartite_complement( R );
   G.logg(DEBUG, "Size of R:" + to_string(R.size()));
   //   G.logg(DEBUG, "Edges in G': " + to_string(G.E.size()));

   for (auto i = R.begin();
	i != R.end();
	++i) {
      (*i)->in_S = true;
   }
   
   // for (auto i = G.V.begin();
   // 	i != G.V.end();
   // 	++i) {
   //    i->in_A = false;
   //    i->in_C = false;
   // }

   // G.logg(DEBUG, "Starting bipartite...");
   // vector< pedge > B;
   // G.bipartite( B );
   // G.logg(DEBUG, "Size of G':" + to_string(B.size() + R.size()));

   G.logg(INFO, "Solution size:" + to_string(R.size() + sizeS));

   return true; //R.size() + sizeS;
}


#endif   
