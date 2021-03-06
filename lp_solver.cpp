#ifndef LP_SOLVER_CPP
#define LP_SOLVER_CPP

#include <ilcplex/ilocplex.h>
#include <unordered_set>
#include "mygraph.cpp"

using namespace mygraph;

ILOSTLBEGIN

bool triangle_valid( pangle& t ) {
   return !( t->e1->in_S || t->e2->in_S || t->e3->in_S );
}


/*
 * Solves LP 1.
 * Output value for each edge e \in E is in v_sol
 * G must have a valid triangle-listing
 * Requires eid of each edge to be valid
 */
bool solve_lp ( Graph& G,
		vector< double >& v_sol,
		IloEnv& env,
		IloCplex& cplex,
		//		IloModel& model,
		unsigned nThreads = 18,
		double max_hours = 4.0) {
   G.logg(INFO, "Beginning LP...");
   
   //Triangles should already be generated
   //   IloEnv env;
   IloModel model(env);
   IloNumVarArray var(env);
   IloRangeArray con(env);

   G.logg(DEBUG, "Adding variables to LP...");
   //For each edge, we need a variable
   // 0 \le x_i \le 1
   for (pedge e = G.E.begin();
	e != G.E.end();
	++e ) {
      //    var.add( IloNumVar(env, 0.0, 1.0) );
      var.add( IloNumVar(env) ); //shouldn't need constraints, >= 0 by default
   }

   // logg(DEBUG, "Building objective of LP...");
   // //Add the objective function, which is to
   // //min. the number of edges
   // IloExpr edge_obj( env, 0 );
   // edge_obj.setName( "Edge_objective" );
   // for (unsigned i = 0; i < G.E.size(); ++i) {
   //   edge_obj += var[i];
   // }

   // model.add( IloMinimize( env, edge_obj ) );

   // logg(DEBUG, "Adding triangle constraints to LP...");
   // //Depends on valid edge ids
   // IloRangeArray c(env);
   // for ( pangle t = G.T.begin();
   // 	t != G.T.end();
   // 	++t ) {
   //   c.add( var[ t->e1->eid ] + var[ t->e2->eid ] + var[ t->e3->eid ] >= 1.0 ); //each triangle must be broken
   // }

   // model.add(c);
   G.logg(DEBUG, "Adding constraints...");
  
   IloRangeArray c(env);
   IloObjective obj = IloMinimize( env );
   for (pangle t = G.T.begin();
	t != G.T.end();
	++t ) {
      if (triangle_valid( t ))
	 c.add(IloRange( env, 1.0, IloInfinity ));
   }

   G.logg(DEBUG, "Setting objective coefficients...");
   unsigned eid = 0;
   for (pedge e = G.E.begin();
	e != G.E.end();
	++e ) {
      obj.setLinearCoef( var[ eid ], 1.0 );
      ++eid;
   }

   G.logg(DEBUG, "Setting constraint coefficients...");
   unsigned tid = 0;
   for (pangle t = G.T.begin();
	t != G.T.end();
	++t ) {
      if (triangle_valid( t )) {
	 c[tid].setLinearCoef(var[ t->e1->eid ], 1.0);
	 c[tid].setLinearCoef(var[ t->e2->eid ], 1.0);
	 c[tid].setLinearCoef(var[ t->e3->eid ], 1.0);
	 ++tid;
      }
   }

   model.add(obj);
   model.add(c);

   G.logg(DEBUG, "Building cplex from model...");
   //   IloCplex cplex(model);
   cplex.extract(model);
   cplex.setParam(IloCplex::IntParam::Threads, nThreads);
   cplex.setParam(IloCplex::ClockType, 2); //Tells cplex to use wall-clock time
   double hour = 60.0 * 60;
   cplex.setParam(IloCplex::TiLim, max_hours * hour );
   cplex.setOut(env.getNullStream());
   G.logg(DEBUG, "Solving LP...");
   cplex.solve();

   eid = 0;

   v_sol.clear();
   v_sol.reserve( G.E.size() );
   for (pedge e = G.E.begin();
        e != G.E.end();
        ++e ) {
      v_sol.push_back( cplex.getValue( var[eid] ) );
      ++eid;
   }

   G.logg(INFO, "Solution finished, objective value: " + to_string(cplex.getObjValue()));
   if (cplex.getCplexStatus() == IloCplex::Optimal)
      return true;
   else
      return false;

   //   cplex.clearModel();
   cplex.end();
   obj.end();
   c.end();
   var.end();
   con.end();
   model.end();
   env.end();
}

/*
 * Solves IP 1 optimally.
 * Output value for each edge e \in E is in v_sol
 * G must have a valid triangle-listing
 * Requires eid of each edge to be valid
 */
unsigned solve_ip ( Graph& G,
		    vector< double >& v_sol,
		    double max_hours,
		    unsigned nThreads = 18) {
   G.logg(INFO, "Beginning IP with maximum of " + to_string( nThreads ) + " threads...");
   
   //Triangles should already be generated
   IloEnv env;
   IloModel model(env);
   IloNumVarArray var(env);//, G.E.size(), 0, 1, ILOBOOL );
   IloRangeArray con(env);

   G.logg(DEBUG, "Adding variables to IP...");
   //   For each edge, we need a variable
   //   x_i \in {0,1}
   for (pedge e = G.E.begin();
   	e != G.E.end();
   	++e ) {
      //         var.add( IloNumVar(env, 0.0, 1.0) );
      var.add( IloNumVar(env,0.0,1.0, ILOINT )); 
   }

   G.logg(DEBUG, "Building objective of LP...");
   //Add the objective function, which is to
   //min. the number of edges
   IloExpr edge_obj( env, 0 );
   edge_obj.setName( "Edge_objective" );
   for (unsigned i = 0; i < G.E.size(); ++i) {
      edge_obj += var[i];
   }

   model.add( IloMinimize( env, edge_obj ) );

   G.logg(DEBUG, "Adding triangle constraints to LP...");
   //Depends on valid edge ids
   IloRangeArray c(env);
   for ( pangle t = G.T.begin();
	 t != G.T.end();
	 ++t ) {
      c.add( var[ t->e1->eid ] + var[ t->e2->eid ] + var[ t->e3->eid ] >= 1.0 ); //each triangle must be broken
   }
   model.add(c);
   // G.logg(DEBUG, "Adding constraints...");
  
   // IloRangeArray c(env);
   // IloObjective obj = IloMinimize( env );
   // for (pangle t = G.T.begin();
   // 	t != G.T.end();
   // 	++t ) {
   //    if (triangle_valid( t ))
   // 	 c.add(IloRange( env, 1.0, IloInfinity ));
   // }

   // G.logg(DEBUG, "Setting objective coefficients...");
   // unsigned eid = 0;
   // for (pedge e = G.E.begin();
   // 	e != G.E.end();
   // 	++e ) {
   //    obj.setLinearCoef( var[ eid ], 1.0 );
   //    ++eid;
   // }

   // G.logg(DEBUG, "Setting constraint coefficients...");
   // unsigned tid = 0;
   // for (pangle t = G.T.begin();
   // 	t != G.T.end();
   // 	++t ) {
   //    if (triangle_valid( t )) {
   // 	 c[tid].setLinearCoef(var[ t->e1->eid ], 1.0);
   // 	 c[tid].setLinearCoef(var[ t->e2->eid ], 1.0);
   // 	 c[tid].setLinearCoef(var[ t->e3->eid ], 1.0);
   // 	 ++tid;
   //    }
   // }

   //   model.add(obj);
   //   model.add(c);

   G.logg(DEBUG, "Building cplex from model...");
   IloCplex cplex(model);
   cplex.setParam(IloCplex::IntParam::Threads, nThreads);
   cplex.setParam(IloCplex::IntParam::ParallelMode, -1 ); //Opportunistic
   //   cplex.setParam(IloCplex::IntParam::Threads, 1);
   //Tells cplex to use wall-clock time
   cplex.setParam(IloCplex::ClockType, 2); 
   double hour = 60.0 * 60;
   cplex.setParam(IloCplex::TiLim, max_hours * hour );
   cplex.setOut(env.getNullStream());
   //   cplex.setParam(IloCplex::EpAGap, 1.0 / (2 * G.E.size()));
   G.logg(DEBUG, "Solving IP...");
   cplex.solve();

   unsigned   eid = 0;

   v_sol.clear();
   v_sol.reserve( G.E.size() );
   for (pedge e = G.E.begin();
        e != G.E.end();
        ++e ) {
      v_sol.push_back( cplex.getValue( var[eid] ) );
      if ( cplex.getValue(var[eid]) == 1 )
	 e->in_S = true;
      ++eid;

   }

   G.logg(INFO, "Solution finished, objective value: " + to_string(cplex.getObjValue()));
   if (cplex.getCplexStatus() == IloCplex::Optimal)
      return static_cast<unsigned> ( cplex.getObjValue() );
   else
      return 0;
}

/*
 * Solves IP 1 optimally.
 * Output value for each edge e \in E is in v_sol
 * G must have a valid triangle-listing
 * Requires eid of each edge to be valid
 */
unsigned solve_ip2 ( Graph& G,
		vector< double >& v_sol,
		unsigned nThreads = 18) {
   G.logg(INFO, "Beginning IP with maximum of " + to_string( nThreads ) + " threads...");
   
   //Triangles should already be generated
   IloEnv env;
   IloModel model(env);
   IloNumVarArray var(env); // G.E.size(), 0, 1, ILOBOOL );
   IloRangeArray con(env);

   G.logg(DEBUG, "Adding variables to IP...");
   //   For each edge, we need a variable
   //   x_i \in {0,1}
     for (pedge e = G.E.begin();
   	e != G.E.end();
   	++e ) {
	//         var.add( IloNumVar(env, 0.0, 1.0) );
	var.add( IloNumVar(env,0.0,1.0, ILOINT )); //shouldn't need constraints, >= 0 by default
     }

     G.logg(DEBUG, "Building objective of LP...");
   //Add the objective function, which is to
   //min. the number of edges
   IloExpr edge_obj( env, 0 );
   edge_obj.setName( "Edge_objective" );
   for (unsigned i = 0; i < G.E.size(); ++i) {
     edge_obj += var[i];
   }

   model.add( IloMinimize( env, edge_obj ) );

   G.logg(DEBUG, "Adding triangle constraints to LP...");
   //Depends on valid edge ids
   IloRangeArray c(env);
   for ( pangle t = G.T.begin();
   	t != G.T.end();
   	++t ) {
     c.add( var[ t->e1->eid ] + var[ t->e2->eid ] + var[ t->e3->eid ] >= 1.0 ); //each triangle must be broken
   }

   model.add(c);
   // G.logg(DEBUG, "Adding constraints...");
  
   // IloRangeArray c(env);
   // IloObjective obj = IloMinimize( env );
   // for (pangle t = G.T.begin();
   // 	t != G.T.end();
   // 	++t ) {
   //    if (triangle_valid( t ))
   // 	 c.add(IloRange( env, 1.0, IloInfinity ));
   // }

   // G.logg(DEBUG, "Setting objective coefficients...");
   // unsigned eid = 0;
   // for (pedge e = G.E.begin();
   // 	e != G.E.end();
   // 	++e ) {
   //    obj.setLinearCoef( var[ eid ], 1.0 );
   //    ++eid;
   // }

   // G.logg(DEBUG, "Setting constraint coefficients...");
   // unsigned tid = 0;
   // for (pangle t = G.T.begin();
   // 	t != G.T.end();
   // 	++t ) {
   //    if (triangle_valid( t )) {
   // 	 c[tid].setLinearCoef(var[ t->e1->eid ], 1.0);
   // 	 c[tid].setLinearCoef(var[ t->e2->eid ], 1.0);
   // 	 c[tid].setLinearCoef(var[ t->e3->eid ], 1.0);
   // 	 ++tid;
   //    }
   // }

   //   model.add(obj);
   //   model.add(c);

   G.logg(DEBUG, "Building cplex from model...");
   IloCplex cplex(model);
   cplex.setParam(IloCplex::IntParam::Threads, nThreads);
   cplex.setParam(IloCplex::IntParam::Threads, 1);
   cplex.setParam(IloCplex::ClockType, 1); //Tells cplex to use cpu time
   double hour = 60.0 * 60;
   cplex.setParam(IloCplex::TiLim, 10 * hour );
   cplex.setOut(env.getNullStream());
   //   cplex.setParam(IloCplex::EpAGap, 1.0 / (2 * G.E.size()));
   G.logg(DEBUG, "Solving IP...");
   cplex.solve();

   unsigned   eid = 0;

   v_sol.clear();
   v_sol.reserve( G.E.size() );
   for (pedge e = G.E.begin();
        e != G.E.end();
        ++e ) {
      v_sol.push_back( cplex.getValue( var[eid] ) );
      if ( cplex.getValue(var[eid]) == 1 )
	 e->in_S = true;
      ++eid;

   }

   G.logg(INFO, "Solution finished, objective value: " + to_string(cplex.getObjValue()));
   if (cplex.getCplexStatus() == IloCplex::Optimal)
      return static_cast<unsigned> ( cplex.getObjValue() );
   else
      return 0;
}


bool is_triangle_redundant( pedge& e ) {
   bool r = true;
   for (auto i = e->Delta.begin();
	i != e->Delta.end();
	++i ) {
      if ( triangle_valid(*i) ) {
	 r = false;
	 break;
      }
   }

   return r;
}

/*
 * Kortsarz et al.
 * 2-approximation
 *
 */
unsigned kortsarz( Graph& G, unsigned nThreads = 18, double max_hours = 4.0 ) {
   G.logg(INFO, "Starting Kortsarz et al....");
   vector< double > lpsol;

   bool b_solve;
   unsigned eid;
   unsigned sizeS = 0;
   IloEnv env;
   IloCplex cplex( env );
   do {
      bool status = solve_lp( G, lpsol, env, cplex, nThreads, max_hours );
      if (status == false)
	 return 0;
      
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
	 if ( is_triangle_redundant( e ) ) {
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

   return 0; //R.size() + sizeS;
}

/*
 * TARL algorithm
 */
unsigned tarl( Graph& G, unsigned nThreads = 18, double max_hours = 4.0 ) {
   vector< double > lpsol;
   IloEnv env;
   IloCplex cplex( env );

   bool status = solve_lp( G, lpsol, env, cplex, nThreads, max_hours );
   if (status == false)
      return 0;
   
   unsigned eid = 0;
   vector< pedge > W;

   for (pedge e = G.E.begin();
	e != G.E.end();
	++e ) {
      if ( lpsol[eid] < 1.0 / 5 ) {
	 e-> in_W = true;
	 W.push_back( e );
	 //	 lpsol[eid] = 0.0;
	    
	 if ( !(e-> in_S) ) {

	    //ADD other edges of triangles containing e
	    
	    for ( auto it = e->Delta.begin(); it != e->Delta.end(); ++it ) {
	       if (triangle_valid(*it)) {
		  pedge eToAdd;
		  //put in edge from this triangle with max weight
		  pangle& T = *it;
		  if ( lpsol[ T->e1->eid ] > lpsol[ T->e2->eid ] )
		     eToAdd = T->e1;
		  else
		     eToAdd = T->e2;

		  if ( lpsol[ T->e3->eid ] > lpsol[ eToAdd->eid ] )
		     eToAdd = T->e3;

		  eToAdd->in_S = true;
	       }
	    }
	 }
      }
     
      ++eid;
   }

   G.logg(DEBUG, "Starting bipartite (complement)...");
   vector< pedge > R;
   G.bipartite_complement( R );
   G.logg(DEBUG, "Size of R:" + to_string(R.size()));

   for (auto i = R.begin();
	i != R.end();
	++i) {
      (*i)->in_S = true;
   }
   
   return G.countS();
}


#endif
