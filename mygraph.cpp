#ifndef MYGRAPH_CPP
#define MYGRAPH_CPP
#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <map>
#include <unordered_set>
#include <thread>
#include <iomanip>
#include <cmath>
#include "logger.cpp"


using namespace std;

template<typename T>
void print_vector( vector< T >& v, ostream& os = cout ) {
  for (size_t i = 0; i < v.size(); ++i) {
    os << v[i] << ' ';
  }
  os << endl;
}

template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(
					  const std::vector<T>& vec,
					  Compare& compare)
{
  std::vector<std::size_t> p(vec.size());
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(),
	    [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
  return p;
}

template <typename T>
void apply_permutation_in_place(
				std::vector<T>& vec,
				const std::vector<std::size_t>& p)
{
  std::vector<bool> done(vec.size());
  for (std::size_t i = 0; i < vec.size(); ++i)
    {
      if (done[i])
	{
	  continue;
	}
      done[i] = true;
      std::size_t prev_j = i;
      std::size_t j = p[i];
      while (i != j)
	{
	  std::swap(vec[prev_j], vec[j]);
	  done[j] = true;
	  prev_j = j;
	  j = p[j];
	}
    }
}



namespace mygraph {

  typedef unsigned node_id;
  bool mycompare( const node_id& a, const node_id& b ) {
    return a < b;
  }  
  class Edge;
  class Triangle;

  typedef list< Edge >::iterator pedge; //this is how we will access edges
  typedef list< Triangle >::iterator pangle;
  
  class Triangle {
  public:
    pedge e1;
    pedge e2;
    pedge e3;

    Triangle() {
	
    }
     
    Triangle( pedge in1, pedge in2, pedge in3) {
      e1 = in1;
      e2 = in2;
      e3 = in3;
    }

    Triangle( const Triangle& rhs ) {
      e1 = rhs.e1;
      e2 = rhs.e2;
      e3 = rhs.e3;
    }

    Triangle& operator=( const Triangle& rhs ) {
      e1 = rhs.e1;
      e2 = rhs.e2;
      e3 = rhs.e3;

      return *this;
    }
  };
  
  bool operator==( const Triangle& lhs, const Triangle& rhs ) {
    bool b1 = (lhs.e1 == rhs.e1) || (lhs.e1 == rhs.e2) || (lhs.e1 == rhs.e3);
    bool b2 = (lhs.e2 == rhs.e1) || (lhs.e2 == rhs.e2) || (lhs.e2 == rhs.e3);
    bool b3 = (lhs.e3 == rhs.e1) || (lhs.e3 == rhs.e2) || (lhs.e3 == rhs.e3);

    return b1 && b2 && b3;
  }

  bool operator!=( const Triangle& lhs, const Triangle& rhs ) {
    return !(lhs == rhs);
  }
   
  class Edge {
  public:
    node_id from; //node id's
    node_id to;
     
    node_id eid;
    list< pangle > Delta; //the triangles containing this edge
    bool in_S; //whether this edge is selected for membership in S
    bool in_W; //whether this edge is selected for membership in W
    pangle T_e; //If this edge e is included in S, this triangle is h(e) (see paper)

    void print( ostream& os ) {
      os << "(" << from << "," << to << ")";
    }
      
    //Returns an incident vertex not equal to the input id
    node_id other( node_id one ) {
      if (from != one)
	return from;
      else
	return to;
    }
    
    bool incident( node_id id ) {
      return ( (from == id) || (to == id) );
    }

    Edge( node_id in_from, node_id in_to ) {
      in_S = false;
      in_W = false;
      from = in_from;
      to = in_to;
    }
  };
  
  /*
   * Vertex class
   */
  class Node {
  public:
    //    node_id id; Can allow id to be index in vector
    //for dynamic graph
    list< pedge > neighbors; //modified adjacency list
    list< node_id > nei_ids;
      
    //for static graph
    vector< pedge > v_neighbors;
    vector< node_id > v_nei_ids;
      
    bool in_A; // for bipartite procedure
    bool in_C; // for bipartite procedure
    Node() { in_A = false; in_C = false;
    }

    void insert_sort( pedge& e, node_id myId ) {
      node_id other_e = e->other( myId );
      node_id other_nei;
	 
      auto it = neighbors.begin();
      auto it2 = nei_ids.begin();
      while (it != neighbors.end()) {
	other_nei = (*it)->other( myId );
	if (other_e > other_nei) {
	  ++it;
	  ++it2;
	}
	else
	  break;
      }
      neighbors.insert( it, e );
      nei_ids.insert( it2, other_e );
    }
      
    void insert( pedge& e, node_id myId, node_id otherId ) {
      //	 neighbors.push_back( e );
      v_neighbors.push_back( e );
      //	 cerr << "Adding " << otherId << " to " << myId << "'s list" << endl;
      this->v_nei_ids.push_back( otherId );
      //	 print_vector( v_nei_ids, cerr );
    }

    void nid_insert( node_id nei_id ) {
      v_nei_ids.push_back( nei_id );
    }
      
    bool incident( node_id id ) {
      list< pedge >::iterator it1 = neighbors.begin();
      while (it1 != neighbors.end()) {
	if ( (*(*it1)).incident( id ) ) {
	  return true;
	}
	++it1;
      }

      return false;
    }

    bool incident_static( node_id id ) {
      vector< pedge >::iterator it1 = v_neighbors.begin();
      while (it1 != v_neighbors.end()) {
	if ( (*(*it1)).incident( id ) ) {
	  return true;
	}
	++it1;
      }

      return false;
    }

    Node( const Node& in ) {
      list< pedge >::const_iterator it1 = in.neighbors.begin();
      while (it1 != in.neighbors.end()) {
	neighbors.push_back( *it1 );
	++it1;
      }
      in_A = in.in_A;
      in_C = in.in_C;
      v_nei_ids.assign( in.v_nei_ids.begin(), in.v_nei_ids.end() );
      v_neighbors.assign( in.v_neighbors.begin(), in.v_neighbors.end() );
      nei_ids.assign( in.nei_ids.begin(), in.nei_ids.end() );
    }

    Node& operator=( const Node& in ) {
      neighbors.clear();
      list< pedge >::const_iterator it1 = in.neighbors.begin();
      while (it1 != in.neighbors.end()) {
	neighbors.push_back( *it1 );
	++it1;
      }

      in_A = in.in_A;
      in_C = in.in_C;
      v_nei_ids.assign( in.v_nei_ids.begin(), in.v_nei_ids.end() );
      v_neighbors.assign( in.v_neighbors.begin(), in.v_neighbors.end() );
      nei_ids.assign( in.nei_ids.begin(), in.nei_ids.end() );
      return *this;
    }

    
  };

  /*
   * For now, this only represents UNdirected graphs
   */
  
  class Graph {
  public:
    node_id n; //number of nodes
    node_id m; //number of edges
    vector< Node > V; //set of vertices
    list< Edge > E;   //set of edges
    vector< Edge > vE;
    list< Triangle > T;   //triangle list
    list< Triangle > Tsol;   //list of triangles in solution of DART
    Logger logg;
    unsigned sizeS;
    double runningTime;
     double tTriangles;
    double preprocessTime;
    
    void save( ostream& os ) {
	 
    }
      
    Graph() {
      n = 0;
      m = 0;
    }

    Graph( LogType LT, ostream& os, bool echo = false) : logg( LT, os, echo ) {
      n = 0;
      m = 0;
    }

    //     Graph( const Graph& in ) : logg( in.logg.loglevel, in.logg.of ) {
	
    //     }

    /*
     * Prepares data structure for efficient insertion, deletion of edges
     * By using std::list for neighboring edges in each node,
     * instead of vector
     */
    void init_dynamic( bool b_delete_static_info = false ) {
      for (size_t i = 0; i < V.size(); ++i ) {
	Node& node = V[i];
	node.neighbors.assign( node.v_neighbors.begin(), node.v_neighbors.end() );
	node.nei_ids.assign( node.v_nei_ids.begin(), node.v_nei_ids.end() );
	if (b_delete_static_info) {
	  node.v_neighbors.clear();
	  node.v_nei_ids.clear();
	}
      }
    }

    /*
     * Prepares data structure for efficient insertion, deletion of edges
     * By using std::list for neighboring edges in each node,
     * instead of vector
     */
    void init_static( ) {
      for (node_id i = 0; i < V.size(); ++i) {
	Node& node = V[i];
	node.v_neighbors.assign( node.neighbors.begin(), node.neighbors.end() );
	node.v_nei_ids.clear();
	for (auto it1 = node.neighbors.begin();
	     it1 != node.neighbors.end();
	     ++it1) {
	  node.v_nei_ids.push_back( (*it1)->other( i ) );
	}
      }
    }
      
    void clear_edges() {
      for (auto i = E.begin();
	   i != E.end();
	   ++i) {
	i->in_S = false;
	i->in_W = false;
      }
    }

     void clear_graph() {
	n = 0;
	m = 0;
	V.clear();
	E.clear();
	vE.clear();
	T.clear();
	Tsol.clear();
	runningTime = 0.0;
	preprocessTime = 0.0;
	sizeS = 0;
     }
     
    void init_empty_graph() {
      Node tmp;
      V.clear();
      V.assign(n, tmp );
      //	 for (unsigned i = 0; i < n; ++i) {
      //	    V.push_back( tmp );
      //
      //	 }
    }

    void print_graph( ostream& os ) {
      os << V.size() << ' ' << E.size() << endl;
      for (size_t i = 0; i < V.size(); ++i) {
	os << i << endl;
	for (auto e = V[i].neighbors.begin();
	     e != V[i].neighbors.end();
	     ++e ) {
	  (*e)->print(os);
	  os << ' ';
	}
	os << endl;
	for (unsigned j = 0; j < V[i].v_nei_ids.size(); ++j) {
	  os << V[i].v_nei_ids[j] << ' ';
	}
	os << endl;
      }
    }
      
    void bipartite( vector< pedge>& B ) {
      unsigned nC;
      unsigned nA;
      node_id this_node = 0;
      for (auto i = V.begin();
	   i != V.end();
	   ++i) {
	nC = 0;
	nA = 0;
	for (auto e = i->v_neighbors.begin();
	     e != i->v_neighbors.end();
	     ++e) {
	  if ( !( (*e)->in_S || (*e)->in_W ) ) { //Ignore edges in S or W
	    node_id other = (*e)->other( this_node );
	    if (V[other].in_A)
	      ++nA;
	    if (V[other].in_C)
	      ++nC;
	  }
	}

	if (nA >= nC) {
	  i->in_C = true;
	} else {
	  i->in_A = true;
	}
	
	++this_node;
      }

      //return the edges that cross the partition
      B.clear();
      for (pedge e = E.begin();
	   e != E.end();
	   ++e ) {
	if ( !( e->in_S || e->in_W ) ) { //Ignore edges in S or W
	  if ( V[ e->from ].in_A ) {
	    if ( V[ e->to ].in_C ) {
	      //e crosses
	      B.push_back( e );
	    }
	  } else {
	    if ( V[ e->to ].in_A ) {
	      //e crosses
	      B.push_back( e );
	    }
	  }
	    
	}
      }
    }

    /*
     * Returns complement of 'bipartite'
     */
    void bipartite_complement( vector< pedge>& R ) {
      unsigned nC;
      unsigned nA;
      node_id this_node = 0;
      for (auto i = V.begin();
	   i != V.end();
	   ++i) {
	nC = 0;
	nA = 0;
	for (auto e = i->v_neighbors.begin();
	     e != i->v_neighbors.end();
	     ++e) {
	  if ( !( (*e)->in_S || (*e)->in_W ) ) { //Ignore edges in S or W
	    node_id other = (*e)->other( this_node );
	    if (V[other].in_A)
	      ++nA;
	    if (V[other].in_C)
	      ++nC;
	  }

	}

	if (nA >= nC) {
	  i->in_C = true;
	} else {
	  i->in_A = true;
	}
	
	++this_node;
      }

      //return the edges that cross the partition
      R.clear();
      for (pedge e = E.begin();
	   e != E.end();
	   ++e ) {
	if ( !( e->in_S || e->in_W ) ) { //Ignore edges in S or W
	  if ( V[ e->from ].in_A ) {
	    if ( V[ e->to ].in_A ) {
	      //e crosses
	      R.push_back( e );
	    }
	  } else {
	    if ( V[ e->to ].in_C ) {
	      //e crosses
	      R.push_back( e );
	    }
	  }
	}
      }
    }

    void clear_triangles() {
      T.clear();
      for (auto it = E.begin();
	   it != E.end();
	   ++it) {
	it->Delta.clear();
      }
    }

    bool find_disjoint_triangle( pedge it1, Triangle& T_e, vector< pedge >& unprunedEdges ) {
      node_id& from = it1->from;
      node_id& to = it1->to;
      Node& From = V[ from ];
      Node& To = V[ to ];
      //Neighbor lists are in sorted order by other's id
      auto it2 = From.neighbors.begin();
      auto it3 = To.neighbors.begin();
      auto fend = From.neighbors.end();
      auto tend = To.neighbors.end();
      if (it2 == fend)
	return false;
      if (it3 == tend)
	return false;
      while (1) {
	if ( (*it2)->other( from ) < (*it3)->other( to ) ) {
	  ++it2;
	  if (it2 == fend)
	    break;
	} else {
	  if ( (*it2)->other( from ) > (*it3)->other( to ) ) {
	    ++it3;
	    if (it3 == tend)
	      break;
	  } else {
	    //found a triangle
	    //if ( !( ( (*it2)->in_S || (*it2)->in_W) || ( (*it3)->in_S || (*it3)->in_W ) ) ) {
	    if ((*it2)->in_W) {
	      (*it2)->in_W = false; //unprune
	      (*it2)->in_S = true; //unprune
	      unprunedEdges.push_back( *it2 );
	    } else {
	      if ((*it3)->in_W) {
		(*it3)->in_W = false; //unprune
		(*it3)->in_S = true; //unprune
		unprunedEdges.push_back( *it3 );
	      } else {
		if ( !( ( (*it2)->in_S ) || ( (*it3)->in_S) ) ) {
		  //this triangle is disjoint from S
		  Triangle t1( it1, *it2, *it3 );
		  T_e = t1;
		  return true;
		}
	      }
	    }
		  
	    ++it2;
	    ++it3;
	    if (it2 == fend || it3 == tend)
	      break;
	  }
	}
      }
	 
      //If we make it here,
      //No disjoint triangles from S containing it1 were found
      return false;
    }

    void add_triangle(node_id& v, size_t& v_s, size_t& v_t, size_t& s_t,
		      node_id& s, node_id& t ) {
      //get edges (v,s), (v,t), (s,t)
      pedge& vs = V[ v ].v_neighbors[ v_s ];
      pedge& vt = V[ v ].v_neighbors[ v_t ];
      pedge& st = V[ s ].v_neighbors[ s_t ];

      Triangle T_out;
      T_out.e1 = vs;
      T_out.e2 = vt;
      T_out.e3 = st;

      T.push_front(T_out);
      pangle tt = T.begin();
      vs->Delta.push_back( tt );
      vt->Delta.push_back( tt );
      st->Delta.push_back( tt );
    }
      
    void list_triangles() {
       clock_t t_start = clock();
      vector < vector< node_id > > A( n, vector< node_id >() );
      vector < vector< size_t > > I( n, vector< size_t >() ); //Indices of edges
      size_t count = 0;
	 
      for (node_id s = 0; s < n; ++s ) {
	for (size_t j = 0; j < V[s].v_nei_ids.size(); ++j) {
	  node_id& t = V[s].v_nei_ids[ j ];
	  if (s < t) {
	    vector< node_id >& A_s = A[s];
	    vector< node_id >& A_t = A[t];
	    vector< size_t >& I_s = I[s];
	    vector< size_t >& I_t = I[t];
	    auto it1 = A_s.begin();
	    auto it2 = A_t.begin();
	    if (it1 == A_s.end() || it2 == A_t.end() ) {
	      A_t.push_back( s );
	      I_t.push_back( j );
	      continue;
	    }
	    auto it3 = I_s.begin();
	    auto it4 = I_t.begin();
	    while (1) {
	      if (*it1 < *it2) {
		++it1;
		++it3;
		if (it1 == A_s.end()) {
		  break;
		}
	      } else {
		if (*it2 < *it1) {
		  ++it2;
		  ++it4;
		  if (it2 == A_t.end()) {
		    break;
		  }
		} else {
		  //found a triangle
		  add_triangle( *it1, *it3, *it4, j, s, t);
		  ++count;
		  ++it1; ++it3;
		  ++it2; ++it4;
		  if (it1 == A_s.end() || it2 == A_t.end() )
		    break;
		}
	      }
	    }
	    A_t.push_back( s );
	    I_t.push_back( j );
	  }
	}
      }
      tTriangles = double (clock() - t_start) / CLOCKS_PER_SEC;
      logg( INFO, to_string(count) + " triangles found.");
    }

    void check_triangles(){
      for (auto t = T.begin();
	   t != T.end();
	   ++t){
	if ( (t->e1 == t->e2) || (t->e2 == t->e3) || (t->e3 == t->e1)) {
	  logg( ERROR, "Incorrect triangle: ("
		+ to_string( t->e1->from ) + ","
		+ to_string( t->e1->to ) + ");"
		+ "(" + to_string( t->e2->from ) + ","
		+ to_string( t->e2->to ) + ");"
		+ "(" + to_string( t->e3->from ) + ","
		+ to_string( t->e3->to ) + ");" );
	}
      }
    }
     
    /*
     * Primal-dual algorithm
     */
    unsigned primal_dual() {
      logg(INFO, "Beginning primal-dual...");
      runningTime = 0.0;
      clock_t t_start = clock();
      unsigned size_S = 0;

      for (pangle it = T.begin();
	   it != T.end();
	   ++it ) {
	if (triangle_disjoint( it )) {
	  triangle_add( it );
	  size_S += 3;
	}
      }
      this->runningTime = double (clock() - t_start) / CLOCKS_PER_SEC;
      logg(INFO, "primal-dual finished, size of S: " + to_string(size_S));
      this->sizeS = size_S;
      this->runningTime += tTriangles;
      return size_S;
    }

    void integrated_prune( pedge& it1 ) {
      map< node_id, pedge > MFrom;
      Node From, To;
      From = V[ it1->from ];
      To = V[ it1->to ];
      //Go through From's edges and remember their other vertex
      MFrom.clear();
      for (list<pedge>::iterator it2 = From.neighbors.begin();
	   it2 != From.neighbors.end();
	   ++it2 ) {
	node_id other = (*it2)->other( it1->from );
	MFrom[ other ] = *it2; //this is an iterator of E
      }
      //now go through to's edges
      for (list<pedge>::iterator it2 = To.neighbors.begin();
	   it2 != To.neighbors.end();
	   ++it2 ) {
	node_id other = (*it2)->other( it1->to );
	map< node_id, pedge>::iterator it3 = MFrom.find( other );
	if (it3 != MFrom.end()) {
	  //we have found a triangle
	  if ( !((*it2)->in_S || it3->second->in_S) ) {
	    //this triangle needs edge it1 to be in the solution S
	    it1->in_S = true;
	    return;
	  }
	}
      }

      //If we make it here, it means we didn't find a triangle
      //that required it1 to be in S
      it1->in_S = false;
      it1->in_W = true;
      return;
    }

    /*
     * Dart-base integrated
     * means pruning step is integrated with dart-base
     *
     */ 
    // void dart_base_integrated() {
    // 	 for (auto e = E.begin(); e != E.end(); ++e) {
    // 	    Triangle T_e;
    // 	    if (!e->in_S) {
    // 	       if (find_disjoint_triangle( e, T_e )) {
    // 		  //Add T_e to S
    // 		  T_e.e1->in_S = true;
    // 		  T_e.e2->in_S = true;
    // 		  T_e.e3->in_S = true;
    // 		  Tsol.push_front( T_e );
    // 		  T_e.e1->T_e = Tsol.begin();
    // 		  T_e.e2->T_e = Tsol.begin();
    // 		  T_e.e3->T_e = Tsol.begin();

    // 		  integrated_prune( e );
    // 		  //e.eid < others.eid, so they will be pruned later
    // 	       }
    // 	    } else {
    // 	       //e was already added earlier. See if we can prune it
    // 	       integrated_prune( e );
    // 	    }
    // 	 }

    // }
      
    // /*
    //  * DART-BASE (no triangle listing)
    //  */
    // unsigned dart_base_free() {
    // 	 unsigned size_S = 0;

    // 	 for (auto e = E.begin(); e != E.end(); ++e) {
    // 	    Triangle T_e;
    // 	    if (!e->in_S) {
    // 	       if (find_disjoint_triangle( e, T_e )) {
    // 		  //Add T_e to S
    // 		  T_e.e1->in_S = true;
    // 		  T_e.e2->in_S = true;
    // 		  T_e.e3->in_S = true;
    // 		  size_S += 3;
    // 		  Tsol.push_front( T_e );
    // 		  T_e.e1->T_e = Tsol.begin();
    // 		  T_e.e2->T_e = Tsol.begin();
    // 		  T_e.e3->T_e = Tsol.begin();
    // 	       }
    // 	    }
    // 	 }

    // 	 logg(INFO, "DART_BASE_FREE finished, size of S: " + to_string(size_S));

    // 	 return size_S;
    // }


    bool prune_triangle(
			pedge& vs,
			pedge& vt,
			pedge& st
			) {
      if ((vs->in_S || vt->in_S)) {
	//(s,t) can be pruned , if it is in S
	if (st->in_S)
	  return true;
	else
	  return false;
      }

      //this triangle requires (s,t) to remain in S
      return false;
    }

      
    /*
     * prune_triangle
     * trying to prune (s,t), which we know is in S
     * returns true if can be pruned
     * false otherwise
     */
      
    bool prune_triangle(
			node_id& v, node_id& s, size_t& s_v,
			node_id& t, size_t& t_v,
			pedge& st
			) {
      //	 cerr << "(v,s): " << v << ' ' << v_s << ' ' << V[v].v_neighbors.size() << endl;
      //	 cerr << "(v,t): " << v << ' ' << v_t << ' ' << V[v].v_neighbors.size() << endl;

      pedge& vs = V[ s ].v_neighbors[ s_v ];
      pedge& vt = V[ t ].v_neighbors[ t_v ];
      //	 pedge& st = V[ s ].v_neighbors[ s_t ];
	 
      if ((vs->in_S || vt->in_S)) {
	//(s,t) can be pruned , if it is in S
	if (st->in_S)
	  return true;
	else
	  return false;
      }

      // //Correct overly aggressive pruning
      // if ( !(st->in_S ) ) {
      //    //vs or vt must have been pruned. Put one back in
      //    if (vs->in_W) {
      //       vs->in_W = false;
      //       vs->in_S = true;
      //       return false;
      //    }

      //    if (vt->in_W) {
      //       vt->in_W = false;
      //       vt->in_S = true;
      //       return false;
      //    }
      // }

      //this triangle requires (s,t) to remain in S
      return false;
    }
      

    bool free_triangle(node_id& v, size_t& v_s, size_t& v_t, size_t& s_t,
		       node_id& s, node_id& t ) {
      pedge& vs = V[ v ].v_neighbors[ v_s ];
      pedge& vt = V[ v ].v_neighbors[ v_t ];
      pedge& st = V[ s ].v_neighbors[ s_t ];
	 
      if (!(vs->in_S || vt->in_S || st->in_S)) {
	//triangle is disjoint
	Triangle t1( vs, vt, st );
	Tsol.push_front( t1 );
	vs->T_e = Tsol.begin();
	vt->T_e = Tsol.begin();
	st->T_e = Tsol.begin();
	vs->in_S = true;
	vt->in_S = true;
	st->in_S = true;
	return true;
      }

      return false;
    }

    bool dynamic_add_edge( node_id from, node_id to, pedge& e_added ) {
      if (from == to)
	return false;
      if (V[from].incident( to ) ) {
	return false;
      }
      //Add edge to graph
      Edge ee( from, to );
      E.push_front( ee );
      e_added = E.begin();
      V[from].insert_sort( e_added, from );
      V[to].insert_sort( e_added, to );
      return true;   
    }
      
    /*
     * Update dart solution upon new edge
     * Returns time to update solution
     */
    double dart_add_edge( node_id from, node_id to, pedge& e ) {
      clock_t t_start = clock();
      //before proceeding, locally unprune the solution
      //	 local_unprune( from, to );
      //Find triangles containing e
      //If disjoint one found, add it to solution, and terminate
      Triangle T_e;
      vector< pedge > unprunedEdges;
      if ( find_disjoint_triangle( e, T_e, unprunedEdges ) ) {
	//Add T_e to S
	T_e.e1->in_S = true;
	T_e.e2->in_S = true;
	T_e.e3->in_S = true;
	Tsol.push_front( T_e );
	T_e.e1->T_e = Tsol.begin();
	T_e.e2->T_e = Tsol.begin();
	T_e.e3->T_e = Tsol.begin();
	unprunedEdges.push_back( T_e.e1 );
	unprunedEdges.push_back( T_e.e2 );
	unprunedEdges.push_back( T_e.e3 );
      }

      prune( unprunedEdges );
      return double (clock() - t_start) / CLOCKS_PER_SEC;
    }

    /*
     * Update dart solution upon edge removal
     * Returns time to update solution
     */
    double dart_remove_edge( pedge& e, bool expensivePruning = false ) {
      clock_t t_start = clock();
      if (e->in_S || e->in_W) {
	pangle T_e = e->T_e;
	//get edges f,g, where T_e = efg
	pedge f; pedge g;
	if (T_e->e1 == e) {
	  f = T_e->e2;
	  g = T_e->e3;
	}
	if (T_e->e2 == e) {
	  f = T_e->e1;
	  g = T_e->e3;
	}
	if (T_e->e3 == e) {
	  f = T_e->e1;
	  g = T_e->e2;
	}

	if (expensivePruning) {
	  node_id from = e->from;
	  node_id to = e->to;
	  //Delete e from the graph structure
	  remove_edge( e );
	  //Need to delete T_e, as well
	  Tsol.erase( T_e );

	  //Removing e might enable some edges to be pruned
	  //Attempt to do so, but this step could be expensive
	  local_reprune( from, to );
	} else {
	  //Simply delete the edge
	  remove_edge( e );
	  //Need to delete T_e, as well
	  Tsol.erase( T_e );
	}
	    
	//Remove f,g from solution
	g->in_S = false;
	g->in_W = false;
	f->in_S = false;
	f->in_W = false;
	vector< pedge > unprunedEdges;
	Triangle T_f;
	if ( find_disjoint_triangle( f, T_f, unprunedEdges ) ) {
	  //Add T_f to S
	  T_f.e1->in_S = true;
	  T_f.e2->in_S = true;
	  T_f.e3->in_S = true;
	  Tsol.push_front( T_f );
	  T_f.e1->T_e = Tsol.begin();
	  T_f.e2->T_e = Tsol.begin();
	  T_f.e3->T_e = Tsol.begin();
	  unprunedEdges.push_back( T_f.e1 );
	  unprunedEdges.push_back( T_f.e2 );
	  unprunedEdges.push_back( T_f.e3 );
	}
	Triangle T_g;
	if ( find_disjoint_triangle( g, T_g, unprunedEdges ) ) {
	  //Add T_g to S
	  T_g.e1->in_S = true;
	  T_g.e2->in_S = true;
	  T_g.e3->in_S = true;
	  Tsol.push_front( T_g );
	  T_g.e1->T_e = Tsol.begin();
	  T_g.e2->T_e = Tsol.begin();
	  T_g.e3->T_e = Tsol.begin();
	  unprunedEdges.push_back( T_g.e1 );
	  unprunedEdges.push_back( T_g.e2 );
	  unprunedEdges.push_back( T_g.e3 );
	}

	prune( unprunedEdges );
      } else {
	//simply remove edge from graph
	remove_edge( e );
      }
	 
      return double (clock() - t_start) / CLOCKS_PER_SEC;
    }
      
    /*
     * DART-BASE 
     */
    unsigned dart_base() {
      clock_t t_start = clock();
      vector < vector< node_id > > A( n, vector< node_id >() );
      vector < vector< size_t > > I( n, vector< size_t >() ); //Indices of edges
      unsigned countS = 0;
	 
      for (node_id s = 0; s < n; ++s ) {
	for (size_t j = 0; j < V[s].v_nei_ids.size(); ++j) {
	  node_id& t = V[s].v_nei_ids[ j ];
	  if (s < t) {
	    vector< node_id >& A_s = A[s];
	    vector< node_id >& A_t = A[t];
	    vector< size_t >& I_s = I[s];
	    vector< size_t >& I_t = I[t];
	    auto it1 = A_s.begin();
	    auto it2 = A_t.begin();
	    if (it1 == A_s.end() || it2 == A_t.end() ) {
	      A_t.push_back( s );
	      I_t.push_back( j );
	      continue;
	    }
	    auto it3 = I_s.begin();
	    auto it4 = I_t.begin();
	    while (1) {
	      if (*it1 < *it2) {
		++it1;
		++it3;
		if (it1 == A_s.end()) {
		  break;
		}
	      } else {
		if (*it2 < *it1) {
		  ++it2;
		  ++it4;
		  if (it2 == A_t.end()) {
		    break;
		  }
		} else {
		  //found a triangle
		  if (free_triangle( *it1, *it3, *it4, j, s, t))
		    countS += 3;
			   
		  ++it1; ++it3;
		  ++it2; ++it4;
		  if (it1 == A_s.end() || it2 == A_t.end() )
		    break;
		}
	      }
	    }
	    A_t.push_back( s );
	    I_t.push_back( j );
	  }
	}
      }
      
      free_prune();
      
      runningTime = double (clock() - t_start) / CLOCKS_PER_SEC;
      logg(INFO, "DART_BASE finished, size of S: " + to_string(countS));

      return countS;
    }


    void free_prune_process( pedge& it1 ) {
      if (!(it1->in_S))
	return;
	 
      node_id& from = it1->from;
      node_id& to = it1->to;
      Node& From = V[ from ];
      Node& To = V[ to ];

      //Neighbor lists are in sorted order by other's id
      auto it2 = From.neighbors.begin();
      auto it3 = To.neighbors.begin();
      //we only want to look at greater edges to avoid
      //recomputing triangles
      auto fend = From.neighbors.end();
      auto tend = To.neighbors.end();
      if (it2 == fend) {
	it1->in_S = false;
	it1->in_W = true;
	return;
      }
      if (it3 == tend) {
	it1->in_S = false;
	it1->in_W = true;
	return;
      }
      while (1) {
	if ( (*it2)->other( from ) < (*it3)->other( to ) ) {
	  ++it2;
	  if (it2 == fend)
	    break;
	} else {
	  if ( (*it2)->other( from ) > (*it3)->other( to ) ) {
	    ++it3;
	    if (it3 == tend)
	      break;
	  } else {
	    //found a triangle
	    if ( !((*it2)->in_S || (*it3)->in_S) ) {
	      //this triangle needs edge it1 to be in the solution S
	      it1->in_S = true;
	      it1->in_W = false;
	      return;
	    }

	    ++it2;
	    ++it3;
	    if (it2 == fend || it3 == tend)
	      break;
	  }
	}
      }
      //prune this edge
      it1->in_S = false;
      it1->in_W = true;
      return;
    }

    void free_prune_process_faster( pedge& st ) {
      if (!(st->in_S)) //only want to prune edges in S
	return;
      node_id& s = st->from;
      node_id& t = st->to;
      list< node_id >& N_s = V[s].nei_ids;
      list< node_id >& N_t = V[t].nei_ids;
      auto it1 = N_s.begin();
      auto it2 = N_t.begin();
      auto it3 = V[s].neighbors.begin();
      auto it4 = V[t].neighbors.begin();
      if (it1 == N_s.end() || it2 == N_t.end() ) {
	st->in_S = false;
	st->in_W = true;
	return;
      }
      while (1) {
	if (*it1 < *it2) {
	  ++it1;
	  ++it3;
	  if (it1 == N_s.end()) {
	    break;
	  }
	} else {
	  if (*it2 < *it1) {
	    ++it2;
	    ++it4;
	    if (it2 == N_t.end()) {
	      break;
	    }
	  } else {
	    //found a triangle
	    if (!prune_triangle( *it3, *it4 , st)) {
	      st->in_S = true;
	      st->in_W = false;
	      return;
	    }
	    ++it1; ++it3;
	    ++it2; ++it4;
	    if (it1 == N_s.end() || it2 == N_t.end() )
	      break;
	  }
	}
      }

      //prune this edge
      st->in_S = false;
      st->in_W = true;
      return;
    }

    void local_reprune( node_id& from, node_id& to ) {
      //	 node_id& from = it->from;
      //	 node_id& to = it->to;
      list< pedge >& N_from = V[ from ].neighbors;
      list< pedge >& N_to = V[ to ].neighbors;
      //	 free_prune_process_faster( it );
	 
      for ( auto e = N_from.begin();
	    e != N_from.end();
	    ++e ) {
	free_prune_process( *e );
      }

      for ( auto e = N_to.begin();
	    e != N_to.end();
	    ++e ) {
	free_prune_process( *e );
      }
    }

    void prune( vector< pedge >& unprunedEdges ) {
      for (auto it = unprunedEdges.begin();
	   it != unprunedEdges.end();
	   ++it) {
	free_prune_process( *it );
      }
    }
      
    void local_unprune( node_id from, node_id to ) {
      list< pedge >& N_from = V[ from ].neighbors;
      list< pedge >& N_to = V[ to ].neighbors;

      for ( auto e = N_from.begin();
	    e != N_from.end();
	    ++e ) {
	if ((*e)->in_W) {
	  //e has been pruned. Unprune it
	  (*e)->in_S = true;
	  (*e)->in_W = false;
	}
      }

      for ( auto e = N_to.begin();
	    e != N_to.end();
	    ++e ) {
	if ((*e)->in_W) {
	  //e has been pruned. Unprune it
	  (*e)->in_S = true;
	  (*e)->in_W = false;
	}
      }
    }
      
    void free_prune() {
      logg(DEBUG, "Starting free_prune()");
      clock_t t_start = clock();
      size_t prune_count = 0;
      for (auto st = E.begin();
	   st != E.end();
	   ++st ) {
	if (!(st->in_S)) //only want to prune edges in S
	  continue;
	bool st_prunable = true; //whether we can prune (s,t)
	node_id& s = st->from;
	node_id& t = st->to;
	vector< node_id >& N_s = V[s].v_nei_ids;
	vector< node_id >& N_t = V[t].v_nei_ids;
	size_t N_s_index = 0;
	size_t N_t_index = 0;
	auto it1 = N_s.begin();
	auto it2 = N_t.begin();
	if (it1 == N_s.end() || it2 == N_t.end() ) {
	  continue;
	}
	while (1) {
	  if (*it1 < *it2) {
	    ++it1;
	    ++N_s_index;
	    if (it1 == N_s.end()) {
	      break;
	    }
	  } else {
	    if (*it2 < *it1) {
	      ++it2;
	      ++N_t_index;
	      if (it2 == N_t.end()) {
		break;
	      }
	    } else {
	      //found a triangle
	      ++prune_count;
	      if (!prune_triangle( *it1, s, N_s_index, t, N_t_index, st)) {
		st_prunable = false;
		break;
	      }


		     
	      ++it1; ++N_s_index;
	      ++it2; ++N_t_index;
	      if (it1 == N_s.end() || it2 == N_t.end() )
		break;
	    }
	  }
	}

	if (st_prunable) {
	  st->in_S = false;
	  st->in_W = true;
	}
      }

      runningTime += double (clock() - t_start) / CLOCKS_PER_SEC;
      logg( DEBUG, "free-prune examined " + to_string(prune_count) + " triangles");
      return;
    }

    /*
     * DART-BASE
     */
    unsigned dart_base_rand() {
      logg(INFO, "Beginning DART_BASE_RAND...");

      logg(INFO, "Creating triangle vector...");
      vector< pangle > T_vec;
      T_vec.reserve( T.size() );
      for (auto t = T.begin();
	   t != T.end();
	   ++t) {
	T_vec.push_back( t );
      }

      unsigned size_S = 10000000;
      for (unsigned j = 0; j < 10000; ++j) {
	clear_edges();
	unsigned tmp_size_S = 0;
	random_shuffle ( T_vec.begin(), T_vec.end() );
	for (unsigned  i = 0;
	     i < T_vec.size();
	     ++i ) {
	  if (triangle_disjoint( T_vec[i] )) {
	    triangle_add( T_vec[i] );
	    tmp_size_S += 3;
	  }
	}
	if (tmp_size_S < size_S)
	  size_S = tmp_size_S;
	//	 logg(INFO, "tmp_size_S " + to_string(tmp_size_S));
      }
      logg(INFO, "DART_BASE_RAND finished, size of S: " + to_string(size_S));

      return size_S;
    }

    /*
     * DART-BASE-HEU
     */
    unsigned dart_base_heu() {
      logg(INFO, "Beginning DART_BASE_HEU...");
      unsigned size_S = 0;

      for (pangle it = T.begin();
	   it != T.end();
	   ++it ) {
	if (triangle_disjoint( it )) {
	  triangle_add( it );
	  size_S += 3;
	}
      }

      logg(INFO, "DART_BASE finished, size of S: " + to_string(size_S));
      logg(INFO, "Pruning solution...");
      unsigned pruned = better_prune_S();
      //      unsigned pruned = prune_S();
      logg(INFO, to_string(pruned) + " edges pruned");
      logg(INFO, "Solution size: " + to_string(size_S - pruned));
					      
      return (size_S - pruned);
    }

    bool triangle_valid( pangle& t ) {
      return !( t->e1->in_S || t->e2->in_S || t->e3->in_S );
    }

    unsigned countTvalid() {
      unsigned count = 0;
      for (auto i = T.begin();
	   i != T.end();
	   ++i) {
	
	if (triangle_valid(i))
	  ++count;
      }
      return count;
    }

    unsigned countS() {
      unsigned count = 0;
      for (auto i = E.begin();
	   i != E.end();
	   ++i) {
	
	if (i->in_S)
	  ++count;
      }
      sizeS = count;
      return count;
    }

    void info_check( ostream& os ) {
      os << "GRAPH INFO: "
	 << V.size() << ' ' << E.size() << ' ' << T.size() << endl;
      os << "Number of triangles valid: ";
      unsigned count = 0;
      for (auto t = T.begin();
	   t != T.end();
	   ++t ) {
	if (triangle_valid( t ))
	  ++count;
      }
      os << count << endl;
    }
     
    unsigned prune_S() {
      unsigned n_pruned = 0;
      for (auto e = E.begin();
	   e != E.end();
	   ++e) {
	if (e->in_S) {
	  e->in_S = false;
	  if (!( e->T_e->e1->in_S || e->T_e->e2->in_S || e->T_e->e3->in_S)) {
	    e->in_S = true;
	    continue;
	  }
	  for (auto t = e->Delta.begin();
	       t != e->Delta.end();
	       ++t) {
	    //	 if (triangle_valid( *t ))
	    //		    ++n_triangles;

		 
		 
	    if ( (*t) != e->T_e ) {
	      e->in_S = true;
	      break;
	    }
	  }
	  if (!e->in_S)
	    ++n_pruned;

	   
	}
      }
      return n_pruned;
    }

    bool triangle_counts( Triangle& T, pedge& e ) {
      if (T == *(e->T_e))
	return false;
      if (T.e1 != e && T.e1->in_S)
	return false;
      if (T.e2 != e && T.e2->in_S)
	return false;
      if (T.e3 != e && T.e3->in_S)
	return false;
	
      return true;
    }      
    bool triangle_counts( pangle& t, pedge& e ) {
      if (t == e->T_e)
	return false;
      if (t->e1 != e && t->e1->in_S)
	return false;
      if (t->e2 != e && t->e2->in_S)
	return false;
      if (t->e3 != e && t->e3->in_S)
	return false;
	
      return true;
    }

    /*
     * Ensures all edges in solution are in disjoint triangles
     */
    bool ensure_validity() {
      unsigned countIncidentS;
      for (auto e = E.begin();
	   e != E.end();
	   ++e ) {
	if (e->in_S) {
	  countIncidentS = 1;
	  node_id& from = e->from;
	  node_id& to = e->to;
	  auto f = V[from].v_neighbors.begin();
	  auto g = V[to].v_neighbors.begin();
	  auto it1 = V[from].v_nei_ids.begin();
	  auto it2 = V[to].v_nei_ids.begin();
	  vector< node_id >& A_s = V[from].v_nei_ids;
	  vector< node_id >& A_t = V[to].v_nei_ids;
	  while (1) {
	    if (*it1 < *it2) {
	      ++it1; ++f;
	      if (it1 == A_s.end()) {
		break;
	      }
	    } else {
	      if (*it2 < *it1) {
		++it2; ++g;
		if (it2 == A_t.end()) {
		  break;
		}
	      } else {
		if ((*f)->in_S)
		  ++countIncidentS;
		if ((*g)->in_S)
		  ++countIncidentS;
		++it1; ++f;
		++it2; ++g;
		if (it1 == A_s.end() || it2 == A_t.end() )
		  break;
	      }
	    }
	  }
	  if (countIncidentS > 3)
	    return false;
	}
      }

      return true;
    }
      
    bool ensure_feasibility() {
      for (auto t = T.begin();
	   t != T.end();
	   ++t) {
	if ( triangle_valid( t ) )
	  return false;
      }
      return true;
    }

    size_t count_infeasible() {
      size_t count= 0;
      for (auto t = T.begin();
	   t != T.end();
	   ++t) {
	if ( triangle_valid( t ) )
	  ++count;
      }
      return count;
    }
     
    unsigned better_prune_S() {
      unsigned n_pruned = 0;
      for (auto e = E.begin();
	   e != E.end();
	   ++e) {
	if (e->in_S) {
	  e->in_S = false;
	  if (!( e->T_e->e1->in_S || e->T_e->e2->in_S || e->T_e->e3->in_S)) {
	    e->in_S = true;
	    continue;
	  }
	  for (auto t = e->Delta.begin();
	       t != e->Delta.end();
	       ++t) {
	    //	 if (triangle_valid( *t ))
	    //		    ++n_triangles;
	    if (triangle_counts( *t, e )) {
	      e->in_S = true;
	      break;
	    }
	  }
	  if (!e->in_S)
	    ++n_pruned;

	   
	}
      }
      return n_pruned;
    }
  
    /*Adds all edges of triangle
     * t to S
     */
    void triangle_add( pangle& t ) {
      t -> e1 -> in_S = true;
      t -> e1 -> T_e = t;
      t -> e2 -> in_S = true;
      t -> e2 -> T_e = t;
      t -> e3 -> in_S = true;
      t -> e3 -> T_e = t;
    }
    
    /*
     * Checks to see if S hits a triangle
     */
    bool triangle_disjoint( const pangle& t ) {
      if ( t -> e1 -> in_S )
	return false;

      if ( t -> e2 -> in_S )
	return false;

      if ( t -> e3 -> in_S )
	return false;

      return true;
    }

    /* 
     * Removes tt from e's triangle list
     */
    void update_triangle_list(
			      pangle& tt,
			      pedge& e ) {
      for (auto i = e->Delta.begin();
	   i != e->Delta.end();
	   ++i ) {
	if (*i == tt) {
	  //Erase this triangle, tt
	  e->Delta.erase( i );
	  break;
	}
      }
    }
     
    /*
     * Input is an iterator to an element
     * in the edge list
     */
    pedge remove_edge( pedge& e_to_remove ) {
      //remove from its incident nodes' adjacency list
      node_id from = e_to_remove->from;
      node_id to = e_to_remove->to;
      auto it2 = V[ from ].neighbors.begin();
      while (*it2 != e_to_remove)
	++it2;
      V[ from ].neighbors.erase( it2 );
      it2 = V[ to ].neighbors.begin();
      while (*it2 != e_to_remove)
	++it2;
      V[ to ].neighbors.erase( it2 );

      //Finally, discard the edge itself
      return E.erase( e_to_remove );
    }

    //Only works with undirected edges atm
    void add_edge( node_id from, node_id to ) {
      if (from == to) { //don't allow self-loops
	logg(WARN, "Self-loop (" + to_string( from ) + "," + to_string( to ) + "), not added to graph");
	return;
      }
      if ( from >= n || to >= n ) {
	//	logg(INFO, "Node index out of bounds: ignoring edge (" + to_string( from ) + "," + to_string( to ) + "), since n = " + to_string( n ) );

	//augment graph to have enough nodes
	n = from + 1;
	if (to + 1 > n)
	  n = to + 1;

	Node tmp;
	while (V.size() < n) {
	  V.push_back( tmp );
	}

      }
      
      Node& From = V[ from ];
      Node& To = V[ to ];

      //check if edge already exists
      if ( From.incident_static( to ) ) {
	logg(WARN, "Not adding edge (" + to_string( from ) + "," + to_string( to ) + "), since it already exists");
	return;
      }

      //Need to add a new edge to the graph
      //	 logg(DEBUG, "Adding edge (" + to_string( from ) + "," + to_string( to ) + ")");
      Edge edge( from, to );
      edge.eid = m;

      E.push_back( edge );
      pedge ee = --E.end();
      From.insert( ee, from, to );
      To.insert( ee, to, from );
	 
      //	 From.nid_insert( to );
      //	 To.nid_insert( from );
      //	 From.v_neighbors.push_back( --E.end() );
      //	 To.v_neighbors.push_back( --E.end() );
      ++m; //number of edges has increased
    }

    bool verify_graph() {
      for (node_id i = 0; i < n; ++i) {
	auto it1 = V[i].neighbors.begin();
	auto it2 = V[i].v_nei_ids.begin();
	auto it3 = V[i].v_neighbors.begin();
	if (V[i].neighbors.size() != V[i].v_nei_ids.size()) {
	  logg( ERROR, "Vertex " + to_string(i) );
	  logg( ERROR, "neighbors.size " + to_string( V[i].neighbors.size() )
		+ " v_nei_ids.size " + to_string( V[i].v_nei_ids.size())
		+ " v_neighbors.size " + to_string( V[i].v_neighbors.size()));
	  return false;
	}
	if (V[i].neighbors.size() != V[i].v_neighbors.size()) {
	  logg( ERROR, "neighbors.size " + to_string( V[i].neighbors.size() )
		+ " v_neighbors.size " + to_string( V[i].v_neighbors.size() ));

	  return false;
	}
	while (it1 != V[i].neighbors.end()) {
	  if (*it1 != *it3)
	    return false;
	  if ( (*it1)->other( i ) != *it2 )
	    return false;

	  ++it1; ++it2; ++it3;
	}
      }

      return true;
    }

    void read_edge_list_bin( string fname ) {
      ifstream ifile ( fname.c_str(), ios::in | ios::binary );
      unsigned m; //number of edges in list
      unsigned n;
      ifile.read( (char*) &n, sizeof( unsigned ) );
      ifile.read( (char*) &m, sizeof( unsigned ) );
      unsigned* fromto_arr = new unsigned [2 * m];
	 
      ifile.read( (char *) fromto_arr,  2 * m *sizeof( unsigned ) );

      logg(INFO, "File input finished. Constructing graph..." );
      this->n = n;
      init_empty_graph();
      for (unsigned i = 0; i < m; ++i) {
	unsigned from = fromto_arr[ 2 * i ];
	unsigned to = fromto_arr[ 2 * i + 1 ];

	add_edge( from, to );
      }

      delete [] fromto_arr;

      // unsigned from,to;
      // for (unsigned i = 0; i < m; ++i) {
      //    ifile.read( (char *) &from,  sizeof( unsigned ) );
      //    ifile.read( (char *) &to,  sizeof( unsigned ) );

      //    add_edge( from, to );
      // }
      // ifile.close();

	 
      logg(INFO, "Graph constructed: n = " + to_string(n) + ", m = " + to_string(m));

      logg(INFO, "Sorting neighbor ids..." );
      clock_t t_start = clock();
      for ( auto v = V.begin(); v != V.end(); ++v ) {
	auto p = sort_permutation( v->v_nei_ids, mycompare );
	//				       [](const node_id& a, const node_id& b){ return a < b; } );
	apply_permutation_in_place( v->v_nei_ids, p );
	//	    sort( v->v_nei_ids.begin(), v->v_nei_ids.end() );
	apply_permutation_in_place( v->v_neighbors, p );
	//	    v->neighbors.assign( v->v_neighbors.begin(),
	//				 v->v_neighbors.end() );
      }
      preprocessTime = double (clock() - t_start) / CLOCKS_PER_SEC;
    }
      
    void read_edge_list( string fname ) {
      ifstream ifile ( fname.c_str() ); 
      string line;
      unsigned linenum = 0;
      istringstream iss;
      node_id from, to; //for each edge

      logg(INFO, "Reading graph from edge list: " + fname);
      
      while ( getline( ifile, line ) ) {
	if ( line[0]  != '#' && line[0] != ' ') {
	  if (line.size() == 0)
	    break;
	  
	  iss.clear();
	  iss.str( line );
	  //need to add an edge 
	  iss >> from;
	  iss >> to;
	  
	  add_edge( from, to );
	  //	       print_vector( V[from].v_nei_ids, cerr );
	  //	       print_vector( V[to].v_nei_ids, cerr );
	  ++linenum;
	}
      }

      //Graph has been constructed
      ifile.close();

      logg(INFO, "Graph constructed: n = " + to_string(n) + ", m = " + to_string(m));

      logg(INFO, "Sorting neighbor ids..." );
      clock_t t_start = clock();
      for ( auto v = V.begin(); v != V.end(); ++v ) {
	auto p = sort_permutation( v->v_nei_ids, mycompare );
	//				       [](const node_id& a, const node_id& b){ return a < b; } );
	apply_permutation_in_place( v->v_nei_ids, p );
	//	    sort( v->v_nei_ids.begin(), v->v_nei_ids.end() );
	apply_permutation_in_place( v->v_neighbors, p );
	//	    v->neighbors.assign( v->v_neighbors.begin(),
	//				 v->v_neighbors.end() );
      }
      preprocessTime = double (clock() - t_start) / CLOCKS_PER_SEC;
      // logg(INFO, "Copying edges..." );
      // vE.assign(E.begin(),E.end());
      // for (size_t i = 0; i < V.size(); ++i) {
      //    logg(DEBUG, to_string(i) + " neighbors: ");
      //    for (unsigned j= 0; j < V[i].v_nei_ids.size(); ++j) {
      //       logg(DEBUG, to_string( V[i].v_nei_ids[j] ) );
      //    }
      // }
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
     
  };

  class smTriangle {
  public:
    unsigned n1;
    unsigned n2;
    unsigned n3;

    smTriangle() {

    }
      
    smTriangle( node_id in1, node_id in2, node_id in3 ) {
      n1 = in1;
      n2 = in2;
      n3 = in3;
    }

    smTriangle ( const smTriangle& rhs ) {
      n1 = rhs.n1;
      n2 = rhs.n2;
      n3 = rhs.n3;
    }
  };

  class smEdge {
  public:
    node_id from;
    node_id to;

    smEdge() {

    }
    
    smEdge( const smEdge& rhs ) {
      from = rhs.from;
      to = rhs.to;
    }
  };

  struct smEdge_hash {
    node_id operator()(const smEdge& e) const {
      return e.from + e.to;
    }
  };

  bool operator==( const smEdge& lhs, const smEdge& rhs ) {
    //      return (((lhs.from == rhs.from) && (lhs.to == rhs.to)) || ((lhs.from == rhs.to) && (lhs.to == rhs.from)));
    return ( (lhs.from == rhs.from) && (lhs.to == rhs.to) );
  }

  bool operator<( const smEdge& lhs, const smEdge& rhs ) {
    if (lhs.from < rhs.from)
      return true;
    else {
      if (rhs.from < lhs.from) {
	return false;
      } else {
	return (lhs.to < rhs.to);
      }
    }
  }
   
  // class smallGraph {
  // public:
  //   vector< vector< node_id > > adjList;
  //   unsigned n;
  //   Logger logg;
  //   //      map< smEdge, smTriangle > S; 
  //   unordered_set< smEdge, smEdge_hash > S;
  //   vector< smTriangle > T; //Triangle listing 
      
  //   bool inS( smEdge& st ) {
  //     return (S.find( st ) != S.end());
  //   }
      
  //   smallGraph() {
  //     n = 0;
  //   }

  //   void init_empty_graph() {
  //     vector< unsigned > emptyList;
  //     adjList.assign(n, emptyList);
  //   }

  //   void add_edge( unsigned from, unsigned to ) {
  //     adjList[ from ].push_back( to );
  //     adjList[ to ].push_back( from );
  //   }
      
  //   void read_edge_list_bin( string fname ) {
  //     ifstream ifile ( fname.c_str(), ios::in | ios::binary );
  //     unsigned m; //number of edges in list
  //     unsigned n;
  //     ifile.read( (char*) &n, sizeof( unsigned ) );
  //     ifile.read( (char*) &m, sizeof( unsigned ) );
  //     this->n = n;
  //     init_empty_graph();


  //     // unsigned* fromto_arr = new unsigned [2 * m];
	 
  //     // ifile.read( (char *) fromto_arr,  2 * m *sizeof( unsigned ) );

  //     // logg(INFO, "File input finished. Constructing graph..." );

  //     // for (unsigned i = 0; i < m; ++i) {
  //     //    unsigned from = fromto_arr[ 2 * i ];
  //     //    unsigned to = fromto_arr[ 2 * i + 1 ];

  //     //    add_edge( from, to );
  //     // }

  //     // delete [] fromto_arr;

  //     unsigned from,to;
  //     for (unsigned i = 0; i < m; ++i) {
  // 	ifile.read( (char *) &from,  sizeof( unsigned ) );
  // 	ifile.read( (char *) &to,  sizeof( unsigned ) );

  // 	add_edge( from, to );
  //     }
  //     ifile.close();
	 
  //     logg(INFO, "Graph constructed: n = " + to_string(n) + ", m = " + to_string(m));

  //     logg(INFO, "Sorting neighbor ids..." );
  //     for (unsigned i = 0; i < n; ++i) {
  // 	sort( adjList[i].begin(), adjList[i].end() );
  //     }
  //   }

  //   bool free_triangle( node_id& j, node_id& s, node_id& t ) {
  //     //know s < t. Need to know relation between other two pairs
  //     smEdge st;
  //     st.from = s;
  //     st.to = t;
  //     smEdge js;
  //     smEdge jt;
  //     if (j < s) {
  // 	js.from = j;
  // 	js.to = s;
  // 	jt.from = j;
  // 	jt.to = t;
  //     } else {
  // 	js.from = s;
  // 	js.to = j;
  // 	if (j < t) {
  // 	  jt.from = j;
  // 	  jt.to = t;
  // 	} else {
  // 	  jt.from = t;
  // 	  jt.to = j;
  // 	}
  //     }
	    
  //     if ( ! (inS( js ) || inS(jt) || inS(st) ) ) {
  // 	//triangle is disjoint
  // 	//	    smTriangle t1;
  // 	//	    t1.n1 = j;
  // 	//	    t1.n2 = s;
  // 	//	    t1.n3 = t;
  // 	//	    S.insert( pair< smEdge, smTriangle >( js, t1 ) );
  // 	//	    S.insert( pair< smEdge, smTriangle >( jt, t1 ) );
  // 	//	    S.insert( pair< smEdge, smTriangle >( st, t1 ) );
  // 	S.insert( js );
  // 	S.insert( jt );
  // 	S.insert( st );
  // 	return true;
  //     }

  //     return false;
  //   }
      
  //   /*
  //    * DART-BASE (no triangle listing)
  //    */
  //   unsigned dart_base_free() {
  //     clock_t t_start = clock();
  //     vector < vector< node_id > > A( n, vector< node_id >() );
  //     unsigned countS = 0;
	 
  //     for (node_id s = 0; s < n; ++s ) {
  // 	for (size_t j = 0; j < adjList[s].size(); ++j) {
  // 	  node_id& t = adjList[s][ j ];
  // 	  if (s < t) {
  // 	    vector< node_id >& A_s = A[s];
  // 	    vector< node_id >& A_t = A[t];
  // 	    auto it1 = A_s.begin();
  // 	    auto it2 = A_t.begin();
  // 	    if (it1 == A_s.end() || it2 == A_t.end() ) {
  // 	      A_t.push_back( s );
  // 	      continue;
  // 	    }
  // 	    while (1) {
  // 	      if (*it1 < *it2) {
  // 		++it1;
  // 		if (it1 == A_s.end()) {
  // 		  break;
  // 		}
  // 	      } else {
  // 		if (*it2 < *it1) {
  // 		  ++it2;
  // 		  if (it2 == A_t.end()) {
  // 		    break;
  // 		  }
  // 		} else {
  // 		  //found a triangle
  // 		  if (free_triangle( *it1, s, t)) {
  // 		    countS += 3;
			      
  // 		  }
			   
  // 		  ++it1; 
  // 		  ++it2; 
  // 		  if (it1 == A_s.end() || it2 == A_t.end() )
  // 		    break;
  // 		}
  // 	      }
  // 	    }
  // 	    A_t.push_back( s );
  // 	  }
  // 	}
  //     }

  //     runningTime = double (clock() - t_start) / CLOCKS_PER_SEC;
  //     logg(INFO, "DART_BASE_FREE finished, size of S: " + to_string(countS));

  //     return countS;
  //   }

  //   /*
  //    * DART-BASE (no triangle listing)
  //    */
  //   unsigned dart_base_free_smaller() {
  //     //	 vector < vector< node_id > > A( n, vector< node_id >() );
  //     unsigned countS = 0;
	 
  //     for (node_id s = 0; s < n; ++s ) {
  // 	for (size_t j = 0; j < adjList[s].size(); ++j) {
  // 	  node_id& t = adjList[s][ j ];
  // 	  if (s < t) {
  // 	    vector< node_id >& A_s = adjList[s];
  // 	    vector< node_id >& A_t = adjList[t];
  // 	    auto it1 = A_s.begin();
  // 	    auto it2 = A_t.begin();
  // 	    if (it1 == A_s.end() || it2 == A_t.end() ) {
  // 	      continue;
  // 	    }
  // 	    while (1) {
  // 	      if (*it1 < *it2) {
  // 		++it1;
  // 		if (it1 == A_s.end()) {
  // 		  break;
  // 		}
  // 	      } else {
  // 		if (*it2 < *it1) {
  // 		  ++it2;
  // 		  if (it2 == A_t.end()) {
  // 		    break;
  // 		  }
  // 		} else {
  // 		  //found a triangle
  // 		  if (free_triangle( *it1, s, t)) {
  // 		    countS += 3;
			      
  // 		  }
			   
  // 		  ++it1; 
  // 		  ++it2; 
  // 		  if (it1 == A_s.end() || it2 == A_t.end() )
  // 		    break;
  // 		}
  // 	      }
  // 	    }
  // 	  }
  // 	}
  //     }
  //     logg(INFO, "DART_BASE_FREE finished, size of S: " + to_string(countS));

  //     return countS;
  //   }

  //   void list_triangles() {
  //     vector < vector< node_id > > A( n, vector< node_id >() );
  //     vector < vector< size_t > > I( n, vector< size_t >() ); //Indices of edges
  //     size_t count = 0;
	 
  //     for (node_id s = 0; s < n; ++s ) {
  // 	for (size_t j = 0; j < adjList[s].size(); ++j) {
  // 	  node_id& t = adjList[s][ j ];
  // 	  if (s < t) {
  // 	    vector< node_id >& A_s = A[s];
  // 	    vector< node_id >& A_t = A[t];
  // 	    vector< size_t >& I_s = I[s];
  // 	    vector< size_t >& I_t = I[t];
  // 	    auto it1 = A_s.begin();
  // 	    auto it2 = A_t.begin();
  // 	    if (it1 == A_s.end() || it2 == A_t.end() ) {
  // 	      A_t.push_back( s );
  // 	      I_t.push_back( j );
  // 	      continue;
  // 	    }
  // 	    auto it3 = I_s.begin();
  // 	    auto it4 = I_t.begin();
  // 	    while (1) {
  // 	      if (*it1 < *it2) {
  // 		++it1;
  // 		++it3;
  // 		if (it1 == A_s.end()) {
  // 		  break;
  // 		}
  // 	      } else {
  // 		if (*it2 < *it1) {
  // 		  ++it2;
  // 		  ++it4;
  // 		  if (it2 == A_t.end()) {
  // 		    break;
  // 		  }
  // 		} else {
  // 		  //found a triangle
  // 		  //add_triangle( *it1, *it3, *it4, j, s, t);
  // 		  //			   smTriangle tt( *it1, s, t );
  // 		  //			   T.push_back(tt);
  // 		  ++count;
  // 		  ++it1; ++it3;
  // 		  ++it2; ++it4;
  // 		  if (it1 == A_s.end() || it2 == A_t.end() )
  // 		    break;
  // 		}
  // 	      }
  // 	    }
  // 	    A_t.push_back( s );
  // 	    I_t.push_back( j );
  // 	  }
  // 	}
  //     }
  //     logg( INFO, to_string(count) + " triangles found.");
	       
  //   }

  // };
  

  //   uint32_t three = 3;
  //   uint32_t one = 1;
  uint32_t bitMask = ~(3 << 30);
  uint32_t bitS = (1 << 31);
  uint32_t bitW = (1 << 30);
   
  class tinyEdge {
  public:
    //last 30 bits are target node_id
    //first bit is inS, second bit is inW
    uint32_t target;
    //    uint32_t matePairLoc; //Location of mate pair edge
     
    node_id getId() const {
      return target & bitMask;
    }

    bool inS() {
      return (target >> 31);
    }

    bool inW() {
      return (target >> 30) & 1; //
    }

    void setS() {
      target = target | bitS;
    }

    void setW() {
      target = target | bitW;
    }

    void unsetS() {
      target = target & (~bitS);
    }

    void unsetW() {
      target = target & (~bitW);
    }

    tinyEdge() {
      target = 0;
    }

    //tinyEdge( const node_id& nid ) {
    //	 target = nid; //inS = inW = false;
    //      }
    tinyEdge( node_id nid ) {
      target = nid; //inS = inW = false;
    }

    tinyEdge( node_id nid, uint32_t mp ) {
      target = nid;
      //      matePairLoc = mp;
    }
     
    tinyEdge( const tinyEdge& rhs ) {
      target = rhs.target;
      //      matePairLoc = rhs.matePairLoc;
    }
  };

  /*
   * Only works if inS and inW are 0
   * Faster than operator<
   */
  bool tinyEdgeCompare( const tinyEdge& a, const tinyEdge& b ) {
    return a.target < b.target;
  }

  /*
   * Works regardless of status of front bits
   */
  bool operator<( const tinyEdge& a, const tinyEdge& b ) {
    return a.getId() < b.getId();
  }

  /*
   * n1, the third vertex of this triangle,
   * will be the vertex storing this triangle
   */
  class tinyTriangle {
  public:
    uint32_t n2;
    uint32_t n3;

    tinyTriangle() {}
    
    tinyTriangle( uint32_t n2_in, uint32_t n3_in ) {
      n2 = n2_in;
      n3 = n3_in;
    }

    bool contains( uint32_t& test, uint32_t& other ) {
      if ( n2 == test ) {
	other = n3;
	return true;
      }
      if (n3 == test ) {
	other = n2;
	return true;
      }

      return false;
    }

    tinyTriangle& operator= (const tinyTriangle& rhs) {
      n2 = rhs.n2;
      n3 = rhs.n3;

      return *this;
    }
  };
   
  //Node class
  class tinyNode {
  public:
    vector< tinyEdge > neis;
    vector< tinyTriangle > solutionTriangles;

    tinyNode () { }
     
    tinyNode ( const tinyNode& rhs ) {
      neis.assign( rhs.neis.begin(), rhs.neis.end() );
      solutionTriangles.assign( rhs.solutionTriangles.begin(), rhs.solutionTriangles.end() );
    }
  };
   
  class tinyGraph {
  public:
    vector< tinyNode > adjList;
    unsigned n;
    unsigned m;
    Logger logg;
    vector< smTriangle > T_sol;
    double preprocessTime;
    
    void init_empty_graph() {
      tinyNode emptyNode;
      adjList.assign(n, emptyNode);
    }

    void add_edge_immediate( unsigned from, unsigned to ) {
      tinyEdge FT( to, adjList[to].neis.size() );
      tinyEdge TF( from, adjList[from].neis.size() );
      adjList[ from ].neis.push_back( FT );
      adjList[ to ].neis.push_back( TF );
    }

    bool add_edge_half( node_id from, node_id to, vector< tinyEdge >::iterator& edgeAdded ) {
      if (from == to)
	return false;
      
      vector< tinyEdge >& v1 = adjList[ from ].neis;

      auto it = v1.begin();
      while (it != v1.end()) {
	if (it->getId() >= to)
	  break;

	++it;
      }

      tinyEdge newEdge;
      newEdge.target = to;
      
      if (it != v1.end()) {
	if (it->getId() == to)
	  return false;
	//The element should be inserted
	edgeAdded = v1.insert( it, newEdge ); //O( max_deg )
	return true;
      } 

      edgeAdded = v1.insert( it, newEdge );
      return true;
    }

    bool remove_edge_half( node_id from, node_id to ) {
      if (from == to)
	return false;
      
      vector< tinyEdge >& v1 = adjList[ from ].neis;

      auto it = v1.begin();
      while (it != v1.end()) {
	if (it->getId() >= to)
	  break;

	++it;
      }

      if (it != v1.end()) {
	if (it->getId() == to) {
	  v1.erase( it );
	  return true;
	}
      } 

      return false;
    }

    /*
     * Removes self-loops and multi-edges
     * Assumes sorted adjacency list in each tinyNode
     */
    void simplify() {
      for (unsigned i = 0; i < n; ++i) {
	bool b_continue;
	do {
	  b_continue = false;
	  auto it1 = adjList[i].neis.begin();
	  auto it2 = it1 + 1;
	  while (it2 != adjList[i].neis.end()) {
	    if (it1->getId() == it2->getId()) {
	      //remove multi-edge
	      b_continue = true;
	      adjList[i].neis.erase( it1 );
	      break;
	    }
	    if (it1->getId() == i ) {
	      //remove loop
	      adjList[i].neis.erase( it1 );
	      b_continue = true;
	      break;
	    }
	    if (it2->getId() == i ) {
	      //remove loop
	      adjList[i].neis.erase( it2 );
	      b_continue = true;
	      break;
	    }
		 
	    ++it2; ++it1;
	  }
	} while (b_continue);
      }
    }

    void read_adj_list_bin( string fname ) {
      ifstream ifile ( fname.c_str(), ios::in | ios::binary );
      unsigned n;
      ifile.read( (char*) &n, sizeof( unsigned ) );

      this->n = n;

      init_empty_graph();
      unsigned ss;
      tinyEdge temp;
      unsigned nei_id;
      clock_t t_start = clock();
      for ( unsigned i = 0; i < n; ++i ) {

	ifile.read( (char*) &ss, sizeof( unsigned ) );

	adjList[i].neis.assign( ss, temp );
	for (unsigned j = 0; j < ss; ++j) {
	  ifile.read( (char*) &nei_id, sizeof( unsigned ) );
	  adjList[i].neis[j].target = nei_id;
	}
      }
      preprocessTime = double (clock() - t_start) / CLOCKS_PER_SEC;
    }

    /*
     * Returns pre-processing time for graph
     */
    double read_edge_list_bin( string fname ) {
      ifstream ifile ( fname.c_str(), ios::in | ios::binary );
      unsigned m; //number of edges in list
      unsigned n;
      ifile.read( (char*) &n, sizeof( unsigned ) );
      ifile.read( (char*) &m, sizeof( unsigned ) );
      this->n = n;
      this->m = m;
      init_empty_graph();

      unsigned from,to;
      for (unsigned i = 0; i < m; ++i) {
	ifile.read( (char *) &from,  sizeof( unsigned ) );
	ifile.read( (char *) &to,  sizeof( unsigned ) );

	add_edge_immediate( from, to );
      }
      ifile.close();
	 
      //      logg(INFO, "Sorting neighbor lists..." );
      clock_t t_start = clock();
      for (unsigned i = 0; i < n; ++i) {
	//	    adjList[i].sort( tinyEdgeCompare );
	sort( adjList[i].neis.begin(), adjList[i].neis.end(), tinyEdgeCompare );
	//update location of mate pairs
	//	for (unsigned j = 0; j < adjList[i].neis.size(); ++j) {
	//	  uint32_t& target = adjList[i].neis[ j ].target;
	//	  uint32_t& mp = adjList[i].neis[ j ].matePairLoc;
	//	  (adjList[ target  ].neis[ mp ]).matePairLoc = j;
	//}
      }
      preprocessTime = double (clock() - t_start) / CLOCKS_PER_SEC;
      return double (clock() - t_start) / CLOCKS_PER_SEC;
    }

    void setSInList( vector< tinyEdge >& l, node_id& v ) {
      auto it = l.begin();
	 
      while (it->getId() != v)
	++it;

      it->setS();
    }

    void unpruneInList( vector< tinyEdge >& l, node_id& v ) {
      auto it = l.begin();
	 
      while (it->getId() != v)
	++it;

      it->setS();
      it->unsetW();
    }

    vector< tinyEdge >::iterator clearInList( vector< tinyEdge >& l, node_id& v ) {
      auto it = l.begin();
	 
      while (it->getId() != v)
	++it;

      it->unsetS();
      it->unsetW();

      return it;
    }
    
    void pruneSInList( vector< tinyEdge >& l, node_id& v ) {
      auto it = l.begin();
	 
      while (it->getId() != v)
	++it;

      it->unsetS();
      it->setW();
    }
      
    bool free_triangle( tinyEdge& st,
			tinyEdge& sv,
			tinyEdge& tv,
			node_id& s ) {
      if (sv.inS() || tv.inS() )
	return false;
	 
      //triangle is disjoint
      node_id& t = st.target;
      node_id& v = sv.target;
      vector< tinyEdge >& At = adjList[ t ].neis;
      //	 vector< tinyEdge >& As = adjList[ s ].neis;
      vector< tinyEdge >& Av = adjList[ v ].neis;
      setSInList( At, s );
      setSInList( Av, s );
      setSInList( Av, t );
      //      At[ st.matePairLoc ].setS();
      //      Av[ sv.matePairLoc ].setS();
      //      Av[ tv.matePairLoc ].setS();
	 
      st.setS();
      sv.setS();
      tv.setS();

      return true;
    }

    bool prune_triangle( tinyEdge& st,
			 tinyEdge& sv,
			 tinyEdge& tv,
			 node_id& s ) {
      if (sv.inS() || tv.inS() )
	return true;
      else
	return false;
    }

    bool find_disjoint_triangle( uint32_t s,
				 tinyEdge& st,
				 vector< tinyEdge >::iterator& sv,
				 vector< tinyEdge >::iterator& tv,
				 vector< smEdge >& unprunedEdges) {
      if (st.inS())
	return false;

      vector< tinyEdge >& A_s = adjList[s].neis;
      vector< tinyEdge >& A_t = adjList[ st.getId() ].neis; 
      auto it1 = A_s.begin();
      auto it2 = A_t.begin();
      if (it1 == A_s.end() || it2 == A_t.end() ) {
	return false;
      }
      while (1) {
	if (*it1 < *it2) {
	  ++it1;
	  if (it1 == A_s.end()) {
	    break;
	  }
	} else {
	  if (*it2 < *it1) {
	    ++it2;
	    if (it2 == A_t.end()) {
	      break;
	    }
	  } else {
	    //found a triangle
	    //st, *it1 = sv, *it2 = tv
	    if ((*it1).inW()) {
	      //unprune this edge
	      (*it1).unsetW();
	      unpruneInList( adjList[ (*it1).target ].neis, s );
	      smEdge e_tmp;
	      e_tmp.from = s;
	      e_tmp.to = (*it1).target;
	      unprunedEdges.push_back( e_tmp );		     
	      (*it1).setS();
	    } else {
	      if ((*it2).inW()) {
		//unprune this edge
		(*it2).unsetW();
		node_id t = st.getId();
		unpruneInList( adjList[ (*it2).target ].neis, t );
		smEdge e_tmp;
		e_tmp.from = t;
		e_tmp.to = (*it2).target;
		unprunedEdges.push_back( e_tmp );		     
		(*it2).setS();
	      } else {
		if ( !( (*it1).inS() || (*it2).inS() ) ) {
		  //this triangle is disjoint from S
		  sv = it1;
		  tv = it2;
		  return true;
		}
	      }
	    }
			
	    ++it1; 
	    ++it2; 
	    if (it1 == A_s.end() || it2 == A_t.end() )
	      break;
	  }
	}
      }
	 
      return false;
    }

    tinyEdge* findEdgeInList( vector< tinyEdge >& v1, node_id target ) {
      for (unsigned i = 0; i < v1.size(); ++i) {
	if (v1[i].getId() == target)
	  return &(v1[i]);
      }

      return 0; //Edge not found
    }

    vector< tinyEdge >::iterator findEdgeInList( node_id source, node_id target ) {
      vector< tinyEdge >& v1 = adjList[source].neis;
      for (auto it = v1.begin();
	   it != v1.end();
	   ++it ) {
	if (it->getId() == target)
	  return it;
      }

      return v1.end(); //Edge not found
    }
    
    void prune( smEdge& ee ) {
      node_id& s = ee.from;
      node_id& t = ee.to;
      tinyEdge& st = *(findEdgeInList( adjList[ s ].neis, t ));
      if (!st.inS())
	return;
      bool prunable = true;
      vector< tinyEdge >& A_s = adjList[s].neis;
      vector< tinyEdge >& A_t = adjList[t].neis; 
      auto it1 = A_s.begin();
      auto it2 = A_t.begin();
      if (it1 == A_s.end() || it2 == A_t.end() ) {
	return;
      }
      while (1) {
	if (*it1 < *it2) {
	  ++it1;
	  if (it1 == A_s.end()) {
	    break;
	  }
	} else {
	  if (*it2 < *it1) {
	    ++it2;
	    if (it2 == A_t.end()) {
	      break;
	    }
	  } else {
	    //found a triangle
	    if (prune_triangle( st, *it1, *it2, s)) {

	    } else {
	      //st must remain in S
	      prunable = false;
	      break;
	    }
			
	    ++it1; 
	    ++it2; 
	    if (it1 == A_s.end() || it2 == A_t.end() )
	      break;
	  }
	}
      }
      //prune if we can
      if (prunable) {
	st.unsetS();
	st.setW();
	pruneSInList( A_t, s );
	//	    A_t[ st.matePairLoc ].unsetS();
	//	    A_t[ st.matePairLoc ].setW();
      }
      
    }
    
    void prune( vector< smEdge >& unprunedEdges ) {
      for (auto it = unprunedEdges.begin();
	   it != unprunedEdges.end();
	   ++it ) {
	prune( *it );
      }
    }

    bool incident( node_id& s, node_id& t ) {
      if ( findEdgeInList( adjList[ s ].neis, t ) != 0 )
	return true;
      else
	return false;
    }
    
    double dart_add_edge( node_id s, vector< tinyEdge >::iterator st ) {
      clock_t t_start = clock();
      
      vector<tinyEdge>::iterator sv; //want references so edges in graph can be modified
      vector<tinyEdge>::iterator tv;

      vector< smEdge > unprunedEdges;
      if ( find_disjoint_triangle( s, *st, sv, tv, unprunedEdges )) {
	//Add this triangle to S
	node_id& t = st->target;
	node_id& v = sv->target;

	smEdge e_tmp;
	e_tmp.from = s;
	e_tmp.to = t;
	unprunedEdges.push_back( e_tmp );
	e_tmp.to = v;
	unprunedEdges.push_back( e_tmp );
	e_tmp.from = t;
	unprunedEdges.push_back( e_tmp );
	
	tinyTriangle sT( t, v );
	tinyTriangle tT( s, v );
	tinyTriangle vT( s, t );
	adjList[ s ].solutionTriangles.push_back (sT);
	adjList[ t ].solutionTriangles.push_back (tT);
	adjList[ v ].solutionTriangles.push_back (vT);
	    
	vector< tinyEdge >& At = adjList[ t ].neis;
	vector< tinyEdge >& Av = adjList[ v ].neis;

	setSInList( At, s );
	setSInList( Av, s );
	setSInList( Av, t );
	    
	st->setS();
	sv->setS();
	tv->setS();

    
      }

      prune( unprunedEdges );
      return double (clock() - t_start) / CLOCKS_PER_SEC;
    }

    bool remove_edge( node_id s, node_id t ) {
      if ( remove_edge_half( s , t ) ) {
	--m;
	return remove_edge_half(t, s );
      }

      return false;
    }

    node_id eraseSolTriangle( node_id& s, node_id& t ) {
      node_id v;

      for (auto it = adjList[s].solutionTriangles.begin();
	   it != adjList[s].solutionTriangles.end();
	   ++it) {
	if (it->contains( t, v )) {
	  adjList[s].solutionTriangles.erase( it );
	  return v;
	  break;
	}
      }
      
      return bitS; //As failure value, no node id can be that big
    }
    
    double dart_remove_edge( node_id s, vector< tinyEdge >::iterator st ) {
      clock_t t_start = clock();
      node_id t = st->getId();
      if (st->inS() || st->inW() ) {
	//remove the solution triangles
	node_id v = eraseSolTriangle( s, t );
	eraseSolTriangle( t, s );
	eraseSolTriangle( s, v );
	//remove (s,t)
	remove_edge( s, t );

	//Remove sv, tv from solution
	vector< tinyEdge >::iterator sv = clearInList( adjList[ s ].neis, v );
	clearInList( adjList[ v ].neis, s );
	vector< tinyEdge >::iterator tv = clearInList( adjList[ t ].neis, v );
	clearInList( adjList[ v ].neis, t );

	dart_add_edge( s, sv );
	dart_add_edge( t, tv );
	
      } else {
	remove_edge( s, t );
      }

      return double (clock() - t_start) / CLOCKS_PER_SEC;
    }
      
    /*
     * DART-BASE (no triangle listing)
     */
    void dart_base() {
      for (node_id s = 0; s < n; ++s ) {
	for (auto it0 = adjList[s].neis.begin();
	     it0 != adjList[s].neis.end();
	     ++it0 ) {
	  tinyEdge& st = *it0;
	  if (st.inS())
	    continue;
	  if (!(s < st.target))
	    continue;
	  vector< tinyEdge >& A_s = adjList[s].neis;
	  vector< tinyEdge >& A_t = adjList[ st.target].neis; //know not in S or W
	  auto it1 = A_s.begin();
	  auto it2 = A_t.begin();
	  if (it1 == A_s.end() || it2 == A_t.end() ) {
	    continue;
	  }
	  while (1) {
	    if (*it1 < *it2) {
	      ++it1;
	      if (it1 == A_s.end()) {
		break;
	      }
	    } else {
	      if (*it2 < *it1) {
		++it2;
		if (it2 == A_t.end()) {
		  break;
		}
	      } else {
		//found a triangle
		if (free_triangle( st, *it1, *it2, s)) {
		  unsigned t = st.getId();
		  unsigned v = (*it1).getId();
			      
		  //smTriangle tt( s, st.getId(), (*it1).getId() );
		  //   T_sol.push_back( tt );
		  tinyTriangle sT( t, v );
		  tinyTriangle tT( s, v );
		  tinyTriangle vT( s, t );
		  adjList[ s ].solutionTriangles.push_back (sT);
		  adjList[ t ].solutionTriangles.push_back (tT);
		  adjList[ v ].solutionTriangles.push_back (vT);
			      
		  break;
		}
			
		++it1; 
		++it2; 
		if (it1 == A_s.end() || it2 == A_t.end() )
		  break;
	      }
	    }
	  }
	}
      }

      free_prune();
    }

    void free_prune() {
      bool prunable;
      for (node_id s = 0; s < n; ++s ) {
	for (auto it0 = adjList[s].neis.begin();
	     it0 != adjList[s].neis.end();
	     ++it0 ) {
	  tinyEdge& st = *it0;
	  node_id t = st.getId();
	  if (!(s < t))
	    continue;
	  if (!st.inS())
	    continue;
	  prunable = true;
	  vector< tinyEdge >& A_s = adjList[s].neis;
	  vector< tinyEdge >& A_t = adjList[t].neis; 
	  auto it1 = A_s.begin();
	  auto it2 = A_t.begin();
	  if (it1 == A_s.end() || it2 == A_t.end() ) {
	    continue;
	  }
	  while (1) {
	    if (*it1 < *it2) {
	      ++it1;
	      if (it1 == A_s.end()) {
		break;
	      }
	    } else {
	      if (*it2 < *it1) {
		++it2;
		if (it2 == A_t.end()) {
		  break;
		}
	      } else {
		//found a triangle
		if (prune_triangle( st, *it1, *it2, s)) {

		} else {
		  //st must remain in S
		  prunable = false;
		  break;
		}
			
		++it1; 
		++it2; 
		if (it1 == A_s.end() || it2 == A_t.end() )
		  break;
	      }
	    }
	  }
	  //prune if we can
	  if (prunable) {
	    st.unsetS();
	    st.setW();
	    pruneSInList( A_t, s );
	    //	    A_t[ st.matePairLoc ].unsetS();
	    //	    A_t[ st.matePairLoc ].setW();
	  }
	}
      }
    }

    unsigned countS() {
      unsigned count = 0;
      for (unsigned s = 0; s < adjList.size(); ++s) {
	for (auto it2 = adjList[s].neis.begin();
	     it2 != (adjList[s]).neis.end();
	     ++it2 ) {
	  node_id t = (*it2).getId();
	  if (s < t) {
	    if ((*it2).inS()) {
	      ++count;
	    }
	  }
	}
      }
      return count;
    }
      
    // bool check_validity() {
    // 	 unsigned vertInCommon;
    // 	 for (auto it1 = T_sol.begin();
    // 	      it1 != T_sol.end();
    // 	      ++it1 ) {
    // 	    for (auto it2 = it1;
    // 		 it2 != T_sol.end();
    // 		 ++it2 ) {
    // 	       vertInCommon = 0;
    // 	       smTriangle& t1 = *it1;
    // 	       smTriangle& t2 = *it2;
    // 	       if (t1.n1 == t2.n1)
    // 		  ++vertInCommon;
    // 	       if (t1.n1 == t2.n2)
    // 		  ++vertInCommon;
    // 	       if (t1.n1 == t2.n3)
    // 		  ++vertInCommon;

    // 	       if (t1.n2 == t2.n1)
    // 		  ++vertInCommon;
    // 	       if (t1.n2 == t2.n2)
    // 		  ++vertInCommon;
    // 	       if (t1.n2 == t2.n3)
    // 		  ++vertInCommon;

    // 	       if (t1.n3 == t2.n1)
    // 		  ++vertInCommon;
    // 	       if (t1.n3 == t2.n2)
    // 		  ++vertInCommon;
    // 	       if (t1.n3 == t2.n3)
    // 		  ++vertInCommon;

    // 	       if (vertInCommon > 1) {
    // 		  logg(ERROR, "Solution triangles share an edge.");
    // 		  cerr << t1.n1 << ' ' << t1.n2 << ' ' << t1.n3 << endl;
    // 		  cerr << t2.n1 << ' ' << t2.n2 << ' ' << t2.n3 << endl;
		  
    // 		  return false;
    // 	       }
    // 	    }
    // 	 }
    // 	 return true;
    // }
      
  };

   
  void worker_list_triangles(Graph& G,
			     pedge start, pedge offend,
			     list< Triangle >& worker_triangles,
			     size_t worker_m ) {

    map< node_id, pedge > MFrom;
    Node From, To;
    size_t progress = 0;
    size_t percent = worker_m / 100;
    //      cout << "\r                       \r0% done";
    for (list< Edge >::iterator it1 = start;
	 it1 != offend;
	 ++it1 ) {
      if (percent != 0) {
	if (progress % percent == 0) {
	  cout << "\r                       \r" << progress * 100 / worker_m << "% done";
	  cout.flush();
	}
      }
      ++progress;
      From = G.V[ it1->from ];
      To = G.V[ it1->to ];
      //Go through From's edges and remember their other vertex
      MFrom.clear();
      for (list<pedge>::iterator it2 = From.neighbors.begin();
	   it2 != From.neighbors.end();
	   ++it2 ) {
	//we only want to look at greater edges to avoid
	//recomputing triangles
	if ( (*it2)->eid > (*it1).eid ) { 
	  node_id other = (*it2)->other( it1->from );
	  MFrom[ other ] = *it2; //this is an iterator of E
	}
      }
      //now go through to's edges
      for (list<pedge>::iterator it2 = To.neighbors.begin();
	   it2 != To.neighbors.end();
	   ++it2 ) {
	if ( (*it2)->eid > (*it1).eid ) { 
	  node_id other = (*it2)->other( it1->to );
	  map< node_id, pedge>::iterator it3 = MFrom.find( other );
	  if (it3 != MFrom.end()) {
	    //we have found a triangle
	    // logg( DEBUG, "Found triangle: ("
	    // 	    + to_string( it1->from ) + ","
	    // 	    + to_string( it1->to ) + ");"
	    // 	    + "(" + to_string( (*it2)->from ) + ","
	    // 	    + to_string( (*it2)->to ) + ");"
	    // 	    + "(" + to_string( it3->second->from ) + ","
	    // 	    + to_string( it3->second->to ) + ");" );
	    Triangle t1( it1, *it2, it3->second );
	    worker_triangles.push_front(t1);
	    //	      list< Triangle >::iterator itt = T.begin();
	    //	      (it1)->Delta.push_back( itt);
	    //	      (*it2)->Delta.push_back( itt );
	    //	      (it3->second)->Delta.push_back( itt );
	  }
	}
      }
    }

    //	logg( INFO, to_string(T.size()) + " triangles found.");
    cout << "\r                       \r100% done" << endl;	
  }

   
  void list_triangles_multi( Graph& G,
			     unsigned nThreads = 1 ) {
    pedge start = G.E.begin();
    pedge end = start;
    size_t m_thread = G.E.size() / nThreads;
    list<Triangle> empty_list;
    vector< list< Triangle > > worker_results( nThreads, empty_list );
    G.logg(DEBUG, "m_thread: " + to_string(m_thread));
    thread* worker_threads = new thread[ nThreads ];
    for (unsigned i = 0; i < nThreads; ++i) {
      start = end;
      if (i == (nThreads - 1)) {
	worker_threads[i] = thread( worker_list_triangles,
				    ref( G ),
				    start, G.E.end(),
				    ref ( worker_results[ i ] ),
				    m_thread );
      } else {
	for (size_t i = 0; i < m_thread; ++i)
	  ++end;
	worker_threads[i] = thread( worker_list_triangles,
				    ref( G ),
				    start, end,
				    ref ( worker_results[ i ] ),
				    m_thread );
      }
    }
	
    for (unsigned i = 0; i < nThreads; ++i ) {
      worker_threads[i].join();
    }

    delete [] worker_threads;

    //splice the triangles into G's Triangle list
    for (unsigned i = 0; i < nThreads; ++i ) {
      G.T.splice( G.T.end(), worker_results[ i ] );
    }

    //Build Delta lists:
    for (auto it = G.T.begin();
	 it != G.T.end();
	 ++it) {
      it->e1->Delta.push_back( it );
      it->e2->Delta.push_back( it );
      it->e3->Delta.push_back( it );
    }
  }

   class resultsHandler {
   public:
      map< string, vector< double > > data;
      
      void add( string name, double val ) {
	 data[ name ].push_back( val );
      }

      void set( string name, double val ) {
	 data[ name ].assign(1, val);
      }
      
      void print( ostream& os, bool printStdDev = false ) {
	 //Print names
	 os << '#';
	 unsigned index = 1;
	 for (auto it = data.begin();
	      it != data.end();
	      ++it ) {
	    os << setw(20);
		 
	    if (it->second.size() == 1) {
	       os << to_string( index ) + it->first;
	       ++index;
	    } else {
	       if (printStdDev) {
		  os << (to_string(index) + it->first + "_mean");
		  ++index;
		  os << (to_string(index) + it->first + "_stddev");
		  ++index;
	       } else {
		  os << to_string( index ) + it->first;
		  ++index;
	       }
	    }
	 }
	 os << endl;
	 os << fixed;
	 double mean;
	 double stddev;
	 for (auto it = data.begin();
	      it != data.end();
	      ++it ) {
	    //	    os << setprecision( 3 );
	    if ((it->second).size() > 1) {
	       //compute mean
	       mean = 0;
	       for (unsigned i = 0; i < (it->second).size(); ++i) {
		  mean += (it->second)[i];
	       }
	       mean /= (it->second).size();

	       os << setw(20) << mean;
	       if (printStdDev) {
		  //compute stddev
		  stddev = 0;
		  for (unsigned i = 0; i < (it->second).size(); ++i) {
		     stddev += ((it->second)[i] - mean) * ((it->second)[i] - mean);
		  }
		  stddev /= ((it->second).size() - 1);
		  stddev = sqrt( stddev );
		  os << setw(20) << stddev;
	       }
	    } else {
	       os << setw(20) << (it->second).front();
	    }
	 }
	 os << endl;
      }
   };
   
}



#endif
