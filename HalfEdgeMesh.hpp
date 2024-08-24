#ifndef HH_INDEXED_HALF_EDGE_MESH_DATA_STRUCTURE
#define HH_INDEXED_HALF_EDGE_MESH_DATA_STRUCTURE


#include "Mesh.hpp"
//#include "utilities.hpp"

#include <vector>
#include <unordered_set>
#include <unordered_map>

#include <set>
//#include <stack>

#include <stdexcept>
#include <fstream>

#define HE_INVALID_INDEX -1
#define HE_ISOLATED_INDEX -2	// in _vert_to_hedge it indicates that the vertex is isolated

// ====================================== INFO =============================================== //

/*
	* Index based half edge data structure
	* A half edge stores only indexes for the head/target vertex, next half edge in the loop and begin half edge of the loop
	* Twin/Opposite half edges are next to each other in memory, so if a half edge has index i, its edge id is i/2
	* There is garbage collection, so when you remove a half edge or a vertex, it is only marked as removed
	* Stores face and border loops, but the faces and borders doesnt have indexes, the index they have is the index of the half edge that indicates the beginning of the loop
	* Most methods work with vertex indexes and half edge indexes, even if the method does something
	  on an edge, face or border, this is because for most face and border operations you need any half edge
	  that lies in the same loop and for the edge operations i dont want the overhead of constantly
	  converting from half edge index to edge index
	* _vert_to_hedge or id_hedge(vert_id) gives you the outgoing half-edge of one of its adjacent edges, but if it has no adjacent edges, then HE_ISOLATED_INDEX is returned
	* A vertex is marked as removed when _vert_to_hedge[vert_id] == HE_INVALID_INDEX and its index vi is in _garbage_vertices
	* A half edge is marked as removed when its head_id, next_id and begin_id are set to HE_INVALID_INDEX and its index is in _garbage_edges
*/


// ==================================== REFERENCES =========================================== //

/*
	http://15462.courses.cs.cmu.edu/fall2016/article/14
	http://15462.courses.cs.cmu.edu/fall2016/article/15
	https://www.reddit.com/r/rust/comments/t0ja80/tips_for_writing_a_halfedge_mesh_data_structure/
	https://cs184.eecs.berkeley.edu/sp20/article/17/an-introduction-to-half-idEdge-dat
	https://cs184.eecs.berkeley.edu/sp20/lecture/8/meshes-and-geometry-processing

	https://observablehq.com/@esperanc/half-idEdge-data-structure
	https://cs418.cs.illinois.edu/website/text/halfedge.html
	Grid Subdivision - https://cs418.cs.illinois.edu/website/text/make-geom.html

	https://student.cs.uwaterloo.ca/~cs779/Gallery/Winter2018/anietoro/doc/
	https://www.reddit.com/r/Cplusplus/comments/quf16o/whats_the_most_c_way_of_a_triangle_tetrahedral/
	https://www.reddit.com/r/rust/comments/t0ja80/tips_for_writing_a_halfedge_mesh_data_structure/
	https://ir.library.louisville.edu/cgi/viewcontent.cgi?article=2093&context=etd
	https://ohiostate.pressbooks.pub/app/uploads/sites/45/2020/01/GM-Euler.pdf
	https://graphics.stanford.edu/courses/cs268HM_INVALID_INDEX3-spring/Notes/halfedge.pdf

	https://graphics.stanford.edu/courses/cs368-00-spring/TA/manuals/CGAL/ref-manual2/Topological_map/Chapter_tpm.html
	https://github.com/strncat/subdivision-algorithms
	https://student.cs.uwaterloo.ca/~cs779/Gallery/Winter2020/zrlu/

	https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/model/euler-op.html - Euler Operators
	https://digitalfirst.bfwpub.com/math_applet/asset/1/ - Eulerizing a Graph
	https://pages.mtu.edu/~shene/COURSES/cs3621/LAB/surface/swept.html - Sweeping Surfaces

	https://www.merlin.uzh.ch/contributionDocument/download/14550 - Edge Collapse and Vertex Split
	http://elib.mi.sanu.ac.rs/files/journals/publ/86/n080p023.pdf - Topology Preserving Edge Contraction

*/



struct HalfEdge {
	int head_id;	// head/target/to vertex index
	int next_id;	// next half edge index
	int begin_id;	// the index of the half edge at the beggining of the face or border loop of the current half edge

	HalfEdge(const int& head_id = HE_INVALID_INDEX, const int& next_id = HE_INVALID_INDEX
		, const int& begin_id = HE_INVALID_INDEX) : head_id(head_id), next_id(next_id), begin_id(begin_id) {}

};


template <typename V>
class HEMesh {
protected:

	std::vector<V> _vertices{};
	std::vector<int> _vert_to_hedge{};	// gives the outgoint half edge of the vertex, so id_tail(_vert_to_hedge[vi]) == vi
	std::vector<HalfEdge> _hedges{};

	std::unordered_set<int> _faces{};   // the indexes of the half edges that are the beginnig of a face loop
	std::unordered_set<int> _borders{}; // the indexes of the half edges that are the beginnig of a border loop

	// Garbage collection/recycling

	std::set<int> _garbage_vertices{};	// stores vertex ids that are removed

	// stores even half-edge ids, an id in the set indicates that both the half-edge
	//  and its twin/opposite half-edge are tagged/marked as removed
	std::set<int> _garbage_edges{};

	//bool _are_verts_shrinked;
	//bool _are_hedges_shrinked;

	// frees the memory space occupied by a vertex that was marked as removed
	// swaps the given vertex with the last one in the vector and then pops the back of vertex vector
	// keeps the vertex vector dense, by moving the last vector at the place of the given one
	void free_vert(const int& vert_id);

	// frees the memory space occupied by a half-edge that was marked as removed
	// basically it moves the last edge in the hedge array to the place of the given one and then pops the last two hedges
	void free_edge(const int& hedge_id);


	bool is_in_range_vert(const int& vert_id) const;
	bool is_in_range_hedge(const int& hedge_id) const;
	//bool is_in_range_edge(const int& edge_id) const;

	void check_vert_id(const int& vert_id) const;
	void check_hedge_id(const int& hedge_id) const;
	//void check_edge_id(const int& edge_id) const;

	// These set_begin_id methods only set the begin_ids of the half-edges in the loop of
	//  the given half-edge , it does not update or modify the faces and borders sets !!!
	void set_begin_id(const int& hedge_id); // sets the begin_id of the given half-edge to the given half-edge
	void set_begin_id(const int& hedge_id, const int& begin_id); // sets the begin half-edge index in the loop of the given half-edge to the given begin half-edge index

	// circulates the given half-edge around its head vertex and sets all
	//  the half-edges head vertex to the new given head vertex index
	void set_vert_id(const int& hedge_id, const int& new_head_vert_id);

	// circulates the two half edges around their head vertex until it finds a common face 
	//  or border loop, but prioritizes face loops, so if it finds both a common face
	//  and a common border, it will choose the face, returns the begin half-edge index of the loop
	int find_common_face_around_head(int& he1, int& he2);

	// circulates the half-edge around its head vertex, until it finds a border, 
	//  if no border is found returns HE_INVALID_INDEX and the given half_edge stays the same
	int find_border_around_head(int& hedge_id);

	int id_ltwin(const int& hedge_id) const; // gives the smaller/even half-edge index of the edge of the given half-edge 
	int id_rtwin(const int& hedge_id) const; // gives the bigger/odd half-edge index of the edge of the given half-edge 


	//int new_vert();	// you dont need that, since it is the same as add_vert
	int new_edge(); // creates a new INVALID_INDEX initialized pair of half-edges, returns the half-edge index of the edge
	int new_edge(int tail_id, int head_id);

	// if the HEMesh is in a valid state, you should use id_begin(he1) == id_begin(he2),
	//  but during modification of the begin_ids of the half-edges, you can rely on this method
	//  which checks if they are in the same loop, based on their next_id, so it 
	//  starts to loop from he1 and searches for he2
	bool are_same_loop(const int& hedge_id1, const int& hedge_id2) const;

	void make_face(const int& hedge_id);	// makes the face or border loop of the given half-edge a face
	void make_border(const int& hedge_id);	// makes the face or border loop of the given half-edge a border

	void remove_dub_edges(int vert_id); // removes edges around the vertex untill all edges point to a unique head vertex

	void invalidate_vert(const int& vert_id);
	void invalidate_edge(const int& hedge_id);

	V get_face_center(const int& hedge_id) const;

	friend class vert_iter;
	friend class hedge_iter;
	friend class edge_iter;
	friend class face_iter;

public:

	void check_validity() const;

	class vert_iter;
	class hedge_iter;
	class edge_iter;
	class face_iter;

	HEMesh() {}
	HEMesh(const Mesh<V>& m);
	operator Mesh<V>() const;

	// Supports .obj and .ply(ascii) file formats, if the file format is other an error is thrown
//	HEMesh(const char* filePath);

//	void load_from_obj(const char* file_path);
//	void save_to_obj(const char* file_path);
	//void load_from_ply(const char* file_path);
	//void save_to_ply(const char* file_path);

	void clear();
	void shrink_verts(); // shrinks the memory of the vertices
	void shrink_edges(); // shrinks the memory of the half edges and edges
	void shrink();		 // shrinks both the vertices and edges

	bool empty() const;
	bool are_verts_shrinked() const;
	bool are_edges_shrinked() const;
	bool is_shrinked() const;

	int verts_size() const;
	int hedges_size() const;
	int edges_size() const;
	int faces_size() const;
	int borders_size() const;

	int verts_max_size() const;	// gives the count of all vertices including removed/invalid ones
	int hedges_max_size() const; // gives the count of all half-edges including removed/invalid ones
	int edges_max_size() const;	// gives the count of all edges including removed_invalid ones

	void reserve_verts(const int& capacity);
	void reserve_hedges(const int& capacity);
	void reserve_edges(const int& capacity);
	//void reserve_faces(const int& capacity);
	//void reserve_borders(const int& capacity);


	const int& id_head(const int& hedge_id) const;
	const int& id_tail(const int& hedge_id) const;
	const int& id_next(const int& hedge_id) const;
	int id_prev(const int& hedge_id) const;
	int id_twin(const int& hedge_id) const;
	const int& id_begin(const int& hedge_id) const;
	int id_begin(const int& vert_id1, const int& vert_id2) const; // gives the begin half edge of the common face between the two vertices, returns HE_INVALID_INDEX if not found
	const int& id_hedge(const int& tail_vert_id) const; // gives the outgoing half-edge of the given vertex, such that id_tail(id_hedge(vert_id)) == vert_id
	int id_hedge(const int& tail_id, const int& head_id) const;	// gives the id of the half edge with the given tail and head vertices, returns HE_ISOLATED_INDEX if not found
	//int id_edge(const int& hedge_id) const;

	// Consider changing id_hedge(vert_id) to these more specific in and out hedge methods
	//  then you can make the id_hedge method to convert an edge index to a half-edge index
	// And if you choose not to do that, then be consistent with that when converting
	//  a vertex index to an incident half-edge index, it always gives the out-going half-edge index
//	int id_hedge_out(const int& vert_id) const; // gives the outgoing half-edge of the given vertex, such that id_tail(id_hedge_out(vert_id)) == vert_id
//	int id_hedge_in(const int& vert_id) const;  // gives the ingoing half-edge to the given vertex, such that id_head(id_hedge_in(vert_id)) == vert_id


	bool is_removed_vert(const int& vert_id) const;
	bool is_removed_hedge(const int& hedge_id) const;
	//bool is_removed_edge(const int& edge_id) const;
	bool is_begin_hedge(const int& hedge_id) const;	// is the given half edge the beginning of a face or border loop

	bool is_border_hedge(const int& hedge_id) const;
	bool is_border_edge(const int& hedge_id) const;
	bool is_border_vert(const int& vert_id) const;
	bool is_border_loop(const int& begin_id) const;
	bool is_n_poly(const int& hedge_id, const int& n) const;	// gives the size of the face or border loop of the given half edge // its cheaper to use this instead of get_n_poly(hedge_id) == n

	int get_n_poly(const int& hedge_id) const;	// gives the size of the face or border loop that the half edge is in
	int get_valence(const int& vert_id) const;	// gives the number of adjacent edges of the vertex

	bool is_extremal_vert(const int& vert_id) const;	// is the vertex adjacent to only one edge
	bool is_isolated_vert(const int& vert_id) const;	// true if the vertex doesnt have any adjacent edges
	bool is_isolated_edge(const int& hedge_id) const;
	//bool is_isolated_loop(const int& hedge_id) const;
	bool is_isolated_face(const int& hedge_id) const;	// works both for face and border loops

	// Keep in mind that for the user a vert, edge or face may appear isolated
	//  when drawn on the screen or converted to a mesh,
	//  a vert, edge or face actually may only have adjacent edges which are border from both sides
	//  so consider implementing methods for this kind of situation
	//  (but i have no idea how it is most appropriate to name them)

	// In fact, the HalfEdgeMesh may be composed only out of borders
	//  and no faces (then it just represents a graph),
	//  and when you convert to a mesh it will have a bunch of vertices
	//  , but no triangle faces at all
	// So either make the converting to a mesh to add only vertices, which have an adjacent face
	//  or constraint the HEMesh such that there are no edges which are border from both sides
	//  but this is very restrictive, because a user may want to use the HEMesh only for that


	// ============================ BASIC OPERATIONS =================================== //

	void shift_begin(const int begin_hedge_id, const int how_much = 1);	// shifts the face or border loop begin half-edge

	void fill_holes();	// fills all holes of the surface of the mesh (all border loops are made face loops)
	void triangulate(const int& hedge_id);	// triangulates the face loop of the given half edge
	void triangulate();		// triangulates the whole mesh
	void remove_internal_edges(const int& hedge_id); // removes all internal edges (if it has any) in the face of the given half-edge that may appear after some operations
	void remove_internal_edges();	// removes internal edges of all the face loops in the mesh

	void remove_degrade_faces();	// removes all faces with size 2 (faces that only have 2 edges)
	void remove_isolated_verts();	// clears the mesh of all vertices that doesnt share any edge

	void swap_next(const int& hedge_id1, const int& hedge_id2); // swaps the next indexes of the two half-edges and rearanges the face or border begin hedge ids properly
	//	void swap_next(const int& hedge_id1, const int& hedge_id2, const bool bFace);
	//	void swap_next_same_loop(const int& hedge_id1, const int& hedge_id2);
	//	void swap_next_diff_loop(const int& hedge_id1, const int& hedge_id2);

	// splits the face or border loop adding a new edge at the head vertices
	//  of the given half-edges, it returns the id of the newly created face loop
	//  , if you want the index of the new half-edge, it is id_next(back_hedge_id)
	//  and its twin is id_next(front_hedge_id)
	int split_face_at(int back_hedge_id, int front_hedge_id);	

	//void make_isolated(const int& vert_id);	// removes all adjacent edges of the given vertex
	void remove_edges(const int& vert_id);	// removes all adjacent edges of the given vertex



	int add_vert(const V& v); // returns the index of the new vertex
	int add_edge(const int& tail_id, const V& head);	// returns the index of the new half_edge which head vertex is the id of the new added head vertex
	int add_edge(const int& tail_id, const int& head_id, const bool bMakeFace = true);	// if bMakeFace is true, it makes a face if the two given vertices had shared a common border loop

	// Adds a half-edge, such that its next hedge becomes the next of front_hedge (before the operation)
	//  and the next of its twin becomes the next of back_hedge (before the operation)
	int add_edge_at(const int& back_hedge_id, const int& front_hedge_id, const bool bMakeFace = true);
	int add_edge_at(const int& back_hedge_id, const V& head);

	int add_face(const std::vector<int>& vert_indices);	// returns the half-edge index that indicates the beginning of the face loop
	int add_face_at(const std::vector<int>& hedge_indices, const bool bMakeFace = true);

	void remove_vert(const int& vert_id);
	void remove_edge(const int& hedge_id);
	void remove_edge(const int& vertA, const int& vertB);
	void remove_face(const int& hedge_id);

	// reconnects the half-edge, such that its head vertex points to the head of the front_hedge
	// after the operation: hedge_id == id_twin(id_next(front_hedge_id)
	void move_edge(const int& hedge_id, const int& front_hedge_id);

	// basically it is a faster version of doing 
	//  remove_edge(hedge_id) and then using add_edge_at(back_hedge_id, front_hedge_id)
	// after the operation: hedge_id == id_twin(id_next(front_hedge_id)) and id_twin(hedge_id) == id_twin(id_next(back_hedge_id))
	void move_edge(const int& hedge_id, const int& back_hedge_id, const int& front_hedge_id);


	// ========================== HIGH-LEVEL OPERATIONS ================================ //

	void flip_edge(const int& hedge_id);

	// adds a new vertex between the edge, the given half-edge's head vertex
	//  is the new one, after the operations and its next is the new edge half-edge
	//  it returns the id of the new created vertex
	int refine_edge(const int& hedge_id, const V& v);
	int refine_edge(const int& hedge_id, const float h = 0.5f);

	// given two half-edges, it adds a new edge between them, returns a new half-edge
	//  index (from the created edge) and its head vertex is the new vertex
	int split_vert_to_edge(const int& hedgeA, const int& hedgeB, const V& v);

	// given two half-edges, it adds a new edge and two triangle faces between them
	//  , returns a new half-edge index (from the created edge) 
	//  and its head vertex is the new vertex
	int split_vert_to_faces(const int& hedgeA, const int& hedgeB, const V& v);

	int split_edge(const int& hedge_id, const V& v);
	int split_edge(const int& hedge_id, const float h = 0.5f);

	int clip_corner(const int& hedge_id);	// clips the corner that the hedge points to into a new triangle face

	// Collapses the edge into a vertex, returns the index of the vertex
	int collapse_edge(const int& hedge_id, const V& v);
	int collapse_edge(const int& hedge_id, const float h = 0.5f);

	// Collapses the face into a vertex, returns the index of the vertex
	int collapse_face(const int& hedge_id, const V& v);
	int collapse_face(const int& hedge_id);	// the resulting vertex is the midpoint/center of the face

	// Bevel operations, they return the half-edge index at the begining of the new face 
	int bevel_vert(const int& vert_id, const float h = 0.5f);
	int bevel_edge(const int& hedge_id, const float h = 0.5f);
	int bevel_face(const int& hedge_id, const float h = 0.5f);

	// cuts an edge along its length and creates a new 2-sided degrade
	//  border and connects it with other adjacent borders
	// , it returns the id of the border that is created
	int cut_edge(const int& hedge_id);


	// With Iterators
	/*
	void swap_next(const hedge_iter& hedge1, const hedge_iter& hedge2); // swaps the next indexes of the two hedges and rearanges the face or border begin hedge ids properly

	//vert_iter add_vert(const V& v);
	int add_edge(const vert_iter& tail, const V& head);
	int add_edge(const vert_iter& tail, const vert_iter& head);
	int add_edge_at(const hedge_iter& back_hedge, const hedge_iter& front_hedge);
	int add_face(const std::vector<vert_iter>& indices);

	void remove_vert(const vert_iter& vert);
	void remove_edge(const hedge_iter& hedge);
	void remove_face(const hedge_iter& hedge);

	void move_edge(const hedge_iter& hedge, const hedge_iter& new_next);
	void move_edge(const hedge_iter& hedge, const hedge_iter& new_prev, const hedge_iter& new_next);
	*/


	// =================================== ITERATORS ========================================= //


	  // Convertion from index to iterators
	vert_iter vert(const int& vert_id) const;
	hedge_iter hedge(const int& hedge_id) const;
	hedge_iter hedge(const int& tail_id, const int& head_id) const;
	edge_iter edge(const int& edge_id) const;
	 face_iter face(const int& hedge_id) const;


	// Forward iterators
	vert_iter begin_verts() const;
	hedge_iter begin_hedges() const;
	edge_iter begin_edges() const;
	face_iter begin_faces() const;
	face_iter begin_borders() const;

	vert_iter end_verts() const;
	hedge_iter end_hedges() const;
	edge_iter end_edges() const;
	face_iter end_faces() const;
	face_iter end_borders() const;


	// Reversed iterators (if you want these, you should implement reversed versions
	//  of the vert and hedge iterators, but for now lets keep it simple)
	/*vert_iter rbegin_verts() const;
	hedge_iter rbegin_hedges() const;
	edge_iter rbegin_edges() const;

	vert_iter rend_verts() const;
	hedge_iter rend_hedges() const;
	edge_iter rend_edges() const;*/


	// Iterables
	Iterable<vert_iter> verts() const;
	Iterable<hedge_iter> hedges() const;
	Iterable<edge_iter> edges() const;
	Iterable<face_iter> faces() const;
	Iterable<face_iter> borders() const;

	//Iterable<vert_iter> rvertices() const;
	//Iterable<vert_iter> rhedges() const;
	//Iterable<vert_iter> redges() const;

	// @todo: Hash for the iterators


	class vert_iter {
	protected:

		HEMesh<V>* hm;
		int id;

		// Only HEMesh and hedge_iter will use the private constructor of this class
		friend class HEMesh<V>;
		friend class HEMesh<V>::hedge_iter;

		vert_iter(const int& vert_id, HEMesh<V>* const& hmm) : id(vert_id), hm(hmm) {}

	public:

		vert_iter() : id(HE_INVALID_INDEX), hm(nullptr) {}
		
	//	vert_iter& operator=(const int& id) { this->id = id; return *this; }

		operator const int& () const { return id; }
		//operator int() const { return id; }
		const int& index() const { return id; }

		bool is_removed() const { return hm->is_removed_vert(id); }
		bool is_border() const { return hm->is_border_vert(id); }
		bool is_isolated() const { return hm->is_isolated_vert(id); }
		bool is_extremal() const { return hm->is_extremal_vert(id); }
		//bool degree() const { return hm->get_valence(id); }
		int valence() const { return hm->get_valence(id); }

		V& operator*() { return hm->_vertices[id]; }
		V* operator->() { return &(hm->_vertices[id]); }
		const V& operator*() const { return hm->_vertices[id]; }
		const V* operator->() const { return &(hm->_vertices[id]); }

		vert_iter& operator++() { do { ++id; } while (id < hm->_vertices.size() && hm->is_removed_vert(id)); return *this; }	// Prefix increment
		vert_iter operator++(int) { vert_iter temp = *this; ++(*this); return temp; }	// Postfix increment
		vert_iter& operator--() { do { --id; } while (id >= 0 && hm->is_removed_vert(id)); return *this; }	// Prefix decrement
		vert_iter operator--(int) { vert_iter temp = *this; --(*this); return temp; }	// Postfix decrement

		friend bool operator== (const vert_iter& a, const vert_iter& b) { return a.id == b.id; };
		friend bool operator!= (const vert_iter& a, const vert_iter& b) { return a.id != b.id; };

		// Gives the half edge that points out of the vertex, so this vertex is the source/tail of it
		hedge_iter hedge() const { return hedge_iter(hm->_vert_to_hedge[id], hm); }
	};


	class hedge_iter {
	protected:

		HEMesh<V>* hm;
		int id;

		// Thats because all those will use the private constructor of hedge_iter
		friend class HEMesh<V>;
		friend class HEMesh<V>::edge_iter;
		friend class HEMesh<V>::face_iter;
		friend class HEMesh<V>::vert_iter;

		hedge_iter(const int& hedge_id, HEMesh<V>* const& hmm) : id(hedge_id), hm(hmm) {}

	public:

		hedge_iter() : id(HE_INVALID_INDEX), hm(nullptr) {}

	//	hedge_iter& operator=(const int& id) { this->id = id; return *this; }

		//operator edge_iter() const { return edge_iter(id / 2, hm); }
	//	operator edge_iter() const { return edge_iter(id >> 1, hm); }
		operator edge_iter() const { return edge_iter(id, hm); }

		operator const int& () const { return id; }
		//operator int() const { return id; }
		const int& index() const { return id; }

		bool is_removed() const { return hm->is_removed_hedge(id); }
		bool is_border() const { return hm->is_border_hedge(id); }
		bool is_isolated() const { return hm->is_isolated_edge(id); } // true if it is adjacent only to its twin half edge and no other half edges

		hedge_iter next() const { return hedge_iter(hm->_hedges[id].next_id, hm); }
		hedge_iter twin() const { return hedge_iter(hm->id_twin(id), hm); }
		hedge_iter prev() const { return hedge_iter(hm->id_prev(id), hm); }
		//edge_iter edge() const { return edge_iter(id / 2, hm); }
	//	edge_iter edge() const { return edge_iter(id >> 1, hm); }
		edge_iter edge() const { return edge_iter(id, hm); }
		vert_iter head() const { return vert_iter(hm->_hedges[id].head_id, hm); }
		//vert_iter tail() const { return vert_iter(hm->_hedges[hm->id_twin(id)].head_id, hm); }
		vert_iter tail() const { return vert_iter(hm->_hedges[id ^ 1].head_id, hm); }
		face_iter face() const { return face_iter(id, hm); }


		hedge_iter& operator*() { return *this; }
		hedge_iter* operator->() { return this; }
		const hedge_iter& operator*() const { return *this; }
		const hedge_iter* operator->() const { return this; }

		hedge_iter& operator++() { do { ++id; } while (id < hm->_hedges.size() && hm->is_removed_hedge(id)); return *this; } // Prefix increment
		hedge_iter operator++(int) { hedge_iter temp = *this; ++(*this); return temp; } // Postfix increment
		hedge_iter& operator--() { do { --id; } while (id >= 0 && hm->is_removed_hedge(id)); return *this; } // Prefix decrement
		hedge_iter operator--(int) { hedge_iter temp = *this; --(*this); return temp; } // Postfix decrement

		friend bool operator== (const hedge_iter& a, const hedge_iter& b) { return a.id == b.id; };
		friend bool operator!= (const hedge_iter& a, const hedge_iter& b) { return a.id != b.id; };

	};


	class edge_iter {
	protected:

		HEMesh<V>* hm;
		int hedge_id;

		friend class HEMesh<V>;
		friend class HEMesh<V>::hedge_iter;

		//	edge_iter(const int& edge_id, HEMesh<V>* const& hmm) : hedge_id(edge_id << 1), hm(hmm) {}
		edge_iter(const int& hedge_id, HEMesh<V>* const& hmm)
			: hedge_id((hedge_id >> 1) << 1), hm(hmm) {}

	public:

		//edge_iter() : id(HE_INVALID_INDEX), hm(nullptr) {}
		edge_iter() : hedge_id(HE_INVALID_INDEX), hm(nullptr) {}

	//	edge_iter& operator=(const int& id) { this->id = (id >> 1) << 1; return *this; }

		//operator const int& () const { return id; }
		//operator int() const { return id; }
		//const int& index() const { return id; }
		
	//	operator const int& () const { return hedge_id; }
		//operator int() const { return hedge_id; }
		int index() const { return hedge_id >> 1; }


		//hedge_iter hedge() const { return hedge_iter(id << 1, hm); }
		//operator hedge_iter() const { return hedge_iter(id << 1, hm); }
		hedge_iter hedge() const { return hedge_iter(hedge_id, hm); }
		operator hedge_iter() const { return hedge_iter(hedge_id, hm); }

		//bool is_removed() const { return hm->is_removed_edge(id); }
		//bool is_border() const { return hm->is_border_edge(id << 1); }
		//bool is_isolated() const { return hm->is_isolated_edge(id << 1); }
		bool is_removed() const { return hm->is_removed_hedge(hedge_id); }
		bool is_border() const { return hm->is_border_edge(hedge_id); }
		bool is_isolated() const { return hm->is_isolated_edge(hedge_id); }

		edge_iter& operator*() { return *this; }
		edge_iter* operator->() { return this; }
		const edge_iter& operator*() const { return *this; }
		const edge_iter* operator->() const { return this; }

		//edge_iter& operator++() { do { ++id; } while (id < hm->edges_size() && hm->is_removed_edge(id)); return *this; } // Prefix increment
		edge_iter& operator++() { do { hedge_id += 2; } while (hedge_id < hm->_hedges.size() && hm->is_removed_hedge(hedge_id)); return *this; } // Prefix increment
		edge_iter operator++(int) { edge_iter temp = *this; ++(*this); return temp; } // Postfix increment

		//edge_iter& operator--() { do { --id; } while (id >= 0 && hm->is_removed_edge(id)); return *this; } // Prefix increment
		edge_iter& operator--() { do { hedge_id -= 2; } while (hedge_id >= 0 && hm->is_removed_hedge(hedge_id)); return *this; } // Prefix decrement
		edge_iter operator--(int) { edge_iter temp = *this; --(*this); return temp; } // Postfix decrement

		friend bool operator== (const edge_iter& a, const edge_iter& b) { return a.hedge_id == b.hedge_id; };
		friend bool operator!= (const edge_iter& a, const edge_iter& b) { return a.hedge_id != b.hedge_id; };
	};


	// Wrapper around the unordered_set<int> iterator
	class face_iter {
	protected:
		using SetIter = typename std::unordered_set<int>::iterator;

		HEMesh<V>* hm;
		SetIter iter;

		friend class HEMesh<V>;
		friend class HEMesh<V>::hedge_iter;		// because of hedge_iter::face()

		face_iter(const SetIter& it, HEMesh<V>* const& hmm) : iter(it), hm(hmm) {}

		face_iter(const int& hedgeInFaceID, HEMesh<V>* const& hmm) : hm(hmm) {
			const int& fi = hm->_hedges[hedgeInFaceID].begin_id;
			auto it = hm->_faces.find(fi);
			iter = it != hm->_faces.end() ? it : hm->_borders.find(fi);
		}

	public:

		face_iter() = delete;

		operator const int& () const { return *iter; }
		//operator int() const { return *iter; }
		const int& index() const { return *iter; }	// gives the begin half-edge index !

		operator hedge_iter() const { return hedge(); }
		hedge_iter hedge() const { return hedge_iter(*iter, hm); }

		bool is_removed() const { return hm->is_begin_hedge(*iter); }
		bool is_border() const { return hm->is_border_loop(*iter); }
		bool is_isolated() const { return hm->is_isolated_face(*iter); }
		int size() const { return hm->get_n_poly(*iter); }
		bool is_size(const int& n) const { return hm->is_n_poly(*iter, n); }


		const face_iter& operator*() const { return *this; }
		const face_iter* operator->() const { return this; }
		face_iter& operator*() { return *this; }
		face_iter* operator->() { return this; }

		face_iter& operator++() { ++iter; return *this; }	// Prefix increment
		face_iter operator++(int) { face_iter temp = *this; ++(*this); return temp; }		// Postfix increment

		friend bool operator== (const face_iter& a, const face_iter& b) { return a.iter == b.iter; };
		friend bool operator!= (const face_iter& a, const face_iter& b) { return a.iter != b.iter; };

	};

};


// ================================== IMPLEMENTATION ======================================== //


template<typename V>
inline void HEMesh<V>::free_vert(const int& vert_id) {
	if (!is_in_range_vert(vert_id))
		throw std::out_of_range("Cannot free space of a vertex which index is out of range !");

	if (!is_removed_vert(vert_id))
		throw std::invalid_argument("The vertex should be removed, to free its memory space");

	int ov = (int)_vertices.size() - 1;
	if (vert_id != ov) {
		_vertices[vert_id] = _vertices[ov];

		_vert_to_hedge[vert_id] = _vert_to_hedge[ov];

		if (id_hedge(ov) != -1)
			set_vert_id(id_twin(id_hedge(ov)), vert_id); // should i make set_vert_id to set _vert_to_hedge too ???

		//_vert_to_hedge[vert_id] = _vert_to_hedge[ov];
	}

	_vertices.pop_back();
	_vert_to_hedge.pop_back();
}

template<typename V>
inline void HEMesh<V>::free_edge(const int& hedge_id) {
	if (!is_in_range_hedge(hedge_id))
		throw std::out_of_range("Cannot free space of a half-edge that is out of range !");

	if (!is_removed_hedge(hedge_id))
		throw std::invalid_argument("The edge should be removed, to free its memory space");

	//if (hedge_id >= (int)_hedges.size() - 2) {
	if (hedge_id + 2 >= (int)_hedges.size()) {
		_hedges.pop_back();
		_hedges.pop_back();
		return;
	}

	int t = id_ltwin(hedge_id);
	int c = t + 1;

	int cm = int(_hedges.size()) - 1;
	int tm = int(_hedges.size()) - 2;

	// update the previous of the last hedges to point to the new location 
	int pcm = id_prev(cm);
	int ptm = id_prev(tm);
	_hedges[pcm].next_id = c;
	_hedges[ptm].next_id = t;

	// replase the erased edge with the last edge
	_hedges[c] = _hedges.back();
	_hedges.pop_back();
	_hedges[t] = _hedges.back();
	_hedges.pop_back();

	_vert_to_hedge[id_head(c)] = t;
	_vert_to_hedge[id_head(t)] = c;

	// update the faces, if it happens that the moved hedges are the beginning of a face loop

	if (id_begin(c) == cm) {

		if (_borders.find(cm) == _borders.end()) {
			_faces.erase(cm);
			_faces.insert(c);
		}
		else {
			_borders.erase(cm);
			_borders.insert(c);
		}

		set_begin_id(c);
	}

	if (id_begin(t) == tm) {

		if (_borders.find(tm) == _borders.end()) {
			_faces.erase(tm);
			_faces.insert(t);
		}
		else {
			_borders.erase(tm);
			_borders.insert(t);
		}

		set_begin_id(t);
	}
}

template<typename V>
inline bool HEMesh<V>::is_in_range_vert(const int& vert_id) const {
	return 0 <= vert_id && vert_id < (int)_vertices.size();
}

template<typename V>
inline bool HEMesh<V>::is_in_range_hedge(const int& hedge_id) const {
	return 0 <= hedge_id && hedge_id < (int)_hedges.size();
}

/*template<typename V>
inline bool HEMesh<V>::is_in_range_edge(const int& edge_id) const {
	//int hedge_id = edge_id * 2;
	int hedge_id = edge_id << 1;
	return 0 <= hedge_id && hedge_id < _hedges.size());
}*/

template<typename V>
inline void HEMesh<V>::check_vert_id(const int& vert_id) const {
	if (!is_in_range_vert(vert_id))
		throw std::out_of_range("Vertex index out of range");
	if (is_removed_vert(vert_id))
		throw std::invalid_argument("Removed vertex");
}

template<typename V>
inline void HEMesh<V>::check_hedge_id(const int& hedge_id) const {
	if (!is_in_range_hedge(hedge_id))
		throw std::out_of_range("Half edge index out of range");
	if (is_removed_hedge(hedge_id))
		throw std::invalid_argument("Removed half edge");
}

/*template<typename V>
inline void HEMesh<V>::check_edge_id(const int& edge_id) const {
	if (!is_in_range_edge(edge_id))
		throw std::out_of_range("Edge index out of range");
	if (is_removed_edge(edge_id))
		throw std::invalid_argument("Removed edge");
}*/


template<typename V>
inline void HEMesh<V>::set_begin_id(const int& hedge_id) {
	int e = hedge_id;
	do {
		_hedges[e].begin_id = hedge_id;

		e = id_next(e);
	} while (e != hedge_id);
}

template<typename V>
inline void HEMesh<V>::set_begin_id(const int& hedge_id, const int& begin_id) {
	int e = hedge_id;
	do {
		_hedges[e].begin_id = begin_id;

		e = id_next(e);
	} while (e != hedge_id);
}

template<typename V>
inline void HEMesh<V>::set_vert_id(const int& hedge_id, const int& new_head_vert_id) {
	int beg = hedge_id;	// beggining
	int e = beg;
	do {
		_hedges[e].head_id = new_head_vert_id;
		e = id_twin(id_next(e));
	} while (e != beg);

	//for (int e = hedge_id; id_head(e) != new_head_vert_id; e = id_twin(id_next(e)))
	//	_hedges[e].head_id = new_head_vert_id;

//	_vert_to_hedge[new_head_vert_id] = id_next(hedge_id);	// should this method also update the _vert_to_hedge vector ??
}

template<typename V>
inline int HEMesh<V>::find_common_face_around_head(int& he1, int& he2) {

	int hedge1 = he1;
	int hedge2 = he2;

	int b = HE_INVALID_INDEX;
	int be1 = he1;
	int be2 = he2;

	// brute force
	do {
		do {
			const int& f = id_begin(he1);

			if (f == id_begin(he2)) {
				if (_borders.find(f) == _borders.end())
					return f;

				b = f;
				be1 = he1;
				be2 = he2;
			}

			he1 = id_twin(id_next(he1));	// rotates around the head vertex
		} while (he1 != hedge1);

		he2 = id_twin(id_next(he2));	// rotates around the head vertex
	} while (he2 != hedge2);

	he1 = be1;
	he2 = be2;

	return b;
}

template<typename V>
inline int HEMesh<V>::find_border_around_head(int& hedge_id) {
	int begin_id = hedge_id;

	do {
		const int& b = id_begin(hedge_id);
		if (is_border_loop(b))
			return b;

		hedge_id = id_twin(id_next(hedge_id));
	} while (hedge_id != begin_id);

	return HE_INVALID_INDEX;
}

template<typename V>
inline int HEMesh<V>::id_ltwin(const int& hedge_id) const {
	//return (hedge_id / 2) * 2;
	return (hedge_id >> 1) << 1;
}

template<typename V>
inline int HEMesh<V>::id_rtwin(const int& hedge_id) const {
	return ((hedge_id >> 1) << 1) + 1;
}



template<typename V>
inline int HEMesh<V>::new_edge() {
	int ne;		// new half-edge index

	if (_garbage_edges.empty()) {
		ne = _hedges.size();
		_hedges.push_back(HalfEdge());
		_hedges.push_back(HalfEdge());
	}
	else {
		auto front_iter = _garbage_edges.begin();
		ne = *front_iter;

		//if (ne % 2 != 0 || _hedges[ne] != HalfEdge() || _hedges[ne + 1] != HalfEdge())
		//	throw std::exception("There is an edge that is garbaged, but not invalidated");

		_garbage_edges.erase(front_iter);
	}

	return ne;
}

template<typename V>
inline int HEMesh<V>::new_edge(int tail_id, int head_id) {
	int ne;		// new half-edge index
	int net;	// new twin half-edge index

	if (_garbage_edges.empty()) {
		ne = _hedges.size();
		net = ne + 1;

		_hedges.push_back(HalfEdge(head_id, net, ne));
		_hedges.push_back(HalfEdge(tail_id, ne, ne));
	}
	else {
		auto front_iter = _garbage_edges.begin();
		ne = *front_iter;
		net = ne + 1;
		_garbage_edges.erase(front_iter);

		_hedges[ne] = HalfEdge(head_id, net, ne);
		_hedges[net] = HalfEdge(tail_id, ne, ne);
	}

	return ne;
}




template<typename V>
inline void HEMesh<V>::check_validity() const {
	//std::cout << "======================================================\n";

	for (const int& f : _faces) {

		int e = f;
		do {
			if (is_removed_hedge(e)) {
				std::cout << "Invalid hedge " << e << " with begin of face " << f << "\n";
				break;
			}
			else {
				if (id_begin(e) != f)
					std::cout << "Hedge " << e << " points to " << id_begin(e) << " rather than its face begin " << f << "\n";

				if (e != f && _faces.find(e) != _faces.end())
					std::cout << "Dublicated face hedge " << e << " in face " << f << "\n";

				if (e != f && _borders.find(e) != _borders.end())
					std::cout << "Dublicated border hedge " << e << " in face " << f << "\n";
			}

			e = id_next(e);
		} while (e != f);
	}

	for (const int& b : _borders) {

		int e = b;
		do {
			if (is_removed_hedge(e)) {
				std::cout << "Invalid hedge " << e << " at border " << b << "\n";
				break;
			}
			else {
				if (id_begin(e) != b)
					std::cout << "Hedge " << e << " points to begin " << id_begin(e) << " rather than its border begin " << b << "\n";

				if (e != b && _faces.find(e) != _faces.end())
					std::cout << "Dublicated face hedge " << e << " in border " << b << "\n";

				if (e != b && _borders.find(e) != _borders.end())
					std::cout << "Dublicated border hedge " << e << " in border " << b << "\n";
			}

			e = id_next(e);
		} while (e != b);
	}


	for (int vi = 0; vi < _vert_to_hedge.size(); ++vi) {
		int e = _vert_to_hedge[vi];

		if (e == HE_INVALID_INDEX || e == HE_ISOLATED_INDEX)
			continue;

		if (is_removed_hedge(e))
			std::cout << "Tail vert " << vi << " points to removed hedge " << e << "\n";
		else if (id_tail(e) != vi)
			std::cout << "Tail vert " << vi << " points to hedge " << e << ", but its tail is " << id_tail(e) << "\n";

	
		e = id_twin(e);
		int beg = e;
		do {
			if (id_head(e) != vi)
				std::cout << "Vert " << vi << " has an incident ingoing hedge " << e << ", but the hedge doesnt point to the same head vertex\n";

			e = id_twin(id_next(e));
		} while (e != beg);
	}

	// Check if it points to a valid face or border index
	for (int e = 0; e < _hedges.size(); ++e) {
		if (is_removed_hedge(e)) {
		//if (_hedges[e] == HalfEdge())
		
			if (_garbage_edges.find(id_ltwin(e)) == _garbage_edges.end())
				std::cout << "Hedge " << e << " is invalidated, but is not in the garbage list of removed hedges\n";
			
			continue;
		}

		if (id_head(id_twin(e)) != id_head(id_prev(e)))
			//std::cout << "Hedge " << e << "'s tail vertex is inconsistent\n";
			std::cout << "The twin and prev hedges of hedge " << e << " , should have the same head vertex\n";

		const int& f = id_begin(e);
		if (_faces.find(f) == _faces.end() && _borders.find(f) == _borders.end())
			std::cout << "Hedge " << e << " points to non-existing face or border begin hedge  " << f << "\n";

		if (e == id_next(e))
			std::cout << "Infinite loop hedge " << e << "'s next hedge is itself\n";

		if (id_head(e) == id_tail(e))
			std::cout << "Degrade edge " << e << "(its head and tail vertices are the same)\n";

		if (id_next(id_next(e)) == e)
			std::cout << "Degrade face with begin " << id_begin(e) << "\n";

		if (id_hedge(id_head(e)) == HE_ISOLATED_INDEX)
			std::cout << "Hedge " << e << " points to a head vertex " << id_head(e) << " , which is marked as isolated\n";
		else if (id_hedge(id_head(e)) == HE_INVALID_INDEX)
			std::cout << "Hedge " << e << " points to a head vertex " << id_head(e) << " , which is marked as removed\n";
		else {

			int head_hedge = id_twin(id_hedge(id_head(e)));
			int ee = head_hedge;
			bool isValid = false;

			do {
				if (ee == e)
					isValid = true;

				ee = id_twin(id_next(ee));
			} while (ee != head_hedge);

			if (!isValid)
				std::cout << "Hedge " << e << " points to an adjacent vertex from which you cannot find adjacency with the hedge\n";
		}
	}

	for (const int& he : _garbage_edges) {
		if (he % 2 != 0)
			std::cout << "Only even half-edge indexes should be in the garbage set of half-edges !\n";
		
		//if (_hedges[he] != HalfEdge() || _hedges[he + 1] != HalfEdge())
		if (!is_removed_hedge(he) || !is_removed_hedge(he + 1))
			std::cout << "There is an pair of hedges " << he << " and " << he + 1 << " which are garbaged, but not invalidated\n";
	
	}
		
}



template<typename V>
inline HEMesh<V>::HEMesh(const Mesh<V>& m) {
	_vertices = m.vertices;

	_vert_to_hedge.assign(_vertices.size(), HE_ISOLATED_INDEX);
	_faces.reserve(m.indices.size() / 3ULL);
	_hedges.reserve(m.indices.size());
	//_hedges.assign(m.indices.size(), HalfEdge());
	//_hedges.resize(indices.size());

	std::unordered_map<std::pair<int, int>, size_t> vert_pair_to_hedge;

	//uint found_edges = 0U;

	int curr_hedge_id;
	int prev_hedge_id;
	int curr_begin_id;

	for (size_t i = 0U; i < m.indices.size(); i += 3ULL) {
		//_faces.push_back(i);

		for (size_t j = 0; j < 3ULL; ++j) {
			int id0 = m.indices[i + j];				// tail vertex id
			int id1 = m.indices[i + (j + 1ULL) % 3ULL];	// head vertex id

			if (vert_pair_to_hedge.find({ id0, id1 }) != vert_pair_to_hedge.end())
				throw std::invalid_argument("Trying to convert a Non-Consistently-Orientated or Non-2-Manifold mesh to a half edge data structure !");

			auto mp = vert_pair_to_hedge.find({ id1, id0 });
			if (mp == vert_pair_to_hedge.end()) {
				//curr_hedge_id = 2ULL * found_edges;
				//++found_edges;

				curr_hedge_id = _hedges.size();
				_hedges.push_back(HalfEdge());
				_hedges.push_back(HalfEdge());

				vert_pair_to_hedge[{id0, id1}] = curr_hedge_id;
			}
			else {
				curr_hedge_id = mp->second + 1ULL;
				vert_pair_to_hedge.erase(mp); // erase, so later we are left only with the border half edges
			}

			_vert_to_hedge[id0] = curr_hedge_id;

			if (j != 0ULL)
				_hedges[prev_hedge_id].next_id = curr_hedge_id;
			else
				curr_begin_id = curr_hedge_id;

			prev_hedge_id = curr_hedge_id;

			_hedges[curr_hedge_id].begin_id = curr_begin_id;

			//_hedges[curr_hedge_id] = HalfEdge(id1);
			_hedges[curr_hedge_id].head_id = id1;
		}

		_faces.insert(curr_begin_id);
		_hedges[prev_hedge_id].next_id = curr_begin_id;

	}

	// Assign Borders ========================================================== //

	std::unordered_set<size_t> inner_border_hedges;
	for (const auto& it : vert_pair_to_hedge)
		inner_border_hedges.insert(it.second);

	// Traverses the boundary inner HalfEdges, to create the outer HalfEdges
	while (!inner_border_hedges.empty()) {
		int begin = *inner_border_hedges.begin();   // inner edge border begin
		int e = begin;								// inner border edge
		int prev = HE_INVALID_INDEX;                // previous inner border edge

		do {

			inner_border_hedges.erase(inner_border_hedges.find(e));

			// If not first iteration
			//if (e != begin)
			if (prev != HE_INVALID_INDEX) {
				//_hedges[id_twin(e)].next_id = prev;
				_hedges[id_twin(e)].next_id = id_twin(prev);
				_hedges[id_twin(e)].head_id = _hedges[prev].head_id;
				_hedges[id_twin(e)].begin_id = curr_begin_id;
			}
			else
				curr_begin_id = id_twin(e);

			prev = e;

			// Find Next Inner Edge
			e = _hedges[e].next_id;
			while (e != begin) {
				//if (id_next(id_twin(e)) == HE_INVALID_INDEX)					// !!!!
				if (_hedges[id_twin(e)].next_id == HE_INVALID_INDEX)
					break;

				e = _hedges[id_twin(e)].next_id;
			}

		} while (e != begin);

		_hedges[id_twin(e)].next_id = id_twin(prev);
		_hedges[id_twin(e)].head_id = _hedges[prev].head_id;
		_hedges[id_twin(e)].begin_id = curr_begin_id;

		_borders.insert(curr_begin_id);
	}
}


// For now it is a simple triangulation, where one vertex is chosen to split the whole face
//  into triangles (the chosen vertex is connected with all other vertices)
template<typename V>
inline HEMesh<V>::operator Mesh<V>() const {

	//	if (!are_verts_shrinked())
	//		throw std::exception("Cannot convert a Half-Edge Mesh to a Triangular Mesh, if the vertices are not shrinked");

	Mesh<V> m;
	m.vertices = _vertices;
	m.indices.reserve(_hedges.size());

	int he[3];

	for (const int& f : _faces) {
		if (is_n_poly(f, 2))
			continue;

		he[0] = f;
		he[1] = id_next(he[0]);
		he[2] = id_next(he[1]);

		do {
			for (int i = 0; i < 3; ++i)
				m.indices.push_back(_hedges[he[i]].head_id);

			he[1] = he[2];
			he[2] = id_next(he[2]);

		} while (he[2] != he[0]);
	}

	return m;
}


template<typename V>
inline void HEMesh<V>::clear() {
	_vertices.clear();
	_hedges.clear();
	_faces.clear();
	_borders.clear();
	_garbage_vertices.clear();
	_garbage_edges.clear();
}


template<typename V>
inline void HEMesh<V>::shrink_verts() {
	for (auto it = _garbage_vertices.rbegin(); it != _garbage_vertices.rend(); ++it)
		free_vert(*it);

	_garbage_vertices.clear();
}

template<typename V>
inline void HEMesh<V>::shrink_edges() {
	for (auto it = _garbage_edges.rbegin(); it != _garbage_edges.rend(); ++it)
		free_edge(*it);

	_garbage_edges.clear();
}

template<typename V>
inline void HEMesh<V>::shrink() {
	shrink_verts();
	shrink_edges();
}


template<typename V>
inline bool HEMesh<V>::empty() const {
	return _faces.empty() && _borders.empty();
	//&& _vertices.size() == _garbage_vertices.size()
	//&& 2ULL * _hedges.size() == _garbage_edges.size()
}


template<typename V>
inline bool HEMesh<V>::are_verts_shrinked() const {
	return _garbage_vertices.empty();
}

template<typename V>
inline bool HEMesh<V>::are_edges_shrinked() const {
	return _garbage_edges.empty();
}

template<typename V>
inline bool HEMesh<V>::is_shrinked() const {
	return _garbage_vertices.empty() && _garbage_edges.empty();
}


template<typename V>
inline int HEMesh<V>::verts_size() const {
	return _vertices.size() - _garbage_vertices.size();
}

template<typename V>
inline int HEMesh<V>::hedges_size() const {
	return _hedges.size() - _garbage_edges.size() * 2ULL;
}

template<typename V>
inline int HEMesh<V>::edges_size() const {
	return _hedges.size() / 2ULL - _garbage_edges.size();
}

template<typename V>
inline int HEMesh<V>::faces_size() const {
	return _faces.size();
}

template<typename V>
inline int HEMesh<V>::borders_size() const {
	return _borders.size();
}


template<typename V>
inline int HEMesh<V>::verts_max_size() const {
	return _vertices.size();
}

template<typename V>
inline int HEMesh<V>::hedges_max_size() const {
	return _hedges.size();
}

template<typename V>
inline int HEMesh<V>::edges_max_size() const {
	return _hedges.size() / 2ULL;
}


template<typename V>
inline void HEMesh<V>::reserve_verts(const int& capacity) {
	_vertices.reserve(capacity);
}

template<typename V>
inline void HEMesh<V>::reserve_hedges(const int& capacity) {
	_hedges.reserve(capacity);
}

template<typename V>
inline void HEMesh<V>::reserve_edges(const int& capacity) {
	_hedges.reserve(2 * capacity);
}


template<typename V>
inline const int& HEMesh<V>::id_head(const int& hedge_id) const {
	return _hedges[hedge_id].head_id;
}

template<typename V>
inline const int& HEMesh<V>::id_tail(const int& hedge_id) const {
	return _hedges[id_twin(hedge_id)].head_id;
}

template<typename V>
inline const int& HEMesh<V>::id_next(const int& hedge_id) const {
	return _hedges[hedge_id].next_id;
}

template<typename V>
inline int HEMesh<V>::id_prev(const int& hedge_id) const {
	int e = hedge_id;

	while (_hedges[e].next_id != hedge_id)
		e = _hedges[e].next_id;

	return e;
}

template<typename V>
inline int HEMesh<V>::id_twin(const int& hedge_id) const {
	//int j = hedge_id % 2;
	//return hedge_id - j + (j + 1) % 2;
	return hedge_id ^ 1;
}

template<typename V>
inline const int& HEMesh<V>::id_begin(const int& hedge_id) const {
	return _hedges[hedge_id].begin_id;
}

template<typename V>
inline int HEMesh<V>::id_begin(const int& vert_id1, const int& vert_id2) const {
	const int& beg1 = _vert_to_hedge[vert_id1];
	const int& beg2 = _vert_to_hedge[vert_id2];

	int he1 = beg1;
	int he2 = beg2;

	do {
		do {

			if (id_begin(he1) == id_begin(he2))
				return id_begin(he1);

			he2 = id_next(id_twin(he2));
		} while (he2 != beg2);

		he1 = id_next(id_twin(he1));
	} while (he1 != beg1);

	return HE_INVALID_INDEX;
}

template<typename V>
inline const int& HEMesh<V>::id_hedge(const int& tail_vert_id) const {
	return _vert_to_hedge[tail_vert_id];
}

template<typename V>
inline int HEMesh<V>::id_hedge(const int& tail_id, const int& head_id) const {
	const int& et = _vert_to_hedge[tail_id];
	const int& eh = _vert_to_hedge[head_id];

	if (et != HE_ISOLATED_INDEX && eh != HE_ISOLATED_INDEX) {
		int e = et;

		do {
			if (id_head(e) == head_id)
				return e;

			e = id_next(id_twin(e));
		} while (e != et);
	}

	return HE_ISOLATED_INDEX;
}

/*template<typename V>
inline int HEMesh<V>::id_edge(const int& hedge_id) const {
	//return hedge_id / 2;
	return hedge_id >> 1;
}*/


template<typename V>
inline bool HEMesh<V>::is_removed_vert(const int& vert_id) const {
	return _vert_to_hedge[vert_id] == HE_INVALID_INDEX;
}

template<typename V>
inline bool HEMesh<V>::is_removed_hedge(const int& hedge_id) const {
	return _hedges[hedge_id].next_id == HE_INVALID_INDEX;
	//return _hedges[hedge_id] == HalfEdge();
}

/*template<typename V>
inline bool HEMesh<V>::is_removed_edge(const int& edge_id) const {
	//return is_removed_hedge(edge_id * 2);
	//return is_removed_hedge(edge_id << 1);
	return _hedges[edge_id << 1].next_id == HE_INVALID_INDEX;
}*/

// @todo:
template<typename V>
inline bool HEMesh<V>::is_begin_hedge(const int& hedge_id) const {
	return _faces.find(hedge_id) != _faces.end()
		|| _borders.find(hedge_id) != _borders.end();
	
	//return id_begin(hedge_id) == hedge_id;

	//return (_faces.find(hedge_id) != _faces.end() 
	//	|| _borders.find(hedge_id) != _borders.end())
	//	&& id_begin(hedge_id) == hedge_id;
}


template<typename V>
inline bool HEMesh<V>::is_border_hedge(const int& hedge_id) const {
	check_hedge_id(hedge_id);
	const int& b = id_begin(hedge_id);
	return _borders.find(b) != _borders.end();
}

template<typename V>
inline bool HEMesh<V>::is_border_edge(const int& hedge_id) const {
	check_hedge_id(hedge_id);
	return is_border_hedge(hedge_id) || is_border_hedge(id_twin(hedge_id));
}

template<typename V>
inline bool HEMesh<V>::is_border_vert(const int& vert_id) const {
	check_vert_id(vert_id);
	const int& begin_id = _vert_to_hedge[vert_id];
	int e = begin_id;

	do {

		//if (is_border_edge(e))
		if (is_border_hedge(e))
			return true;

		e = id_next(e);
	} while (e != begin_id);

	return false;
}

template<typename V>
inline bool HEMesh<V>::is_border_loop(const int& begin_id) const {
	//check_hedge_id(begin_id);
	//is_begin_hedge(begin_id);
	if (!is_begin_hedge(begin_id))
		throw std::invalid_argument("The given half edge id is not the begging of a face or border loop");

	return _borders.find(begin_id) != _borders.end();
}

template<typename V>
inline bool HEMesh<V>::is_n_poly(const int& hedge_id, const int& n) const {
	check_hedge_id(hedge_id);
	int e = hedge_id;

	for (int i = 0; i < n; ++i)
		e = id_next(e);

	return e == hedge_id;
}

template<typename V>
inline int HEMesh<V>::get_n_poly(const int& hedge_id) const {
	check_hedge_id(hedge_id);
	int e = hedge_id;
	int n = 0;

	do {
		++n;
		e = id_next(e);
	} while (e != hedge_id);

	return n;
}

template<typename V>
inline int HEMesh<V>::get_valence(const int& vert_id) const {
	check_vert_id(vert_id);

	const int& beg = _vert_to_hedge[vert_id];
	if (beg = HE_ISOLATED_INDEX)
		return 0;

	int e = beg;
	int n = 0;

	do {
		++n;
		e = id_next(id_twin(e));
	} while (e != beg);

	return n;
}

template<typename V>
inline bool HEMesh<V>::is_extremal_vert(const int& vert_id) const {
	check_vert_id(vert_id);
	const int& he = _vert_to_hedge[vert_id];

	if (he == HE_ISOLATED_INDEX)
		return false;

	return id_next(id_twin(he)) == he;
}

template<typename V>
inline bool HEMesh<V>::is_isolated_vert(const int& vert_id) const {
	check_vert_id(vert_id);
	return _vert_to_hedge[vert_id] == HE_ISOLATED_INDEX;
}

template<typename V>
inline bool HEMesh<V>::is_isolated_edge(const int& hedge_id) const {
	check_hedge_id(hedge_id);
	int t = id_twin(hedge_id);
	return id_next(hedge_id) == t && id_next(t) == hedge_id;
}

template<typename V>
inline bool HEMesh<V>::is_isolated_face(const int& hedge_id) const {
	check_hedge_id(hedge_id);

	const int& f = id_begin(hedge_id);
	int t = id_twin(hedge_id);
	int e = t;

	for (int e = id_next(t); e != t; e = id_next(e))
		if (id_begin(id_twin(e)) != f)
			return false;

	return true;
}

template<typename V>
inline bool HEMesh<V>::are_same_loop(const int& hedge_id1, const int& hedge_id2) const {
	check_hedge_id(hedge_id1);
	check_hedge_id(hedge_id2);
	int e = hedge_id1;

	do {
		if (e == hedge_id2)
			return true;

		e = id_next(e);
	} while (e != hedge_id1);

	return false;
}

template<typename V>
//inline void HEMesh<V>::make_face(const int& begin_hedge_id) {
inline void HEMesh<V>::make_face(const int& hedge_id) {
	//if (!is_begin_hedge(begin_hedge_id)) throw std::invalid_argument("The given half-edge should be the beginning of a face or border loop");
	check_hedge_id(hedge_id);
	int f = id_begin(hedge_id);
	auto iter = _borders.find(f);

	if (iter != _borders.end()) {
		_borders.erase(iter);
		_faces.insert(f);
	}
}

template<typename V>
//inline void HEMesh<V>::make_border(const int& begin_hedge_id) {
inline void HEMesh<V>::make_border(const int& hedge_id) {
	//if (!is_begin_hedge(begin_hedge_id)) throw std::invalid_argument("The given half-edge should be the beginning of a face or border loop");
	check_hedge_id(hedge_id);
	int f = id_begin(hedge_id);
	auto iter = _faces.find(f);

	if (iter != _faces.end()) {
		_faces.erase(iter);
		_borders.insert(f);
	}
}

// @todo: Try to do it without using std::unordered_set
template<typename V>
inline void HEMesh<V>::remove_dub_edges(int vert_id) {
	int beg = _vert_to_hedge[vert_id];

	if (beg == HE_ISOLATED_INDEX)
		return;

	std::unordered_set<int> dub_verts;
	int e = beg;

	// because it may happen that the half-edge, that the vertex points to, gets removed
	int non_removed_hedge = HE_ISOLATED_INDEX;

	do {
		if (is_removed_hedge(e))
			break;

		const int& vi = id_head(e);
		int ne = id_next(id_twin(e));

		if (dub_verts.find(vi) != dub_verts.end()) {

			const int& f = id_begin(id_twin(e));
			if (is_border_loop(f))
				make_face(f);

			remove_edge(e);
		}
		else {
			dub_verts.insert(vi);
			non_removed_hedge = e;
		}
			
		e = ne;
	} while (e != beg);

	_vert_to_hedge[vert_id] = non_removed_hedge;
}

template<typename V>
inline void HEMesh<V>::invalidate_vert(const int& vert_id) {
	_vert_to_hedge[vert_id] = HE_INVALID_INDEX;
	_garbage_vertices.insert(vert_id);
}

template<typename V>
inline void HEMesh<V>::invalidate_edge(const int& hedge_id) {
	int he = id_ltwin(hedge_id);

	_hedges[he] = HalfEdge();
	_hedges[he + 1] = HalfEdge();
	_garbage_edges.insert(he);
}

template<typename V>
inline V HEMesh<V>::get_face_center(const int& hedge_id) const {
	int e = hedge_id;
	int n = 0;
	V center;

	do {
		center += _vertices[id_head(e)];
		++n;

		e = id_next(e);
	} while (e != hedge_id);

	return center / float(n);
}


// ===================================================================================== //


template<typename V>
inline void HEMesh<V>::shift_begin(const int begin_hedge_id, const int how_much) {
	check_hedge_id(begin_hedge_id);

	if (!is_begin_hedge(begin_hedge_id))
		throw std::invalid_argument("The half-edge edge should be the begin half-edge of a face or border loop, in order to shift it");

	int nc = begin_hedge_id;
	for (int i = 0; i < how_much; ++i)
		nc = id_next(nc);

	if (nc == begin_hedge_id)
		return;

	if (_borders.find(begin_hedge_id) == _borders.end())
		_faces.insert(nc);
	else
		_borders.insert(nc);
	
	set_begin_id(nc);

	_faces.erase(begin_hedge_id);
	_borders.erase(begin_hedge_id);
}

template<typename V>
inline void HEMesh<V>::fill_holes() {
	_faces.insert(_borders.begin(), _borders.end());
	_borders.clear();
}


// @todo: Consider using a more robust triangulation method, that handles concave or non-planar polygon faces
template<typename V>
inline void HEMesh<V>::triangulate(const int& hedge_id) {
	check_hedge_id(hedge_id);

	if (is_n_poly(hedge_id, 2))
		//throw std::invalid_argument("The face should be non-degrade");
		return;

	// if it is a triangle, it will terminate at the beginning, so no need to check explicidly
	int e = id_next(id_next(hedge_id));
	while (id_next(e) != hedge_id) {
		int ne = id_next(e);

		split_face_at(hedge_id, e);

		e = ne;
	}
}

template<typename V>
inline void HEMesh<V>::triangulate() {
	//for (const int& f : _faces)
	//	triangulate(f);

	std::unordered_set<int> old_faces = _faces;
	for (const int& f : old_faces)
		triangulate(f);

}


template<typename V>
inline void HEMesh<V>::remove_internal_edges(const int& hedge_id) {
	check_hedge_id(hedge_id);

	int e = hedge_id;

	do {
		int ne = id_next(e);

		if (id_next(id_twin(e)) == e)
			remove_edge(e);

		e = ne;
	} while (e != hedge_id);

}

template<typename V>
inline void HEMesh<V>::remove_internal_edges() {
	//for (const int& f : _faces)
	//	remove_internal_edges(f);

	std::unordered_set<int> old_faces = _faces;
	for (const int& f : old_faces)
		remove_internal_edges(f);

}


template<typename V>
inline void HEMesh<V>::remove_degrade_faces() {
	for (const int& f : _faces)
		if (is_n_poly(f, 2))
			remove_edge(f);
	//remove_edge(id_next(f));
}

template<typename V>
inline void HEMesh<V>::remove_isolated_verts() {
	for (int vi = 0; vi < _vertices.size(); ++vi)
		if (_vert_to_hedge[vi] == HE_ISOLATED_INDEX) {
			invalidate_vert(vi);
		}
}


// ===================================================================================== //


template<typename V>
inline void HEMesh<V>::swap_next(const int& hedge_id1, const int& hedge_id2) {

	if (id_head(hedge_id1) != id_head(hedge_id2))
		throw std::invalid_argument("The two half-edges need to point to the same head vertex to swap their next half-edges");

	if (hedge_id1 == hedge_id2)
		return;

	bool bSameLoop = id_begin(hedge_id1) == id_begin(hedge_id2);

	if (!bSameLoop) {	// join the two loops
		// Erase the old face of hedge2 and replace it with the face of hedge1
		int fa = id_begin(hedge_id1);
		int fb = id_begin(hedge_id2);

		// makes such that if fa or fb is a border, it joins the two faces into a border
		if (_borders.find(fa) == _borders.end())
			std::swap(fa, fb);

		set_begin_id(fb, fa); // if fa and fb were swapped, then heA and heB should have also been swapped !

		_borders.erase(fb);
		_faces.erase(fb);
	}

	std::swap(_hedges[hedge_id1].next_id, _hedges[hedge_id2].next_id);

	if (bSameLoop) {	// split it into two loops
		int f = id_begin(hedge_id1);	// current loop begin
		int nf = are_same_loop(hedge_id1, f) ? hedge_id2 : hedge_id1;	// new loop begin

		if (_borders.find(f) == _borders.end()) _faces.insert(nf);
		else _borders.insert(nf);

		set_begin_id(nf);
	}

}


template<typename V>
inline int HEMesh<V>::split_face_at(int back_hedge_id, int front_hedge_id) {
	check_hedge_id(back_hedge_id);
	check_hedge_id(front_hedge_id);

	int f = id_begin(back_hedge_id);

	if (f != id_begin(front_hedge_id))
		throw std::invalid_argument("The two half-edges should be in the same face or border loop");

	int neA = new_edge(); // new half-edge
	int neB = neA + 1;	  // new twin half-edge

	//	std::cout << "Split Face At\n";

	_hedges[neA] = _hedges[front_hedge_id];
	_hedges[neB] = _hedges[back_hedge_id];

	_hedges[back_hedge_id].next_id = neA;
	_hedges[front_hedge_id].next_id = neB;

	int nf = are_same_loop(back_hedge_id, f) ? front_hedge_id : back_hedge_id;

	//	_faces.insert(nf);
	if (_borders.find(f) == _borders.end()) _faces.insert(nf);
	else _borders.insert(nf);

	set_begin_id(nf);

	return nf;
}


template<typename V>
//inline void HEMesh<V>::make_isolated(const int& vert_id) {
inline void HEMesh<V>::remove_edges(const int& vert_id) {
	check_vert_id(vert_id);

	// Remove all of its adjacent edges

	//while (!is_isolated_vert(vert_id))
	//	remove_edge(id_hedge(vert_id));

	for (int he = id_hedge(vert_id); he != HE_ISOLATED_INDEX; he = id_hedge(vert_id))
		remove_edge(he);
}


template <typename V>
inline int HEMesh<V>::add_vert(const V& v) {
	int vi;

	if (_garbage_vertices.empty()) {
		vi = _vertices.size();

		_vertices.push_back(v);
		_vert_to_hedge.push_back(HE_ISOLATED_INDEX);
	}
	else {
		auto front_iter = _garbage_vertices.begin();
		vi = *front_iter;
		_garbage_vertices.erase(front_iter);

		_vertices[vi] = v;
		_vert_to_hedge[vi] = HE_ISOLATED_INDEX;
	}

	return vi;
}


// THIS NEEDS SOME FIXING !!!
template <typename V>
inline int HEMesh<V>::add_edge(const int& tail_id, const V& head) {
	check_vert_id(tail_id);

	int head_id = add_vert(head);

	//int bFace = false;	// should it make the new edge a face or border
	int heA = _vert_to_hedge[tail_id];

	int ne = new_edge(tail_id, head_id);		// new half edge
	int net = ne + 1;	// new half edge twin
	_borders.insert(ne);

	if (heA != HE_ISOLATED_INDEX) {
		heA = id_twin(heA);	// make it point to the tail vertex

		// if its tail is an extremal vertex (adjacent only to one edge)
		//  , it doesn't need to search for a border edge
		if (id_twin(id_next(heA)) != heA) {
			int b = find_border_around_head(heA);

			if (b == HE_INVALID_INDEX)
				throw std::invalid_argument("Edges can only be added around isolated, extremal or border vertices");
		}

		//if (is_border_hedge(heA)) _borders.insert(ne);
		//else _faces.insert(ne);

		if (!is_border_hedge(heA)) {
			_borders.erase(ne);
			_faces.insert(ne);
		}

		swap_next(heA, net);
	}

	_vert_to_hedge[tail_id] = ne;
	_vert_to_hedge[head_id] = net;

	return ne;
}


// @todo: Try to optimize it, or at least organize the implementation code better
template <typename V>
inline int HEMesh<V>::add_edge(const int& tail_id, const int& head_id, const bool bMakeFace) {
	check_vert_id(tail_id);
	check_vert_id(head_id);

	// Already added edge
	int current_hedge = id_hedge(tail_id, head_id);
	if (current_hedge != HE_ISOLATED_INDEX)
		return current_hedge;

	int heA = _vert_to_hedge[tail_id];
	int heB = _vert_to_hedge[head_id];

	// First check if they share a face, and if so, just split the face with the new edge
	if (heA != HE_ISOLATED_INDEX && heB != HE_ISOLATED_INDEX) {
		int eA = id_twin(heA);
		int eB = id_twin(heB);

		int f = find_common_face_around_head(eA, eB);
		if (f != HE_INVALID_INDEX) {

			int nf = split_face_at(eA, eB);

			if (bMakeFace && is_border_loop(f)) {
				_borders.erase(f);
				_faces.insert(f);
			}

			//	if (!is_border_loop(f) && is_isolated_face(f)) {
			//		_faces.erase(f);
			//		_borders.insert(f);
			//	}

			return id_next(eA);
		}
	}

	int ne = new_edge(tail_id, head_id);  // new half-edge
	int net = ne + 1;			// new half edge twin

	//if (bMakeFace)
	//	_faces.insert(ne);
	//else
		_borders.insert(ne);

	bool isValidOper = true;

	if (heA != HE_ISOLATED_INDEX) {
		heA = id_twin(heA);	// now it points towards the vertex

		if (heA == id_twin(id_next(heA)))	// its extremal vertex
			swap_next(heA, net);
		else {
			int b = find_border_around_head(heA);		// its border vertex
			if (b != HE_INVALID_INDEX)
				swap_next(heA, net);
			else
				isValidOper = false;
		}

	}

	if (heB != HE_ISOLATED_INDEX && isValidOper) {
		heB = id_twin(heB);	// now it points towards the vertex

		// try if not this and the other, and use joinHedges not within the conditionals for better handling

		if (heB == id_twin(id_next(heB))) {	// its extremal vertex
			swap_next(heB, ne);

		}
		else {
			int b = find_border_around_head(heB);		// its border vertex
			if (b != HE_INVALID_INDEX) {
				swap_next(heB, ne);

			}
			else
				isValidOper = false;
		}
	}

	if (isValidOper) {
		_vert_to_hedge[tail_id] = ne;
		_vert_to_hedge[head_id] = net;
		return ne;
	}
	else {
		_borders.erase(ne);
		_faces.erase(ne);
		_hedges.pop_back();
		_hedges.pop_back();
		return HE_INVALID_INDEX;
	}
}


template <typename V>
inline int HEMesh<V>::add_edge_at(const int& back_hedge_id, const int& front_hedge_id, const bool bMakeFace) {
	check_hedge_id(back_hedge_id);
	check_hedge_id(front_hedge_id);

	int vb = id_head(front_hedge_id);
	int va = id_head(back_hedge_id);

	int ne = new_edge(va, vb);  // new half-edge
	int net = ne + 1;			// new half edge twin

	bool bSameLoop = id_begin(back_hedge_id) == id_begin(front_hedge_id);
	_borders.insert(ne);

	/*_faces.insert(ne);
	int f = id_begin(back_hedge_id);
	if (bMakeFace && is_border_loop(f) && f == id_begin(front_hedge_id)) {
		_borders.erase(f);
		_faces.insert(f);
	}*/

	swap_next(back_hedge_id, net);
	swap_next(front_hedge_id, ne);

	const int& f = id_begin(ne);
	if (bMakeFace && bSameLoop) { // && is_border_loop(f)
		_borders.erase(f);
		_faces.insert(f);
	}

	_vert_to_hedge[va] = ne;
	_vert_to_hedge[vb] = net;

	return ne;
}


template <typename V>
inline int HEMesh<V>::add_edge_at(const int& back_hedge_id, const V& head) {
	check_hedge_id(back_hedge_id);
	int vi = add_vert(head);

	int ne = new_edge(id_head(back_hedge_id), vi);	// new edge half-edge that points to head_id
	int net = ne + 1;	// new twin half-edge
	_faces.insert(ne);

	_vert_to_hedge[vi] = net;

	swap_next(back_hedge_id, net);
	return ne;
}


// @todo: Try not to use add_edge and constrain it to not allow (inner-edges after the operation) or (edges/faces that cross the new face)
template <typename V>
inline int HEMesh<V>::add_face(const std::vector<int>& vert_indices) {
	if (vert_indices.size() < 2ULL)
		throw std::invalid_argument("The number of new face indices cannot be smaller than 2");

	bool isEdgeAdded = false;

	for (int i = 0, n = (int)vert_indices.size() - 1ULL; i < n; ++i) {
		check_vert_id(vert_indices[i]);

		if (id_hedge(vert_indices[i], vert_indices[i + 1]) == HE_ISOLATED_INDEX) {
			add_edge(vert_indices[i], vert_indices[i + 1], true);
			isEdgeAdded = true;
		}
	}

	int e = id_hedge(vert_indices.back(), vert_indices.front());
	if (e == HE_ISOLATED_INDEX) {
		e = add_edge(vert_indices.back(), vert_indices.front(), true);
		isEdgeAdded = true;
	}

	const int& f = id_begin(e);
	if (!isEdgeAdded) {
		_borders.erase(f);
		_faces.insert(f);
	}

	return f;
}


// @todo: This method here will brake if there is already added edge in the
//		  where a new edge will be added, so try to handle that
template <typename V>
inline int HEMesh<V>::add_face_at(const std::vector<int>& hedge_indices, const bool bMakeFace) {

	split_face_at(hedge_indices.back(), hedge_indices.front());
	int prev_hedge_id = id_next(hedge_indices.back());
	//	int last_hedge = id_twin(prev_hedge_id);

	//	for (int i = 1, n = hedge_indices.size(); i < n - 1; ++i) {
	for (int i = 1, n = hedge_indices.size(); i < n; ++i) {
		const int& he = hedge_indices[i];
		split_face_at(prev_hedge_id, he);
		prev_hedge_id = id_next(prev_hedge_id);
	}

	//	split_face_at(prev_hedge_id, last_hedge);

	const int& f = id_begin(prev_hedge_id);
	if (bMakeFace && is_border_loop(f)) {
		_borders.erase(f);
		_faces.insert(f);
	}

	return f;
}


template<typename V>
inline void HEMesh<V>::remove_vert(const int& vert_id) {

	//	if (is_removed_vert(vert_id))
	//		return;

	check_vert_id(vert_id);

	// Remove all of its adjacent edges

	//while (!is_isolated_vert(vert_id))
	//	remove_edge(id_hedge(vert_id));

	for (int he = id_hedge(vert_id); he != HE_ISOLATED_INDEX; he = id_hedge(vert_id))
		remove_edge(he);

	// Remove The Vertex
	invalidate_vert(vert_id);
}

// @todo: Try to make this method to conserve begin half-edge indexes as much as possible
template<typename V>
inline void HEMesh<V>::remove_edge(const int& hedge_id) {

	//	if (is_removed_hedge(hedge_id))
	//		return;

	check_hedge_id(hedge_id);

	int c = hedge_id;
	int t = id_twin(hedge_id);
	const int& nc = id_next(hedge_id);
	const int& nt = id_next(t);
	int pc = id_prev(hedge_id);
	int pt = id_prev(t);

	// Vert to Half-Edge arrangement

	if (nt != c) {
		_hedges[pc].next_id = nt;
		_vert_to_hedge[id_head(t)] = nt;
		//_vert_to_hedge[id_tail(nt)] == nt;
	}
	else {
		_vert_to_hedge[id_head(t)] = HE_ISOLATED_INDEX;
	}

	if (nc != t) {
		_hedges[pt].next_id = nc;
		_vert_to_hedge[id_head(c)] = nc;
	}
	else {
		_vert_to_hedge[id_head(c)] = HE_ISOLATED_INDEX;
	}

	// Faces and borders arrangement

	int f = id_begin(c);
	int tf = id_begin(t);

	bool bFace = _borders.find(f) == _borders.end() && _borders.find(tf) == _borders.end();
	_faces.erase(f);
	_borders.erase(f);
	_faces.erase(tf);
	_borders.erase(tf);

	if (nt != c) {
		if (bFace) _faces.insert(nt);
		else _borders.insert(nt);

		set_begin_id(nt);
	}

	if (nc != t) {
		//if (nt != c && f == tf) {
		if (nt == c || (nt != c && f == tf)) {
			if (bFace) _faces.insert(nc);
			else _borders.insert(nc);

			set_begin_id(nc);
		}
	}

	invalidate_edge(c);
}


template<typename V>
inline void HEMesh<V>::remove_edge(const int& vertA, const int& vertB) {
	remove_edge(id_hedge(vertA, vertB));
}


template<typename V>
inline void HEMesh<V>::remove_face(const int& hedge_id) {

	check_hedge_id(hedge_id);

	int f = id_begin(hedge_id);
	int e = f;

	std::vector<int> hedgesToRemove;

	do {
		if (is_border_hedge(id_twin(e)))
			hedgesToRemove.push_back(e);

		e = id_next(e);
	} while (e != f);

	if (hedgesToRemove.empty()) {
		_faces.erase(f);
		_borders.insert(f);
	}
	else
		for (const int& he : hedgesToRemove)
			remove_edge(he);


}

template<typename V>
inline void HEMesh<V>::move_edge(const int& hedge_id, const int& front_hedge_id) {
	check_hedge_id(hedge_id);
	check_hedge_id(front_hedge_id);

	if (hedge_id == id_twin(id_next(front_hedge_id)))
		return; // then they are already connected

	int t = id_twin(hedge_id);
	const int& vi = id_head(hedge_id);

	// Disconnect previous connections, if the front of the half-edge was connected with another edge
	if (id_next(hedge_id) != t) {
		_vert_to_hedge[vi] = id_next(hedge_id);
		swap_next(hedge_id, id_prev(t));
	}
	else
		_vert_to_hedge[vi] = HE_ISOLATED_INDEX;

	_hedges[hedge_id].head_id = id_head(front_hedge_id); // Update the vertices
	swap_next(hedge_id, front_hedge_id); // Connect the hedge to the new front hedge
}

template<typename V>
inline void HEMesh<V>::move_edge(const int& hedge_id, const int& back_hedge_id, const int& front_hedge_id) {
	move_edge(hedge_id, front_hedge_id);
	move_edge(id_twin(hedge_id), back_hedge_id);
}

// ======================================================================================= //


template<typename V>
inline void HEMesh<V>::flip_edge(const int& hedge_id) {
	check_hedge_id(hedge_id);

	int c = hedge_id;			// current half-edge
	int t = id_twin(hedge_id);  // twin half-edge

	int nc = id_next(hedge_id); // next half-edge
	int nt = id_next(t);		// next of twin half-edge

	int pc = id_prev(hedge_id); // previous half-edge
	int pt = id_prev(t);		// previous of twin half-edge

	// Update face or border loops
	if (nc == id_begin(nc))
		shift_begin(nc);
	if (nt == id_begin(nt))
		shift_begin(nt);

	// Update the vert to hedges
	_vert_to_hedge[id_head(c)] = nc;
	_vert_to_hedge[id_head(t)] = nt;

	_hedges[pc].next_id = nt;
	_hedges[pt].next_id = nc;

	_hedges[c] = _hedges[nc];
	_hedges[t] = _hedges[nt];

	_hedges[nc].next_id = t;
	_hedges[nt].next_id = c;

	std::swap(_hedges[nc].begin_id, _hedges[nt].begin_id);

}


template<typename V>
inline int HEMesh<V>::refine_edge(const int& hedge_id, const V& v) {
	check_hedge_id(hedge_id);

	int c = hedge_id;
	int t = id_twin(c);

	int ne = new_edge();
	int net = id_twin(ne);

	int vi = add_vert(v);
	_vert_to_hedge[vi] = ne;
	_vert_to_hedge[id_head(hedge_id)] = net;

	_hedges[id_prev(t)].next_id = net;
	_hedges[ne] = _hedges[c];

	_hedges[c] = HalfEdge(vi, ne, id_begin(ne));
	_hedges[net] = HalfEdge(vi, t, id_begin(t));

	return vi;
}

template<typename V>
inline int HEMesh<V>::refine_edge(const int& hedge_id, const float h) {
	const V& va = _vertices[id_head(hedge_id)];	// head vertex
	const V& vb = _vertices[id_tail(hedge_id)]; // tail vertex

	return refine_edge(hedge_id, h * va + (1.0f - h) * vb);
}


template<typename V>
inline int HEMesh<V>::split_vert_to_edge(const int& hedgeA, const int& hedgeB, const V& v) {
	check_hedge_id(hedgeA);
	check_hedge_id(hedgeB);

	if (id_head(hedgeA) != id_head(hedgeB))
		throw std::invalid_argument("The two half-edges should point to the same head vertex, in order to split the vertex to an edge");

	int ne = new_edge();	// new edge half-edge index
	int net = id_twin(ne);	// new edge twin half-edge index

	int vi = add_vert(v);	// new vertex index
	_vert_to_hedge[vi] = net;
	_vert_to_hedge[id_head(hedgeA)] = ne;

	_hedges[ne] = _hedges[hedgeA];
//	_hedges[ne].head_id = vi;
	_hedges[hedgeA].next_id = ne;

	_hedges[net] = _hedges[hedgeB];
	_hedges[hedgeB].head_id = vi;
	_hedges[hedgeB].next_id = net;

	set_vert_id(ne, vi);	// !!!

	return ne;
}

template<typename V>
inline int HEMesh<V>::split_vert_to_faces(const int& hedgeA, const int& hedgeB, const V& v) {
	
	int ne = split_vert_to_edge(hedgeA, hedgeB, v);
	int net = id_twin(ne);

	split_face_at(ne, id_prev(hedgeA));
	split_face_at(net, id_prev(hedgeB));

	make_face(ne);
	make_face(net);

	return ne;
}

template<typename V>
inline int HEMesh<V>::split_edge(const int& hedge_id, const V& v) {

	int vi = refine_edge(hedge_id, v);

	int t = id_twin(hedge_id);
	int ne = id_next(hedge_id); // new half-edge after the refinement

	if (!is_border_hedge(hedge_id))
		split_face_at(hedge_id, id_next(ne));

	if (!is_border_hedge(t))
		split_face_at(id_twin(ne), id_next(t));

	return vi;
}

template<typename V>
inline int HEMesh<V>::split_edge(const int& hedge_id, const float h) {
	const V& va = _vertices[id_head(hedge_id)];	// head vertex
	const V& vb = _vertices[id_tail(hedge_id)]; // tail vertex

	return split_edge(hedge_id, h * va + (1.0f - h) * vb);
}

template<typename V>
inline int HEMesh<V>::clip_corner(const int& hedge_id) {
	check_hedge_id(hedge_id);

	//if (get_n_poly(hedge_id) <= 3)
	if (is_n_poly(hedge_id, 2) || is_n_poly(hedge_id, 3))
		//throw std::invalid_argument("Cannot clip the corner of a degrade or triangle face");
		return HE_INVALID_INDEX;


	return split_face_at(id_prev(hedge_id), id_next(hedge_id));
}


template <typename V>
inline int HEMesh<V>::collapse_edge(const int& hedge_id, const V& v) {
	check_hedge_id(hedge_id);

	int c = hedge_id;
	int t = id_twin(c);
	const int& nc = id_next(c);
	const int& nt = id_next(t);

	int vt = id_head(t);	// the vertex that we collapse the edge into
	int vc = id_head(c);	// the deleted/removed vertex after the collapse
	_vertices[vt] = v;

	// Checking if it has adjacent faces (that are not the edge itself) and if so, are those faces degrade (2-sided), then it is a degrade face
	// if (!is_isolated_edge(c) && is_n_poly(c, 2))
	if (nc != t && id_next(nc) == c)
		remove_edge(nc);
	if (nt != c && id_next(nt) == t)
		remove_edge(nt);

	// Remove the faces, if isolated or shift the face beging hedges, if it happens that it has adjacent edges
	// and one of the adjacent faces begin half edge is in the collapsed edge
	if (is_isolated_edge(c)) {
		_faces.erase(id_begin(c));
		_borders.erase(id_begin(c));
	}
	else {
		if (id_begin(t) == t) {
			int by = nt != c ? 1 : 2;
			shift_begin(t, by);
		}

		if (id_begin(c) == c) {
			int by = nc != t ? 1 : 2;
			shift_begin(c, by);
		}
	}

	int pc = id_prev(c);
	int pt = id_prev(t);

	// Rearrange the adjacent of the collapsed edge half-edges and update their head vertex 
	if (nc != t) {
		set_vert_id(id_twin(nc), vt);
		_hedges[pt].next_id = nt != c ? nt : nc;
	}
	if (nt != c) {
		_hedges[pc].next_id = nc != t ? nc : nt;
		_vert_to_hedge[vt] = id_next(t);
	}
	else
		_vert_to_hedge[vt] = nc != t ? nc : HE_ISOLATED_INDEX; // !!!

	// Remove dublicate edges
	remove_dub_edges(vt);

	invalidate_edge(c);
	invalidate_vert(vc);

	return vt;
}

template<typename V>
inline int HEMesh<V>::collapse_edge(const int& hedge_id, const float h) {
	const V& va = _vertices[id_head(hedge_id)];	// head vertex
	const V& vb = _vertices[id_tail(hedge_id)]; // tail vertex

	return collapse_edge(hedge_id, h * va + (1.0f - h) * vb);
}


// @todo: FIX! It works most of the time, but on a long border loop, it breaks
template<typename V>
inline int HEMesh<V>::collapse_face(const int& hedge_id, const V& v) {
	check_hedge_id(hedge_id);

	int c = hedge_id;
	int f = id_begin(c);

	// Collect the face half-edges and vertices to remove
	std::vector<int> face_hedges;
	std::vector<int> face_verts;
	int e = c;
	do {
		face_hedges.push_back(e);
		face_verts.push_back(id_head(e));
		e = id_next(e);
	} while (e != c);


	// That is the vertex that we collapse the face into
	//const int& vi = face_verts.front();
	int vi = face_verts.front();
	_vertices[vi] = v;
	//_vert_to_hedge[vi] = id_next(id_twin(id_next(c)));
	_vert_to_hedge[vi] = HE_ISOLATED_INDEX; // half-edge which edge is not on the face

	// Remove face edges
	for (int& he : face_hedges) {
		int temp = id_prev(id_twin(he));

		int oe = id_twin(temp);
		if (oe != id_next(he)) // then it is an edge that is not on the face
			_vert_to_hedge[vi] = oe;
		
		remove_edge(he);
		he = temp;
	}

	f = id_begin(face_hedges.front());	// any begin_id of the hedges in face_hedges
	if (_vert_to_hedge[vi] == HE_ISOLATED_INDEX) { // then the face is isolated, so remove it
		_faces.erase(f);
		_borders.erase(f);
	}

	// Make all the edges that were adjacent to the face to point to the collapse vertex
	for (int i = 1, n = face_hedges.size(); i < n; ++i) {
		const int& he = face_hedges[i];

		invalidate_vert(face_verts[i]);

		if (!is_removed_hedge(he))
			set_vert_id(he, vi);
	}


	std::vector<int> non_removed_hedges;

	for (const int& he : face_hedges)
		if (!is_removed_hedge(he))
			non_removed_hedges.push_back(he);

	face_hedges = std::move(non_removed_hedges);


	// Now there is one star shaped face that we need to separate into 2poly faces
	int temp_hedge = face_hedges.front();
	for (auto it = ++face_hedges.begin(); it != face_hedges.end(); ++it)
		swap_next(temp_hedge, *it);

	// Now we can safely remove any dublicated head_id  edges
	remove_dub_edges(vi);

	return vi;
}


template<typename V>
inline int HEMesh<V>::collapse_face(const int& hedge_id) {
	check_hedge_id(hedge_id);

	V center = get_face_center(hedge_id);
	return collapse_face(hedge_id, center);
}


/*template<typename V>
inline int HEMesh<V>::bevel_vert(const int& vert_id, const float h) {
	check_vert_id(vert_id);

	if (is_isolated_vert(vert_id))
		throw std::invalid_argument("Cannot bevel an isolated vertex");

	V center = _vertices[vert_id];

	// Collect the adjacent outgoing half-edges (one per each adjacent edge)
	int beg = id_twin(id_hedge(vert_id));
	int e = beg;
	std::vector<int> adj_hedges; // adjacent half-edges that point to the vertex
	do {
		adj_hedges.push_back(e);
		e = id_twin(id_next(e));
	} while (e != beg);

	bool bFace = false;

	// Make all adjacent faces to be one single face
	for (const int& he : adj_hedges) {
		_hedges[he].next_id = id_twin(he);

		const int& f = id_begin(he);
		auto fit = _faces.find(f);
		if (fit != _faces.end()) {
			_faces.erase(fit);
			bFace = true;
		}
		else
			_borders.erase(f);
	}

	const int& f = adj_hedges.front();
	if (bFace) _faces.insert(f);
	else _borders.insert(f);
	set_begin_id(f);
	
	// Make all edge (except the first) to point to a new vertex
	for (int i = 1, n = adj_hedges.size(); i < n; ++i) {
		const int& he = adj_hedges[i]; // adjacent half-edge that points to the center vertex
		int t = id_twin(he);
		const V& v = _vertices[id_head(t)];	// an adjacent vertex to the center vertex

		int nv = add_vert(h * center + (1.0f - h) * v);
		_hedges[he].head_id = nv;
		_vert_to_hedge[nv] = t;
	}

	const V& v = _vertices[id_tail(adj_hedges.front())];
	_vertices[vert_id] = h * center + (1.0f - h) * v;
	
	// The orientation of the face is opposite to the direction 
	//  that the adjacent vertex edges are iterated
	std::reverse(adj_hedges.begin(), adj_hedges.end());

	return add_face_at(adj_hedges);
}*/


// @todo: Think about how vertices with only 1 or 2 adjacent edges, should be handled
// This version preserves the faces and borders
// The older version would have made all adjacent faces/borders into
// either a face or a border, but this one makes such that if an
// adjacent face was face/border it stays as a face/border after the beveling
template<typename V>
inline int HEMesh<V>::bevel_vert(const int& vert_id, const float h) {
	check_vert_id(vert_id);

	if (is_isolated_vert(vert_id))
		throw std::invalid_argument("Cannot bevel an isolated vertex");

	// Collect the adjacent outgoing half-edges
	//  (the half-edges, which tail vertex is vert_id)
	int beg = id_hedge(vert_id);
	int e = beg;
	std::vector<int> adj_hedges;
	do {
		adj_hedges.push_back(e);
		e = id_next(id_twin(e));
	} while (e != beg);
	int n = adj_hedges.size();

	// Collect the new vertices and the one that will be reused
	std::vector<int> adj_verts;
	for (int i = 0; i < n - 1; ++i)
		adj_verts.push_back(add_vert(V()));
	adj_verts.push_back(vert_id);

	// Create the new edges, update their connectivity with adj_hedges
	//  , update vertices blended position and update the adj_hedges
	//  to store the new hedges that are for the new face
	for (int i = 0; i < n; ++i) {
		int e = adj_hedges[i]; // points out of vert_id
		int t = id_twin(e); // points towards vert_id

		// Update the new blended vertices
		const V& ov = _vertices[id_head(e)];
		const int& nvi = adj_verts[(i + 1) % n];
		const int& vi = adj_verts[i];
		_vertices[vi] = h * ov + (1.0f - h) * _vertices[vert_id];
		_vert_to_hedge[vi] = e;

		// Create a new edge that will be on the new face
		int ne = new_edge();	// new half-edge
		int net = ne + 1;		// new twin half-edge
		_hedges[ne] = HalfEdge(nvi, id_next(t), id_begin(t)); // the one that is out of the new face
		_hedges[net] = HalfEdge(vi, HE_INVALID_INDEX, HE_INVALID_INDEX); // the one that is in the new face

		// Update the adjacent hedge to now point to the new vertex and new half-edge
		_hedges[t].head_id = vi;
		_hedges[t].next_id = ne;

		adj_hedges[i] = net;	// reuse the adj_edges vector to store the new hedges

	}

	const int& f = adj_hedges.front();
	_faces.insert(f);

	// Now that adj_edges stores the new hedges, we can update them, so they all connect
	//  and form the loop of the new face
	for (int i = 0; i < n; ++i) {

		// reverse the order that we set the next_ids, because the order at which
		//  adjacent hedges are iterated around the vertex
		//  is opposite to the orientation of the new face
		const int& e = adj_hedges[(i + 1) % n];
		const int& ne = adj_hedges[i];

		_hedges[e].begin_id = f;
		_hedges[e].next_id = ne;
	}

	return f;
}

template<typename V>
inline int HEMesh<V>::bevel_edge(const int& hedge_id, const float h) {
	check_hedge_id(hedge_id);

	if (is_isolated_edge(hedge_id))
		throw std::invalid_argument("Cannot bevel an isolated edge");
	
	int t = id_twin(hedge_id);
	int e;
	std::vector<int> front_hedges;
	std::vector<int> back_hedges;

	for (int e = id_next(hedge_id); e != t; e = id_next(id_twin(e)))
		front_hedges.push_back(e);
	for (int e = id_next(t); e != hedge_id; e = id_next(id_twin(e)))
		back_hedges.push_back(e);

	split_face_at(hedge_id, id_prev(hedge_id));
	_vertices.reserve(_vertices.size() + front_hedges.size() + back_hedges.size());

	int ne = hedge_id;
	t = id_next(hedge_id);
	
	for (auto it = front_hedges.begin(); it != front_hedges.end(); ++it) {
		const int& he = *it;
		V cv = _vertices[id_head(ne)];
		const V& ov = _vertices[id_head(he)];

		_vertices[id_head(ne)] = h * ov + (1.0f - h) * cv;

		if (std::next(it) != front_hedges.end())
			split_vert_to_edge(id_twin(he), ne, cv);
	}
	
	ne = t;

	for (auto it = back_hedges.begin(); it != back_hedges.end(); ++it) {
		const int& he = *it;
		V cv = _vertices[id_head(ne)];
		const V& ov = _vertices[id_head(he)];

		_vertices[id_head(ne)] = h * ov + (1.0f - h) * cv;

		if (std::next(it) != back_hedges.end())
			split_vert_to_edge(id_twin(he), ne, cv);
	}

	return id_begin(hedge_id);
}


template<typename V>
inline int HEMesh<V>::bevel_face(const int& hedge_id, const float h) {
	check_hedge_id(hedge_id);

	// Calculate mid-point/center
	V center = get_face_center(hedge_id);
	int e = hedge_id;

	// Create the new edges that connect the old face and the new inner face
	std::vector<int> new_hedges;
	do {
		int ne = id_next(e);

		const V& v = _vertices[id_head(e)];
		new_hedges.push_back(add_edge_at(e, h * center + (1.0f - h) * v));

		e = ne;
	} while (e != hedge_id);
	
	return add_face_at(new_hedges, false);
}


template<typename V>
inline int HEMesh<V>::cut_edge(const int& hedge_id) {
	check_hedge_id(hedge_id);

	if (is_border_edge(hedge_id))
		//throw std::invalid_argument("Cannot cut around a border (where it has already been cut)");
		return HE_INVALID_INDEX;

	split_face_at(hedge_id, id_prev(hedge_id));
	int f = id_begin(hedge_id);
	_faces.erase(f);
	_borders.insert(f);

	int hedge_ids[2] = { hedge_id, id_next(hedge_id) };

	for (int i = 0; i < 2; ++i) {
		const int& he = hedge_ids[i];

		int bh = id_twin(id_next(he)); // border half-edge index  // one step ahead, so find_border_around_head doesnt always stop at he
		int b = find_border_around_head(bh); // border begin half-edge index

		// it will always be a valid index because he is a border half-edge
		if (bh != he) { // && b != HE_INVALID_INDEX 
			int vi = id_head(he);
			_vert_to_hedge[vi] = id_twin(he);

			// split the vertex when two borders are joined
			int nvi = _vertices.size();
			_vertices.push_back(_vertices[vi]);
			_vert_to_hedge.push_back(id_twin(bh));

			swap_next(bh, he); // join the borders
			set_vert_id(bh, nvi);
		}
	}

	//return id_begin(hedge_id);
	return id_begin(hedge_ids[1]);	// better return the other border, if it happens that the other hedge is now separated from hedge_id
}



// ===================================================================================== //


template <typename V>
inline typename HEMesh<V>::vert_iter HEMesh<V>::vert(const int& vert_id) const {
	check_vert_id(vert_id);
	return vert_iter(vert_id, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::hedge_iter HEMesh<V>::hedge(const int& hedge_id) const {
	check_hedge_id(hedge_id);
	return hedge_iter(hedge_id, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::hedge_iter HEMesh<V>::hedge(const int& tail_id, const int& head_id) const {
	check_vert_id(tail_id);
	check_vert_id(head_id);
	return hedge_iter(id_hedge(tail_id, head_id), const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::edge_iter HEMesh<V>::edge(const int& edge_id) const {
	//check_edge_id(edge_id);
	check_hedge_id(edge_id << 1);

	//	return edge_iter(edge_id, const_cast<HEMesh<V>*>(this));
	return edge_iter(edge_id << 1, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::face_iter HEMesh<V>::face(const int& hedge_id) const {
	check_hedge_id(hedge_id);
	return face_iter(hedge_id, const_cast<HEMesh<V>*>(this));
}


template <typename V>
inline typename HEMesh<V>::vert_iter HEMesh<V>::begin_verts() const {
	int id = 0;

	while (id < (int)_vertices.size() && is_removed_vert(id))
		++id;

	return vert_iter(id, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::hedge_iter HEMesh<V>::begin_hedges() const {
	int id = 0;

	while (id < (int)_hedges.size() && is_removed_hedge(id))
		++id;

	return hedge_iter(id, const_cast<HEMesh<V>*>(this));
}

/*template <typename V>
inline typename HEMesh<V>::edge_iter HEMesh<V>::begin_edges() const {
	int id = 0;

	while ((2 * id < _hedges.size()) && is_removed_edge(id))
		++id;

	return edge_iter(id, const_cast<HEMesh<V>*>(this));
}*/

template <typename V>
inline typename HEMesh<V>::edge_iter HEMesh<V>::begin_edges() const {
	int hedge_id = 0;

	while (hedge_id < (int)_hedges.size() && is_removed_hedge(hedge_id))
		hedge_id += 2;

	//	return edge_iter(hedge_id >> 1, const_cast<HEMesh<V>*>(this));
	return edge_iter(hedge_id, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::face_iter HEMesh<V>::begin_faces() const {
	return face_iter(_faces.begin(), const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::face_iter HEMesh<V>::begin_borders() const {
	return face_iter(_borders.begin(), const_cast<HEMesh<V>*>(this));
}


template <typename V>
inline typename HEMesh<V>::vert_iter HEMesh<V>::end_verts() const {
	return vert_iter(_vertices.size(), const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::hedge_iter HEMesh<V>::end_hedges() const {
	return hedge_iter(_hedges.size(), const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::edge_iter HEMesh<V>::end_edges() const {
	//return edge_iter(_hedges.size() / 2, const_cast<HEMesh<V>*>(this));
//	return edge_iter(_hedges.size() >> 1, const_cast<HEMesh<V>*>(this));

	return edge_iter(_hedges.size(), const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::face_iter HEMesh<V>::end_faces() const {
	return face_iter(_faces.end(), const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::face_iter HEMesh<V>::end_borders() const {
	return face_iter(_borders.end(), const_cast<HEMesh<V>*>(this));
}


/*template <typename V>
inline typename HEMesh<V>::vert_iter HEMesh<V>::rbegin_verts() const {
	int id = (int)_vertices.size() - 1;

	while (id >= 0 && is_removed_vert(id))
		--id;

	return vert_iter(id, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::hedge_iter HEMesh<V>::rbegin_hedges() const {
	int id = (int)_hedges.size() - 1;

	while (id >= 0 && is_removed_hedge(id))
		--id;

	return hedge_iter(id, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::edge_iter HEMesh<V>::rbegin_edges() const {
	int id = int(_hedges.size()) / 2 - 1;

	while (id >= 0 && is_removed_edge(id))
		--id;

	return edge_iter(id, const_cast<HEMesh<V>*>(this));
}


template <typename V>
inline typename HEMesh<V>::vert_iter HEMesh<V>::rend_verts() const {
	return vert_iter(-1, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::hedge_iter HEMesh<V>::rend_hedges() const {
	return hedge_iter(-1, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::edge_iter HEMesh<V>::rend_edges() const {
	return edge_iter(-1, const_cast<HEMesh<V>*>(this));
}*/


template <typename V>
inline Iterable<typename HEMesh<V>::vert_iter> HEMesh<V>::verts() const {
	return Iterable<vert_iter>(begin_verts(), end_verts());
}

template <typename V>
inline Iterable<typename HEMesh<V>::hedge_iter> HEMesh<V>::hedges() const {
	return Iterable<hedge_iter>(begin_hedges(), end_hedges());
}

template <typename V>
inline Iterable<typename HEMesh<V>::edge_iter> HEMesh<V>::edges() const {
	return Iterable<edge_iter>(begin_edges(), end_edges());
}

template <typename V>
inline Iterable<typename HEMesh<V>::face_iter> HEMesh<V>::faces() const {
	return Iterable<face_iter>(begin_faces(), end_faces());
}

template <typename V>
inline Iterable<typename HEMesh<V>::face_iter> HEMesh<V>::borders() const {
	return Iterable<face_iter>(begin_borders(), end_borders());
}


/*template <typename V>
inline Iterable<typename HEMesh<V>::vert_iter> HEMesh<V>::rvertices() const {
	return Iterable<face_iter>(rbegin_verts(), rend_verts());
}

template <typename V>
inline Iterable<typename HEMesh<V>::vert_iter> HEMesh<V>::rhedges() const {
	return Iterable<face_iter>(rbegin_hedges(), rend_hedges());
}

template <typename V>
inline Iterable<typename HEMesh<V>::vert_iter> HEMesh<V>::redges() const {
	return Iterable<face_iter>(rbegin_edges(), rend_edges());
}*/


#endif