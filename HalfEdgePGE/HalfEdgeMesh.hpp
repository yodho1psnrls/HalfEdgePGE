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
	* Stores face and border loops, but the faces and borders doesnt have indexes, the index they have is the index of the half edge at the beginning of the loop
	* Most methods work with vertex indexes and half edge indexes, even if the method does something
	  on an edge, face or border, this is because for most face and border operations you need any half edge
	  that lies in the same loop and for the edge operations i dont want the overhead of constantly
	  converting from half edge index to edge index
	* A vertex is marked as removed when _vert_to_hedge[vi] == HE_INVALID_INDEX and its index vi is in _garbage_vertices
	* A half edge is marked as removed when its head_id, next_id and begin_id are set to HE_INVALID_INDEX and its index is in _garbage_hedges
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
	std::set<int> _garbage_hedges{};
	
	//bool _are_verts_shrinked;
	//bool _are_hedges_shrinked;

	bool is_in_range_vert(const int& vert_id) const;
	bool is_in_range_hedge(const int& hedge_id) const;
	//bool is_in_range_edge(const int& edge_id) const;

	void check_vert_id(const int& vert_id) const;
	void check_hedge_id(const int& hedge_id) const;
	//void check_edge_id(const int& edge_id) const;


	void set_begin_id(const int& hedge_id); // sets the begin_id of the given half-edge to the given half-edge
	void set_begin_id(const int& hedge_id, const int& begin_id); // sets the begin half-edge index in the loop of the given half-edge to the given begin half-edge index
	//void shift_begin(const int& begin_id, const int& how_much = 1);

	// circulates the two half edges around their head vertex until it finds a common face 
	//  or border loop, but prioritizes face loops, so if it finds both a common face
	//  and a common border, it will choose the face, returns the begin half-edge index of the loop
	int find_common_face_around_head(int& he1, int& he2);
	
	// circulates the half-edge around its head vertex, until it finds a border, 
	//  if no border is found returns HE_INVALID_INDEX and the given half_edge stays the same
	int find_border_around_head(int& hedge_id);

	// gives the smaller/even half-edge index of the edge of the given half-edge 
	int id_ltwin(const int& hedge_id) const;

	// gives the bigger/odd half-edge index of the edge of the given half-edge 
	int id_rtwin(const int& hedge_id) const;

	friend class vert_iter;
	friend class hedge_iter;
	friend class edge_iter;
	friend class face_iter;

public:

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
	const int& id_hedge(const int& vert_id) const;
	int id_hedge(const int& tail_id, const int& head_id) const;	// gives the id of the half edge with the given tail and head vertices, returns HE_ISOLATED_INDEX if not found
	//int id_edge(const int& hedge_id) const;

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

	void fill_holes();	// all border loops are made face loops
	void triangulate(const int& hedge_id);	// triangulates the face loop of the given half edge
	void triangulate();		// triangulates the whole mesh
	void remove_internal_edges(const int& hedge_id);	// if the face or border loop of the given half edge has any internal edges (that may appear after some operations), it removes them
	void remove_internal_edges();	// removes internal edges of all face or border loops in the mesh

	void remove_degrade_faces() const;	// removes all faces with size 2 (faces that only have 2 edges)
	void remove_isolated_verts() const;	// clears the mesh of all vertices that doesnt share any edge

	//bool do_share_loop(const int& hedge_id1, const int& hedge_id2) const; // true if the two half-edges share the same face or border loop
	bool are_same_loop(const int& hedge_id1, const int& hedge_id2) const; // true if the two half-edges share the same face or border loop


	// Operations

	void swap_next(const int& hedge_id1, const int& hedge_id2); // swaps the next indexes of the two half-edges and rearanges the face or border begin hedge ids properly
	int split_face_at(const int& hedge_id1, const int& hedge_id2);	// splits the face or border loop adding a new edge at the head vertices of the given half-edges

	int add_vert(const V& v);
	int add_edge(const int& tail_id, const V& head);	// returns the index of the new half_edge which head vertex is the id of the new added head vertex
	int add_edge(const int& tail_id, const int& head_id, const bool make_face = false);	// if make_face is true, it makes a face if the two given vertices had shared a common border loop
	int add_edge_at(const int& prev_id, const int& next_id);
	int add_face(const std::vector<int>& indices);

	void remove_vert(const int& vert_id);
	void remove_edge(const int& hedge_id);
	void remove_face(const int& hedge_id);

	void move_edge(const int& hedge_id, const int& new_next_id);
	void move_edge(const int& hedge_id, const int& new_prev_id, const int& new_next_id);


	// With Iterators

	void swap_next(const hedge_iter& hedge1, const hedge_iter& hedge2); // swaps the next indexes of the two hedges and rearanges the face or border begin hedge ids properly

	//vert_iter add_vert(const V& v);
	int add_edge(const vert_iter& tail, const V& head);
	int add_edge(const vert_iter& tail, const vert_iter& head);
	int add_edge_at(const hedge_iter& prev, const hedge_iter& next);
	int add_face(const std::vector<vert_iter>& indices);

	void remove_vert(const vert_iter& vert);
	void remove_edge(const hedge_iter& hedge);
	void remove_face(const hedge_iter& hedge);

	void move_edge(const hedge_iter& hedge, const hedge_iter& new_next);
	void move_edge(const hedge_iter& hedge, const hedge_iter& new_prev, const hedge_iter& new_next);



  // =================================== ITERATORS ========================================= //


	// Convertion from index to iterators
	vert_iter vert(const int& vert_id) const;
	hedge_iter hedge(const int& hedge_id) const;
	hedge_iter hedge(const int& tail_id, const int& head_id) const;
	edge_iter edge(const int& edge_id) const;
	//	face_iter face(const int& hedge_id) const;


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

		//operator const int& () const { return id; }
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
		vert_iter& operator--() { do { --id; } while (id >= 0 && hm->is_removed_vert(id)); return *this; }	// Prefix increment
		vert_iter operator--(int) { vert_iter temp = *this; --(*this); return temp; }	// Postfix increment

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

		//operator edge_iter() const { return edge_iter(id / 2, hm); }
		operator edge_iter() const { return edge_iter(id >> 1, hm); }

		//operator const int& () const { return id; }
		//operator int() const { return id; }
		const int& index() const { return id; }

		bool is_removed() const { return hm->is_removed_hedge(id); }
		bool is_border() const { return hm->is_border_hedge(id); }
		bool is_isolated() const { return hm->is_isolated_edge(id); } // true if it is adjacent only to its twin half edge and no other half edges

		hedge_iter next() const { return hedge_iter(hm->_hedges[id].next_id, hm); }
		hedge_iter twin() const { return hedge_iter(hm->id_twin(id), hm); }
		hedge_iter prev() const { return hedge_iter(hm->id_prev(id), hm); }
		//edge_iter edge() const { return edge_iter(id / 2, hm); }
		edge_iter edge() const { return edge_iter(id >> 1, hm); }
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
		hedge_iter& operator--() { do { --id; } while (id >= 0 && hm->is_removed_hedge(id)); return *this; } // Prefix increment
		hedge_iter operator--(int) { hedge_iter temp = *this; --(*this); return temp; } // Postfix increment

		friend bool operator== (const hedge_iter& a, const hedge_iter& b) { return a.id == b.id; };
		friend bool operator!= (const hedge_iter& a, const hedge_iter& b) { return a.id != b.id; };

	};


	class edge_iter {
	protected:

		HEMesh<V>* hm;
		int hedge_id;

		friend class HEMesh<V>;
		friend class HEMesh<V>::hedge_iter;

		edge_iter(const int& edge_id, HEMesh<V>* const& hmm) : hedge_id(edge_id << 1), hm(hmm) {}

	public:

		//edge_iter() : id(HE_INVALID_INDEX), hm(nullptr) {}
		edge_iter() : hedge_id(HE_INVALID_INDEX), hm(nullptr) {}

		//operator const int& () const { return id; }
		//operator int() const { return id; }
		const int& index() const { return hedge_id >> 1; }

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
		edge_iter& operator++() { do { hedge_id += 2; } while (hedge_id < hm->edges_size() && hm->is_removed_hedge(hedge_id)); return *this; } // Prefix increment
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

		operator hedge_iter() const { return hedge(); }

		hedge_iter hedge() const { return hedge_iter(*iter, hm); }

		bool is_valid() const { return hm->is_begin_hedge(*iter); }
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
	Mesh<V> m;

	m.vertices = _vertices;
	m.indices.reserve(_hedges.size());
	//m.indices.reserve(_faces.size() * size_per_face[i]);

	int he[3];

	for (const int& f : _faces) {

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
	_garbage_hedges.clear();
}


template<typename V>
inline bool HEMesh<V>::empty() const {
	return _faces.empty() && _borders.empty();
		//&& _vertices.size() == _garbage_vertices.size()
		//&& 2ULL * _hedges.size() == _garbage_hedges.size()
}


template<typename V>
inline bool HEMesh<V>::are_verts_shrinked() const {
	return _garbage_vertices.empty();
}

template<typename V>
inline bool HEMesh<V>::are_edges_shrinked() const {
	return _garbage_hedges.empty();
}


template<typename V>
inline int HEMesh<V>::verts_size() const {
	return _vertices.size() - _garbage_vertices.size();
}

template<typename V>
inline int HEMesh<V>::hedges_size() const {
	return _hedges.size() - _garbage_hedges.size() * 2ULL;
}

template<typename V>
inline int HEMesh<V>::edges_size() const {
	return _hedges.size() / 2ULL - _garbage_hedges.size();
}

template<typename V>
inline int HEMesh<V>::faces_size() const{
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
 
			if (id_face(he1) == id_face(he2))
				return id_face(he1);

			he2 = id_next(id_twin(he2));
		} while (he2 != beg2);

		he1 = id_next(id_twin(he1));
	} while (he1 != beg1);

	return HE_INVALID_INDEX;
}

template<typename V>
inline const int& HEMesh<V>::id_hedge(const int& vert_id) const {
	return _vert_to_hedge[vert_id];
}

template<typename V>
inline int HEMesh<V>::id_hedge(const int& tail_id, const int& head_id) const {
	const int& et = _vert_to_hedge[tail_id];
	const int& eh = _vert_to_hedge[head_id];
	int e;

	if (et != HE_ISOLATED_INDEX && eh != HE_ISOLATED_INDEX) {
		e = et;
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
}

/*template<typename V>
inline bool HEMesh<V>::is_removed_edge(const int& edge_id) const {
	//return is_removed_hedge(edge_id * 2);
	//return is_removed_hedge(edge_id << 1);
	return _hedges[edge_id << 1].next_id == HE_INVALID_INDEX;
}*/

template<typename V>
inline bool HEMesh<V>::is_begin_hedge(const int& hedge_id) const {
	return _faces.find(hedge_id) != _faces.end() || _borders.find(hedge_id) != _borders.end();
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
inline void HEMesh<V>::fill_holes() {
	_faces.insert(_borders.begin(), _borders.end());
	_borders.clear();
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

  // ===================================================================================== //


template<typename V>
inline void HEMesh<V>::swap_next(const int& hedge_id1, const int& hedge_id2) {

	if (id_head(hedge_id1) != id_head(hedge_id2))
		throw std::invalid_argument("The two half-edges need to point to the same head vertex to swap their next half-edges");
	
	if (hedge_id1 == hedge_id2)
		return;

	bool bSameLoop = id_begin(hedge_id1) == id_begin(hedge_id2);

	if (!bSameLoop) {	// join the two faces
		// Erase the old face of hedge2 and replace it with the face of hedge1
		int fa = id_begin(hedge_id1);
		int fb = id_begin(hedge_id2);

		// if fa is face, then swap them so fb is face (which will be erased and the border will win)
		if (_borders.find(fa) == _borders.end())
			//if (_borders.find(fb) == _borders.end())
			std::swap(fa, fb);

		//setFaceID(heB, fa); // if fa and fb were swapped, then heA and heB should have also been swapped !
		set_begin_id(fb, fa);

		_borders.erase(fb);
		_faces.erase(fb);
	}

	std::swap(_hedges[hedge_id1].next_id, _hedges[hedge_id2].next_id);

	// YOU SHOULD HANDLE THE CASE WHERE THEY SHARE A FACE, because then you need to create a new face
	// If same face, then split it into two faces
	if (bSameLoop) {
		int f = id_begin(hedge_id1);	// face
		int nf = are_same_loop(hedge_id1, f) ? hedge_id2 : hedge_id1;	// new face

		//_faces.insert(nf);
		if (_borders.find(f) == _borders.end()) _faces.insert(nf);
		else _borders.insert(nf);

		set_begin_id(nf);
	}

}

template<typename V>
inline int HEMesh<V>::split_face_at(const int& hedge_id1, const int& hedge_id2) {
	check_hedge_id(hedge_id1);
	check_hedge_id(hedge_id2);

	int f = id_begin(hedge_id1);

	if (f != id_begin(hedge_id2))
		throw std::invalid_argument("The two half-edges should be in the same face or border loop");

	int neA = _hedges.size();
	int neB = neA + 1;

	_hedges.push_back(_hedges[hedge_id2]);
	_hedges.push_back(_hedges[hedge_id1]);

	_hedges[hedge_id1].next_id = neA;
	_hedges[hedge_id2].next_id = neB;

	int nf = are_same_loop(hedge_id1, f) ? hedge_id2 : hedge_id1;
	if (_borders.find(f) == _borders.end()) _faces.insert(nf);
	else _borders.insert(nf);
	set_begin_id(nf);

	return nf;
}

template <typename V>
inline int HEMesh<V>::add_vert(const V& v) {
	int vi = _vertices.size();
	
	_vertices.push_back(v);
	_vert_to_hedge.push_back(HE_ISOLATED_INDEX);

	return vi;
}


// THIS NEEDS SOME FIXING !!!
template <typename V>
inline int HEMesh<V>::add_edge(const int& tail_id, const V& head) {
	check_vert_id(tail_id);

	int nv = add_vert(head);

	int ne = _hedges.size();	// new half edge
	int net = ne + 1;			// new half edge twin

	_hedges.push_back(HalfEdge(nv, net, net));
	_hedges.push_back(HalfEdge(tail_id, ne, net));
	_borders.insert(net);

	int heA = _vert_to_hedge[tail_id];
	if (heA != HE_ISOLATED_INDEX) {
		heA = id_twin(heA);	// make it point to vertA
		int b = find_border_around_head(heA);

		if (b == HE_INVALID_INDEX)
			throw std::invalid_argument("Edges can only be added around isolated or border vertices");

		if (!is_border_hedge(heA)) {
		//if (b == HE_INVALID_INDEX)
			_borders.erase(net);
			_faces.insert(net);
		}

		swap_next(heA, net);
	}

	return ne;
}


// @todo: Try to optimize it
template <typename V>
inline int HEMesh<V>::add_edge(const int& tail_id, const int& head_id, const bool make_face) {
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

			int ne = split_face_at(eA, eB);

			//if (make_face) {
			if (make_face && is_border_loop(f)) {
				int f = id_begin(ne);
				//int f = id_begin(id_twin(ne));
				_borders.erase(f);
				_faces.insert(f);
			}
			
			return ne;
		}
	}

	int ne = _hedges.size();	// new half edge
	int net = ne + 1;			// new half edge twin
	_hedges.push_back(HalfEdge(head_id, net, net));
	_hedges.push_back(HalfEdge(tail_id, ne, net));
	_borders.insert(net);

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
		_borders.erase(net);
		_hedges.pop_back();
		_hedges.pop_back();
		return HE_INVALID_INDEX;
	}
}


template <typename V>
inline int HEMesh<V>::add_edge_at(const int& prev_id, const int& next_id) {
	check_hedge_id(prev_id);
	check_hedge_id(next_id);

	int ne = _hedges.size();	// new half edge
	int net = ne + 1;			// new half edge twin

	int vb = id_head(next_id);
	int va = id_head(prev_id);

	_hedges.push_back(HalfEdge(vb, net, net));
	_hedges.push_back(HalfEdge(va, ne, net));
	//_borders.insert(net);
	//_faces.insert(net);

	if (!(is_border_hedge(prev_id) || is_border_hedge(next_id)))
		_faces.insert(net);
	else
		_borders.insert(net);


	swap_next(prev_id, net);
	swap_next(next_id, ne);

	_vert_to_hedge[va] = ne;
	_vert_to_hedge[vb] = net;

	return ne;
}


// @todo: Try not to use add_edge and constrain it to not allow (inner-edges after the operation) or (edges/faces that cross the new face)
template <typename V>
inline int HEMesh<V>::add_face(const std::vector<int>& indices) {
	if (indices.size() < 2ULL)
		throw std::invalid_argument("The number of new face indices cannot be smaller than 2");

	bool isEdgeAdded = false;

	for (int i = 0, n = (int)indices.size() - 1ULL; i < n; ++i) {
		check_vert_id(indices[i]);

		if (id_hedge(indices[i], indices[i + 1]) == HE_ISOLATED_INDEX) {
			add_edge(indices[i], indices[i + 1], true);
			isEdgeAdded = true;
		}
	}

	int e = id_hedge(indices.back(), indices.front());
	if (e == HE_ISOLATED_INDEX) {
		e = add_edge(indices.back(), indices.front(), true);
		isEdgeAdded = true;
	}

	if (!isEdgeAdded && is_border_hedge(e)) {
		const int& f = id_begin(e);
		_borders.erase(f);
		_faces.insert(f);
	}

	return id_hedge(indices[0], indices[1]);
}


template<typename V>
inline void HEMesh<V>::remove_vert(const int& vert_id) {
	check_vert_id(vert_id);

	// Remove all of its adjacent edges
	while (!is_isolated_vert(vert_id))
		remove_edge(id_hedge(vert_id));

	// Remove The Vertex
	_vert_to_hedge[vert_id] = HE_INVALID_INDEX;
	_garbage_vertices.insert(vert_id);
}


template<typename V>
inline void HEMesh<V>::remove_edge(const int& hedge_id) {

	if (is_removed_hedge(hedge_id))
		return;

	check_hedge_id(hedge_id);

	int t = id_twin(hedge_id);
	const int& nc = id_next(hedge_id);
	const int& nt = id_next(t);
	int pc = id_prev(hedge_id);
	int pt = id_prev(t);

	// Vert to Half-Edge arrangement

	if (nt != hedge_id) {
		_hedges[pc].next_id = nt;
		_vert_to_hedge[id_head(t)] = nt;
		//_vert_to_hedge[id_tail(nt)] == nt;
	}
	else {
		_vert_to_hedge[id_head(t)] = HE_ISOLATED_INDEX;
	}

	if (nc != t) {
		_hedges[pt].next_id = nc;
		_vert_to_hedge[id_head(hedge_id)] = nc;
	}
	else {
		_vert_to_hedge[id_head(hedge_id)] = HE_ISOLATED_INDEX;
	}

	// Faces and borders arrangement

	int f = id_begin(hedge_id);
	int tf = id_begin(t);

	bool bFace = _borders.find(f) == _borders.end() && _borders.find(tf) == _borders.end();
	_faces.erase(f);
	_borders.erase(f);
	_faces.erase(tf);
	_borders.erase(tf);

	if (nt != hedge_id) {
		if (bFace) _faces.insert(nt);
		else _borders.insert(nt);

		set_begin_id(nt);
	}

	if (nc != t) {
		//if (nt != hedge_id && f == tf) {
		if (nt == hedge_id || (nt != hedge_id && f == tf)) {
			if (bFace) _faces.insert(nc);
			else _borders.insert(nc);

			set_begin_id(nc);
		}
	}

	// Invalidate the edge
	_hedges[hedge_id] = HalfEdge();
	_hedges[t] = HalfEdge();
	_garbage_hedges.insert(id_ltwin(hedge_id));

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
	return edge_iter(edge_id, const_cast<HEMesh<V>*>(this));
}

/*template <typename V>
inline typename HEMesh<V>::face_iter HEMesh<V>::face(const int& hedge_id) const {
	check_hedge_id(hedge_id);
	return face_iter(hedge_id, const_cast<HEMesh<V>*>(this));
}*/


template <typename V>
inline typename HEMesh<V>::vert_iter HEMesh<V>::begin_verts() const {
	int id = 0;
	
	while (is_removed_vert(id) && id < _vertices.size())
		++id;
	
	return vert_iter(id, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::hedge_iter HEMesh<V>::begin_hedges() const {
	int id = 0;
	
	while (is_removed_hedge(id) && id < _hedges.size())
		++id;
	
	return hedge_iter(id, const_cast<HEMesh<V>*>(this));
}

/*template <typename V>
inline typename HEMesh<V>::edge_iter HEMesh<V>::begin_edges() const {
	int id = 0;
	
	while (is_removed_edge(id) && (2 * id < _hedges.size()))
		++id;
	
	return edge_iter(id, const_cast<HEMesh<V>*>(this));
}*/

template <typename V>
inline typename HEMesh<V>::edge_iter HEMesh<V>::begin_edges() const {
	int hedge_id = 0;

	while (is_removed_hedge(hedge_id) && (hedge_id < _hedges.size()))
		hedge_id += 2;

	return edge_iter(hedge_id >> 1, const_cast<HEMesh<V>*>(this));
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
	return edge_iter(_hedges.size() >> 1, const_cast<HEMesh<V>*>(this));
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
	
	while (is_removed_vert(id) && id >= 0)
		--id;
	
	return vert_iter(id, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::hedge_iter HEMesh<V>::rbegin_hedges() const {
	int id = (int)_hedges.size() - 1;
	
	while (is_removed_hedge(id) && id >= 0)
		--id;
	
	return hedge_iter(id, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::edge_iter HEMesh<V>::rbegin_edges() const {
	int id = int(_hedges.size()) / 2 - 1;
	
	while (is_removed_edge(id) && id >= 0)
		--id;
	
	return edge_iter(id, const_cast<HEMesh<V>*>(this));
}


template <typename V>
inline typename HEMesh<V>::vert_iter HEMesh<V>::rend_verts() const {
	return vert_iter(HE_INVALID_INDEX, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::hedge_iter HEMesh<V>::rend_hedges() const {
	return hedge_iter(HE_INVALID_INDEX, const_cast<HEMesh<V>*>(this));
}

template <typename V>
inline typename HEMesh<V>::edge_iter HEMesh<V>::rend_edges() const {
	return edge_iter(HE_INVALID_INDEX, const_cast<HEMesh<V>*>(this));
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