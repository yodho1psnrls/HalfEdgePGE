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
	int begin_id;	// the index of the half edge at the beggining of the face of the current half edge

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
	std::set<int> _garbage_edges{};		// stores edge ids that are removed // i would prefer it to store half edge ids, but if one half edge is marked as removed, storing edge_id will also mark its twin half edge as removed
	//bool _are_verts_shrinked;
	//bool _are_hedges_shrinked;

public:

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

	int vertices_size() const;
	int hedges_size() const;
	int edges_size() const;
	int faces_size() const;
	int borders_size() const;

	void reserve_vertices(const int& capacity);
	void reserve_hedges(const int& capacity);
	void reserve_edges(const int& capacity);
	//void reserve_faces(const int& capacity);
	//void reserve_borders(const int& capacity);


	const int& id_head(const int& hedge_id) const;
	const int& id_tail(const int& hedge_id) const;
	const int& id_next(const int& hedge_id) const;
	int id_prev(const int& hedge_id) const;
	const int& id_twin(const int& hedge_id) const;
	const int& id_begin(const int& hedge_id) const;
	int id_begin(const int& vert_id1, const int& vert_id2) const; // gives the begin half edge of the common face between the two vertices, returns HE_INVALID_INDEX if not found
	const int& id_hedge(const int& vert_id) const;
	int id_hedge(const int& tail_id, const int& head_id) const;	// gives the id of the half edge with the given tail and head vertices, returns HE_ISOLATED_INDEX if not found
	const int& id_edge(const int& hedge_id) const;

	bool is_removed_vert(const int& vert_id) const;
	bool is_removed_hedge(const int& hedge_id) const;
	bool is_removed_edge(const int& hedge_id) const;

	bool is_border_hedge(const int& hedge_id) const;
	bool is_border_edge(const int& hedge_id) const;
	bool is_border_vert(const int& vert_id) const;
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

	bool do_share_loop(const int& hedge_id1, const int& hedge_id2) const;
	//void shift_begin(const int& begin_id, const int& how_much = 1);

	

	int add_vert(const V& v);
	int add_edge(const int& tail_id, const V& head);
	int add_edge(const int& tail_id, const int& head_id);
	int add_edge_at(const int& prev_hedge_id, const int& hext_hedge_id);
	int add_face(const std::vector<int>& face_indices);

	int remove_vert(const int& vert_id);
	int remove_edge(const int& hedge_id);
	int remove_face(const int& hedge_id);

	void move_edge(const int& hedge_id, const int& new_next_id);
	void move_edge(const int& hedge_id, const int& new_prev_id, const int& new_next_id);




};


 // ================================== IMPLEMENTATION ======================================== //


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
			//e->id_twin = new HalfEdge();
			//e->id_twin->id_twin = e;
			// CHECK IF YOU NEED TO 

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
				if (id_next(id_twin(e)) == HE_INVALID_INDEX)					// !!!!
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
	_garbage_edges.clear();
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
inline int HEMesh<V>::vertices_size() const {
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
inline int HEMesh<V>::faces_size() const{
	return _faces.size();
}

template<typename V>
inline int HEMesh<V>::borders_size() const {
	return _borders.size();
}


template<typename V>
inline void HEMesh<V>::reserve_vertices(const int& capacity) {
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
inline const int& HEMesh<V>::id_twin(const int& hedge_id) const {
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

template<typename V>
inline const int& HEMesh<V>::id_edge(const int& hedge_id) const {
	return hedge_id / 2;
}

template<typename V>
inline bool HEMesh<V>::is_removed_vert(const int& vert_id) const {
	return _vert_to_hedge[vert_id] == HE_INVALID_INDEX;
}

template<typename V>
inline bool HEMesh<V>::is_removed_hedge(const int& hedge_id) const {
	return _hedges[hedge_id].next_id == HE_INVALID_INDEX;
}

template<typename V>
inline bool HEMesh<V>::is_removed_edge(const int& hedge_id) const {
	//return is_removed_hedge(hedge_id) && is_removed_hedge(id_twin(hedge_id));
	return is_removed_hedge(hedge_id) || is_removed_hedge(id_twin(hedge_id));
}

template<typename V>
inline bool HEMesh<V>::is_border_hedge(const int& hedge_id) const {
	return _borders.find(hedge_id) != _borders.end();
}

template<typename V>
inline bool HEMesh<V>::is_border_edge(const int& hedge_id) const {
	return is_border_hedge(hedge_id) || is_border_hedge(id_twin(hedge_id));
}

template<typename V>
inline bool HEMesh<V>::is_border_vert(const int& vert_id) const {
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
inline bool HEMesh<V>::is_n_poly(const int& hedge_id, const int& n) const {
	int e = hedge_id;
	
	for (int i = 0; i < n; ++i)
		e = id_next(e);

	return e == hedge_id;
}

template<typename V>
inline int HEMesh<V>::get_n_poly(const int& hedge_id) const {
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
	const int& he = _vert_to_hedge[vert_id];
	
	if (he == HE_ISOLATED_INDEX)
		return false;

	return id_next(id_twin(he)) == he;
}

template<typename V>
inline bool HEMesh<V>::is_isolated_vert(const int& vert_id) const {
	return _vert_to_hedge[vert_id] == HE_ISOLATED_INDEX;
}

template<typename V>
inline bool HEMesh<V>::is_isolated_edge(const int& hedge_id) const {
	int t = id_twin(hedge_id);
	return id_next(hedge_id) == t && id_next(t) == hedge_id;
}

template<typename V>
inline bool HEMesh<V>::is_isolated_face(const int& hedge_id) const {
	
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
inline bool HEMesh<V>::do_share_loop(const int& hedge_id1, const int& hedge_id2) const {
	int e = hedge_id1;

	do {
		if (e == hedge_id2)
			return true;

		e = id_next(e);
	} while (e != hedge_id1);

	return false;
}
