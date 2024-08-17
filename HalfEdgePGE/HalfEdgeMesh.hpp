#include "Mesh.hpp"

#include <vector>
#include <unordered_set>
#include <unordered_map>

#include <set>
//#include <stack>

#include <stdexcept>
#include <fstream>

#define HE_INVALID_INDEX -1
#define HE_ISOLATED_INDEX -2



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
	std::vector<int> _vert_to_hedge{};
	std::vector<HalfEdge> _hedges{};

	std::unordered_set<int> _faces{};   // the indexes of the half edges that are the beginnig of a face
	std::unordered_set<int> _borders{}; // the indexes of the half edges that are the beginnig of a border

	// Garbage collection/recycling
	std::set<int> _garbage_vertices;	// stores vertex ids that are removed
	std::set<int> _garbage_edges;		// stores edge ids that are removed

	

public:

	HEMesh() {}
	HEMesh(const Mesh<V>& mesh);

	// Supports .obj and .ply(ascii) file formats, if the file format is other an error is thrown
	HEMesh(const char* filePath);
	
	void load_from_obj(const char* file_path);
	void save_to_obj(const char* file_path);
	//void load_from_ply(const char* file_path);
	//void save_to_ply(const char* file_path);

	
	const int& id_head(const int& hedge_id) const;
	const int& id_tail(const int& hedge_id) const;
	const int& id_next(const int& hedge_id) const;
	const int& id_prev(const int& hedge_id) const;
	const int& id_twin(const int& hedge_id) const;
	const int& id_begin(const int& hedge_id) const;
	int id_begin(const int& vert_id1, const int& vert_id2) const;	// gives the begin half edge of the common face between the two vertices, returns HM_INVALID_INDEX if not found
	const int& id_hedge(const int& vert_id) const;
	int id_hedge(const int& tail_id, const int& head_id) const;	// gives the id of the half edge with the given tail and head vertices, returns HM_ISOLATED_INDEX if not found
	const int& id_edge(const int& hedge_id) const;



	V& vert(const int& vert_id);
	const V& vert(const int& vert_id) const;


};
