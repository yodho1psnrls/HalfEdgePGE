#ifndef HH_BASE_MESH_CLASS_H
#define HH_BASE_MESH_CLASS_H

#include "utilities.hpp" // for generate_plane
#include <vector>

using uint = unsigned int;
//using uint8 = unsigned char;


// @todo: Polygonal mesh (assumes that all polygons have the same number N of vertices/indices)

// Triangular Mesh // OpenGL Style
template <typename V>
struct Mesh {

	std::vector<V> vertices{};
	std::vector<uint> indices{};

	Mesh() {}

	Mesh(const std::vector<V>& vertices, const std::vector<uint>& indices)
		: vertices(vertices), indices(indices) {}

	
	bool is_valid() const { return indices.size() % 3ULL == 0ULL; }


};



template <typename V>
Mesh<V> generate_plane(const int div, const float screen_percentage = 1.0f) {
	Mesh<V> mesh;

	for (int i = 0; i < div + 1; i++)
		for (int j = 0; j < div + 1; j++) {
			V v;

			v.uv = olc::vf2d(float(i) / div, float(j) / div); // from 0 to 1
			v.pos = 2.0f * v.uv - olc::vf2d(1.0f, 1.0f);	  // from -1 to 1
			v.pos *= screen_percentage;
			//pos = -pos.y
			v.pos = (v.pos + olc::vf2d(1.0f, 1.0f)) * 0.5f;
			v.pos *= olc::vf2d(scrW, scrH);

			mesh.vertices.push_back(v);
		}

	for (int i = 0; i < div; i++) {
		for (int j = 0; j < div; j++) {

			int index = i + j * (div + 1);

			mesh.indices.push_back(index);					// 0--*
			mesh.indices.push_back(index + 1);				// |  |
			mesh.indices.push_back(index + (div + 1));		// 0--0

			mesh.indices.push_back(index + 1);				// *--0
			mesh.indices.push_back(index + (div + 1) + 1);	// |  |
			mesh.indices.push_back(index + (div + 1));		// 0--0
		}
	}

	return mesh;
}


#endif