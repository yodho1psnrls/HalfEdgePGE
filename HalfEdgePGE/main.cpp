#include "HalfEdgeMesh.hpp"
#include "Mesh.hpp"
#include "utilities.hpp"
#include "Vertex.h"
#include <iostream>

const float PINDIST = 8.0f;	// the range that the mouse will select a vertex
const int INIT_HEDGE_ID = 0;

template <typename V>
class Example : public olc::PixelGameEngine {

	Mesh<V> mesh;
	HEMesh<V> he_mesh;

	int selected_vert_id;
	int hedge_id;

	//@todo: Save the state of the mesh after some operations, so you can undo and outdo
	//std::vector<Mesh<V>> prevHistory;
	//std::vector<Mesh<V>> nextHistory;

public:
	Example() { sAppName = "Example"; }

	void drawMesh(const Mesh<V>& mesh, const olc::Pixel& color = olc::WHITE) {
		for (int i = 0U; i < mesh.indices.size(); i += 3U) {
			const V& a = mesh.vertices[mesh.indices[i]];
			const V& b = mesh.vertices[mesh.indices[i + 1U]];
			const V& c = mesh.vertices[mesh.indices[i + 2U]];

			// To better visualize the different triangles, we will color them based on
			// their uv.x coordinate
			float diffuse = (a.uv.x + b.uv.x + c.uv.x) / 3.0f;
			//float diffuse = (b.pos - a.pos).norm().perp().dot((c.pos - a.pos).norm());
			//float diffuse = 0.5f;


			FillTriangle(
				mesh.vertices[mesh.indices[i]].pos,
				mesh.vertices[mesh.indices[i + 1U]].pos,
				mesh.vertices[mesh.indices[i + 2U]].pos,
				color * diffuse);
		}
	}

	void drawMeshWireFrame(const Mesh<V>& mesh, const olc::Pixel& color = olc::WHITE) {
		for (int i = 0U; i < mesh.indices.size(); i += 3U) {
			DrawTriangle(
				mesh.vertices[mesh.indices[i]].pos,
				mesh.vertices[mesh.indices[i + 1U]].pos,
				mesh.vertices[mesh.indices[i + 2U]].pos,
				color);
		}
	}

	void moveVertexWithCursor() {
		olc::vf2d mouse_pos = GetMousePos();

		if (GetMouse(0).bPressed)
			for (int i = 0; i < mesh.vertices.size(); ++i) {
				const V& v = mesh.vertices[i];
				if ((v.pos - mouse_pos).mag() < PINDIST) {
					selected_vert_id = i;
					break;
				}
			}
				

		if (GetMouse(0).bReleased)
			selected_vert_id = -1;

		if (selected_vert_id != -1)
			mesh.vertices[selected_vert_id].pos = mouse_pos;

	}

	void print_info(const Mesh<V>& m) const {
		std::cout << "vertices: " << m.vertices.size() << "\n";
		std::cout << "triangles: " << m.indices.size() / 3ULL << "\n";
	}

	void print_info(const HEMesh<V>& hm) const {
		std::cout << "vertices: " << hm.vertices_size() << "\n";
		std::cout << "edges: " << hm.edges_size() << "\n";
		std::cout << "faces: " << hm.faces_size() << "\n";
		std::cout << "borders: " << hm.borders_size() << "\n";
	}


public:
	bool OnUserCreate() override {
		
		mesh = generate_plane<V>(5, 0.9f);
		print_info(mesh);
		he_mesh = mesh;
		mesh = he_mesh;
		print_info(mesh);

		selected_vert_id = -1;
		hedge_id = INIT_HEDGE_ID;

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override {
		Clear(olc::BLACK);

		drawMesh(mesh);

		moveVertexWithCursor();


		return true;
	}
};

int main() {
	Example<Vertex> demo;

	if (demo.Construct(scrW, scrH, 2, 2, false, true))
		demo.Start();

	return 0;
}