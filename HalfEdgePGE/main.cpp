#include "Mesh.hpp"
#include "utilities.hpp"
#include "Vertex.h"

const float PINDIST = 8.0f;	// the range that the mouse will select a vertex


template <typename V>
class Example : public olc::PixelGameEngine {

	Mesh<V> mesh;

	int selected_vert_id;

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


public:
	bool OnUserCreate() override {
		
		mesh = generate_plane<V>(5, 0.9f);
		
		selected_vert_id = -1;

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