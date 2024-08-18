#include "HalfEdgeMesh.hpp"
#include "Mesh.hpp"
#include "utilities.hpp"
#include "Vertex.h"
#include <iostream>

const float PINDIST = 8.0f;	// the range that the mouse will select a vertex
const int INIT_HEDGE_ID = 0;

template <typename V>
class Example : public olc::PixelGameEngine {

	using vert_iter = typename HEMesh<V>::vert_iter;
	using hedge_iter = typename HEMesh<V>::hedge_iter;
	using edge_iter = typename HEMesh<V>::edge_iter;
	using face_iter = typename HEMesh<V>::face_iter;

	Mesh<V> mesh;
	HEMesh<V> he_mesh;

	int selected_vert_id;
	//int hedge_id;
	hedge_iter hedge;

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

	void drawWireFrame(const Mesh<V>& mesh, const olc::Pixel& color = olc::WHITE) {
		for (int i = 0U; i < mesh.indices.size(); i += 3U) {
			DrawTriangle(
				mesh.vertices[mesh.indices[i]].pos,
				mesh.vertices[mesh.indices[i + 1U]].pos,
				mesh.vertices[mesh.indices[i + 2U]].pos,
				color);
		}
	}

	void drawWireFrame(const HEMesh<V>& hm, const olc::Pixel& color = olc::WHITE) {
		hedge_iter e;	// curr half edge
		hedge_iter ne;	// next half edge

		for (const face_iter& f : hm.faces()) {
			e = *f;
			ne = e.next();

			do {

				//if (hm.idTwin(e))
				DrawLine(e.head()->pos, ne.head()->pos, color);

				e = ne;
				ne = ne.next();
			} while (e != f);
		}
	}

	void drawHalfArrow(const olc::vf2d& from, const olc::vf2d& to, const olc::Pixel& color = olc::WHITE) {
		float sep = 6.0f;
		olc::vf2d n = (to - from).norm();

		for (float i = 0.0f; i < 3.0f; i += 1.0f) {
			olc::vf2d dif = n.perp() * (0.25f * sep + i);
			DrawLine(from + dif, to + dif, color);
			DrawLine(to + dif, to - n * sep * 1.25f + n.perp() * sep * 0.75f + dif, color);
		}

		//DrawLine(to, to - n * sep + n.perp() * sep, color);
	}

	//void drawHalfEdge(const HalfEdgeMesh<vf2d>::HalfEdge* edge, const olc::Pixel& color = olc::WHITE) {
	void drawHalfEdge(const hedge_iter& e, const olc::Pixel& color = olc::WHITE) {
		if (e.index() == -1)
			return;

		drawHalfArrow(e.tail()->pos, e.head()->pos, color);

		// Draw the head vertex index
		DrawStringDecal(e.head()->pos, std::to_string(e.head().index()), olc::GREEN);

		// Draw the half edge index
		DrawStringDecal(0.5f * (e.tail()->pos + e.head()->pos) + olc::vf2d(0.0f, 4.0f)
			, std::to_string(e.index()), olc::MAGENTA);

	}

	void drawBorders(const HEMesh<V>& hm, const olc::Pixel& color = olc::RED) {

		for (const face_iter& face : hm.borders()) {
			hedge_iter begin = face;
			hedge_iter e = begin;

			do {
				//drawHalfEdge(e, color);

				const olc::vf2d& from = e.tail()->pos;
				const olc::vf2d& to = e.head()->pos;
				DrawLine(from, to, color);

				olc::vf2d n = (to - from).norm().perp();
				DrawLine(from + n, to + n, color);


				e = e.next();
			} while (e != begin);
		}
	}


	void moveVertexWithCursor() {
		olc::vf2d mouse_pos = GetMousePos();

		if (GetMouse(0).bPressed)
			/*for (int i = 0; i < mesh.vertices.size(); ++i) {
				const V& v = mesh.vertices[i];
				if ((v.pos - mouse_pos).mag() < PINDIST) {
					selected_vert_id = i;
					break;
				}
			}*/
			for (auto vit = he_mesh.begin_verts(); vit != he_mesh.end_verts(); ++vit)
				if ((vit->pos - mouse_pos).mag() < PINDIST) {
					selected_vert_id = vit.index();
					break;
				}
				

		if (GetMouse(0).bReleased)
			selected_vert_id = -1;

		if (selected_vert_id != -1)
			//mesh.vertices[selected_vert_id].pos = mouse_pos;
			he_mesh.vert(selected_vert_id)->pos = mouse_pos;

	}

	void keyInput() {
		if (GetKey(olc::I).bReleased) print_info(he_mesh);
		if (GetKey(olc::G).bReleased) mesh = he_mesh;
		if (GetKey(olc::N).bReleased) hedge = hedge.next();
		if (GetKey(olc::P).bReleased) hedge = hedge.prev();
		if (GetKey(olc::T).bReleased) hedge = hedge.twin();
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
		he_mesh = mesh;
		print_info(he_mesh);

		selected_vert_id = -1;
		//hedge_id = INIT_HEDGE_ID;
		hedge = he_mesh.hedge(INIT_HEDGE_ID);

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override {
		Clear(olc::BLACK);

		drawMesh(mesh);
		drawWireFrame(he_mesh);
		drawBorders(he_mesh);

		drawHalfEdge(hedge, olc::DARK_CYAN);
		
		
		moveVertexWithCursor();
		keyInput();

		return true;
	}
};

int main() {
	Example<Vertex> demo;

	if (demo.Construct(scrW, scrH, 2, 2, false, true))
		demo.Start();

	return 0;
}