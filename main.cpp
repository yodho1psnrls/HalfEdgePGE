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

	std::vector<std::pair<HEMesh<V>, int>> prev_history;
	std::vector<std::pair<HEMesh<V>, int>> next_history;

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
		//	float diffuse = 0.5f;


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

		/*for (const face_iter& f : hm.faces()) {
			e = *f;
			ne = e.next();

			do {

				//if (hm.idTwin(e))
				DrawLine(e.head()->pos, ne.head()->pos, color);

				e = ne;
				ne = ne.next();
			} while (e != f);
		}*/

		/*for (const auto& edge : hm.edges()) {
			hedge_iter he = edge;
			DrawLine(he.head()->pos, he.tail()->pos, color);
		}*/

for (hedge_iter half_edge : hm.edges())
DrawLine(half_edge.head()->pos, half_edge.tail()->pos, color);

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
		if (GetKey(olc::V).bReleased) he_mesh.check_validity();

		if (GetKey(olc::Z).bReleased && !prev_history.empty()) {
			next_history.push_back({ he_mesh, hedge.index() });

			auto pairr = prev_history.back();
			he_mesh = pairr.first;
			hedge = he_mesh.hedge(pairr.second);

			prev_history.pop_back();
		}

		if (GetKey(olc::X).bReleased && !next_history.empty()) {
			prev_history.push_back({ he_mesh, hedge.index() });

			auto pairr = next_history.back();
			he_mesh = pairr.first;
			hedge = he_mesh.hedge(pairr.second);

			next_history.pop_back();
		}


		if (GetKey(olc::S).bReleased) {
			he_mesh.remove_isolated_verts();
			he_mesh.shrink();
		}
		
		if (GetKey(olc::E).bReleased) {
			int given_hedge_id;
			std::cin >> given_hedge_id;
			hedge = he_mesh.hedge(given_hedge_id);
		}

		// Draws the begin half-edge of the loop and indicates is it a face(GREEN) or border(RED)
		if (GetKey(olc::SPACE).bHeld) {
			drawHalfEdge(hedge.face(), olc::WHITE);

			if (hedge.is_border())
				FillRect({ 0, 0 }, { 15, 15 }, olc::RED);
			else
				FillRect({ 0, 0 }, { 15, 15 }, olc::GREEN);
		}
		

		if (GetKey(olc::R).bReleased) {
			next_history.clear();
			prev_history.push_back({ he_mesh, hedge.index() });

			auto temp = hedge;

			if (hedge.next().twin() == hedge)
				temp = hedge.next().next();
			else
				temp = hedge.next();

			he_mesh.remove_edge(hedge.index());

			hedge = temp;

			he_mesh.check_validity();
		}

		if (GetKey(olc::F).bReleased) {
			next_history.clear();
			prev_history.push_back({ he_mesh, hedge.index() });
			
			//he_mesh.triangulate();
			//auto temp = hedge.prev();
			//he_mesh.remove_vert(hedge.head());
			
			//he_mesh.shift_begin(hedge.face(), 13);
			//he_mesh.flip_edge(hedge.index());
			//he_mesh.refine_edge(hedge);

			V v = *hedge.head();
			v.pos += olc::vf2d(-10.0f, -17.0f);
			//he_mesh.split_vert_to_edge(hedge, hedge.next().twin(), v);
			he_mesh.split_vert_to_faces(hedge, hedge.next().twin(), v);
			

		//	hedge = temp;
			he_mesh.check_validity();
		}

	}

	void print_info(const Mesh<V>& m) const {
		std::cout << "================\n";
		std::cout << "vertices: " << m.vertices.size() << "\n";
		std::cout << "triangles: " << m.indices.size() / 3ULL << "\n";
		std::cout << "\n";
	}

	void print_info(const HEMesh<V>& hm) const {
		std::cout << "================\n";
		std::cout << "vertices: " << hm.verts_size() << "\n";
	//	std::cout << "max_verts: " << hm.verts_max_size() << "\n";
	//	std::cout << "\n";
		std::cout << "edges: " << hm.edges_size() << "\n";
	//	std::cout << "max_edges: " << hm.edges_max_size() << "\n";
	//	std::cout << "\n";
		std::cout << "faces: " << hm.faces_size() << "\n";
		std::cout << "borders: " << hm.borders_size() << "\n";
		std::cout << "\n";
	}


	HEMesh<V> test_mesh1() {
		using namespace olc;

		HEMesh<V> hm;

		hm.add_vert(vf2d(-0.5f, -0.5f));
		hm.add_vert(vf2d(-0.5f, 0.5f));
		hm.add_vert(vf2d(0.5f, 0.5f));
		hm.add_vert(vf2d(0.5f, -0.5f));

		//hm.add_edge(0, 1, true);
		//hm.add_edge(1, 2, true);
		//hm.add_edge(2, 3, true);
		//hm.add_edge(3, 0, true);

		int nf = hm.add_face({ 0, 1, 2, 3 });
		//std::cout << nf << "\n";

		//hm.add_edge_at(0, 4);
		//hm.add_edge_at(2, 6);


		for (int i = 0; i < 1000; ++i) {
			hm.remove_edge(hm.id_hedge(0, 1));
			hm.add_edge(0, 1, true);
		}


		for (V& v : hm.verts()) {
			v.pos = to_screen(v.pos);
			v.uv += vf2d(0.5f, 0.5f);
		}


		return hm;
	}

	HEMesh<V> test_mesh2() {
		using namespace olc;

		HEMesh<V> hm;

		hm.add_vert(vf2d(-0.5f, -0.5f));

		hm.add_edge(0, vf2d(-0.5f, 0.5f));
		hm.add_edge(1, vf2d(0.5f, 0.5f));
		hm.add_edge(2, vf2d(0.5f, -0.5f));

		hm.add_edge(3, 0);
		//hm.split_face_at(1, 4);
		//hm.add_edge_at(1, 4);

		std::vector<int> test_face_hedges;
		hedge = hm.hedge(0);
		auto e = hedge;

		do {
			auto ne = e.next();
			
			//test_face_hedges.push_back(hm.add_edge_at(e.index(), 0.5f * e.head()->pos));
			test_face_hedges.push_back(hm.add_edge_at(e.index(), 0.5f * *e.head()));

			e = ne;
		} while (e != hedge);

		hm.add_face_at(test_face_hedges);


		//int vi = hm.add_vert(vf2d(0.0f, 0.0f));
		//hm.add_edge(0, vi);
		//hm.add_edge(vi, 1);

		//hm.add_edge(vi, 1);
		//hm.add_edge(vi, 2);
		//hm.add_edge(vi, 3);


		/*int he;
		he = hm.add_edge(0, vf2d(-0.5f, 0.5f));
		he = hm.add_edge(hm.hedge(he).head().index(), vf2d(0.5f, 0.5f));
		he = hm.add_edge(hm.hedge(he).head().index(), vf2d(0.5f, -0.5f));
		he = hm.add_edge(hm.hedge(he).head().index(), 0);*/


		for (V& v : hm.verts()) {
			v.pos = to_screen(v.pos);
			v.uv += vf2d(0.5f, 0.5f);
		}

		return hm;
	}

	void test_mesh3() {
		mesh = generate_plane<V>(5, 0.9f);

		for (const V& v : mesh.vertices)
			he_mesh.add_vert(v);

		for (int i = 0; i < mesh.indices.size(); i += 3) {
			for (int j = 0; j < 3; ++j) {
				const int& vi0 = mesh.indices[i + j];
				const int& vi1 = mesh.indices[(i + j + 1) % 3];

				he_mesh.add_edge(vi0, vi1);
			}
		}

	}


public:
	bool OnUserCreate() override {
		// Try removing some faces of the grid and see if most of the times
		//  the begin index of the loop is preserved, if the begin index is not removed

		mesh = generate_plane<V>(3, 0.9f);
		he_mesh = mesh;

		//he_mesh = test_mesh2();
		//mesh = he_mesh;

		// test_mesh3();

		he_mesh.check_validity();
		print_info(he_mesh);
		selected_vert_id = -1;
		hedge = he_mesh.hedge(INIT_HEDGE_ID);

		return true;
	}

	bool OnUserUpdate(float fElapsedTime) override {
		Clear(olc::BLACK);

		mesh = he_mesh;

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