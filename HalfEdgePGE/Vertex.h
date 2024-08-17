#pragma once
#include "olcPixelGameEngine.h"

struct Vertex {

	using vf2d = olc::vf2d;

	vf2d pos;
	vf2d uv;

	//vf2d tan;
	//vf2d bitan;

	Vertex() : pos(0.0f, 0.0f), uv(0.0f, 0.0f) {}
	Vertex(const vf2d& Pos) : pos(Pos), uv(0.0f, 0.0f) {}

	//operator vf2d& () { return pos; }
	//operator const vf2d& () const { return pos; }
	//operator vi2d () const { return pos; }


	// Component/Attribute wise operations on vertices

	Vertex& operator+=(const Vertex& other) {
		pos += other.pos;
		uv += other.uv;

		return *this;
	}

	Vertex& operator-=(const Vertex& other) {
		pos -= other.pos;
		uv -= other.uv;

		return *this;
	}

	Vertex& operator*=(const float& scalar) {
		pos *= scalar;
		uv *= scalar;

		return *this;
	}

	Vertex& operator/=(const float& scalar) {
		pos /= scalar;
		uv /= scalar;

		return *this;
	}


	Vertex operator+(const Vertex& rhs) const {
		Vertex v(*this);
		return v += rhs;
	}

	Vertex operator-(const Vertex& rhs) const {
		Vertex v(*this);
		return v -= rhs;
	}

	Vertex operator*(const float& scalar) const {
		Vertex v(*this);
		return v *= scalar;
	}

	Vertex operator/(const float& scalar) const {
		Vertex v(*this);
		return v /= scalar;
	}

	/*friend Vertex operator*(const float& scalar, const Vertex& vert) {	// static	// friend // constexpr
		Vertex v(vert);
		return v *= scalar;
	}*/

};

//Vertex operator*(const float& scalar, Vertex vert) { return vert *= scalar; }
Vertex operator*(const float& scalar, const Vertex& vert) {
	Vertex v(vert);
	return v *= scalar;
}