#pragma once
#include "olcPixelGameEngine.h"
#include <cmath>

//const float scrW = 640.0f;
//const float scrH = 420.0f;
const float scrW = 380.0f;
const float scrH = 340.0f;

typedef unsigned int uint;

// Isnt this std::fmodf(x, 1.0f) ??
inline float fract(const float x) {
	return x - (x >= 0.0f ? floor(x) : ceil(x));
	//return x - floor(x);
}

// from x(-1, 1) and y(-1, 1) to x(0, scrW) and y(0, scrH), and reverts Y axis 
inline olc::vf2d to_screen(olc::vf2d v) {
	v.y = -v.y;
	v = (v + olc::vf2d(1.0f, 1.0f)) * 0.5f; // to (0, 1)
	return v * olc::vf2d(scrW, scrH);
}

// template specialization for std::hash to work with std::pair
namespace std {
	template<typename T, typename Y> struct std::hash<std::pair<T, Y>> {
		size_t operator()(const std::pair<T, Y>& p) const {
			size_t h1 = std::hash<T>()(p.first);
			size_t h2 = std::hash<Y>()(p.second);
			return h1 ^ (h2 << 1);
		}
	};
}

template <typename IterType>	// Iterator Type
class Iterable {
protected:

	IterType begin_iter;
	IterType end_iter;

public:

	Iterable(const IterType& begin_iterator, const IterType& end_iterator)
		: begin_iter(begin_iterator), end_iter(end_iterator) {}

	//IterType begin() const { return begin_iter; }
	//IterType end() const { return end_iter; }

	const IterType& begin() const { return begin_iter; }
	const IterType& end() const { return end_iter; }

};