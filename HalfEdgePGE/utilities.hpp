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