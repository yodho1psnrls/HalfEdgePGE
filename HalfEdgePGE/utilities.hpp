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