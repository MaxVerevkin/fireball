#include "structs.h"

#include <math.h>

// Dot product
double vec2d_t::operator*(const vec2d_t &vec) const {
    return x*vec.x + y*vec.y;
}

vec2d_t vec2d_t::operator*(double a) const {
    return vec2d_t {x*a, y*a};
}

vec2d_t vec2d_t::operator-(const vec2d_t &vec) const {
    return vec2d_t {x-vec.x, y-vec.y};
}
vec2d_t vec2d_t::operator+(const vec2d_t &vec) const {
    return vec2d_t {x+vec.x, y+vec.y};
}

double vec2d_t::length() const {
    return sqrt(x*x+y*y);
}


// Dot product
double vec3d_t::operator*(const vec3d_t &vec) const {
    return x*vec.x + y*vec.y + z*vec.z;
}

vec3d_t vec3d_t::operator*(double a) const {
    return vec3d_t {x*a, y*a, z*a};
}

vec3d_t vec3d_t::operator-(const vec3d_t &vec) const {
    return vec3d_t {x-vec.x, y-vec.y, z-vec.z};
}
vec3d_t vec3d_t::operator+(const vec3d_t &vec) const {
    return vec3d_t {x+vec.x, y+vec.y, z+vec.z};
}

double vec3d_t::length() const {
    return sqrt(x*x+y*y+z*z);
}
