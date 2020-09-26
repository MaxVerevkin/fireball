#include "structs.h"
#include <cmath>


/////////////////
/// 2D vector ///
/////////////////

vec2d_t vec2d_t::operator*(double a) const {
    return {x*a, y*a};
}
vec2d_t vec2d_t::operator+(const vec2d_t &vec) const {
    return {x+vec.x, y+vec.y};
}
vec2d_t vec2d_t::operator-(const vec2d_t &vec) const {
    return {x-vec.x, y-vec.y};
}

double vec2d_t::length() const {
    return sqrt(x*x+y*y);
}


/////////////////
/// 3D vector ///
/////////////////

// Dot product
double vec3d_t::operator*(const vec3d_t &vec) const {
    return x*vec.x + y*vec.y + z*vec.z;
}

vec3d_t vec3d_t::operator*(double a) const {
    return {x*a, y*a, z*a};
}
vec3d_t vec3d_t::operator/(double a) const {
    return {x/a, y/a, z/a};
}

vec3d_t vec3d_t::operator-(const vec3d_t &vec) const {
    return {x-vec.x, y-vec.y, z-vec.z};
}
vec3d_t vec3d_t::operator+(const vec3d_t &vec) const {
    return {x+vec.x, y+vec.y, z+vec.z};
}

double vec3d_t::length() const {
    return sqrt(x*x+y*y+z*z);
}

vec3d_t vec3d_t::normalized() const {
    return *this/length();
}

vec2d_t vec3d_t::to2d() const {
    return {x, y};
}


/////////////////
/// Dataset   ///
/////////////////

data_set_t::data_set_t(int n) {
    z0 = new double[n];
    h0 = new double[n];
    zb = new double[n];
    hb = new double[n];
    a =  new double[n];
    t =  new double[n];
}
