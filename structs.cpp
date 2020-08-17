#include "structs.h"
#include <cmath>


// Dot product
double vec3d_t::operator*(const vec3d_t &vec) const {
    return x*vec.x + y*vec.y + z*vec.z;
}

vec3d_t vec3d_t::operator*(double a) const {
    return vec3d_t {x*a, y*a, z*a};
}
vec3d_t vec3d_t::operator/(double a) const {
    return vec3d_t {x/a, y/a, z/a};
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

vec3d_t vec3d_t::normalized() const {
    return *this/length();
}


data_set_t::data_set_t(int n) {
    z0 = new double[n];
    h0 = new double[n];
    zb = new double[n];
    hb = new double[n];
    a =  new double[n];
    t =  new double[n];
}
