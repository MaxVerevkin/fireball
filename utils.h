#ifndef UTILS_H
#define UTILS_H


#include <math.h>

#include "structs.h"
#include "simd.h"

#define PI 3.141592653589793
#define EARTH_R 6371000.

// Translate geographical location to XYZ.
vec3d_t geo_to_xyz(const vec3d_t &pos);

// Calculate delta of two angles.
inline double angle_delta(double a1, double a2) {
    //return PI - fmod((PI - a2 + a1), 2*PI);
    double delta = fmod(a2 - a1, 2*PI);
    double retval = abs(delta);

    //if (retval > PI)
        //retval = 2*PI - retval;
    int big_retval = retval > PI;
    retval *= big_retval * -2 + 1;
    retval += 2*PI*big_retval;

    int sign = ((delta >= 0 && delta <= PI) || (delta >= -2*PI && delta <= -PI)) * 2 - 1;
    return retval * sign;
}
__m128d angle_delta_sq_pd(double *addr1, double *addr2);

double calc_height(vec3d_t p, vec2d_t p0, double h);

// Calculate tje length of arc on a sphere
double arc_len(double h1, double z1, double h2, double z2);
inline double arc_len(const vec2d_t &a, const vec2d_t &b) {
    return arc_len(a.x, a.y, b.x, b.y);
}

// Calculate the azimuth given XYZ coordinates, North and East vectors.
double azimuth(const vec3d_t &point, const vec3d_t &normal, const vec2d_t &observer, double r);

// Calculate the disent angle for the begining of the path.
double desent_angle(double h, double z, double h0, double z0);
inline double desent_angle(const vec2d_t &start, const vec2d_t &end) {
    return desent_angle(start.x, start.y, end.x, end.y);
}

// Calculate normal
vec3d_t normal_vec(double lat, double lon);
inline vec3d_t normal_vec(const vec2d_t &vec) {
    return normal_vec(vec.x, vec.y);
}

// Genetare vector pointing to Notrh
vec3d_t north_vec(const vec2d_t &ob);
// Genetare vector pointing to East
vec3d_t east_vec(double lon);

// Translate global vector to local
vec3d_t global_to_local(vec3d_t vec, double lat, double lon);


#endif
