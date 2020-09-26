#ifndef UTILS_H
#define UTILS_H


#include <math.h>
#include "structs.h"
#include "simd.h"


/*
 * Constants
 */
#define PI 3.141592653589793
#define EARTH_R 6371000.


/*
 * Translate geographical location to XYZ and back.
 */
vec3d_t geo_to_xyz(const vec3d_t &geo);
vec3d_t xyz_to_geo(const vec3d_t &xyz);



/*
 * Calculates height of flash above the sea, given
 * geolocation of observer, geolocation of flash
 * and altitude angle and back.
 */
double altitude_to_height(vec3d_t p, vec2d_t p0, double h);
double height_to_altitude(vec3d_t p, vec3d_t p0);


/*
 * Calculate delta of two angles
 */
inline double angle_delta(double a1, double a2) {
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
inline __m128d angle_delta_sq_pd(double *addr1, double *addr2) {
    double t1 = angle_delta(addr1[0], addr2[0]);
    double t2 = angle_delta(addr1[1], addr2[1]);
    t1 *= t1;
    t2 *= t2;
    return _mm_set_pd(t1, t2);
}


/*
 * Calculate the length of arc on a sphere
 */
inline double arc_len(const vec2d_t &a, const vec2d_t &b) {
    double cos_l = sin(a.x)*sin(b.x) + cos(a.x)*cos(b.x)*cos(a.y-b.y);
    return acos(cos_l);
}


/*
 * Calculate the disent angle for two points
 */
double desent_angle(const vec2d_t &start, const vec2d_t &end);


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
