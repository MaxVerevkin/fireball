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
 * Macros
 */
#define DEG(x) (x*180./PI)
#define RAD(x) (x*PI/180.)


/*
 * Returns 1 if observation is acceptable.
 */
double is_observation_correct(const vec3d_t &flash,
        const vec3d_t &traj,
        const vec3d_t &observer,
        const vec3d_t &observation);


/*
 * Translate point from polar to xyz and back
 */
vec3d_t polar_to_xyz(const vec2d_t &polar);
vec2d_t xyz_to_polar(const vec3d_t &xyz);


/*
 * Translate geographical location to XYZ and back
 */
vec3d_t geo_to_xyz(const vec2d_t &geo);
vec3d_t geo_to_xyz(const vec3d_t &geo);
vec3d_t xyz_to_geo(const vec3d_t &xyz);


/* 
 * Translate global vector to local and back
 */
vec3d_t global_to_local(const vec3d_t &vec, const vec2d_t pos_geo);
vec3d_t local_to_global(const vec3d_t &vec, const vec2d_t pos_geo);



/*
 * Calculate the height of flash above the sea, given
 * geo-location of observer, geo-location of flash
 * and altitude angle. And vise versa.
 */
double altitude_to_height(vec3d_t p, vec2d_t p0, double h);
double height_to_altitude(vec3d_t p, vec3d_t p0);


/*
 * Calculate the absolute difference of 2 angles.
 */
inline double angle_delta(double a1, double a2) {
    //return atan2(sin(a2-a1), cos(a2-a1));

    //double a = a2 - a1;
    //return mod((a + PI), 2*PI) - PI;

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
double descent_angle(const vec2d_t &start, const vec2d_t &end);


// Calculate normal
vec3d_t normal_vec(const vec2d_t &p);

// Genetare vector pointing to Notrh
vec3d_t north_vec(const vec2d_t &ob);
// Genetare vector pointing to East
vec3d_t east_vec(double lon);



#endif
