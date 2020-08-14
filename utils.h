#ifndef UTILS_H
#define UTILS_H


#include <math.h>

#include "structs.h"
#include "simd.h"

#define PI 3.141592653589793
#define EARTH_R 6371000.

// Linearly interpolate between a0 and a1
#define lerp(a0, a1, w) ((a0) + (w)*((a1) - (a0))) 

// Translate geographical location to XYZ.
vec3d_t geo_to_xyz(double lat, double lon, double z);
vec3d_t geo_to_xyz(const vec3d_t &pos);

// Calculate delta of two angles.
double angle_delta(double a1, double a2);
__m128d angle_delta_sq_pd(double *addr1, double *addr2);

// Calculate the azimuth given XY coordinates.
double azimuth(const vec3d_t &point, const vec3d_t &north, const vec3d_t east);

// Calculate the disent angle for the begining of the path.
double desent_angle(double h, double z, double h0, double z0);


#endif
