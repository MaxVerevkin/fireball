#ifndef UTILS_H
#define UTILS_H


#include <math.h>

#include "structs.h"
#include "simd.h"

#define PI 3.141592653589793
#define EARTH_R 6371000.

// Translate geographical location to XYZ.
vec3d_t geo_to_xyz(double lat, double lon, double z);
vec3d_t geo_to_xyz(const vec3d_t &pos);

// Calculate delta of two angles.
double angle_delta(double a1, double a2);
__m128d angle_delta_sq_pd(double *addr1, double *addr2);

// Calculate the azimuth given XY coordinates.
double azimuth(const vec3d_t &observer, const vec3d_t &point, const vec3d_t normal, double ob_lat, double ob_lon, double r);

// Calculate the disent angle for the begining of the path.
double desent_angle(double h, double z, double h0, double z0);

// Calculate normal
vec3d_t normal_vec(double lat, double lon);

// Genetare vector pointing to Notrh
vec3d_t north_vec(double lat, double lon);
// Genetare vector pointing to East
vec3d_t east_vec(double lat, double lon);

// Translate global vector to local
vec3d_t global_to_local(vec3d_t vec, double lat, double lon);


#endif
