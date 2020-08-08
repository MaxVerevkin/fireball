#ifndef UTILS_H
#define UTILS_H


#include <math.h>

#include "structs.h"
#include "simd.h"

#define PI 3.141592653589793
#define EARTH_R 6371000.

// Linearly interpolate between a0 and a1
#define lerp(a0, a1, w) ((a0) + (w)*((a1) - (a0))) 

// Translate geographical location to XY.
#define geo_to_xy(lat, lon) (vec3d_t {(EARTH_R * lon * cos(lat)), (EARTH_R * lat), 0})

// Return minimum value
double min_8d(double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8);

// Calculate delta of two angles.
double angle_delta(double a1, double a2);
__m128d angle_delta_sq_pd(double *addr1, double *addr2);

// Calculate the azimuth given XY coordinates.
double azimuth(double x, double y);

// Calculre the altitude for a 3D point (x,y,z).
double altitude(const vec3d_t &point);

// Calculate the disent angle for the begining of the path.
double desent_angle(double h, double z, double h0, double z0);


#endif
