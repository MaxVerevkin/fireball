#include "utils.h"
#include <math.h>


// Return minimum value
double min_8d(double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8) {
    if (x1 < x2 && x1 < x3 && x1 < x4 && x1 < x5 && x1 < x6 && x1 < x7 && x1 < x8)
        return x1;
    if (x2 < x3 && x2 < x4 && x2 < x5 && x2 < x6 && x2 < x7 && x2 < x8)
        return x2;
    if (x3 < x4 && x3 < x5 && x3 < x6 && x3 < x7 && x3 < x8)
        return x3;
    if (x4 < x5 && x4 < x6 && x4 < x7 && x4 < x8)
        return x4;
    if (x5 < x6 && x5 < x7 && x5 < x8)
        return x5;
    if (x6 < x7 && x6 < x8)
        return x6;
    if (x7 < x8)
        return x7;
    return x8;
}

inline __m128d abs_pd(__m128d x) {
    static const __m128d __minus_zero = _mm_set1_pd(-0.);
    return _mm_andnot_pd(__minus_zero, x);
}


// Calculate delta of two angles.
double angle_delta(double a1, double a2) {
    double delta = a2 - a1;
    double abs_delta = abs(delta);
    if (abs_delta > PI)
        return 2*PI - abs_delta;
    return delta;
}
__m128d angle_delta_sq_pd(double *addr1, double *addr2) {
    static const __m128d __pi = _mm_set1_pd(PI);
    static const __m128d __2pi = _mm_set1_pd(2*PI);

    __m128d a1 = _mm_load_pd(addr1);
    __m128d a2 = _mm_load_pd(addr2);

    a1 = _mm_sub_pd(a1, a2); // a1 = delta
    a2 = abs_pd(a1);         // a2 = abs(delta)

    __m128d a3 = _mm_sub_pd(__2pi, a2); // a3 = 360 - abs(delta)
    __m128d a4 = _mm_sub_pd(__pi, a2);  // a4 = 180 - abs(delta)

    // if (abs(delta) > 180)
    //     a1 = 360 - abs(delta);
    a1 = _mm_blendv_pd(a1, a3, a4);

    return _mm_mul_pd(a1, a1); // return a1^2
}


// Calculate the azimuth given XY coordinates.
double azimuth(double x, double y) {
    // Handle y = 0
    if (y == 0) {
        if (x >= 0)
            return PI / 2.f;
        if (x < 0)
            return 3.f / 2.f * PI;
    }

    // Calculate azimuth
    double z = atan(x / y);

    // I
    if (x >= 0 and y > 0)
        return z;
    // II
    if (x < 0 and y > 0)
        return z + 2*PI;
    // III and IV
    return z + PI;
}


// Calculre the altitude for a 3D point (x,y,z).
double altitude(const vec3d_t &point) {
    if (point.x == 0 and point.y == 0)
        return PI/2;
    return atan(point.z / sqrt(point.x*point.x + point.y*point.y));
}


// Calculate the disent angle for the begining of the path.
double desent_angle(double h, double z, double h0, double z0) {
    // Delta z
    double dz = z0 - z;

    // Compute l, its sin and cos
    double cos_l = sin(h0)*sin(h) + cos(h0)*cos(h)*cos(dz);
    double l = acos(cos_l);
    // Undefined angle (zero path)
    if (l == 0)
        return 0;
    double sin_l = sin(l);

    // Compute angle and its cos
    double sin_a = cos(h0) * sin(dz) / sin_l;
    double a = asin(sin_a);
    double cos_a = (sin(h0) - sin(h)*cos_l) / (cos(h) * sin_l);

    // Decide
    if (sin_a >= 0 and cos_a >= 0)
        return a;
    if (sin_a < 0 and cos_a >= 0)
        return a + 2*PI;
    return PI - a;
}
