#include "utils.h"
#include <math.h>


// Translate geographical location to XYZ.
vec3d_t geo_to_xyz(double lat, double lon, double z) {
    vec3d_t pos;
    double r = EARTH_R + z;
    double xy = r * cos(lat);
    pos.x = xy * cos(lon);
    pos.y = xy * sin(lon);
    pos.z = r * sin(lat);
    return pos;
}
vec3d_t geo_to_xyz(const vec3d_t &pos) {
    return geo_to_xyz(pos.x, pos.y, pos.z);
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

    // Abs delta
    __m128d delta = abs_pd(a1 - a2);

    a1 = __2pi - delta; // a1 = 360 - delta
    a2 = __pi  - delta; // a2 = 180 - delta

    // if (abs(delta) > 180)
    //     a1 = 360 - abs(delta);
    delta = _mm_blendv_pd(delta, a1, a2);

    return delta * delta; // return a1^2
}

// Calculate the length of arc on a sphere
double arc_len(double h1, double z1, double h2, double z2) {
    double cos_l = sin(h1)*sin(h2) + cos(h1)*cos(h2)*cos(z1-z2);
    return acos(cos_l);
}

// Calculate the azimuth given XYZ coordinates, North and East vectors.
double azimuth(const vec3d_t &point, const vec3d_t normal, double ob_lat, double ob_lon, double r) {
    double inv_point_len = 1. / point.length();

    vec3d_t north = north_vec(ob_lat, ob_lon);
    vec3d_t east = east_vec(ob_lat, ob_lon);

    double zc = acos(point * north * inv_point_len / north.length());
    double cos_ze = point * east * inv_point_len / east.length();
    return cos_ze >= 0 ? zc : 2*PI - zc;
}

// Calculate the disent angle for the begining of the path.
double desent_angle(double h, double z, double h0, double z0) {
    // Delta z
    double dz = angle_delta(z, z0);

    // Compute l
    double cos_l = sin(h0)*sin(h) + cos(h0)*cos(h)*cos(dz);
    double sin_l = sqrt(1 - cos_l*cos_l);

    // Compute angle.
    double a = acos((sin(h0) - sin(h)*cos_l) / (cos(h)*sin_l));

    // Undefined angle
    if (isnan(a))
        return 0;

    return dz < 0 ? 2*PI - a : a;
}

// Calculate normal
vec3d_t normal_vec(double lat, double lon) {
    vec3d_t normal;
    normal.x = cos(lat) * cos(lon);
    normal.y = cos(lat) * sin(lon);
    normal.z = sin(lat);
    return normal;
}

// Genetare vector pointing to Notrh
vec3d_t north_vec(double lat, double lon) {
    return {-sin(lat)*cos(lon), -sin(lat)*sin(lon), cos(lat)};
}

// Genetare vector pointing to East
vec3d_t east_vec(double lat, double lon) {
    return {-sin(lon), cos(lon), 0};
}

// Translate global vector to local
vec3d_t global_to_local(vec3d_t vec, double lat, double lon) {
    vec3d_t north = north_vec(lat, lon);
    vec3d_t east = east_vec(lat, lon);
    vec3d_t normal = normal_vec(lat, lon);

    vec3d_t ex = {east.x, north.x, normal.x};
    vec3d_t ey = {east.y, north.y, normal.y};
    vec3d_t ez = {east.z, north.z, normal.z};

    return ex*vec.x + ey*vec.y + ez*vec.z;
}
