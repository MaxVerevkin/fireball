#include "utils.h"
#include <math.h>


// Translate geographical location to XYZ.
vec3d_t geo_to_xyz(const vec3d_t &geo) {
    vec3d_t pos;
    double r = EARTH_R + geo.z;
    double xy = r * cos(geo.x);
    pos.x = xy * cos(geo.y);
    pos.y = xy * sin(geo.y);
    pos.z = r * sin(geo.y);
    return pos;
}


inline __m128d abs_pd(__m128d x) {
    static const __m128d __minus_zero = _mm_set1_pd(-0.);
    return _mm_andnot_pd(__minus_zero, x);
}


// Calculate delta of two angles.
__m128d angle_delta_sq_pd(double *addr1, double *addr2) {
    double t1 = angle_delta(addr1[0], addr2[0]);
    double t2 = angle_delta(addr1[1], addr2[1]);
    t1 *= t1;
    t2 *= t2;
    return _mm_set_pd(t1, t2);

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

double calc_height(vec3d_t p, vec2d_t p0, double h) {
    double cos_l = sin(p.x)*sin(p0.x) + cos(p.x)*cos(p0.x)*cos(p.y-p0.y);
    double sin_l = sqrt(1 - cos_l*cos_l);
    return (1/(cos_l - tan(h) * sin_l) - 1) * (EARTH_R+p.z) + p.z;
}

// Calculate the length of arc on a sphere
double arc_len(double h1, double z1, double h2, double z2) {
    double cos_l = sin(h1)*sin(h2) + cos(h1)*cos(h2)*cos(z1-z2);
    return acos(cos_l);
}

// Calculate the azimuth given XYZ coordinates, North and East vectors.
//double azimuth(const vec3d_t &point, const vec3d_t normal, double ob_lat, double ob_lon, double r) {
double azimuth(const vec3d_t &point, const vec3d_t &normal, const vec2d_t &observer, double r) {
    double inv_point_len = 1. / point.length();

    vec3d_t north = north_vec(observer);
    vec3d_t east = east_vec(observer.y);

    double zc = acos(point * north * inv_point_len / north.length());
    double cos_ze = point * east * inv_point_len / east.length();
    return cos_ze >= 0 ? zc : 2*PI - zc;
}

// Calculate the disent angle for the begining of the path.
double desent_angle(double h1, double z1, double h2, double z2) {
    // Delta z
    double dz = angle_delta(z1, z2);

    // Compute l
    double cos_l = sin(h1)*sin(h2) + cos(h1)*cos(h2)*cos(dz);
    double sin_l = sqrt(1 - cos_l*cos_l);

    // Compute angle.
    double a = acos((sin(h2) - sin(h1)*cos_l) / (cos(h1)*sin_l));

    // Undefined angle
    if (isnan(a))
        return 0;

    int left = dz < 0;
    return 2*PI*left - (left*2 - 1) * a;
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
vec3d_t north_vec(const vec2d_t &ob) {
    vec3d_t north;
    north.x = -sin(ob.x)*cos(ob.y);
    north.y = -sin(ob.x)*sin(ob.y);
    north.z = cos(ob.x);
    return north;
}

// Genetare vector pointing to East
vec3d_t east_vec(double lon) {
    return {-sin(lon), cos(lon), 0};
}

// Translate global vector to local
vec3d_t global_to_local(vec3d_t vec, double lat, double lon) {
    vec3d_t north = north_vec({lat, lon});
    vec3d_t east = east_vec(lon);
    vec3d_t normal = normal_vec(lat, lon);

    vec3d_t ex = {east.x, north.x, normal.x};
    vec3d_t ey = {east.y, north.y, normal.y};
    vec3d_t ez = {east.z, north.z, normal.z};

    return ex*vec.x + ey*vec.y + ez*vec.z;
}
