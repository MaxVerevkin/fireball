#ifndef STRUCTURES_H
#define STRUCTURES_H


/*
 * Vector of 2 doubles and some opertors.
 */
struct vec2d_t {
    double x, y;

    vec2d_t operator*(double a) const;
    vec2d_t operator+(const vec2d_t &vec) const;
    vec2d_t operator-(const vec2d_t &vec) const;

    double length() const;
};


/*
 * Vector of 3 doubles and some opertors.
 */
struct vec3d_t {
    double x, y, z;

    double operator*(const vec3d_t &vec) const;

    vec3d_t operator*(double a) const;
    vec3d_t operator/(double a) const;

    vec3d_t operator-(const vec3d_t &vec) const;
    vec3d_t operator+(const vec3d_t &vec) const;

    double length() const;
    vec3d_t normalized() const;

    vec2d_t to2d() const;
};


/*
 * Represents data given by observer.
 */
struct data_set_t {
    alignas(16) double *z0;
    alignas(16) double *h0;
    alignas(16) double *zb;
    alignas(16) double *hb;
    alignas(16) double *a;
    alignas(16) double *t;

    data_set_t(int n);
};

#endif
