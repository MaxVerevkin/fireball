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
 * Vector of 2 2D vectors.
 */
struct vec2d2_t {
    vec2d_t v1, v2;
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

    vec3d_t xprod(const vec3d_t &vec) const;

    double length() const;
    vec3d_t normalized() const;

    vec2d_t to2d() const;
};

/*
 * Vector of 2 2D vectors.
 */
struct line3d_t {
    vec3d_t start, end;

    vec3d_t vec() const;
};


/*
 * Represents data given by observer.
 */
struct data_set_t {
    double *z0;
    double *h0;
    double *zb;
    double *hb;
    double *a;
    double *t;

    data_set_t(int n);
};

#endif
