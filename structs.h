#ifndef STRUCTURES_H
#define STRUCTURES_H


struct vec3d_t {
    double x, y, z;

    double operator*(const vec3d_t &vec) const;

    vec3d_t operator*(double a) const;
    vec3d_t operator/(double a) const;

    vec3d_t operator-(const vec3d_t &vec) const;
    vec3d_t operator+(const vec3d_t &vec) const;

    double length() const;
    vec3d_t normalized() const;
};


#endif
