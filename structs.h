#ifndef STRUCTURES_H
#define STRUCTURES_H


struct vec3d_t {
    double x, y, z;

    vec3d_t operator*(double a) const;
    double operator*(const vec3d_t &vec) const;
    vec3d_t operator-(const vec3d_t &vec) const;
    vec3d_t operator+(const vec3d_t &vec) const;

    double length() const;
};


#endif
