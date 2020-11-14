#include "utils.h"
#include <math.h>


/*
 * Returns 1 if observation is acceptable.
 */
double is_observation_correct(const vec3d_t &flash,
        const vec3d_t &traj,
        const vec3d_t &observer,
        const vec3d_t &observation) {

    vec3d_t flash_rel = flash - observer;
    vec3d_t perpendic = flash_rel.xprod(traj);
    perpendic = traj.xprod(perpendic);

    return perpendic * observation > 0;
}


/*
 * Translate point from polar to xyz and back
 */
vec3d_t polar_to_xyz(const vec2d_t &polar) {
    vec3d_t retval;
    double xy = cos(polar.x);
    retval.x = xy * sin(polar.y);
    retval.y = xy * cos(polar.y);
    retval.z = sin(polar.x);
    return retval;
}
vec2d_t xyz_to_polar(const vec3d_t &xyz) {
    vec2d_t retval;
    retval.x = atan(xyz.z / xyz.to2d().length());
    retval.y = atan2(xyz.x, xyz.y);
    if (retval.y < 0)
        retval.y += 2*PI;
    return retval;
}


/*
 * Translate geographical location to XYZ and back.
 */
vec3d_t geo_to_xyz(const vec2d_t &geo) {
    vec3d_t pos;
    double xy = cos(geo.x);
    pos.x = xy * cos(geo.y);
    pos.y = xy * sin(geo.y);
    pos.z = sin(geo.x);
    return pos;
}
vec3d_t geo_to_xyz(const vec3d_t &geo) {
    vec3d_t pos;
    double r = EARTH_R + geo.z;
    double xy = r * cos(geo.x);
    pos.x = xy * cos(geo.y);
    pos.y = xy * sin(geo.y);
    pos.z = r * sin(geo.x);
    return pos;
}
vec3d_t xyz_to_geo(const vec3d_t &xyz) {
    vec3d_t pos;
    pos.z = xyz.length() - EARTH_R;
    pos.x = atan(xyz.z / xyz.to2d().length());
    pos.y = atan2(xyz.y, xyz.x);
    return pos;
}


/*
 * Translate geographical location to XYZ and back
 */
vec3d_t global_to_local(const vec3d_t &vec, const vec2d_t pos_geo) {
    vec3d_t north = north_vec(pos_geo);
    vec3d_t east = east_vec(pos_geo.y);
    vec3d_t normal = normal_vec(pos_geo);

    vec3d_t ex = {east.x, north.x, normal.x};
    vec3d_t ey = {east.y, north.y, normal.y};
    vec3d_t ez = {east.z, north.z, normal.z};

    return ex*vec.x + ey*vec.y + ez*vec.z;
}
vec3d_t local_to_global(const vec3d_t &vec, const vec2d_t pos_geo) {
    vec3d_t north = north_vec(pos_geo);
    vec3d_t east = east_vec(pos_geo.y);
    vec3d_t normal = normal_vec(pos_geo);

    return east*vec.x + north*vec.y + normal*vec.z;
}


/*
 * Calculates height of flash above the sea, given
 * geolocation of observer, geolocation of flash
 * and altitude angle and back.
 */
double altitude_to_height(vec3d_t p, vec2d_t p0, double h) {
    double cos_l = sin(p.x)*sin(p0.x) + cos(p.x)*cos(p0.x)*cos(p.y-p0.y);
    double sin_l = sqrt(1 - cos_l*cos_l);
    
    return (1/(cos_l - tan(h) * sin_l) - 1) * (EARTH_R+p.z) + p.z;
}
double height_to_altitude(vec3d_t p, vec3d_t p0) {
    double cos_l = sin(p.x)*sin(p0.x) + cos(p.x)*cos(p0.x)*cos(p.y-p0.y);
    double sin_l = sqrt(1 - cos_l*cos_l);

    return atan((cos_l - (EARTH_R+p.z)/(EARTH_R+p0.z)) / sin_l);
}



/*
 * Calculate the disent angle for two points
 */
double desent_angle(const vec2d_t &start, const vec2d_t &end) {
    // Delta z
    double dz = angle_delta(start.y, end.y);

    // Compute l
    double cos_l = sin(start.x)*sin(end.x) + cos(start.x)*cos(end.x)*cos(dz);
    double sin_l = sqrt(1 - cos_l*cos_l);

    // Compute angle.
    double a = acos((sin(end.x) - sin(start.x)*cos_l) / (cos(start.x)*sin_l));

    if (dz == 0) {
        if (start.x < end.x)
            return 0;
        return PI;
    }

    // Undefined angle
    if (isnan(a))
        return -1;


    if (dz > 0)
        return a;
    return 2*PI - a;

    //int left = dz < 0;
    //return 2*PI*left - (left*2 - 1) * a;
}


// Calculate normal
vec3d_t normal_vec(const vec2d_t &p) {
    vec3d_t normal;
    normal.x = cos(p.x) * cos(p.y);
    normal.y = cos(p.x) * sin(p.y);
    normal.z = sin(p.x);
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

