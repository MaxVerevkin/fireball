#include "data.h"
#include <cstdio>

#include "utils.h"
#include "simd.h"

#include "hyperparams.h"

data_t::data_t() {

    // Prepreocess
    for (int i = 0; i < data_N; i++) {

        // data_Ne is a sum of all ob_e
        data_Ne += ob_e[i];

        // Height relative to earth center
        r[i] = ob_height[i] + EARTH_R;

        // Initialize K=1
        k_z0[i] = 1;
        k_h0[i] = 1;
        k_zb[i] = 1;
        k_hb[i] = 1;
        k_a[i] = 1;

        // Translate
        ob_lat[i] *= PI / 180;
        ob_lon[i] *= PI / 180;
        ob_z0[i] *= PI / 180;
        ob_h0[i] *= PI / 180;
        ob_zb[i] *= PI / 180;
        ob_hb[i] *= PI / 180;
        ob_a[i] *= PI / 180;
        
        // Calculate in 3D position
        ob_pos[i] = geo_to_xyz(ob_lat[i], ob_lon[i], ob_height[i]);

        // Calculate normal
        normal[i].x = cos(ob_lat[i]) * cos(ob_lon[i]);
        normal[i].y = cos(ob_lat[i]) * sin(ob_lon[i]);
        normal[i].z = sin(ob_lat[i]);
    }
}


/*
 * Sets K=0 to all data which square-error is
 * max_error_k times greater than mean square-error
 */
void data_t::eliminate_inconsistent_flash_data(const vec3d_t &flash) {
    k_count = 0;
    for (int i = 0; i < data_N; i++) {
        k_z0[i] = 1;
        k_h0[i] = 1;
    }

    double mean_error = rate_flash_pos(flash, ex_data) / (data_Ne * 2);
    double max_error = mean_error * MAX_ERROR;

    for (int i = 0; i < data_N; i++) {
        k_z0[i] = !(pow(angle_delta(ex_data.z0[i], ob_z0[i]), 2) > max_error);
        k_h0[i] = !(pow(angle_delta(ex_data.h0[i], ob_h0[i]), 2) > max_error);

        k_count += !k_z0[i] * ob_e[i];
        k_count += !k_h0[i] * ob_e[i];
    }
}
//void data_t::eliminate_inconsistent_traj_data(const vec3d_t &flash, const vec3d_t params) {
    //k_count = 0;
    //for (int i = 0; i < data_N; i++) {
        //k_zb[i] = 1;
        //k_hb[i] = 1;
        //k_a[i] = 1;
    //}

    //double mean_error = rate_flash_traj(flash, params, ex_data) / (data_Ne * 3);
    //double max_error = mean_error * MAX_ERROR;

    //for (int i = 0; i < data_N; i++) {
        //k_zb[i] = !(pow(angle_delta(ex_data.zb[i], ob_zb[i]), 2) > max_error);
        //k_hb[i] = !(pow(angle_delta(ex_data.hb[i], ob_hb[i]), 2) > max_error);
        //k_a[i] = !(pow(angle_delta(ex_data.a[i], ob_a[i]), 2) > max_error);

        //k_count += !k_zb[i] * ob_e[i];
        //k_count += !k_hb[i] * ob_e[i];
        //k_count += !k_a[i] * ob_e[i];
    //}
//}


/*
 * Returns square-error of given answer
 */
double data_t::rate_flash_pos(const vec3d_t &flash, processed_answer &dest) {
    process_flash_pos(flash, dest);

    __m128d error = _mm_setzero_pd();
    __m128d err;
    __m128d K;

    for (int i = 0; i+1 < data_N; i+=2) {
        err = angle_delta_sq_pd(ob_z0 + i, dest.z0 + i); // err = sq_delta(ob_z0, dest.z0)
        K = _mm_load_pd(k_z0 + i);                       // K = k_z0
        err = _mm_mul_pd(err, K);                        // err *= K
        K = _mm_load_pd(ob_e + i);                       // K = ob_e
        err = _mm_mul_pd(err, K);                        // err *= K
        error = _mm_add_pd(error, err);                  // error += err

        err = angle_delta_sq_pd(ob_h0 + i, dest.h0 + i);
        K = _mm_load_pd(k_h0 + i);
        err = _mm_mul_pd(err, K);
        K = _mm_load_pd(ob_e + i);
        err = _mm_mul_pd(err, K);
        error = _mm_add_pd(error, err);
    }

    double x[2];
    _mm_storeu_pd(x, error);
    return x[0] + x[1];
}
//double data_t::rate_flash_traj(const vec3d_t &flash, const vec3d_t &params, processed_answer &dest) {
    //process_flash_traj(flash, params, dest);

    //__m128d error = _mm_setzero_pd();
    //__m128d err;
    //__m128d K;

    //for (int i = 0; i+1 < data_N; i+=2) {
        //err = angle_delta_sq_pd(ob_zb + i, dest.zb + i);
        //K = _mm_load_pd(k_zb + i);
        //err = _mm_mul_pd(err, K);
        //K = _mm_load_pd(ob_e + i);
        //err = _mm_mul_pd(err, K);
        //error = _mm_add_pd(error, err);

        //err = angle_delta_sq_pd(ob_hb + i, dest.hb + i);
        //K = _mm_load_pd(k_hb + i);
        //err = _mm_mul_pd(err, K);
        //K = _mm_load_pd(ob_e + i);
        //err = _mm_mul_pd(err, K);
        //error = _mm_add_pd(error, err);

        //err = angle_delta_sq_pd(ob_a + i, dest.a + i);
        //K = _mm_load_pd(k_a + i);
        //err = _mm_mul_pd(err, K);
        //K = _mm_load_pd(ob_e + i);
        //err = _mm_mul_pd(err, K);
        //error = _mm_add_pd(error, err);
    //}

    //double x[2];
    //_mm_storeu_pd(x, error);
    //return x[0] + x[1];
//}


/*
 * Return sigma (standard deviation)
 * of a given answer.
 */
//vec3d_t data_t::sigma_flash_pos(const vec3d_t &flash) {
    //double sigma_x = 0;
    //double sigma_y = 0;
    //double sigma_z = 0;
    //for (int i = 0; i < data_N; i++) {
        //double x_rel = flash.x - ob_pos[i].x;
        //double y_rel = flash.y - ob_pos[i].y;
        //double z_rel = flash.z - ob_pos[i].z;
        //double tan_z0 = tan(ob_z0[i]);
        //sigma_x += pow(y_rel * tan_z0 - x_rel, 2) * k_z0[i] * ob_e[i];
        //sigma_y += pow(x_rel / tan_z0 - y_rel, 2) * k_z0[i] * ob_e[i];
        //sigma_z += pow(sqrt(x_rel*x_rel + y_rel*y_rel) * tan(ob_h0[i]) - z_rel, 2) * k_h0[i] * ob_e[i];
    //}

    //sigma_x /= data_N * (data_Ne - 1);
    //sigma_y /= data_N * (data_Ne - 1);
    //sigma_z /= data_N * (data_Ne - 1);

    //return vec3d_t {sqrt(sigma_x), sqrt(sigma_y), sqrt(sigma_z)};
//}
//vec3d_t data_t::sigma_flash_traj(const vec3d_t &flash, const vec3d_t &params, processed_answer &dest) {
    //process_flash_traj(flash, params, dest);

    //double sigma_x = 0;
    //double sigma_y = 0;
    //double sigma_z = 0;
    //for (int i = 0; i < data_N; i++) {
        //double x_rel = flash.x - ob_pos[i].x + params.x * dest.t[i];
        //double y_rel = flash.y - ob_pos[i].y + params.y * dest.t[i];
        //double z_rel = flash.z - ob_pos[i].z + params.z * dest.t[i];
        //double tan_zb = tan(ob_zb[i]);
        //sigma_x += pow((y_rel * tan_zb - x_rel) / dest.t[i], 2) * k_zb[i] * ob_e[i];
        //sigma_y += pow((x_rel / tan_zb - y_rel) / dest.t[i], 2) * k_zb[i] * ob_e[i];
        //sigma_z += pow((sqrt(x_rel*x_rel + y_rel*y_rel) * tan(ob_hb[i]) - z_rel) / dest.t[i], 2) * k_hb[i] * ob_e[i];
    //}

    //sigma_x /= data_N * (data_N - 1);
    //sigma_y /= data_N * (data_N - 1);
    //sigma_z /= data_N * (data_N - 1);

    //return vec3d_t {sqrt(sigma_x), sqrt(sigma_y), sqrt(sigma_z)};
//}


/*
 * Fills in the 'expected' values.
 */
void data_t::process_flash_pos(const vec3d_t &flash, processed_answer &dest) {
    for (int i = 0; i < data_N; i++) {
        vec3d_t flash_rel = flash - ob_pos[i];

        double d = flash * normal[i] - r[i];
        double l = flash_rel.length();

        vec3d_t flash_proj_rel = flash - normal[i]*d - ob_pos[i];
        vec3d_t north = vec3d_t {-ob_pos[i].x, -ob_pos[i].y, -ob_pos[i].z+r[i]/sin(ob_lat[i])};
        vec3d_t east = vec3d_t {-sin(ob_lon[i]), cos(ob_lon[i]), 0};

        dest.h0[i] = l == 0 ? PI/2 : asin(d/l);
        dest.z0[i] = azimuth(flash_proj_rel, north, east);
        printf("%i: %f %f\n", i, dest.h0[i]/PI*180, dest.z0[i]/PI*180);
        // Relative position of flash
        //vec3d_t rel_flash;
        //rel_flash.x = flash.x - ob_pos[i].x;
        //rel_flash.y = flash.y - ob_pos[i].y;
        //rel_flash.z = flash.z - ob_height[i];

        // Azimuth and altitude of flash
        //dest.z0[i] = azimuth(rel_flash.x, rel_flash.y);
        //dest.h0[i] = altitude(rel_flash);
    }
}
//void data_t::process_flash_traj(const vec3d_t &flash, const vec3d_t &params, processed_answer &dest) {
    //for (int i = 0; i < data_N; i++) {
        //// Relative position of flash
        //vec3d_t rel_flash;
        //rel_flash.x = flash.x - ob_pos[i].x;
        //rel_flash.y = flash.y - ob_pos[i].y;
        //rel_flash.z = flash.z - ob_height[i];

        //// Binary search for 't'
        //double min = T_SEARCH_MIN;
        //double max = T_SEARCH_MAX;
        //dest.t[i] = (min+max) / 2;
        //for (int j = 0; j < T_SEARCH_DEPTH; j++) {
            //double e1 = traj_error_i(rel_flash, params, dest.t[i]-.0001, i, dest);
            //double e2 = traj_error_i(rel_flash, params, dest.t[i]+.0001, i, dest);
            //if (e1 < e2) {
                //max = dest.t[i];
                //dest.t[i] = (min + dest.t[i]) / 2;
            //} else {
                //min = dest.t[i];
                //dest.t[i] = (dest.t[i] + max) / 2;
            //}
        //}
        //// Actually process current answer
        //process_flash_traj_i(rel_flash, params, dest.t[i], i, dest);
    //}
//}

/*
 * Proceese answer for the trajectory for one observer given 't'.
 */
//void data_t::process_flash_traj_i(const vec3d_t &rel_flash, const vec3d_t params, double t, int i, processed_answer &dest) {
    //// Relative begining position
    //vec3d_t rel_begin;
    //rel_begin.x = rel_flash.x + params.x*t;
    //rel_begin.y = rel_flash.y + params.y*t;
    //rel_begin.z = rel_flash.z + params.z*t;
    //// Azimuth and altitude of the begining
    //dest.zb[i] = azimuth(rel_begin.x, rel_begin.y);
    //dest.hb[i] = altitude(rel_begin);
    //// Desent angle
    //dest.a[i] = desent_angle(dest.hb[i], dest.zb[i], dest.h0[i], dest.z0[i]);
//}

/*
 * Returns square-error the trajectory for current answer given 't'
 */
//double data_t::traj_error_i(const vec3d_t &rel_flash, const vec3d_t params, double t, int i, processed_answer &dest) {
    //process_flash_traj_i(rel_flash, params, t, i, dest);
    //// Error
    //double error = 0;
    //error += pow(angle_delta(dest.zb[i], ob_zb[i]), 2) * k_zb[i] * ob_e[i];
    //error += pow(angle_delta(dest.hb[i], ob_hb[i]), 2) * k_hb[i] * ob_e[i];
    //error += pow(angle_delta(dest.a[i], ob_a[i]), 2) * k_a[i] * ob_e[i];
    //return error;
//}

