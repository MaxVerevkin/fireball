#include "data.h"

#include "utils.h"
#include "simd.h"

#include "hyperparams.h"

data_t::data_t() {

    // Prepreocess
    for (int i = 0; i < data_N; i++) {
        // data_Ne is a sum of all ob_e
        data_Ne += ob_e[i];

        // Initialize K=1
        k_z0[i] = 1;
        k_h0[i] = 1;
        k_zb[i] = 1;
        k_hb[i] = 1;
        k_a[i] = 1;

        // Fill in pos_x and pos_y
        pos_2d[i] = geo_to_xy(ob_lat[i]*PI/180, ob_lon[i]*PI/180);

        // Translate
        ob_z0[i] *= PI / 180;
        ob_h0[i] *= PI / 180;
        ob_zb[i] *= PI / 180;
        ob_hb[i] *= PI / 180;
        ob_a[i] *= PI / 180;
    }
}


/*
 * Sets K=0 to all data which square-error is
 * max_error_k times greater than mean square-error
 */
void data_t::eliminate_inconsistent_flash_data(const vec3d_t &pos) {
    k_count = 0;
    for (int i = 0; i < data_N; i++) {
        k_z0[i] = 1;
        k_h0[i] = 1;
    }

    double mean_error = rate_flash_pos(pos, ex_data) / (data_Ne * 2);
    double max_error = mean_error * MAX_ERROR;

    for (int i = 0; i < data_N; i++) {
        k_z0[i] = !(pow(angle_delta(ex_data.z0[i], ob_z0[i]), 2) > max_error);
        k_h0[i] = !(pow(angle_delta(ex_data.h0[i], ob_h0[i]), 2) > max_error);

        k_count += !k_z0[i] * ob_e[i];
        k_count += !k_h0[i] * ob_e[i];
    }
}
void data_t::eliminate_inconsistent_traj_data(const vec3d_t &flash, const vec3d_t params) {
    k_count = 0;
    for (int i = 0; i < data_N; i++) {
        k_zb[i] = 1;
        k_hb[i] = 1;
        k_a[i] = 1;
    }

    double mean_error = rate_flash_traj(flash, params, ex_data) / (data_Ne * 3);
    double max_error = mean_error * MAX_ERROR;

    for (int i = 0; i < data_N; i++) {
        k_zb[i] = !(pow(angle_delta(ex_data.zb[i], ob_zb[i]), 2) > max_error);
        k_hb[i] = !(pow(angle_delta(ex_data.hb[i], ob_hb[i]), 2) > max_error);
        k_a[i] = !(pow(angle_delta(ex_data.a[i], ob_a[i]), 2) > max_error);

        k_count += !k_zb[i] * ob_e[i];
        k_count += !k_hb[i] * ob_e[i];
        k_count += !k_a[i] * ob_e[i];
    }
}


/*
 * Returns square-error of given answer
 */
double data_t::rate_flash_pos(const vec3d_t &pos, processed_answer &dest) {
    process_flash_pos(pos, dest);

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
double data_t::rate_flash_traj(const vec3d_t &flash, const vec3d_t &params, processed_answer &dest) {
    process_flash_traj(flash, params, dest);

    __m128d error = _mm_setzero_pd();
    __m128d err;
    __m128d K;

    for (int i = 0; i+1 < data_N; i+=2) {
        err = angle_delta_sq_pd(ob_zb + i, dest.zb + i);
        K = _mm_load_pd(k_zb + i);
        err = _mm_mul_pd(err, K);
        K = _mm_load_pd(ob_e + i);
        err = _mm_mul_pd(err, K);
        error = _mm_add_pd(error, err);

        err = angle_delta_sq_pd(ob_hb + i, dest.hb + i);
        K = _mm_load_pd(k_hb + i);
        err = _mm_mul_pd(err, K);
        K = _mm_load_pd(ob_e + i);
        err = _mm_mul_pd(err, K);
        error = _mm_add_pd(error, err);

        err = angle_delta_sq_pd(ob_a + i, dest.a + i);
        K = _mm_load_pd(k_a + i);
        err = _mm_mul_pd(err, K);
        K = _mm_load_pd(ob_e + i);
        err = _mm_mul_pd(err, K);
        error = _mm_add_pd(error, err);
    }

    double x[2];
    _mm_storeu_pd(x, error);
    return x[0] + x[1];
}


/*
 * Return sigma (standard deviation)
 * of a given answer.
 */
vec3d_t data_t::sigma_flash_pos(const vec3d_t &pos) {
    double sigma_x = 0;
    double sigma_y = 0;
    double sigma_z = 0;
    for (int i = 0; i < data_N; i++) {
        double x_rel = pos.x - pos_2d[i].x;
        double y_rel = pos.y - pos_2d[i].y;
        double z_rel = pos.z - ob_height[i];
        double tan_z0 = tan(ob_z0[i]);
        sigma_x += pow(y_rel * tan_z0 - x_rel, 2) * k_z0[i] * ob_e[i];
        sigma_y += pow(x_rel / tan_z0 - y_rel, 2) * k_z0[i] * ob_e[i];
        sigma_z += pow(sqrt(x_rel*x_rel + y_rel*y_rel) * tan(ob_h0[i]) - z_rel, 2) * k_h0[i] * ob_e[i];
    }

    sigma_x /= data_N * (data_Ne - 1);
    sigma_y /= data_N * (data_Ne - 1);
    sigma_z /= data_N * (data_Ne - 1);

    return vec3d_t {sqrt(sigma_x), sqrt(sigma_y), sqrt(sigma_z)};
}
vec3d_t data_t::sigma_flash_traj(const vec3d_t &flash, const vec3d_t &params, processed_answer &dest) {
    process_flash_traj(flash, params, dest);

    double sigma_x = 0;
    double sigma_y = 0;
    double sigma_z = 0;
    for (int i = 0; i < data_N; i++) {
        double x_rel = flash.x - pos_2d[i].x + params.x * dest.t[i];
        double y_rel = flash.y - pos_2d[i].y + params.y * dest.t[i];
        double z_rel = flash.z - ob_height[i] + params.z * dest.t[i];
        double tan_zb = tan(ob_zb[i]);
        sigma_x += pow((y_rel * tan_zb - x_rel) / dest.t[i], 2) * k_zb[i] * ob_e[i];
        sigma_y += pow((x_rel / tan_zb - y_rel) / dest.t[i], 2) * k_zb[i] * ob_e[i];
        sigma_z += pow((sqrt(x_rel*x_rel + y_rel*y_rel) * tan(ob_hb[i]) - z_rel) / dest.t[i], 2) * k_hb[i] * ob_e[i];
    }

    sigma_x /= data_N * (data_N - 1);
    sigma_y /= data_N * (data_N - 1);
    sigma_z /= data_N * (data_N - 1);

    return vec3d_t {sqrt(sigma_x), sqrt(sigma_y), sqrt(sigma_z)};
}


/*
 * Fills in the 'expected' values.
 */
void data_t::process_flash_pos(const vec3d_t &pos, processed_answer &dest) {
    for (int i = 0; i < data_N; i++) {
        // Relative position of flash
        vec3d_t rel_flash;
        rel_flash.x = pos.x - pos_2d[i].x;
        rel_flash.y = pos.y - pos_2d[i].y;
        rel_flash.z = pos.z - ob_height[i];

        // Azimuth and altitude of flash
        dest.z0[i] = azimuth(rel_flash.x, rel_flash.y);
        dest.h0[i] = altitude(rel_flash);
    }
}

double data_t::traj_error(const vec3d_t &rel_flash, const vec3d_t params, double t, int i, double best_error, processed_answer &dest) {
    // Relative begining position
    vec3d_t rel_begin;
    rel_begin.x = rel_flash.x + params.x*t;
    rel_begin.y = rel_flash.y + params.y*t;
    rel_begin.z = rel_flash.z + params.z*t;
    // Azimuth and altitude of the begining
    double zbi = azimuth(rel_begin.x, rel_begin.y);
    double hbi = altitude(rel_begin);
    // Desent angle
    double ai = desent_angle(hbi, zbi, dest.h0[i], dest.z0[i]);

    // Error
    double error = 0;
    error += pow(angle_delta(zbi, ob_zb[i]), 2) * k_zb[i] * ob_e[i];
    error += pow(angle_delta(hbi, ob_hb[i]), 2) * k_hb[i] * ob_e[i];
    error += pow(angle_delta(ai, ob_a[i]), 2) * k_a[i] * ob_e[i];

     if (error < best_error) {
         dest.zb[i] = zbi;
         dest.hb[i] = hbi;
         dest.a[i] = ai;
         dest.t[i] = t;
     }


    return error;
}

void data_t::process_flash_traj(const vec3d_t &flash, const vec3d_t &params, processed_answer &dest) {
    for (int i = 0; i < data_N; i++) {
        // Relative position of flash
        vec3d_t rel_flash;
        rel_flash.x = flash.x - pos_2d[i].x;
        rel_flash.y = flash.y - pos_2d[i].y;
        rel_flash.z = flash.z - ob_height[i];

        double min = 1;
        double max = 4;
        double t = 2.5;
        double best_error = traj_error(rel_flash, params, t, i, INFINITY, dest);
        for (int j = 0; j < T_SEARCH_DEPTH; j++) {
            double e1 = traj_error(rel_flash, params, t-.0001, i, best_error, dest);
            double e2 = traj_error(rel_flash, params, t+.0001, i, best_error, dest);
            if (e1 < e2) {
                max = t;
                t = (min + t) / 2;
                best_error = e1;
            } else {
                min = t;
                t = (t + max) / 2;
                best_error = e2;
            }
        }
    }
}
