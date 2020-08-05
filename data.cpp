#include "data.h"

#include "utils.h"
#include "simd.h"


data_t::data_t() {

    // in deg
    ob_lat = new double[data_N] {
        34.87,
        35.13,
        34.72,
        34.73,
        35.17,
        35.36
    };
    // in deg
    ob_lon = new double[data_N] {
        32.92,
        33.30,
        33.04,
        33.69,
        32.78,
        32.89
    };
    // in m
    ob_height = new double[data_N] {
        1200,
        800,
        1100,
        1400,
        700,
        200
    };
    // in deg
    ob_a = new double[data_N] {
        236,
        288,
        241,
        203,
        84,
        227
    };
    // in deg
    ob_zb = new double[data_N] {
        83,
        130,
        78,
        359,
        114,
        124
    };
    // in deg
    ob_hb = new double[data_N] {
        40,
        55,
        48,
        68,
        39,
        37
    };
    // in deg
    ob_z0 = new double[data_N] {
        51,
        52,
        34,
        340,
        83,
        104
    };
    // in deg
    ob_h0 = new double[data_N] {
        14,
        40,
        33,
        14,
        14,
        17
    };


    // Prepreocess
    for (int i = 0; i < data_N; i++) {
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
void data_t::eliminate_inconsistent_flash_data(const vec3d_t &pos, double max_error_k) {
    k_count = 0;

    double mean_error = rate_flash_pos(pos, ex_data) / (data_N * 2);
    double max_error = mean_error * max_error_k;

    for (int i = 0; i < data_N; i++) {
        k_z0[i] = !(pow(angle_delta(ex_data.z0[i], ob_z0[i]), 2) > max_error);
        k_h0[i] = !(pow(angle_delta(ex_data.h0[i], ob_h0[i]), 2) > max_error);

        k_count += !k_z0[i];
        k_count += !k_h0[i];
    }
}
void data_t::eliminate_inconsistent_traj_data(const vec3d_t &flash, const vec3d_t params, double max_error_k) {
    k_count = 0;

    double mean_error = rate_flash_traj(flash, params, ex_data) / (data_N * 3);
    double max_error = mean_error * max_error_k;

    for (int i = 0; i < data_N; i++) {
        k_zb[i] = !(pow(angle_delta(ex_data.zb[i], ob_zb[i]), 2) > max_error);
        k_hb[i] = !(pow(angle_delta(ex_data.hb[i], ob_hb[i]), 2) > max_error);
        k_a[i] = !(pow(angle_delta(ex_data.a[i], ob_a[i]), 2) > max_error);

        k_count += !k_zb[i];
        k_count += !k_hb[i];
        k_count += !k_a[i];
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
        err = angle_delta_sq_pd(ob_z0 + i, dest.z0 + i);
        K = _mm_load_pd(k_z0 + i);
        err = _mm_mul_pd(err, K);
        error = _mm_add_pd(error, err);

        err = angle_delta_sq_pd(ob_h0 + i, dest.h0 + i);
        K = _mm_load_pd(k_h0 + i);
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
        error = _mm_add_pd(error, err);

        err = angle_delta_sq_pd(ob_hb + i, dest.hb + i);
        K = _mm_load_pd(k_hb + i);
        err = _mm_mul_pd(err, K);
        error = _mm_add_pd(error, err);

        err = angle_delta_sq_pd(ob_a + i, dest.a + i);
        K = _mm_load_pd(k_a + i);
        err = _mm_mul_pd(err, K);
        error = _mm_add_pd(error, err);
    }

    double x[2];
    _mm_storeu_pd(x, error);
    return x[0] + x[1];
}


/*
 * Fills in the 'expected' values.
 */
void data_t::process_flash_pos(const vec3d_t &pos, processed_answer &dest) {
    for (int i = 0; i < data_N; i++) {
        // Relative position of flash
        double x0i = pos.x - pos_2d[i].x;
        double y0i = pos.y - pos_2d[i].y;
        double z0i = pos.z - ob_height[i];

        // Azimuth and altitude of flash
        dest.z0[i] = azimuth(x0i, y0i);
        dest.h0[i] = altitude(x0i, y0i, z0i);
    }
}
void data_t::process_flash_traj(const vec3d_t &flash, const vec3d_t &params, processed_answer &dest) {
    for (int i = 0; i < data_N; i++) {
        // FIXME
        double t = 2;

        // Relative position of flash
        double x0i = flash.x - pos_2d[i].x;
        double y0i = flash.y - pos_2d[i].y;
        double z0i = flash.z - ob_height[i];

        // Relative begining position
        double xbi = x0i + params.x*t;
        double ybi = y0i + params.y*t;
        double zbi = z0i + params.z*t;

        // Azimuth and altitude of the begining
        dest.zb[i] = azimuth(xbi, ybi);
        dest.hb[i] = altitude(xbi, ybi, zbi);
        
        // Desent angle
        dest.a[i] = desent_angle(dest.hb[i], dest.zb[i], dest.h0[i], dest.z0[i]);
    }
}
