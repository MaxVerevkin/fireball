#include "data.h"

#include "utils.h"
#include "simd.h"


data_t::data_t() {

    // in deg
    ob_lat = new double[data_N] {
        34.91,
        35.22,
        34.72,
        34.73,
        35.26,
        35.49
    };
    // in deg
    ob_lon = new double[data_N] {
        32.84,
        33.16,  
        32.94,
        33.48,
        32.72,
        32.82,
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
    //ob_a = new double[data_N] {
        //236,
        //288,
        //241,
        //203,
        //84,
        //227
    //};
    //// in deg
    //ob_zb = new double[data_N] {
        //83,
        //130,
        //78,
        //359,
        //114,
        //124
    //};
    //// in deg
    //ob_hb = new double[data_N] {
        //40,
        //55,
        //48,
        //68,
        //39,
        //37
    //};
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

        // Fill in pos_x and pos_y
        pos_2d[i] = geo_to_xy(ob_lat[i]*PI/180, ob_lon[i]*PI/180);
        // Translate
        ob_z0[i] *= PI / 180;
        ob_h0[i] *= PI / 180;
    }
}


/*
 * Sets K=0 to all data which square-error is
 * max_error_k times greater than mean square-error
 */
void data_t::eliminate_inconsistent(const answer_t &answer, double max_error_k) {
    for (int i = 0; i < data_N; i++) {
        k_z0[i] = 1;
        k_h0[i] = 1;
    }
    k_count = 0;

    double mean_error = rate_answer(answer) / (data_N * 5);
    double max_error = mean_error * max_error_k;

    for (int i = 0; i < data_N; i++) {
        k_z0[i] = !(pow(angle_delta(ex_data.z0[i], ob_z0[i]), 2) > max_error);
        k_h0[i] = !(pow(angle_delta(ex_data.h0[i], ob_h0[i]), 2) > max_error);

        k_count += !k_z0[i];
        k_count += !k_h0[i];
    }
}


/*
 * Returns square-error of given answer
 */
double data_t::rate_answer(const answer_t &answer) {
    return rate_answer(answer, ex_data);
}
double data_t::rate_answer(const answer_t &answer, processed_answer &dest) {
    process_answer(answer, dest);

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


/*
 * Fills in the 'expected' values.
 */
void data_t::process_answer(const answer_t &answer) {
    process_answer(answer, ex_data);
}
void data_t::process_answer(const answer_t &answer, processed_answer &dest) {
    for (int i = 0; i < data_N; i++) {
        // Relative position of flash
        double x0i = answer.x0 - pos_2d[i].x;
        double y0i = answer.y0 - pos_2d[i].y;
        double z0i = answer.z0 - ob_height[i];

        // Azimuth and altitude of flash
        dest.z0[i] = azimuth(x0i, y0i);
        dest.h0[i] = altitude(x0i, y0i, z0i);
    }
}
