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
    // in sec
    ob_t = new double[data_N] {
        1.7,
        2.7,
        2.3,
        1.7,
        2.6,
        2.0
    };


    // Prepreocess
    for (int i = 0; i < data_N; i++) {
        k_z0[i] = 1;
        k_h0[i] = 1;
        k_zb[i] = 1;
        k_hb[i] = 1;
        k_a[i]  = 1;

        // Fill in pos_x and pos_y
        pos_2d[i] = geo_to_xy(ob_lat[i]*PI/180, ob_lon[i]*PI/180);
        // Translate
        ob_a[i]  *= PI / 180;
        ob_zb[i] *= PI / 180;
        ob_hb[i] *= PI / 180;
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
        k_zb[i] = 1;
        k_hb[i] = 1;
        k_a[i]  = 1;
    }
    k_count = 0;

    double mean_error = rate_answer(answer) / (data_N * 5);

    //__m128d err;
    //__m128d max_err = _mm_set1_pd(mean_error * max_error_k);
    //__m128d K_total = _mm_set_pd(data_N*5, 0);
    //static const __m128d one = _mm_set1_pd(1);
    //static const __m128d zero = _mm_set1_pd(0);
    //for (int i = 0; i+1 < data_N; i+=2) {

        //err = angle_delta_sq_pd(ob_z0 + i, ex_z0 + i);
        //err = _mm_sub_pd(max_err, err);
        //err = _mm_blendv_pd(one, zero, err);
        //K_total = _mm_sub_pd(K_total, err);
        //_mm_store_pd(k_z0 + i, err);

        //err = angle_delta_sq_pd(ob_h0 + i, ex_h0 + i);
        //err = _mm_sub_pd(max_err, err);
        //err = _mm_blendv_pd(one, zero, err);
        //K_total = _mm_sub_pd(K_total, err);
        //_mm_store_pd(k_h0 + i, err);

        //err = angle_delta_sq_pd(ob_zb + i, ex_zb + i);
        //err = _mm_sub_pd(max_err, err);
        //err = _mm_blendv_pd(one, zero, err);
        //K_total = _mm_sub_pd(K_total, err);
        //_mm_store_pd(k_zb + i, err);

        //err = angle_delta_sq_pd(ob_hb + i, ex_hb + i);
        //err = _mm_sub_pd(max_err, err);
        //err = _mm_blendv_pd(one, zero, err);
        //K_total = _mm_sub_pd(K_total, err);
        //_mm_store_pd(k_hb + i, err);

        //err = angle_delta_sq_pd(ob_a + i, ex_a + i);
        //err = _mm_sub_pd(max_err, err);
        //err = _mm_blendv_pd(one, zero, err);
        //K_total = _mm_sub_pd(K_total, err);
        //_mm_store_pd(k_a + i, err);
    //}
    
    //double x[2];
    //_mm_storeu_pd(x, K_total);
    //k_count = x[0] + x[1];

    double max_error = mean_error * max_error_k;
    for (int i = 0; i < data_N; i++) {
        k_z0[i] = !(pow(angle_delta(ex_z0[i], ob_z0[i]), 2) > max_error);
        k_h0[i] = !(pow(angle_delta(ex_h0[i], ob_h0[i]), 2) > max_error);
        k_zb[i] = !(pow(angle_delta(ex_zb[i], ob_zb[i]), 2) > max_error);
        k_hb[i] = !(pow(angle_delta(ex_hb[i], ob_hb[i]), 2) > max_error);
        k_a[i]  = !(pow(angle_delta(ex_a[i],  ob_a[i]),  2) > max_error);

        k_count += !k_z0[i];
        k_count += !k_h0[i];
        k_count += !k_zb[i];
        k_count += !k_hb[i];
        k_count += !k_a[i];
    }
}


/*
 * Returns square-error of given answer
 */
double data_t::rate_answer(const answer_t &answer) {
    //ex_process(*this, answer);
    //return ex_rate(*this);
    
    process_answer(answer);

    __m128d error = _mm_setzero_pd();
    __m128d err;
    __m128d K;

    for (int i = 0; i+1 < data_N; i+=2) {

        err = angle_delta_sq_pd(ob_z0 + i, ex_z0 + i);
        K = _mm_load_pd(k_z0 + i);
        err = _mm_mul_pd(err, K);
        error = _mm_add_pd(error, err);

        err = angle_delta_sq_pd(ob_h0 + i, ex_h0 + i);
        K = _mm_load_pd(k_h0 + i);
        err = _mm_mul_pd(err, K);
        error = _mm_add_pd(error, err);

        err = angle_delta_sq_pd(ob_zb + i, ex_zb + i);
        K = _mm_load_pd(k_zb + i);
        err = _mm_mul_pd(err, K);
        error = _mm_add_pd(error, err);

        err = angle_delta_sq_pd(ob_hb + i, ex_hb + i);
        K = _mm_load_pd(k_hb + i);
        err = _mm_mul_pd(err, K);
        error = _mm_add_pd(error, err);

        err = angle_delta_sq_pd(ob_a + i, ex_a + i);
        K = _mm_load_pd(k_a + i);
        err = _mm_mul_pd(err, K);
        error = _mm_add_pd(error, err);
    }

    double x[2];
    _mm_storeu_pd(x, error);
    return x[0] + x[1];
}



void data_t::process_answer(const answer_t &answer) {
    //ex_process(*this, answer);
    for (int i = 0; i < data_N; i++) {
        // TODO: non-linear tau
        double t=ob_t[i];

        // Relative position of flash
        double x0i = answer.x0 - pos_2d[i].x;
        double y0i = answer.y0 - pos_2d[i].y;
        double z0i = answer.z0 - ob_height[i];

        // Relative begining position
        double xbi = x0i + answer.kx*t;
        double ybi = y0i + answer.ky*t;
        double zbi = z0i + answer.kz*t;

        // Azimuth and altitude of flash
        ex_z0[i] = azimuth(x0i, y0i);
        ex_h0[i] = altitude(x0i, y0i, z0i);

        // Azimuth and altitude of the begining
        ex_zb[i] = azimuth(xbi, ybi);
        ex_hb[i] = altitude(xbi, ybi, zbi);
        
        // Desent angle
        ex_a[i] = desent_angle(ex_hb[i], ex_zb[i], ex_h0[i], ex_z0[i]);
    }
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
void data_t::process_answer(const answer_t &answer, processed_answer &dest) {
    for (int i = 0; i < data_N; i++) {
        // TODO: non-linear tau
        double t=ob_t[i];

        // Relative position of flash
        double x0i = answer.x0 - pos_2d[i].x;
        double y0i = answer.y0 - pos_2d[i].y;
        double z0i = answer.z0 - ob_height[i];

        // Relative begining position
        double xbi = x0i + answer.kx*t;
        double ybi = y0i + answer.ky*t;
        double zbi = z0i + answer.kz*t;

        // Azimuth and altitude of flash
        dest.z0[i] = azimuth(x0i, y0i);
        dest.h0[i] = altitude(x0i, y0i, z0i);

        // Azimuth and altitude of the begining
        dest.zb[i] = azimuth(xbi, ybi);
        dest.hb[i] = altitude(xbi, ybi, zbi);
        
        // Desent angle
        dest.a[i] = desent_angle(dest.hb[i], dest.zb[i], dest.h0[i], dest.z0[i]);
    }
}

//processed_answer::processed_answer(const data_t &data, const answer_t &answer) {
    //process(data, answer);
//}

//void processed_answer::process(const data_t &data, const answer_t &answer) {
    //for (int i = 0; i < data_N; i++) {
        //// TODO: non-linear tau
        //double t=data.ob_t[i];

        //// Relative position of flash
        //double x0i = answer.x0 - data.pos_2d[i].x;
        //double y0i = answer.y0 - data.pos_2d[i].y;
        //double z0i = answer.z0 - data.ob_height[i];

        //// Relative begining position
        //double xbi = x0i + answer.kx*t;
        //double ybi = y0i + answer.ky*t;
        //double zbi = z0i + answer.kz*t;

        //// Azimuth and altitude of flash
        //z0[i] = azimuth(x0i, y0i);
        //h0[i] = altitude(x0i, y0i, z0i);

        //// Azimuth and altitude of the begining
        //zb[i] = azimuth(xbi, ybi);
        //hb[i] = altitude(xbi, ybi, zbi);
        
        //// Desent angle
        //a[i] = desent_angle(hb[i], zb[i], h0[i], z0[i]);
    //}

//}

/*
 * Returns square-error of given answer
 */
//double processed_answer::rate(const data_t &data) {

    //__m128d error = _mm_setzero_pd();
    //__m128d err;
    //__m128d K;

    //for (int i = 0; i+1 < data_N; i+=2) {

        //err = angle_delta_sq_pd(data.ob_z0 + i, z0 + i);
        //K = _mm_load_pd(data.k_z0 + i);
        //err = _mm_mul_pd(err, K);
        //error = _mm_add_pd(error, err);

        //err = angle_delta_sq_pd(data.ob_h0 + i, h0 + i);
        //K = _mm_load_pd(data.k_h0 + i);
        //err = _mm_mul_pd(err, K);
        //error = _mm_add_pd(error, err);

        //err = angle_delta_sq_pd(data.ob_zb + i, zb + i);
        //K = _mm_load_pd(data.k_zb + i);
        //err = _mm_mul_pd(err, K);
        //error = _mm_add_pd(error, err);

        //err = angle_delta_sq_pd(data.ob_hb + i, hb + i);
        //K = _mm_load_pd(data.k_hb + i);
        //err = _mm_mul_pd(err, K);
        //error = _mm_add_pd(error, err);

        //err = angle_delta_sq_pd(data.ob_a + i, a + i);
        //K = _mm_load_pd(data.k_a + i);
        //err = _mm_mul_pd(err, K);
        //error = _mm_add_pd(error, err);
    //}

    //double x[2];
    //_mm_storeu_pd(x, error);
    //return x[0] + x[1];
//}
