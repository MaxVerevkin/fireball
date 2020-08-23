#include "data.h"

#include "utils.h"
#include "simd.h"
#include <cstdio>

#include "hyperparams.h"


// Paralelisation
#ifdef OP_PARALEL
#include <omp.h>
#endif


data_t::data_t(const char *file) {

    // Open file
    FILE *infile = fopen(file, "r");
    if (!infile) {
        fprintf(stderr, "Failed to open file %s\n", file);
        exit(1);
    }
    fscanf(infile, "%d", &data_N);

    // Allocate memory
    ob_lat =    new double[data_N];
    ob_lon =    new double[data_N];
    ob_height = new double[data_N];
    ob_e =      new double[data_N];
    k_z0 =      new double[data_N];
    k_h0 =      new double[data_N];
    k_zb =      new double[data_N];
    k_hb =      new double[data_N];
    k_a =       new double[data_N];
    r =         new double[data_N];
    ob_pos =    new vec3d_t[data_N];
    normal =    new vec3d_t[data_N];
    ob_data =   new data_set_t(data_N);
    ex_data =   new data_set_t(data_N);

    for (int i = 0; i < data_N; i++) {
        // Read from file
        fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                ob_lat+i, ob_lon+i, ob_height+i,
                ob_data->a+i, ob_data->zb+i, ob_data->hb+i,
                ob_data->z0+i, ob_data->h0+i,
                ob_data->t+i, ob_e+i);

        // data_Ne is a sum of all ob_e
        data_Ne += ob_e[i];

        // Height relative to earth center
        r[i] = ob_height[i] + EARTH_R;

        // Translate
        ob_lat[i] *= PI / 180;
        ob_lon[i] *= PI / 180;
        ob_data->z0[i] *= PI / 180;
        ob_data->h0[i] *= PI / 180;
        ob_data->zb[i] *= PI / 180;
        ob_data->hb[i] *= PI / 180;
        ob_data->a[i] *= PI / 180;
        mean_lat += ob_lat[i] / data_N;
        mean_lon += ob_lon[i] / data_N;
        
        // Calculate global 3D position
        ob_pos[i] = geo_to_xyz(ob_lat[i], ob_lon[i], ob_height[i]);

        // Calculate normal
        normal[i] = normal_vec(ob_lat[i], ob_lon[i]);
    }

    // Set k=1
    reset_k_flash();
    reset_k_traj();

    // Close file
    fclose(infile);
}


/*
 * Reset k-values to '1'
 */
void data_t::reset_k_flash() {
    k_count_flash = 0;
    for (int i = 0; i < data_N; i++) {
        k_z0[i] = ob_data->z0[i] >= 0;
        k_h0[i] = ob_data->h0[i] >= 0;

        k_count_flash += !k_z0[i];
        k_count_flash += !k_h0[i];
    }
}
void data_t::reset_k_traj() {
    k_count_traj = 0;
    for (int i = 0; i < data_N; i++) {
        k_zb[i] = ob_data->zb[i] >= 0;
        k_hb[i] = ob_data->hb[i] >= 0;
        k_a[i]  = ob_data->a[i]  >= 0;

        k_count_traj += !k_zb[i];
        k_count_traj += !k_hb[i];
        k_count_traj += !k_a[i];
    }
}


/*
 * Sets K=0 to all data which square-error is
 * max_error_k times greater than mean square-error
 */
void data_t::eliminate_inconsistent_flash_data(const vec3d_t &flash_geo) {
    reset_k_flash();

    double mean_error = rate_flash_pos(flash_geo) / (data_Ne * 2 - k_count_flash);
    double max_error = mean_error * MAX_ERROR;

    for (int i = 0; i < data_N; i++) {
        if (ob_data->z0[i] >= 0) {
            k_z0[i] = !(pow(angle_delta(ex_data->z0[i], ob_data->z0[i]), 2) > max_error);
            k_count_flash += !k_z0[i] * ob_e[i];
        }
        if (ob_data->h0[i] >= 0) {
            k_h0[i] = !(pow(angle_delta(ex_data->h0[i], ob_data->h0[i]), 2) > max_error);
            k_count_flash += !k_h0[i] * ob_e[i];
        }
    }
}
void data_t::eliminate_inconsistent_traj_data(const vec3d_t &flash_geo, const vec3d_t params) {
    reset_k_traj();

    double mean_error = rate_flash_traj(flash_geo, params) / (data_Ne * 3 - k_count_traj);
    double max_error = mean_error * MAX_ERROR;

    for (int i = 0; i < data_N; i++) {
        if (ob_data->zb[i] >= 0) {
            k_zb[i] = !(pow(angle_delta(ex_data->zb[i], ob_data->zb[i]), 2) > max_error);
            k_count_traj += !k_zb[i] * ob_e[i];
        }
        if (ob_data->hb[i] >= 0) {
            k_hb[i] = !(pow(angle_delta(ex_data->hb[i], ob_data->hb[i]), 2) > max_error);
            k_count_traj += !k_hb[i] * ob_e[i];
        }
        if (ob_data->a[i] >= 0) {
            k_a[i] = !(pow(angle_delta(ex_data->a[i], ob_data->a[i]), 2) > max_error);
            k_count_traj += !k_a[i] * ob_e[i];
        }
    }
}

/*
 * Translate flash trajectory vector to local velocity.
 */
vec3d_t data_t::get_flash_vel(const vec3d_t &flash_geo, const vec3d_t &traj) {
    vec3d_t vel = global_to_local(traj, flash_geo.x, flash_geo.y);

    // Calculate mean 'k'
    double mean_k = 0;
    double n = 0;
    for (int i = 0; i < data_N; i++) {
        if (ob_data->t[i] > 0) {
            mean_k += ex_data->t[i] / ob_data->t[i];
            n++;
        }
    }
    mean_k /= n;

    // Calculate mean error
    double mean_error = 0;
    for (int i = 0; i < data_N; i++)
        if (ob_data->t[i] > 0)
            mean_error += pow(mean_k - ex_data->t[i] / ob_data->t[i], 2);
    mean_error /= n;

    // Relculate 'k', ignoring inconsistent
    double k = 0;
    n = 0;
    for (int i = 0; i < data_N; i++) {
        if (ob_data->t[i] > 0) {
            double ki = ex_data->t[i] / ob_data->t[i];
            if (pow(mean_k-ki ,2) < mean_error*MAX_ERROR) {
                k += ki;
                n++;
            }
        }
    }
    k /= n;

    return vel * (-k);
}


/*
 * Normalize aobserver's 't'
 */
void data_t::normalize_t(const vec3d_t &flash_vel) {
    double vel = flash_vel.length();
    for (int i = 0; i < data_N; i++)
        ex_data->t[i] /= vel;
}


/*
 * Returns square-error of given answer
 */
double data_t::rate_flash_pos(const vec3d_t &flash_geo) {
    process_flash_pos(flash_geo);

    __m128d total_error = _mm_setzero_pd();
    for (int i = 0; i+1 < data_N; i+=2) {

        // h0
        __m128d err1 = _mm_load_pd(k_h0 + i) * _mm_set1_pd(DELTA_H0_W);
        err1 *= angle_delta_sq_pd(ob_data->h0 + i, ex_data->h0 + i);

        // z0
        __m128d err2 = _mm_load_pd(k_z0 + i) * _mm_set1_pd(DELTA_Z0_W);
        err2 *= angle_delta_sq_pd(ob_data->z0 + i, ex_data->z0 + i);

        err1 += err2;
        err1 *= _mm_load_pd(ob_e + i);
        total_error += err1;
    }

    double x[2];
    _mm_storeu_pd(x, total_error);
    return x[0] + x[1];
}
double data_t::rate_flash_traj(const vec3d_t &flash_geo, const vec3d_t &params) {
    process_flash_traj(flash_geo, params);

    __m128d total_error = _mm_setzero_pd();
    for (int i = 0; i+1 < data_N; i+=2) {

        // hb
        __m128d err1 = _mm_load_pd(k_hb + i) * _mm_set1_pd(DELTA_HB_W);
        err1 *= angle_delta_sq_pd(ob_data->hb + i, ex_data->hb + i);

        // zb
        __m128d err2 = _mm_load_pd(k_zb + i) * _mm_set1_pd(DELTA_ZB_W);
        err2 *= angle_delta_sq_pd(ob_data->zb + i, ex_data->zb + i);

        // a
        __m128d err3 = _mm_load_pd(k_a + i) * _mm_set1_pd(DELTA_A_W);
        err3 *= angle_delta_sq_pd(ob_data->a + i, ex_data->a + i);

        err1 += err2;
        err1 += err3;
        err1 *= _mm_load_pd(ob_e + i);
        total_error += err1;
    }

    double x[2];
    _mm_storeu_pd(x, total_error);
    return x[0] + x[1];
}


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
//vec3d_t data_t::sigma_flash_traj(const vec3d_t &flash, const vec3d_t &params) {
    //process_flash_traj(flash, params);

    //double sigma_x = 0;
    //double sigma_y = 0;
    //double sigma_z = 0;
    //for (int i = 0; i < data_N; i++) {
        //double x_rel = flash.x - ob_pos[i].x + params.x * ex_data->t[i];
        //double y_rel = flash.y - ob_pos[i].y + params.y * ex_data->t[i];
        //double z_rel = flash.z - ob_pos[i].z + params.z * ex_data->t[i];
        //double tan_zb = tan(ob_zb[i]);
        //sigma_x += pow((y_rel * tan_zb - x_rel) / ex_data->t[i], 2) * k_zb[i] * ob_e[i];
        //sigma_y += pow((x_rel / tan_zb - y_rel) / ex_data->t[i], 2) * k_zb[i] * ob_e[i];
        //sigma_z += pow((sqrt(x_rel*x_rel + y_rel*y_rel) * tan(ob_hb[i]) - z_rel) / ex_data->t[i], 2) * k_hb[i] * ob_e[i];
    //}

    //sigma_x /= data_N * (data_N - 1);
    //sigma_y /= data_N * (data_N - 1);
    //sigma_z /= data_N * (data_N - 1);

    //return vec3d_t {sqrt(sigma_x), sqrt(sigma_y), sqrt(sigma_z)};
//}


/*
 * Fills in the 'expected' values.
 */
void data_t::process_flash_pos(const vec3d_t &flash_geo) {
    // Translate from geo position to 3D point
    vec3d_t flash = geo_to_xyz(flash_geo);

    #ifdef OP_PARALEL
    #pragma omp parallel for
    #endif
    for (int i = 0; i < data_N; i++) {
        // Relative flash_position
        vec3d_t flash_rel = flash - ob_pos[i];

        // Height and length
        double d = flash * normal[i] - r[i];
        double l = flash_rel.length();

        vec3d_t flash_proj_rel = flash - normal[i]*d - ob_pos[i];

        ex_data->h0[i] = l == 0 ? PI/2 : asin(d/l);
        ex_data->z0[i] = azimuth(ob_pos[i], flash_proj_rel, normal[i], ob_lat[i], ob_lon[i], r[i]);
    }
}
void data_t::process_flash_traj(const vec3d_t &flash_geo, const vec3d_t &params) {
    #ifdef OP_PARALEL
    #pragma omp parallel for
    #endif
    for (int i = 0; i < data_N; i++) {
        // Binary search for 't'
        double min = T_SEARCH_MIN;
        double max = T_SEARCH_MAX;
        for (int j = 0; j < T_SEARCH_DEPTH; j++) {
            ex_data->t[i] = (min+max) / 2;
            double e1 = traj_error_i(flash_geo, params, ex_data->t[i]-(max-min)/100, i);
            double e2 = traj_error_i(flash_geo, params, ex_data->t[i]+(max-min)/100, i);
            if (e1 < e2)
                max = ex_data->t[i];
            else
                min = ex_data->t[i];
        }
        // Actually process current answer
        process_flash_traj_i(flash_geo, params, ex_data->t[i], i);
    }
}

/*
 * Proceese answer for the trajectory for one observer given 't'.
 */
void data_t::process_flash_traj_i(const vec3d_t &flash_geo, const vec3d_t params, double t, int i) {
    // Translate from geo position to 3D point
    vec3d_t flash = geo_to_xyz(flash_geo);

    // Begin position
    vec3d_t begin = flash + params*t;

    // Relative begin position
    vec3d_t begin_rel = begin - ob_pos[i];

    // Height and length
    double d = begin * normal[i] - r[i];
    double l = begin_rel.length();

    vec3d_t begin_proj_rel = begin - normal[i]*d - ob_pos[i];

    ex_data->hb[i] = l == 0 ? PI/2 : asin(d/l);
    ex_data->zb[i] = azimuth(ob_pos[i], begin_proj_rel, normal[i], ob_lat[i], ob_lon[i], r[i]);
    ex_data->a[i] = desent_angle(ex_data->hb[i], ex_data->zb[i], ex_data->h0[i], ex_data->z0[i]);
}

/*
 * Returns square-error of trajectory for current answer given 't'
 */
double data_t::traj_error_i(const vec3d_t &flash_geo, const vec3d_t params, double t, int i) {
    process_flash_traj_i(flash_geo, params, t, i);
    // Error
    double error = 0;
    error += pow(angle_delta(ex_data->zb[i], ob_data->zb[i]), 2) * k_zb[i] * ob_e[i];
    error += pow(angle_delta(ex_data->hb[i], ob_data->hb[i]), 2) * k_hb[i] * ob_e[i];
    error += pow(angle_delta(ex_data->a[i], ob_data->a[i]), 2) * k_a[i] * ob_e[i];
    return error;
}

