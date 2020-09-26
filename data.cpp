#include "data.h"

#include "utils.h"
#include "simd.h"
#include <cstdio>

#include "hyperparams.h"


// Paralelisation
#ifdef OP_PARALEL
#include <omp.h>
#endif


/*
 * Helper functions.
 */
typedef double (*error_func)(data_t*, int);
double z0_error(data_t *data, int i) {
    double err = abs(angle_delta(data->ob_data->z0[i], data->ex_data->z0[i]));
    err *= cos(data->ob_data->h0[i]);
    return err*err * data->k_z0[i] * data->ob_e[i];
}
double *calc_errors(data_t *data, double *k, error_func func, double *total, double *n) {
    *total = 0;
    *n = 0;
    double *errors = new double[data->data_N];
    for (int i = 0; i < data->data_N; i++) {
        if (k[i]) {
            errors[i] = func(data, i);
            *total += errors[i];
            *n += data->ob_e[i];
        }
    }
    return errors;
}

int eliminate(data_t *data, double ne, double total, double *errors, double *k, double max_e, double acc) {
    int k_count = 0;
    bool clear;
    do {
        clear = true;
        for (int i = 0; i < data->data_N; i++) {
            // Skip eliminated
            if (!k[i])
                continue;

            // Mean w/o this observer
            double mean = (total - errors[i]) / (ne - data->ob_e[i]);

            // Eliminate
            if (errors[i] > mean * max_e && errors[i] > acc*data->ob_e[i]) {
                clear = false; // Should repeat
                k[i] = 0;
                k_count += data->ob_e[i];
                total -= errors[i];
                ne -= data->ob_e[i];
            }
        }
    } while (!clear);
    return k_count;
}

data_t::data_t(const char *file) {

    // Open file
    FILE *infile = fopen(file, "r");
    if (!infile) {
        fprintf(stderr, "Failed to open file %s\n", file);
        exit(1);
    }
    fscanf(infile, "%d", &data_N);

    // Allocate memory
    ob_e =       new double[data_N];
    k_z0 =       new double[data_N];
    k_h0 =       new double[data_N];
    k_zb =       new double[data_N];
    k_hb =       new double[data_N];
    k_a =        new double[data_N];
    r =          new double[data_N];
    ob_pos =     new vec3d_t[data_N];
    ob_pos_geo = new vec3d_t[data_N];
    normal =     new vec3d_t[data_N];
    ob_data =    new data_set_t(data_N);
    ex_data =    new data_set_t(data_N);

    // Used for finding mean longitude
    double sin_lon = 0;
    double cos_lon = 0;
    for (int i = 0; i < data_N; i++) {
        // Read from file
        fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                &(ob_pos_geo+i)->x, &(ob_pos_geo+i)->y, &(ob_pos_geo+i)->z,
                ob_data->a+i, ob_data->zb+i, ob_data->hb+i,
                ob_data->z0+i, ob_data->h0+i,
                ob_data->t+i, ob_e+i);

        // data_Ne is a sum of all ob_e
        data_Ne += ob_e[i];

        // Height relative to earth center
        r[i] = ob_pos_geo[i].z + EARTH_R;

        // Translate
        ob_pos_geo[i].x *= PI / 180;
        ob_pos_geo[i].y *= PI / 180;
        ob_data->z0[i] *= PI / 180;
        ob_data->h0[i] *= PI / 180;
        ob_data->zb[i] *= PI / 180;
        ob_data->hb[i] *= PI / 180;
        ob_data->a[i] *= PI / 180;
        
        // Find mean
        sin_lon += sin(ob_pos_geo[i].y);
        cos_lon += cos(ob_pos_geo[i].y);
        mean_lat += ob_pos_geo[i].x;
        
        // Calculate global 3D position
        ob_pos[i] = geo_to_xyz(ob_pos_geo[i]);

        // Calculate normal
        normal[i] = normal_vec(ob_pos_geo[i].to2d());
    }

    // Calculate mean lon and lat
    mean_lon = atan2(sin_lon, cos_lon);
    mean_lat /= data_N;

    // Set k=1
    reset_k_z0();
    reset_k_h0();
    reset_k_traj();

    // Close file
    fclose(infile);
}


/*
 * Reset k-values to '1'
 */
void data_t::reset_k_z0() {
    k_count_z0 = 0;
    for (int i = 0; i < data_N; i++) {
        k_z0[i] = ob_data->z0[i] >= 0;
        k_count_z0 += !k_z0[i] * ob_e[i];
    }
}
void data_t::reset_k_h0() {
    k_count_h0 = 0;
    for (int i = 0; i < data_N; i++) {
        k_h0[i] = ob_data->h0[i] >= 0;
        k_count_h0 += !k_h0[i] * ob_e[i];
    }
}
void data_t::reset_k_traj() {
    k_count_traj = 0;
    for (int i = 0; i < data_N; i++) {
        k_zb[i] = ob_data->zb[i] >= 0;
        k_hb[i] = ob_data->hb[i] >= 0;
        k_a[i]  = ob_data->a[i]  >= 0;

        k_count_traj += !k_zb[i] * ob_e[i];
        k_count_traj += !k_hb[i] * ob_e[i];
        k_count_traj += !k_a[i] * ob_e[i];
    }
}


/*
 * Calculate the height of a flash given it's location.
 */
double data_t::calc_flash_height(const vec2d_t &flash_geo) {
    reset_k_h0();

    // Compute heights
    double *heights = new double[data_N];
    double total_height = 0;
    double count = 0;
    for (int i = 0; i < data_N; i++) {
        if (k_h0[i]) {
            heights[i] = altitude_to_height(ob_pos_geo[i], flash_geo, ob_data->h0[i]) * ob_e[i];
            total_height += heights[i];
            count += ob_e[i];
        }
    }
    // Compute errors
    double *errors = new double[data_N];
    double total_error = 0;
    double mean_height = total_height / count;
    for (int i = 0; i < data_N; i++) {
        if (k_h0[i]) {
            errors[i] = pow(heights[i] - mean_height, 2) * ob_e[i];
            total_error += errors[i];
        }
    }

    //int eliminate(data_t *data, double ne, double total, double *errors, double *k, double max_e, double acc);
    int x = eliminate(this, count, total_error, errors, k_h0, MAX_ERROR_H0, 0);
    delete errors;
    k_count_h0 += x;
    count -= x;
    
    double height = 0;
    for (int i = 0; i < data_N; i++)
        if (k_h0[i])
            height += heights[i] * ob_e[i];
    height /= count;
    delete heights;

    
    process_h0({flash_geo.x, flash_geo.y, height});
    return height;
}



/*
 * Sets K=0 to all data which square-error is
 * max_error_k times greater than mean square-error
 */
void data_t::eliminate_inconsistent_z0(const vec2d_t &flash_geo) {
    reset_k_z0();
    process_z0(flash_geo);

    // Compute errors
    double total_error = 0, n = 0;
    double *errors = calc_errors(this, k_z0, z0_error, &total_error, &n);

    //int eliminate(data_t *data, double ne, double total, double *errors, double *k, double max_e, double acc);
    k_count_z0 += eliminate(this, n, total_error, errors, k_z0, MAX_ERROR_Z, ACCURACY_Z);
    delete errors;
}
void data_t::eliminate_inconsistent_traj_data(const vec3d_t &flash_pos, const vec3d_t params) {
    reset_k_traj();

    double n = (data_Ne * 3) - k_count_traj;
    double total_error = rate_flash_traj(flash_pos, params);

    for (int i = 0; i < data_N; i++) {
        if (k_hb[i]) {
            double error = pow(angle_delta(ex_data->hb[i], ob_data->hb[i]), 2);
            double mean_error = (total_error - error*ob_e[i]) / (n - ob_e[i]);
            k_hb[i] = error < mean_error * MAX_ERROR_H;
            k_count_traj += !k_hb[i] * ob_e[i];
        }
        if (k_zb[i]) {
            double error = pow(angle_delta(ex_data->zb[i], ob_data->zb[i]), 2);
            double mean_error = (total_error - error*ob_e[i]) / (n - ob_e[i]);
            k_zb[i] = error < mean_error * MAX_ERROR_Z;
            k_count_traj += !k_zb[i] * ob_e[i];
        }
        if (k_a[i]) {
            double error = pow(angle_delta(ex_data->a[i], ob_data->a[i]), 2);
            double mean_error = (total_error - error*ob_e[i]) / (n - ob_e[i]);
            k_a[i] = error < mean_error * MAX_ERROR_A;
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
    //return vel * (-mean_k);

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
            //if (pow(mean_k-ki ,2) < mean_error*MAX_ERROR_T) {
                k += ki;
                n++;
            //}
        }
    }
    k /= n;

    return vel * (-k);
}


/*
 * Normalize observer's 't'
 */
void data_t::normalize_t(double velocity) {
    for (int i = 0; i < data_N; i++) {
        ex_data->t[i] /= velocity;
    }
}


/*
 * Returns square-error of given answer
 */
double data_t::rate_z0(const vec2d_t &flash_geo) {
    process_z0(flash_geo);

    double error = 0;
    for (int i = 0; i < data_N; i++)
        error += z0_error(this, i);
    return error;
}
double data_t::rate_flash_traj(const vec3d_t &flash_pos, const vec3d_t &params) {
    process_flash_traj(flash_pos, params);

    __m128d hb_error = _mm_setzero_pd();
    __m128d zb_error = _mm_setzero_pd();
    __m128d a_error = _mm_setzero_pd();
    for (int i = 0; i+1 < data_N; i+=2) {
        __m128d ei = _mm_load_pd(ob_e + i);

        // h0
        hb_error += angle_delta_sq_pd(ob_data->hb + i, ex_data->hb + i)
            * _mm_load_pd(k_hb + i)
            * ei;

        // z0
        zb_error += angle_delta_sq_pd(ob_data->zb + i, ex_data->zb + i)
            * _mm_load_pd(k_zb + i)
            * ei;

        // z0
        a_error += angle_delta_sq_pd(ob_data->a + i, ex_data->a + i)
            * _mm_load_pd(k_a + i)
            * ei;
    }

    double x[2];
    _mm_storeu_pd(x, hb_error + zb_error + a_error);
    return x[0] + x[1];
}


/*
 * Return sigma (standard deviation)
 * of a given answer.
 */
vec3d_t data_t::sigma_flash_pos(const vec3d_t &flash) {
    double sigma_lat = 0;
    double sigma_lon = 0;
    double sigma_z = 0;
    double lat_n = 0;
    double lon_n = 0;
    double z_n = 0;

    for (int i = 0; i < data_N; i++) {
        double z0 = ob_data->z0[i];
        double d_lat = angle_delta(ob_pos_geo[i].x, flash.x);
        double d_lon = angle_delta(ob_pos_geo[i].y, flash.y);

        bool valid_lat = (d_lon > 0 && z0 < PI) || (d_lon < 0 && z0 > PI);
        bool valid_lon = (d_lat > 0 && z0 < PI/2 && z0 > 1.5*PI) || (d_lat < 0 && z0 > PI/2 && z0 < 1.5*PI);

        // Latitude
        if (k_z0[i] && valid_lat) {
            int sign = (z0>PI)*(-2)+1;
            double cos_alpha = sin(z0) * sin(d_lon) * sin(ob_pos_geo[i].x) - cos(z0) * cos(d_lon);
            double sin_alpha = sqrt(1 - cos_alpha*cos_alpha);
            double sin_lat = (cos(z0) + cos(d_lon)*cos_alpha) / (sin(d_lon)*sin_alpha);
            double lat = sign * asin(sin_lat);
            if (isnan(lat))
                lat = sin_lat * PI * .5;
            lat_n += ob_e[i];
            sigma_lat += pow(flash.x - lat, 2) *  ob_e[i];
        }

        // Longitude
        if (k_z0[i] && valid_lon) {
            double gamma = abs(asin(cos(ob_pos_geo[i].x) * sin(z0) / cos(flash.x)));
            if (flash.x > ob_pos_geo[i].x)
                gamma = PI - gamma;
            double lon = 2*atan(1. / (tan((ob_pos_geo[i].x - flash.x)/2) * tan((z0 - gamma)/2))) + ob_pos_geo[i].y;
            lon_n += ob_e[i];
            sigma_lon += pow(flash.y - lon, 2) *  ob_e[i];
        }

        // Height
        if (k_h0[i]) {
            double z = altitude_to_height(ob_pos_geo[i], flash.to2d(), ob_data->h0[i]);
            sigma_z += pow(flash.z - z, 2) *  ob_e[i];
            z_n += ob_e[i];
        }
    }


    sigma_lat /= lat_n * (lat_n-1);
    sigma_lon /= lon_n * (lon_n-1);
    sigma_z /= z_n * (z_n-1);


    return vec3d_t {sqrt(sigma_lat), sqrt(sigma_lon), sqrt(sigma_z)};
}


/*
 * Fills in the 'expected' values.
 */
void data_t::process_z0(const vec2d_t &flash_geo) {
    #ifdef OP_PARALEL
    #pragma omp parallel for
    #endif
    for (int i = 0; i < data_N; i++)
        ex_data->z0[i] = desent_angle(ob_pos_geo[i].to2d(), flash_geo);
}
void data_t::process_h0(const vec3d_t &flash_geo) {
    #ifdef OP_PARALEL
    #pragma omp parallel for
    #endif
    for (int i = 0; i < data_N; i++)
        ex_data->h0[i] = height_to_altitude(ob_pos_geo[i], flash_geo);
}
void data_t::process_flash_traj(const vec3d_t &flash_pos, const vec3d_t &params) {
    #ifdef OP_PARALEL
    #pragma omp parallel for
    #endif
    for (int i = 0; i < data_N; i++) {
        // Binary search for 't'
        double min = T_SEARCH_MIN;
        double max = T_SEARCH_MAX;
        for (int j = 0; j < T_SEARCH_DEPTH; j++) {
            ex_data->t[i] = (min+max) / 2;
            double e1 = traj_error_i(flash_pos, params, ex_data->t[i]-(max-min)/100, i);
            double e2 = traj_error_i(flash_pos, params, ex_data->t[i]+(max-min)/100, i);
            if (e1 < e2)
                max = ex_data->t[i];
            else
                min = ex_data->t[i];
        }
        // Actually process current answer
        process_flash_traj_i(flash_pos, params, ex_data->t[i], i);
    }
}

/*
 * Proceese answer for the trajectory for one observer given 't'.
 */
void data_t::process_flash_traj_i(const vec3d_t &flash_pos, const vec3d_t params, double t, int i) {
    // Begin position
    vec3d_t begin = flash_pos + params*t;
    vec3d_t begin_geo = xyz_to_geo(begin);

    ex_data->hb[i] = height_to_altitude(ob_pos_geo[i], begin_geo);
    ex_data->zb[i] = desent_angle(ob_pos_geo[i].to2d(), begin_geo.to2d());
    ex_data->a[i] = desent_angle({ex_data->hb[i], ex_data->zb[i]}, {ex_data->h0[i], ex_data->z0[i]});
}

/*
 * Returns square-error of trajectory for current answer given 't'
 */
double data_t::traj_error_i(const vec3d_t &flash_pos, const vec3d_t params, double t, int i) {
    process_flash_traj_i(flash_pos, params, t, i);
    // Error
    double error = 0;
    error += pow(angle_delta(ex_data->zb[i], ob_data->zb[i]), 2) * k_zb[i];
    error += pow(angle_delta(ex_data->hb[i], ob_data->hb[i]), 2) * k_hb[i];
    error += pow(angle_delta(ex_data->a[i], ob_data->a[i]), 2) * k_a[i];
    return error * ob_e[i];
}

