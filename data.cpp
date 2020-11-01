#include "data.h"

#include "utils.h"
#include "simd.h"
#include <cstdio>

#include "hyperparams.h"



/*
 * Init
 */
data_t::data_t(const char *file) {

    // Open file
    FILE *infile = fopen(file, "r");
    if (!infile) {
        fprintf(stderr, "Failed to open file %s\n", file);
        exit(1);
    }
    fscanf(infile, "%d", &data_N);

    // Allocate memory
    traj_error_start = new double[data_N];
    traj_error_end =   new double[data_N];
    ob_e =             new double[data_N];
    k_z0 =             new double[data_N];
    k_h0 =             new double[data_N];
    k_traj_start =     new double[data_N];
    k_traj_start_tmp = new double[data_N];
    k_traj_end =       new double[data_N];
    k_traj_end_tmp =   new double[data_N];
    k_traj_a =         new double[data_N];
    ob_pos =           new vec3d_t[data_N];
    ob_pos_geo =       new vec3d_t[data_N];
    ob_data =          new data_set_t(data_N);
    ex_data =          new data_set_t(data_N);

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
    for (int i = 0; i < data_N; i++)
        k_z0[i] = ob_data->z0[i] >= 0;
}
void data_t::reset_k_h0() {
    for (int i = 0; i < data_N; i++)
        k_h0[i] = ob_data->h0[i] >= 0;
}
void data_t::reset_k_traj() {
    for (int i = 0; i < data_N; i++) {
        k_traj_start[i] = ob_data->zb[i] >= 0 && ob_data->hb >= 0;
        k_traj_end[i] = ob_data->z0[i] >= 0 && ob_data->h0 >= 0;
        k_traj_a[i] = ob_data->a[i] >= 0;
    }
}


/*
 * Sets K=0 to all data which square-error is
 * max_error_k times greater than mean square-error
 */
int data_t::eliminate(double ne, double total, double *errors, double *k, double max_e, double acc) {
    //printf("%f\n", ne);
    int k_count = 0;
    bool clear;
    do {
        clear = true;
        for (int i = 0; i < data_N; i++) {
            // Skip eliminated
            if (!k[i])
                continue;

            // Mean w/o this observer
            double mean = (total - errors[i]) / (ne - ob_e[i]);

            // Eliminate
            if (errors[i] > mean * max_e && errors[i] > acc*ob_e[i]) {
                clear = false; // Should repeat
                k[i] = 0;
                k_count += ob_e[i];
                total -= errors[i];
                ne -= ob_e[i];
            }
        }
    } while (!clear);
    return k_count;
}
void data_t::eliminate_inconsistent_z0(const vec2d_t &flash_geo) {
    reset_k_z0();
    process_z0(flash_geo);

    // Compute errors
    double total_error = 0, n = 0;
    double *errors = new double[data_N];
    for (int i = 0; i < data_N; i++) {
        double k_i = k_z0[i] * ob_e[i];
        errors[i] = pow(angle_delta(ob_data->z0[i], ex_data->z0[i]), 2) * k_i;
        total_error += errors[i] * k_z0[i];
        n += k_i;
    }

    //int eliminate(double ne, double total, double *errors, double *k, double max_e, double acc);
    eliminate(n, total_error, errors, k_z0, MAX_ERROR_Z0, ACCURACY_Z);
    delete errors;
}
void data_t::eliminate_inconsistent_traj_data(const line3d_t &traj) {
    reset_k_traj();
    process_traj(traj);

    double start = 0, start_n = 0;
    double end = 0, end_n = 0;
    double a = 0, a_n = 0;
    double *a_errors = new double[data_N];
    for (int i = 0; i < data_N; i++) {
        // Update temporal k's
        k_traj_start[i] *= k_traj_start_tmp[i];
        k_traj_end[i] *= k_traj_end_tmp[i];

        // Shortcuts
        double e = ob_e[i];

        // Counts
        start_n += e * k_traj_start[i];
        end_n += e * k_traj_end[i];
        a_n += e * k_traj_a[i];

        // Errors
        start += traj_error_start[i] * k_traj_start[i];
        end += traj_error_end[i] * k_traj_end[i];
        a_errors[i] = pow(ob_data->a[i] - ex_data->a[i], 2) * e * k_traj_a[i];
        a += a_errors[i];
    }
    //int eliminate(double ne, double total, double *errors, double *k, double max_e, double acc);
    eliminate(start_n, start, traj_error_start, k_traj_start, MAX_ERROR_START, 0);
    eliminate(end_n, end, traj_error_end, k_traj_end, MAX_ERROR_END, 0);
    eliminate(a_n, a, a_errors, k_traj_a, MAX_ERROR_A, ACCURACY_A);
    delete a_errors;
}


/*
 * Returns square-error of given answer
 */
double data_t::rate_z0(const vec2d_t &flash_geo) {
    process_z0(flash_geo);

    double error = 0;
    for (int i = 0; i < data_N; i++)
        error += pow(angle_delta(ex_data->z0[i], ob_data->z0[i]), 2)
            * k_z0[i] * ob_e[i];
    return error;
}
double data_t::rate_zb(const vec2d_t &flash_geo) {
    process_zb(flash_geo);

    double error = 0;
    for (int i = 0; i < data_N; i++)
        error += pow(angle_delta(ex_data->zb[i], ob_data->zb[i]), 2)
            * (ob_data->zb[i] >= 0) * ob_e[i];
    return error;
}
double data_t::rate_traj(const line3d_t &traj) {
    process_traj(traj);

    double error = 0, count = 0;
    for (int i = 0; i < data_N; i++) {
        // Shortcuts
        double e = ob_e[i];
        double k_start_i = k_traj_start[i] * k_traj_start_tmp[i];
        double k_end_i = k_traj_end[i] * k_traj_end_tmp[i];
        count += e*(k_traj_a[i] + k_start_i + k_end_i);

        // Sum up errors
        error += pow(ob_data->a[i] - ex_data->a[i], 2) * k_traj_a[i] * e;
        error += traj_error_start[i] * k_start_i;
        error += traj_error_end[i] * k_end_i;
    }
    return error / count;
}


/*
 * Calculate the height of a flash given it's location.
 */
vec2d_t data_t::calc_flash_height(const vec2d_t &flash_2d) {
    reset_k_h0();

    // Compute height
    double *heights = new double[data_N];
    double height = 0;
    double count = 0;
    for (int i = 0; i < data_N; i++) {
        if (k_h0[i]) {
            heights[i] = altitude_to_height(ob_pos_geo[i], flash_2d, ob_data->h0[i]);
            if (heights[i] > 120000. || heights[i] < 0)
                k_h0[i] = 0;
            else {
                height += heights[i] * ob_e[i];
                count += ob_e[i];
            }
        }
    }
    height /= count;

    // Compute errors
    double *errors = new double[data_N];
    double total_error = 0;
    for (int i = 0; i < data_N; i++) {
        errors[i] = pow(heights[i] - height, 2) * ob_e[i] * k_h0[i];
        total_error += errors[i];
    }

    // Eliminate inconsistent
    count -= eliminate(count, total_error, errors, k_h0, MAX_ERROR_H0, ACCURACY_H);
    
    // Recompute height
    height = 0;
    for (int i = 0; i < data_N; i++)
        height += heights[i] * ob_e[i] * k_h0[i];
    height /= count;

    
    // Sigma
    double d_z = 0;
    for (int i = 0; i < data_N; i++)
        d_z += pow(height - heights[i], 2) * ob_e[i] * k_h0[i];
    d_z /= count - 1;

    // Free memory
    delete errors;
    delete heights;

    // Return {height, sigma}
    return {height, sqrt(d_z)};
}




/*
 * Return sigma (standard deviation)
 * of a given answer.
 */
vec2d_t data_t::sigma_z0(const vec2d_t &flash) {
    double d_lat = 0;
    double d_lon = 0;
    double count = 0;
    
    vec3d_t flash_xyz = geo_to_xyz(flash);

    for (int i = 0; i < data_N; i++) {
        if (!k_z0[i])
            continue;

        double z0 = ob_data->z0[i];
        double e = ob_e[i];

        vec3d_t k {sin(z0), cos(z0), 0};
        k = east_vec(ob_pos_geo[i].y)*k.x + north_vec(ob_pos_geo[i].to2d())*k.y;

        vec3d_t n = ob_pos[i].xprod(k).normalized();


        double d = flash_xyz * n;
        vec3d_t pp = xyz_to_geo(flash_xyz - n*d);

        d_lat += pow(flash.x - pp.x, 2) * e;
        d_lon += pow(flash.y - pp.y, 2) * e;
        count += e;
    }

    vec2d_t retval;
    retval.x = sqrt(d_lat / (count-1));
    retval.y = sqrt(d_lon / (count-1));
    return retval;

}
vec2d_t data_t::sigma_zb(const vec2d_t &flash) {
    double d_lat = 0;
    double d_lon = 0;
    double count = 0;
    
    vec3d_t flash_xyz = geo_to_xyz(flash);

    for (int i = 0; i < data_N; i++) {
        if (!k_z0[i])
            continue;

        double z0 = ob_data->zb[i];
        double e = ob_e[i];

        vec3d_t k {sin(z0), cos(z0), 0};
        k = east_vec(ob_pos_geo[i].y)*k.x + north_vec(ob_pos_geo[i].to2d())*k.y;

        vec3d_t n = ob_pos[i].xprod(k).normalized();


        double d = flash_xyz * n;
        vec3d_t pp = xyz_to_geo(flash_xyz - n*d);
        //printf("%f %f\n", pp.x/PI*180, pp.y/PI*180);

        d_lat += pow(flash.x - pp.x, 2) * e;
        d_lon += pow(flash.y - pp.y, 2) * e;
        count += e;
    }

    vec2d_t retval;
    retval.x = sqrt(d_lat / (count-1));
    retval.y = sqrt(d_lon / (count-1));
    return retval;

}
double data_t::sigma_h0(const vec3d_t &flash) {
    vec2d_t flash_2d = flash.to2d();
    double d_z = 0;
    double count = 0;

    for (int i = 0; i < data_N; i++) {
        if (!k_h0[i])
            continue;

        double h0 = ob_data->h0[i];
        double e = ob_e[i];

        double z = altitude_to_height(ob_pos_geo[i], flash_2d, h0);

        d_z += pow(flash.z - z, 2) * e;
        count += e;
    }

    return sqrt(d_z / (count * (count-1)));
}


/*
 * Fills in the 'expected' values.
 */
void data_t::process_z0(const vec2d_t &flash_geo) {
    for (int i = 0; i < data_N; i++)
        ex_data->z0[i] = desent_angle(ob_pos_geo[i].to2d(), flash_geo);
}
void data_t::process_zb(const vec2d_t &flash_geo) {
    for (int i = 0; i < data_N; i++)
        ex_data->zb[i] = desent_angle(ob_pos_geo[i].to2d(), flash_geo);
}
void data_t::process_traj(const line3d_t &traj_line) {
    vec3d_t flash = geo_to_xyz(traj_line.end);
    vec3d_t traj = traj_line.vec();
    for (int i = 0; i < data_N; i++) {
        // Shortcuts
        double &z0 = ob_data->z0[i];
        double &h0 = ob_data->h0[i];
        double &zb = ob_data->zb[i];
        double &hb = ob_data->hb[i];
        vec2d_t observer = ob_pos_geo[i].to2d();
        vec3d_t flash_rel = flash - ob_pos[i];

        // Transalte start/end points to polar
        vec2d_t start_pol = {
            height_to_altitude(ob_pos[i], traj_line.start),
            desent_angle(observer, traj_line.start.to2d())
        };
        vec2d_t end_pol = {
            height_to_altitude(ob_pos[i], traj_line.end),
            desent_angle(observer, traj_line.end.to2d())
        };

        // Compute the global vectors of observations
        vec3d_t k_start = polar_to_xyz({hb, zb});
        vec3d_t k_end = polar_to_xyz({h0, z0});
        k_start = local_to_global(k_start, observer);
        k_end =   local_to_global(k_end, observer);

        // Construct a plane by the observer and trajectory vector
        vec3d_t plane = flash_rel.xprod(traj).normalized();

        // Get errors
        traj_error_start[i] = ob_e[i] * pow(asin(plane*k_start),2);
        traj_error_end[i]   = ob_e[i] * pow(asin(plane*k_end),2);

        // Project observations to that plane and normalize
        k_start = (k_start - plane*(plane * k_start)).normalized();
        k_end = (k_end - plane*(plane * k_end)).normalized();

        // Check if observations are correct
        k_traj_start_tmp[i] = is_observation_correct(flash, traj, ob_pos[i], k_start);
        k_traj_end_tmp[i] =   is_observation_correct(flash, traj, ob_pos[i], k_end);

        // Translate the observation back to local
        k_start = global_to_local(k_start, observer);
        k_end = global_to_local(k_end, observer);

        // Conver them to azimuths and altitudes
        vec2d_t start = xyz_to_polar(k_start);
        vec2d_t end = xyz_to_polar(k_end);
        ex_data->hb[i] = start.x;
        ex_data->zb[i] = start.y;
        ex_data->h0[i] = end.x;
        ex_data->z0[i] = end.y;

        // Do not accept negative altitudes
        if (start.x < 0)
            k_traj_start_tmp[i] = 0;
        if (end.x < 0)
            k_traj_end_tmp[i] = 0;

        // Compute desent angle
        if (k_traj_start[i] * k_traj_start_tmp[i] == 0)
            start = start_pol;
        if (k_traj_end[i] * k_traj_end_tmp[i] == 0)
            end = end_pol;
        ex_data->a[i] = desent_angle(start, end);
    }
}
