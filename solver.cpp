#include <stdio.h>
#include <random>
#include <omp.h>
#include "data.h"

#include "hyperparams.h"


/////////////////
/// Functions ///
/////////////////

/*
 * Run one epoch of guessing.
 */
std::random_device r;
std::uniform_real_distribution<> dist(-3, 3);
void run_random_epoch_on_flash(data_t &data, vec3d_t &best_pos, double *best_error) {
    vec3d_t sigma = data.sigma_flash_pos(best_pos);
    #pragma omp parallel
    {
        std::mt19937 e2(r()*(omp_get_thread_num()+878));

        processed_answer proc_ans;
        vec3d_t local_best_pos;
        double local_best_error = INFINITY;

        vec3d_t pos;
        for (int i = 0; i < TRIES_N; i++) {
            pos.x = best_pos.x + sigma.x * dist(e2);
            pos.y = best_pos.y + sigma.y * dist(e2);
            pos.z = best_pos.z + sigma.z * dist(e2);
            
            double error = data.rate_flash_pos(pos, proc_ans);
            if (error < local_best_error) {
                local_best_error = error;
                local_best_pos = pos;
            }
        }


        #pragma omp critical
        if (local_best_error < *best_error) {
            *best_error = local_best_error;
            best_pos = local_best_pos;
        }
    }
}
void run_random_epoch_on_traj(data_t &data, const vec3d_t &flash, vec3d_t &best_traj, double *best_error) {
    vec3d_t sigma = data.sigma_flash_traj(flash, best_traj, data.ex_data);
    #pragma omp parallel
    {
        std::mt19937 e2(r()*(omp_get_thread_num()+878));

        processed_answer proc_ans;
        vec3d_t local_best_traj;
        double local_best_error = INFINITY;

        vec3d_t traj;
        for (int i = 0; i < TRIES_N; i++) {
            traj.x = best_traj.x + sigma.x * dist(e2);
            traj.y = best_traj.y + sigma.y * dist(e2);
            traj.z = best_traj.z + sigma.z * dist(e2);
            
            double error = data.rate_flash_traj(flash, traj, proc_ans);
            if (error < local_best_error) {
                local_best_error = error;
                local_best_traj = traj;
            }
        }


        #pragma omp critical
        if (local_best_error < *best_error) {
            *best_error = local_best_error;
            best_traj = local_best_traj;
        }
    }
}

double min_8d(double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8) {
    if (x1 < x2 && x1 < x3 && x1 < x4 && x1 < x5 && x1 < x6 && x1 < x7 && x1 < x8)
        return x1;
    if (x2 < x3 && x2 < x4 && x2 < x5 && x2 < x6 && x2 < x7 && x2 < x8)
        return x2;
    if (x3 < x4 && x3 < x5 && x3 < x6 && x3 < x7 && x3 < x8)
        return x3;
    if (x4 < x5 && x4 < x6 && x4 < x7 && x4 < x8)
        return x4;
    if (x5 < x6 && x5 < x7 && x5 < x8)
        return x5;
    if (x6 < x7 && x6 < x8)
        return x6;
    if (x7 < x8)
        return x7;
    return x8;
}


/*
 * main
 */
int main() {
    data_t data;
    printf("Data is initialized\n");

    
    ///////////////////////////
    /// Find flash position ///
    ///////////////////////////

    vec3d_t flash_pos = {x0_default, y0_default, z0_default};
    double flash_error = INFINITY;
    for (int i = 0; i+4 < RANDOM_FLASH_N; i+=5) {
        run_random_epoch_on_flash(data, flash_pos, &flash_error);
        run_random_epoch_on_flash(data, flash_pos, &flash_error);
        run_random_epoch_on_flash(data, flash_pos, &flash_error);
        run_random_epoch_on_flash(data, flash_pos, &flash_error);
        run_random_epoch_on_flash(data, flash_pos, &flash_error);
        data.eliminate_inconsistent_flash_data(flash_pos);
    }

    printf("\nSummary on finding flash position:\n");
    printf("    Total square-error (rad): %#9.6f\n", flash_error);
    printf("    Mean square-error  (rad): %#9.6f\n", flash_error / (data.data_Ne * 2 - data.k_count));
    printf("    Standard error     (deg): %#9.6f\n", sqrt(flash_error / (data.data_Ne * 2 - data.k_count))/PI*180);


    /////////////////////////////
    /// Find flash trajectory ///
    /////////////////////////////

    vec3d_t flash_traj;
    double traj_error = INFINITY;

    for (int j = 0; j < 2; j++) {
        vec3d_t min_val {kx_min, ky_min, kz_min};
        vec3d_t max_val {kx_max, ky_max, kz_max};
        flash_traj.x = (min_val.x + max_val.x) / 2;
        flash_traj.y = (min_val.y + max_val.y) / 2;
        flash_traj.z = (min_val.z + max_val.z) / 2;
        for (int i = 0; i < 20; i++) {
            double e1 = data.rate_flash_traj(flash_pos, {flash_traj.x+10,flash_traj.y+10,flash_traj.z+10}, data.ex_data);
            double e2 = data.rate_flash_traj(flash_pos, {flash_traj.x+10,flash_traj.y+10,flash_traj.z-10}, data.ex_data);
            double e3 = data.rate_flash_traj(flash_pos, {flash_traj.x+10,flash_traj.y-10,flash_traj.z+10}, data.ex_data);
            double e4 = data.rate_flash_traj(flash_pos, {flash_traj.x+10,flash_traj.y-10,flash_traj.z-10}, data.ex_data);
            double e5 = data.rate_flash_traj(flash_pos, {flash_traj.x-10,flash_traj.y+10,flash_traj.z+10}, data.ex_data);
            double e6 = data.rate_flash_traj(flash_pos, {flash_traj.x-10,flash_traj.y+10,flash_traj.z-10}, data.ex_data);
            double e7 = data.rate_flash_traj(flash_pos, {flash_traj.x-10,flash_traj.y-10,flash_traj.z+10}, data.ex_data);
            double e8 = data.rate_flash_traj(flash_pos, {flash_traj.x-10,flash_traj.y-10,flash_traj.z-10}, data.ex_data);
            double min = min_8d(e1, e2, e3, e4, e5, e6, e7, e8);

            if (min == e1 || min == e2 || min == e3 || min == e4) {
                min_val.x = flash_traj.x;
                flash_traj.x = (flash_traj.x + max_val.x) / 2;
            } else {
                max_val.x = flash_traj.x;
                flash_traj.x = (flash_traj.x + min_val.x) / 2;
            }

            if (min == e1 || min == e2 || min == e5 || min == e6) {
                min_val.y = flash_traj.y;
                flash_traj.y = (flash_traj.y + max_val.y) / 2;
            } else {
                max_val.y = flash_traj.y;
                flash_traj.y = (flash_traj.y + min_val.y) / 2;
            }

            if (min == e1 || min == e3 || min == e5 || min == e7) {
                min_val.z = flash_traj.z;
                flash_traj.z = (flash_traj.z + max_val.z) / 2;
            } else {
                max_val.z = flash_traj.z;
                flash_traj.z = (flash_traj.z + min_val.z) / 2;
            }
        }
        data.eliminate_inconsistent_traj_data(flash_pos, flash_traj);
        traj_error = data.rate_flash_traj(flash_pos, flash_traj, data.ex_data);
    }


    printf("\nSummary on finding flash trajectory:\n");
    printf("    Total square-error (rad): %#9.6f\n", traj_error);
    printf("    Mean square-error  (rad): %#9.6f\n", traj_error / (data.data_Ne * 3 - data.k_count));
    printf("    Standard error     (deg): %#9.6f\n", sqrt(traj_error / (data.data_Ne * 3 - data.k_count))/PI*180);



    // Print answer
    vec3d_t flash_pos_sigma = data.sigma_flash_pos(flash_pos);
    vec3d_t flash_traj_sigma = data.sigma_flash_traj(flash_pos, flash_traj, data.ex_data);
    printf("\nAnswer: %9.0f±%4.0f %9.0f±%4.0f %9.0f±%4.0f %9.0f±%4.0f %9.0f±%4.0f %9.0f±%4.0f\n",
            flash_pos.x, flash_pos_sigma.x,
            flash_pos.y, flash_pos_sigma.y,
            flash_pos.z, flash_pos_sigma.z,
            flash_traj.x, flash_traj_sigma.x,
            flash_traj.y, flash_traj_sigma.y,
            flash_traj.z, flash_traj_sigma.z);

    // Print processed answer for each observer
    printf("\n                  i     z0        h0        zb        hb         a         T   \n");
    for (int i = 0; i < data_N; i++) {
        printf("Processed Answer (%i): %#7.3f | %7.3f | %7.3f | %7.3f | %7.3f | %7.3f\n", i+1,
                data.ex_data.z0[i]/PI*180, data.ex_data.h0[i]/PI*180,
                data.ex_data.zb[i]/PI*180, data.ex_data.hb[i]/PI*180,
                data.ex_data.a[i]/PI*180, data.ex_data.t[i]);
    }

    // Print ignored data
    printf("\n");
    for (int i = 0; i < data_N; i++) {
        if (!data.k_z0[i])
            printf("Ignore: 'azimuth end'    for observer %i\n", i+1);
        if (!data.k_h0[i])
            printf("Ignore: 'altitude end'   for observer %i\n", i+1);
        if (!data.k_zb[i])
            printf("Ignore: 'azimuth begin'  for observer %i\n", i+1);
        if (!data.k_hb[i])
            printf("Ignore: 'altitude begin' for observer %i\n", i+1);
        if (!data.k_a[i])
            printf("Ignore: 'desent_angle'   for observer %i\n", i+1);
    }

    // Exit
    printf("\n");
    return 0;
}
