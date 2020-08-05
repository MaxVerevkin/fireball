#include <stdio.h>
#include <random>
#include <omp.h>

#include "data.h"


///////////////////////
/// Hyperparameters ///
///////////////////////

#define RANDOM_FLASH_N 20
#define RANDOM_TRAJ_N 100
#define TRIES_N 4000
#define MAX_ERROR 6

#define FLASH_AMPL_MUL (1. / 1.45) // from 1 to inf
#define FLASH_ALPHA_MEAN 1        // from 0 to 1
#define TRAJ_AMPL_MUL (1. / 1.03) // from 1 to inf
#define TRAJ_ALPHA_MEAN 1        // from 0 to 1


///////////////////////////////////////////
/// Default range for random estimation ///
///////////////////////////////////////////

double x0_r_mean = 4000000.;
double x0_r_ampl = 500000.;

double y0_r_mean = 4000000.;
double y0_r_ampl = 500000.;

double z0_r_mean = 50000.;
double z0_r_ampl = 50000.;

double kx_r_mean = 0.;
double kx_r_ampl = 72000.;

double ky_r_mean = 0.;
double ky_r_ampl = 72000.;

double kz_r_mean = 36000.;
double kz_r_ampl = 36000.;


/////////////////
/// Functions ///
/////////////////

/*
 * Update range of the answer, given current best answer.
 *
 * New range will satisfy those constraints:
 *      1) Amplitude of the range will be multiplied by RANDOM_AMPL_MUL.
 *      2) New mean will be in linearly interpolated
 *      between old mean and current answer, i.e.:
 *          new_mean=lerp(mean,answer,RANDOM_ALPHA_MEAN)
 */
void update_flash_range(const vec3d_t &pos) {
    // Calculate new mean
    x0_r_mean = lerp(x0_r_mean, pos.x, FLASH_ALPHA_MEAN);
    y0_r_mean = lerp(y0_r_mean, pos.y, FLASH_ALPHA_MEAN);
    z0_r_mean = lerp(z0_r_mean, pos.z, FLASH_ALPHA_MEAN);

    // Calculate new amplitude
    x0_r_ampl *= FLASH_AMPL_MUL;
    y0_r_ampl *= FLASH_AMPL_MUL;
    z0_r_ampl *= FLASH_AMPL_MUL;
}
void update_traj_range(const vec3d_t &traj) {
    // Calculate new mean
    kx_r_mean = lerp(kx_r_mean, traj.x, TRAJ_ALPHA_MEAN);
    ky_r_mean = lerp(ky_r_mean, traj.y, TRAJ_ALPHA_MEAN);
    kz_r_mean = lerp(kz_r_mean, traj.z, TRAJ_ALPHA_MEAN);

    // Calculate new amplitude
    kx_r_ampl *= TRAJ_AMPL_MUL;
    ky_r_ampl *= TRAJ_AMPL_MUL;
    kz_r_ampl *= TRAJ_AMPL_MUL;
}



/*
 * Run one epoch of guessing.
 */
std::random_device r;
std::uniform_real_distribution<> dist(-1, 1);
void run_random_epoch_on_flash(data_t &data, vec3d_t &best_pos, double *best_error) {
    #pragma omp parallel
    {
        std::mt19937 e2(r()*(omp_get_thread_num()+878));

        processed_answer proc_ans;
        vec3d_t local_best_pos;
        double local_best_error = INFINITY;

        vec3d_t pos;
        for (int i = 0; i < TRIES_N; i++) {
            pos.x = x0_r_mean + x0_r_ampl * dist(e2);
            pos.y = y0_r_mean + y0_r_ampl * dist(e2);
            pos.z = z0_r_mean + z0_r_ampl * dist(e2);
            
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
    #pragma omp parallel
    {
        std::mt19937 e2(r()*(omp_get_thread_num()+878));

        processed_answer proc_ans;
        vec3d_t local_best_traj;
        double local_best_error = INFINITY;

        vec3d_t traj;
        for (int i = 0; i < TRIES_N; i++) {
            traj.x = kx_r_mean + kx_r_ampl * dist(e2);
            traj.y = ky_r_mean + ky_r_ampl * dist(e2);
            traj.z = kz_r_mean + kz_r_ampl * dist(e2);
            
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


/*
 * main
 */
int main() {
    data_t data;
    printf("Data is initialized\n");

    vec3d_t flash_pos;
    vec3d_t flash_traj;
    double flash_error = INFINITY;
    double traj_error = INFINITY;


    ///////////////////////////
    /// Find flash position ///
    ///////////////////////////

    for (int i = 0; i+4 < RANDOM_FLASH_N; i+=5) {
        run_random_epoch_on_flash(data, flash_pos, &flash_error);
        update_flash_range(flash_pos);
        run_random_epoch_on_flash(data, flash_pos, &flash_error);
        update_flash_range(flash_pos);
        run_random_epoch_on_flash(data, flash_pos, &flash_error);
        update_flash_range(flash_pos);
        run_random_epoch_on_flash(data, flash_pos, &flash_error);
        update_flash_range(flash_pos);
        run_random_epoch_on_flash(data, flash_pos, &flash_error);
        update_flash_range(flash_pos);
        data.eliminate_inconsistent_flash_data(flash_pos, MAX_ERROR);
    }

    printf("\nSummary on finding flash position:\n");
    printf("    Total square-error (rad): %#9.6f\n", flash_error);
    printf("    Mean square-error  (rad): %#9.6f\n", flash_error / (data_N * 2 - data.k_count));
    printf("    Standard error     (deg): %#9.6f\n", sqrt(flash_error / (data_N * 2 - data.k_count))/PI*180);


    /////////////////////////////
    /// Find flash trajectory ///
    /////////////////////////////

    for (int i = 0; i+4 < RANDOM_TRAJ_N; i+=5) {
        run_random_epoch_on_traj(data, flash_pos, flash_traj, &traj_error);
        update_traj_range(flash_traj);
        run_random_epoch_on_traj(data, flash_pos, flash_traj, &traj_error);
        update_traj_range(flash_traj);
        run_random_epoch_on_traj(data, flash_pos, flash_traj, &traj_error);
        update_traj_range(flash_traj);
        run_random_epoch_on_traj(data, flash_pos, flash_traj, &traj_error);
        update_traj_range(flash_traj);
        run_random_epoch_on_traj(data, flash_pos, flash_traj, &traj_error);
        update_traj_range(flash_traj);
        data.eliminate_inconsistent_traj_data(flash_pos, flash_traj, MAX_ERROR);
    }

    printf("\nSummary on finding flash trajectory:\n");
    printf("    Total square-error (rad): %#9.6f\n", traj_error);
    printf("    Mean square-error  (rad): %#9.6f\n", traj_error / (data_N * 3 - data.k_count));
    printf("    Standard error     (deg): %#9.6f\n", sqrt(traj_error / (data_N * 3 - data.k_count))/PI*180);



    // Print answer
    printf("\nAnswer: %#9.0f %#9.0f %#9.0f %#9.0f %#9.0f %#9.0f\n",
            flash_pos.x, flash_pos.y, flash_pos.z, flash_traj.x, flash_traj.y, flash_traj.z);

    // Print processed answer for each observer
    printf("\n");
    //data.process_flash_traj(flash_pos, flash_traj, data.ex_data);
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
