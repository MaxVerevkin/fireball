#include <stdio.h>
#include <random>
#include <omp.h>

#include "data.h"


///////////////////////
/// Hyperparameters ///
///////////////////////

#define RANDOM_N 20
#define TRIES_N 10000
#define MAX_ERROR 11

#define RANDOM_AMPL_MUL (1. / 1.45) // from 1 to inf
#define RANDOM_ALPHA_MEAN 1        // from 0 to 1


///////////////////////////////////////////
/// Default range for random estimation ///
///////////////////////////////////////////

double x0_r_mean = 4000000.;
double x0_r_ampl = 1000000.;

double y0_r_mean = 4000000.;
double y0_r_ampl = 1000000.;

double z0_r_mean = 50000.;
double z0_r_ampl = 50000.;


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
void update_range(const vec3d_t &answer) {
    // Calculate new mean
    x0_r_mean = lerp(x0_r_mean, answer.x, RANDOM_ALPHA_MEAN);
    y0_r_mean = lerp(y0_r_mean, answer.y, RANDOM_ALPHA_MEAN);
    z0_r_mean = lerp(z0_r_mean, answer.z, RANDOM_ALPHA_MEAN);

    // Calculate new amplitude
    x0_r_ampl *= RANDOM_AMPL_MUL;
    y0_r_ampl *= RANDOM_AMPL_MUL;
    z0_r_ampl *= RANDOM_AMPL_MUL;
}



/*
 * Run one epoch of guessing.
 */
std::random_device r;
std::uniform_real_distribution<> dist(-1, 1);
void run_random_epoch(data_t &data, vec3d_t &best_answer, double *best_error) {
    std::mt19937 _e2(r());
    #pragma omp parallel
    {
        std::mt19937 e2(_e2()*(omp_get_thread_num()+878));

        processed_answer proc_ans;
        vec3d_t local_best_answer;
        double local_best_error = INFINITY;
        vec3d_t answer;
        double error;


        //#pragma omp for
        for (int i = 0; i < TRIES_N; i++) {
            answer.x = x0_r_mean + x0_r_ampl * dist(e2);
            answer.y = y0_r_mean + y0_r_ampl * dist(e2);
            answer.z = z0_r_mean + z0_r_ampl * dist(e2);
            
            error = data.rate_flash_pos(answer, proc_ans);

            if (error < local_best_error) {
                local_best_error = error;
                local_best_answer = answer;
            }
        }


        #pragma omp critical
        if (local_best_error < *best_error) {
            *best_error = local_best_error;
            best_answer = local_best_answer;
        }
    }
}


/*
 * main
 */
int main() {
    data_t data;
    printf("Data is initialized\n");

    vec3d_t answer;
    double error = INFINITY;

    for (int i = 0; i+4 < RANDOM_N; i+=5) {
        run_random_epoch(data, answer, &error);
        update_range(answer);
        run_random_epoch(data, answer, &error);
        update_range(answer);
        run_random_epoch(data, answer, &error);
        update_range(answer);
        run_random_epoch(data, answer, &error);
        update_range(answer);
        run_random_epoch(data, answer, &error);
        update_range(answer);
        data.eliminate_inconsistent(answer, MAX_ERROR);
    }

    error = data.rate_flash_pos(answer);
    printf("\nSummary:\n");
    printf("    Total square-error (rad): %#9.6f\n", error);
    printf("    Mean square-error  (rad): %#9.6f\n", error / (data_N * 5 - data.k_count));
    printf("    Standard error     (deg): %#9.6f\n", sqrt(error / (data_N * 5 - data.k_count))/PI*180);


    // Print answer
    printf("\nAnswer: %#9.0f %#9.0f %#9.0f\n",
            answer.x, answer.y, answer.z);

    // Print processed answer for each observer (include ignored)
    printf("\n");
    data.process_answer(answer);
    for (int i = 0; i < data_N; i++) {
        printf("Processed Answer (%i): %#7.3f | %7.3f\n", i+1,
                data.ex_data.z0[i]/PI*180, data.ex_data.h0[i]/PI*180);
    }

    // Print ignored data
    printf("\n");
    for (int i = 0; i < data_N; i++) {
        if (!data.k_z0[i])
            printf("Ignore: 'azimuth end'    for observer %i\n", i+1);
        if (!data.k_h0[i])
            printf("Ignore: 'altitude end'   for observer %i\n", i+1);
    }

    // Exit
    printf("\n");
    return 0;
}
