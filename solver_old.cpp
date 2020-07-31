#include <stdio.h>
#include <random>
#include <time.h>
#include <omp.h>

#include "utils.h"
#include "data.h"


///////////////////////
/// Hyperparameters ///
///////////////////////

#define random_epoch_N 20
#define deriv_epoch_N 2000
#define tries_N 10000
#define MAX_ERROR 3.5
#define DERIV_W 10.

#define alpha_ampl (1. / 1.4) // from 1 to inf
#define alpha_mean 1          // from 0 to 1


///////////////////////////////////////////
/// Default range for random estimation ///
///////////////////////////////////////////

double x0_r_mean = 4000000.;
double x0_r_ampl = 2000000.;

double y0_r_mean = 4000000.;
double y0_r_ampl = 2000000.;

double z0_r_mean = 0.;
double z0_r_ampl = 2000000.;

double kx_r_mean = 0.;
double kx_r_ampl = 200000.;

double ky_r_mean = 0.;
double ky_r_ampl = 200000.;

double kz_r_mean = 100000.;
double kz_r_ampl = 100000.;


/////////////////
/// Functions ///
/////////////////

/*
 * Compute derivative vector given
 * given current answer:
 *      d error / d answer
 *
 * Note:
 *      This function returns derivative
 *      multiplied by DERIV_W, because the
 *      derivatives are very small anyway.
 */
answer_t compute_deriv(data_t &data, answer_t answer) {

    // Current error
    double error = data.rate_answer(answer);
    
    // Derivative vector
    answer_t deriv;

    answer.x0 += DERIV_W;
    deriv.x0 = (data.rate_answer(answer) - error);
    answer.x0 -= DERIV_W;

    answer.y0 += DERIV_W;
    deriv.y0 = (data.rate_answer(answer) - error);
    answer.y0 -= DERIV_W;

    answer.z0 += DERIV_W;
    deriv.z0 = (data.rate_answer(answer) - error);
    answer.z0 -= DERIV_W;

    answer.kx += DERIV_W;
    deriv.kx = (data.rate_answer(answer) - error);
    answer.kx -= DERIV_W;

    answer.ky += DERIV_W;
    deriv.ky = (data.rate_answer(answer) - error);
    answer.ky -= DERIV_W;

    answer.kz += DERIV_W;
    deriv.kz = (data.rate_answer(answer) - error);
    answer.kz -= DERIV_W;

    return deriv;
}


/*
 * Update range of the answer, given current best answer.
 *
 * New range will satisfy those constraints:
 *      1) Amplitude of the range will be multiplied by alpha.
 *      2) New mean will be in linearly interpolated
 *      between old mean and current answer, i.e.:
 *          new_mean=lerp(mean,answer,alpha_mean)
 */
void update_range(const answer_t &answer) {
    // Calculate new mean
    x0_r_mean = lerp(x0_r_mean, answer.x0, alpha_mean);
    y0_r_mean = lerp(y0_r_mean, answer.y0, alpha_mean);
    z0_r_mean = lerp(z0_r_mean, answer.z0, alpha_mean);
    kx_r_mean = lerp(kx_r_mean, answer.kx, alpha_mean);
    ky_r_mean = lerp(ky_r_mean, answer.ky, alpha_mean);
    kz_r_mean = lerp(kz_r_mean, answer.kz, alpha_mean);

    // Calculate new amplitude
    x0_r_ampl *= alpha_ampl;
    y0_r_ampl *= alpha_ampl;
    z0_r_ampl *= alpha_ampl;
    kx_r_ampl *= alpha_ampl;
    ky_r_ampl *= alpha_ampl;
    kz_r_ampl *= alpha_ampl;
}



/*
 * Run one epoch of guessing.
 */
std::random_device r;
std::uniform_real_distribution<> dist(-1, 1);
void run_random_epoch(data_t &data, answer_t &best_answer, double *best_error) {
    std::mt19937 _e2(r());
    //#pragma omp parallel
    {
        long int ddd = _e2()*(omp_get_thread_num()+878);
        //printf("ran: %ld\n", ddd);
        std::mt19937 e2(ddd);

        answer_t local_best_answer;
        double local_best_error = INFINITY;
        answer_t answer;
        double error;

        for (int i = 0; i < tries_N; i++) {
            answer.x0 = x0_r_mean + x0_r_ampl * dist(e2);
            answer.y0 = y0_r_mean + y0_r_ampl * dist(e2);
            answer.z0 = z0_r_mean + z0_r_ampl * dist(e2);
            answer.kx = kx_r_mean + kx_r_ampl * dist(e2);
            answer.ky = ky_r_mean + ky_r_ampl * dist(e2);
            answer.kz = kz_r_mean + kz_r_ampl * dist(e2);
    
            error = data.rate_answer(answer);
            //if (error < local_best_error) {
                //local_best_error = error;
                //local_best_answer = answer;
            //}
            if (error < *best_error) {
                *best_error = error;
                best_answer = answer;
            }
        }


        //#pragma omp critical
        //if (local_best_error < *best_error) {
            //*best_error = local_best_error;
            //best_answer = local_best_answer;
        //}
    }
}



/*
 * Run one epoch of estimation.
 */
void run_deriv_epoch(data_t &data, answer_t &answer) {
    //for (int i = 0; i < 10; i++) {
    answer_t answer_deriv = compute_deriv(data, answer);
    double x = 1000000000;
    double y = 100000000;
    //printf("%f, %f, %f, %f, %f\ni", answer_deriv.x0*x, answer_deriv.y0*x, answer_deriv.kx*y, answer_deriv.ky*y, answer_deriv.kz * 6 * 100000);
    answer.x0 -= answer_deriv.x0 * x;
    answer.y0 -= answer_deriv.y0 * x;
    answer.z0 -= answer_deriv.z0 * x;
    answer.kx -= answer_deriv.kx * y;
    answer.ky -= answer_deriv.ky * y;
    answer.kz -= answer_deriv.kz * 6 * 100000;
    //}
}


int main() {
    data_t data;
    printf("Data is initialized\n");

    answer_t answer;
    double error = INFINITY;

    
    for (int i = 0; i+9 < random_epoch_N; i+=10) {
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

    error = data.rate_answer(answer);
    printf("\nSummary of randomness-based estimation:\n");
    printf("    Total square-error (rad): %#9.6f\n", error);
    printf("    Mean square-error  (rad): %#9.6f\n", error / (data_N * 5 - data.k_count));
    printf("    Standard error     (deg): %#9.6f\n", sqrt(error / (data_N * 5 - data.k_count))/PI*180);


    for (int i = 0; i+9 < deriv_epoch_N; i+=10) {
        run_deriv_epoch(data, answer);
        run_deriv_epoch(data, answer);
        run_deriv_epoch(data, answer);
        run_deriv_epoch(data, answer);
        run_deriv_epoch(data, answer);
        run_deriv_epoch(data, answer);
        run_deriv_epoch(data, answer);
        run_deriv_epoch(data, answer);
        run_deriv_epoch(data, answer);
        run_deriv_epoch(data, answer);

        data.eliminate_inconsistent(answer, MAX_ERROR);
    }

    // Print error-summary
    error = data.rate_answer(answer);
    printf("\nSummary of derivative-based estimation:\n");
    printf("    Total square-error (rad): %#9.6f\n", error);
    printf("    Mean square-error  (rad): %#9.6f\n", error / (data_N * 5 - data.k_count));
    printf("    Standard error     (deg): %#9.6f\n", sqrt(error / (data_N * 5 - data.k_count))/PI*180);

    // Print answer
    printf("\nAnswer: %#9.0f %#9.0f %#9.0f %#9.0f %#9.0f %#9.0f \n",
            answer.x0, answer.y0, answer.z0,
            answer.kx, answer.ky, answer.kz);

    // Print processed answer for each observer (include ignored)
    printf("\n");
    data.process_answer(answer);
    for (int i = 0; i < data_N; i++) {
        printf("Processed Answer (%i): %#7.3f | %7.3f | %#7.3f | %#7.3f | %#7.3f\n", i+1,
                data.ex_z0[i]/PI*180, data.ex_h0[i]/PI*180,
                data.ex_zb[i]/PI*180, data.ex_hb[i]/PI*180 ,
                data.ex_a[i]/PI*180);
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
            printf("Ignore: 'desent angle'   for observer %i\n", i+1);
    }

    // Exit
    printf("\n");
    return 0;
}
