#ifndef DATA_H
#define DATA_H

# define data_N 6 // Should be divisible be 2

#include "utils.h"
#include "structs.h"


/*
 * Observer's expected data.
 */
struct processed_answer {
    alignas(16) double z0[data_N];
    alignas(16) double h0[data_N];
    alignas(16) double zb[data_N];
    alignas(16) double hb[data_N];
    alignas(16) double a[data_N];
    alignas(16) double t[data_N];
};


/*
 * This struct represents observer's data,
 * as well as some helper functions, that
 * operate on that data.
 */
struct data_t {

    // Observer's data
    alignas(16) double *ob_lat;
    alignas(16) double *ob_lon;
    alignas(16) double *ob_height;
    alignas(16) double *ob_z0;
    alignas(16) double *ob_h0;
    alignas(16) double *ob_zb;
    alignas(16) double *ob_hb;
    alignas(16) double *ob_a;

    // Observer's expected data given some answer
    processed_answer ex_data;

    // Additional variables
    int k_count = 0;
    alignas(16) double k_z0[data_N];
    alignas(16) double k_h0[data_N];
    alignas(16) double k_zb[data_N];
    alignas(16) double k_hb[data_N];
    alignas(16) double k_a[data_N];
    alignas(16) vec2d_t pos_2d[data_N];

    /*
     * Initialize the data
     * and prepreocess it
     */
    data_t();

    /*
     * Sets K=0 to all data which square-error is
     * max_error_k times greater than mean square-error
     */
    void eliminate_inconsistent_flash_data(const vec3d_t &pos, double max_error_k);
    void eliminate_inconsistent_traj_data(const vec3d_t &flash, const vec3d_t params, double max_error_k);


    /*
     * Return square-error of a given answer.
     */
    double rate_flash_pos(const vec3d_t &pos, processed_answer &dest);
    double rate_flash_traj(const vec3d_t &flash, const vec3d_t &params, processed_answer &dest);
    vec3d_t sigma_flash_pos(const vec3d_t &pos);


    /*
     * Fills in the 'expected' values.
     */
    void process_flash_pos(const vec3d_t &pos, processed_answer &dest);
    void process_flash_traj(const vec3d_t &flash, const vec3d_t &params, processed_answer &dest);
};



#endif
