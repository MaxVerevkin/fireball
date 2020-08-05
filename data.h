#ifndef DATA_H
#define DATA_H

# define DATA_N 6 // Should be divisible be 2
# define data_N 6 // Should be divisible be 2

#include "utils.h"
#include "structs.h"


/*
 * Observer's expected data.
 */
struct processed_answer {    
    alignas(16) double z0[data_N];
    alignas(16) double h0[data_N];
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

    // Observer's expected data given some answer
    processed_answer ex_data;

    // Additional variables
    int k_count = 0;
    alignas(16) double k_z0[data_N];
    alignas(16) double k_h0[data_N];
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
    void eliminate_inconsistent(const vec3d_t &answer, double max_error);


    /*
     * Return square-error of a given answer.
     */
    double rate_flash_pos(const vec3d_t &pos, processed_answer &dest);
    inline double rate_flash_pos(const vec3d_t &pos) {
        return rate_flash_pos(pos, ex_data);
    }


    /*
     * Fills in the 'expected' values.
     */
    void process_answer(const vec3d_t &answer, processed_answer &dest);
    inline void process_answer(const vec3d_t &answer) {
        return process_answer(answer, ex_data);
    }
};



#endif
