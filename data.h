#ifndef DATA_H
#define DATA_H

# define DATA_N 6 // Should be divisible be 2
# define data_N 6 // Should be divisible be 2

#include "utils.h"
#include "structs.h"

struct processed_answer;


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
    alignas(16) double *ob_a;
    alignas(16) double *ob_zb;
    alignas(16) double *ob_hb;
    alignas(16) double *ob_z0;
    alignas(16) double *ob_h0;
    alignas(16) double *ob_t;

    // Observer's expected data given some answer
    //processed_answer *expected_data;
    alignas(16) double ex_a[data_N];
    alignas(16) double ex_zb[data_N];
    alignas(16) double ex_hb[data_N];
    alignas(16) double ex_z0[data_N];
    alignas(16) double ex_h0[data_N];

    // Additional variables
    int k_count = 0;
    alignas(16) double k_a[data_N];
    alignas(16) double k_zb[data_N];
    alignas(16) double k_hb[data_N];
    alignas(16) double k_z0[data_N];
    alignas(16) double k_h0[data_N];
    alignas(16) vector2f_t pos_2d[data_N];

    /*
     * Initialize the data
     * and prepreocess it
     */
    data_t();

    /*
     * Sets K=0 to all data which square-error is
     * max_error_k times greater than mean square-error
     */
    void eliminate_inconsistent(const answer_t &answer, double max_error);


    /*
     * Returns square-error of a given answer.
     */
    double rate_answer(const answer_t &answer);
    double rate_answer(const answer_t &answer, processed_answer &dest);


    /*
     * Fills in the 'expected' values.
     */
    void process_answer(const answer_t &answer);
    void process_answer(const answer_t &answer, processed_answer &dest);
};


struct processed_answer {
    
    // Observer's expected data
    alignas(16) double a[data_N];
    alignas(16) double zb[data_N];
    alignas(16) double hb[data_N];
    alignas(16) double z0[data_N];
    alignas(16) double h0[data_N];

    /*
     * Just call process().
     */
    //processed_answer(const data_t &data, const answer_t &answer);

    /*
     * Fill in the 'expected' values.
     */
    //void process(const data_t &data, const answer_t &answer);

    /*
     * Returns square-error of given answer.
     */
    //double rate(const data_t &data);
};

#endif
