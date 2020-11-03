#ifndef DATA_H
#define DATA_H

#include "utils.h"
#include "structs.h"


/*
 * This struct represents observer's data,
 * as well as some helper functions, that
 * operate on that data.
 */
struct data_t {

    // Functions' defenitions
    typedef double (data_t::*rate_2d_fn)(const vec2d_t&);
    typedef vec2d_t (data_t::*sigma_2d_fn)(const vec2d_t&);

    // Observers' number
    int data_N;

    // Observer's data
    vec3d_t *ob_pos;     // Position in global coorinates (x,y,z)
    vec3d_t *ob_pos_geo; // Position in geographical coordinates (lat,lon,z)

    data_set_t *ob_data; // Data given by observer
    data_set_t *ex_data; // Data produced by answer

    // Error of trajecroty
    double *traj_error_start;
    double *traj_error_end;
    double *traj_accept_start;
    double *traj_accept_end;

    // Additional variables
    double mean_lat;
    double mean_lon;
    // Used to ignore inconsistent data
    double *k_z0;
    double *k_h0;
    double *k_traj_start;
    double *k_traj_end;
    double *k_traj_a;

    /*
     * Initialize the data
     * and prepreocess it
     */
    data_t(const char *file);

    /*
     * Reset k-values to '1'
     */
    void reset_k_z0();
    void reset_k_h0();
    void reset_k_traj();

    /*
     * Sets K=0 for all data which square-error is
     * max_error_k times greater than mean square-error
     */
    void eliminate(double count, double total, double *errors, double *k, double max_e, double acc); // Abstracttial function
    void eliminate_inconsistent_z0(const vec2d_t &flash_geo);
    void eliminate_inconsistent_traj_data(const line3d_t &traj);
    

    /*
     * Return square-error of a given answer.
     */
    double rate_z0(const vec2d_t &flash_geo);
    double rate_zb(const vec2d_t &flash_geo);
    double rate_traj(const line3d_t &traj);


    /*
     * Calculate the height of a flash given it's location.
     */
    vec2d_t calc_flash_height(const vec2d_t &flash_2d);

    /*
     * Return sigma (standard deviation)
     * of a given answer.
     */
    vec2d_t sigma_z0(const vec2d_t &flash);
    vec2d_t sigma_zb(const vec2d_t &flash);


    /*
     * Fills in the 'expected' values.
     */
    void process_z0(const vec2d_t &flash_geo);
    void process_zb(const vec2d_t &flash_geo);
    void process_traj(const line3d_t &traj_line);


    /*
     * Calculate the magnitude of velocity given trajectory.
     */
    double calc_speed(const line3d_t &traj_line);
};



#endif
