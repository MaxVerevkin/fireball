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

    int data_N;

    // Observer's data
    alignas(16) double *ob_lat;    // Latitude
    alignas(16) double *ob_lon;    // Longitude
    alignas(16) double *ob_height; // Height
    alignas(16) double *ob_e;      // Experience

    data_set_t *ob_data; // Data given by observer
    data_set_t *ex_data; // Data produced by answer

    // Additional variables
    int k_count_z0 = 0;
    int k_count_h0 = 0;
    int k_count_traj = 0;
    int data_Ne = 0;
    double mean_lat;
    double mean_lon;
    alignas(16) vec3d_t *ob_pos;
    alignas(16) vec3d_t *normal;
    alignas(16) double *r;
    // Used to ignore inconsistent data
    alignas(16) double *k_z0;
    alignas(16) double *k_h0;
    alignas(16) double *k_zb;
    alignas(16) double *k_hb;
    alignas(16) double *k_a;

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
     * Sets K=0 to all data which square-error is
     * max_error_k times greater than mean square-error
     */
    void eliminate_inconsistent_z0(const vec2d_t &flash_geo);
    void eliminate_inconsistent_h0(const vec3d_t &flash_geo);
    void eliminate_inconsistent_traj_data(const vec3d_t &flash_pos, const vec3d_t params);

    /*
     * Translate flash trajectory vector to local velocity.
     */
    vec3d_t get_flash_vel(const vec3d_t &flash_geo, const vec3d_t &traj);

    /*
     * Normalize aobserver's 't'
     */
    void normalize_t(const vec3d_t &flash_vel);


    /*
     * Return square-error of a given answer.
     */
    double rate_z0(const vec2d_t &flash_geo);
    double rate_h0(const vec3d_t &flash_geo);
    double rate_flash_traj(const vec3d_t &flash_pos, const vec3d_t &params);

    /*
     * Return sigma (standard deviation)
     * of a given answer.
     */
    //vec3d_t sigma_flash_pos(const vec3d_t &flash_geo);
    //vec3d_t sigma_flash_traj(const vec3d_t &flash_pos, const vec3d_t &params);


    /*
     * Fills in the 'expected' values.
     */
    void process_z0(const vec2d_t &flash_geo);
    void process_h0(const vec3d_t &flash_geo);
    void process_flash_traj(const vec3d_t &flash_pos, const vec3d_t &params);

    /*
     * Proceese answer for the trajectory for one observer given 't'.
     */
    void process_flash_traj_i(const vec3d_t &flash_pos, const vec3d_t params, double t, int i);

    /*
     * Returns square-error the trajectory for current answer given 't'
     */
    double traj_error_i(const vec3d_t &flash_pos, const vec3d_t params, double t, int i);
};



#endif
