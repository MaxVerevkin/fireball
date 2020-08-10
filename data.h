#ifndef DATA_H
#define DATA_H

#include "data_values.h"

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
    alignas(16) double ob_lat[data_N] =    {OB_LAT};    // Latitude
    alignas(16) double ob_lon[data_N] =    {OB_LON};    // Longitude
    alignas(16) double ob_height[data_N] = {OB_HEIGHT}; // Height
    alignas(16) double ob_z0[data_N] =     {OB_Z0};     // Azimuth end
    alignas(16) double ob_h0[data_N] =     {OB_H0};     // Altitude end
    alignas(16) double ob_zb[data_N] =     {OB_ZB};     // Azimuth begin
    alignas(16) double ob_hb[data_N] =     {OB_HB};     // Altitude begin
    alignas(16) double ob_a[data_N] =      {OB_A};      // Desent angle
    alignas(16) double ob_e[data_N] =      {OB_E};      // Experience

    // Observer's expected data given some answer
    processed_answer ex_data;

    // Additional variables
    int k_count = 0;
    int data_Ne = 0;
    alignas(16) vec3d_t ob_pos[data_N];
    alignas(16) vec3d_t normal[data_N];
    alignas(16) double r[data_N];
    // Used to ignore inconsistent data
    alignas(16) double k_z0[data_N];
    alignas(16) double k_h0[data_N];
    alignas(16) double k_zb[data_N];
    alignas(16) double k_hb[data_N];
    alignas(16) double k_a[data_N];

    /*
     * Initialize the data
     * and prepreocess it
     */
    data_t();

    /*
     * Sets K=0 to all data which square-error is
     * max_error_k times greater than mean square-error
     */
    void eliminate_inconsistent_flash_data(const vec3d_t &flash_geo);
    void eliminate_inconsistent_traj_data(const vec3d_t &flash_geo, const vec3d_t params);


    /*
     * Return square-error of a given answer.
     */
    double rate_flash_pos(const vec3d_t &flash_geo, processed_answer &dest);
    double rate_flash_traj(const vec3d_t &flash_geo, const vec3d_t &params, processed_answer &dest);

    /*
     * Return sigma (standard deviation)
     * of a given answer.
     */
    //vec3d_t sigma_flash_pos(const vec3d_t &flash_geo);
    //vec3d_t sigma_flash_traj(const vec3d_t &flash_geo, const vec3d_t &params, processed_answer &dest);


    /*
     * Fills in the 'expected' values.
     */
    void process_flash_pos(const vec3d_t &flash_geo, processed_answer &dest);
    void process_flash_traj(const vec3d_t &flash_geo, const vec3d_t &params, processed_answer &dest);

    /*
     * Proceese answer for the trajectory for one observer given 't'.
     */
    void process_flash_traj_i(const vec3d_t &rel_flash, const vec3d_t params, double t, int i, processed_answer &dest);

    /*
     * Returns square-error the trajectory for current answer given 't'
     */
    double traj_error_i(const vec3d_t &rel_flash, const vec3d_t params, double t, int i, processed_answer &dest);
};



#endif
