#include <cstdio>
#include "data.h"
#include "hyperparams.h"

#include <random>

std::default_random_engine e(time(0));

/*
 * Run binary-tree-like search.
 */
vec2d_t btree_flash_2d_search(data_t &data) {

    // For the search, latitude and longitude
    // are defined in rangle of 20 degrees.
    //
    // This 20 degrees are calculated with:
    //      a = 2 * acos(R/(R+h)),
    // where R - radius of Earth
    //       h - maximum height of flash
    //
    // TODO get latitude into account.

    //vec2d_t min_val = { data.mean_lat - PI/18, data.mean_lon - PI/18};
    //vec2d_t max_val = { data.mean_lat + PI/18, data.mean_lon + PI/18};
    vec2d_t min_val = { data.mean_lat - PI/36, data.mean_lon - PI/36};
    vec2d_t max_val = { data.mean_lat + PI/36, data.mean_lon + PI/36};
    vec2d_t pos = (min_val + max_val) * .5;


    std::uniform_real_distribution<double> rand(min_val.x, max_val.x);
    double best_error = INFINITY;

    for (int i = 0; i < 100000; i++) {
        vec2d_t v = {rand(e), rand(e)};

        double error = data.rate_z0(v);
        if (error < best_error) {
            best_error = error;
            pos = v;
        }
    }

    //for (int i = 0; i < FLASH_2D_SEARCH_N; i++) {

        //// Lat
        //for (int j = 0; j < FLASH_2D_SEARCH_DEPTH; j++) {
            //vec2d_t correction = {(max_val.x - min_val.x) * .01, 0};
            //double e1 = data.rate_z0(pos - correction);
            //double e2 = data.rate_z0(pos + correction);

            //if (e1 < e2)
                //max_val.x = pos.x;
            //else
                //min_val.x = pos.x;

            //pos.x = (min_val.x + max_val.x) * .5;
        //}
        //// Lon
        //for (int j = 0; j < FLASH_2D_SEARCH_DEPTH; j++) {
            //vec2d_t correction = {0, (max_val.y - min_val.y) * .01};
            //double e1 = data.rate_z0(pos - correction);
            //double e2 = data.rate_z0(pos + correction);

            //if (e1 < e2)
                //max_val.y = pos.y;
            //else
                //min_val.y = pos.y;
            
            //pos.y = (min_val.y + max_val.y) * .5;
        //}
        
        //// Reset bounds
        //min_val = { data.mean_lat - PI/18, data.mean_lon - PI/18};
        //max_val = { data.mean_lat + PI/18, data.mean_lon + PI/18};
    //}

    return pos;
}


inline vec3d_t gen_vec(double z, double theta) {
    double r = sqrt(1. - z*z);
    return {cos(theta) * r, sin(theta) * r, z};
}
vec3d_t btree_traj_search(data_t &data, const vec3d_t flash_pos) {
    //return {7.5/27.4, 23.5/27.4, 11.9/27.4};


    std::normal_distribution<double> rand(0.,1.);
    double best_error = INFINITY;
    vec3d_t best_ans;

    for (int i = 0; i < 2000; i++) {
        vec3d_t traj = {rand(e), rand(e), rand(e)};
        traj = traj.normalized();

        double error = data.rate_flash_traj(flash_pos, traj);
        if (error < best_error) {
            best_error = error;
            best_ans = traj;
        }
    }
    return best_ans;
}


/*
 * main
 */
int main(int argc, char **argv) {

    /////////////////
    /// Init data ///
    /////////////////

    data_t data(argc == 2 ? argv[1] : "data.txt");
    printf("Data is initialized\n");

    
    ///////////////////////////
    /// Find flash position ///
    ///////////////////////////

    // (lat, lon) coordinates
    vec2d_t flash_geo = btree_flash_2d_search(data);
    for (int i = 0; i < 2; i++) {
        data.eliminate_inconsistent_z0(flash_geo);
        flash_geo = btree_flash_2d_search(data);
    }

    double flash_geo_error = data.rate_z0(flash_geo);
    double n = data.data_Ne - data.k_count_z0;
    printf("\nSummary on finding flash position:\n");
    printf("    Total square-error         : %#9.6f rad\n", flash_geo_error);
    printf("    Standard deviation of meaan: %#9.6f°\n", sqrt(flash_geo_error / n / (n-1))/PI*180);

    // height
    vec3d_t flash_pos_geo = {flash_geo.x, flash_geo.y, data.calc_flash_height(flash_geo)};

    //vec2d_t flash_geo = {33.1/180*PI, 34.3/180*PI};
    //vec3d_t flash_pos_geo = {33.1/180*PI, 34.3/180*PI, 27400};
    //data.eliminate_inconsistent_z0(flash_geo);
    //data.eliminate_inconsistent_h0(flash_pos_geo);

    vec3d_t flash_pos = geo_to_xyz(flash_pos_geo);
    vec3d_t flash_sigma = data.sigma_flash_pos(flash_pos_geo);


    ///////////////////////////////
    ///// Find flash trajectory ///
    ///////////////////////////////

    vec3d_t flash_traj = btree_traj_search(data, flash_pos);
    for (int i = 0; i < 5; i++) {
        data.eliminate_inconsistent_traj_data(flash_pos, flash_traj);
        flash_traj = btree_traj_search(data, flash_pos);
    }
    vec3d_t flash_vel = data.get_flash_vel(flash_pos_geo, flash_traj);


    double traj_error = data.rate_flash_traj(flash_pos, flash_traj);
    n = data.data_Ne * 3 - data.k_count_traj;

    double velocity = flash_vel.length();
    data.normalize_t(velocity);

    printf("\nSummary on finding flash trajectory:\n");
    printf("    Total square-error:          %#9.6f rad\n", traj_error);
    printf("    Standard deviation of meaan: %#9.6f°\n", sqrt(traj_error / n / (n-1))/PI*180);


    // Print answer
    printf("\nAnswer:\n");
    printf("  Location:\n");
    printf("    lat =  %8.4f ± %6.4f(°)\n", flash_pos_geo.x*180/PI, flash_sigma.x*180/PI);
    printf("    lon =  %8.4f ± %6.4f(°)\n", flash_pos_geo.y*180/PI, flash_sigma.y*180/PI);
    printf("    z   =  %8.4f ± %6.4f(km)\n", flash_pos_geo.z/1000, flash_sigma.z/1000);
    printf("  Velocity: %8.3f(km/s)\n", velocity/1000);
    printf("    Local:\n");
    printf("      v_East  =  %8.3f(km/s)\n", flash_vel.x/1000);
    printf("      v_North =  %8.3f(km/s)\n", flash_vel.y/1000);
    printf("      v_z     =  %8.3f(km/s)\n", flash_vel.z/1000);
    printf("    Global:\n");
    printf("      v_x = %8.3f(km/s)\n", -flash_traj.x * velocity / 1000);
    printf("      v_y = %8.3f(km/s)\n", -flash_traj.y * velocity / 1000);
    printf("      v_z = %8.3f(km/s)\n", -flash_traj.z * velocity / 1000);

    // Print processed answer for each observer
    printf("\n                  i     z0°       h0°       zb°       hb°        A°        t(s)    w(°/s)\n");
    for (int i = 0; i < data.data_N; i++) {
        printf("Processed Answer (%i): %#7.3f | %7.3f | %7.3f | %7.3f | %7.3f | %7.3f | ", i+1,
                data.ex_data->z0[i]*180/PI, data.ex_data->h0[i]*180/PI,
                data.ex_data->zb[i]*180/PI, data.ex_data->hb[i]*180/PI,
                data.ex_data->a[i]*180/PI, data.ex_data->t[i]);
        printf("%7.3f\n", arc_len({data.ex_data->h0[i],
                    data.ex_data->z0[i]},
                    {data.ex_data->hb[i],
                    data.ex_data->zb[i]})
                /data.ex_data->t[i]/PI*180);
    }

    // Print ignored data
    printf("\n");
    int n_ignored = 0;
    for (int i = 0; i < data.data_N; i++) {
        if (!data.k_z0[i] && data.ob_data->z0[i] >= 0) {
            printf("Ignore: 'azimuth end'    for observer %i\n", i+1);
            n_ignored++;
        }
        if (!data.k_h0[i] && data.ob_data->h0[i] >= 0) {
            printf("Ignore: 'altitude end'   for observer %i\n", i+1);
            n_ignored++;
        }
        if (!data.k_zb[i] && data.ob_data->zb[i] >= 0) {
            printf("Ignore: 'azimuth begin'  for observer %i\n", i+1);
            n_ignored++;
        }
        if (!data.k_hb[i] && data.ob_data->hb[i] >= 0) {
            printf("Ignore: 'altitude begin' for observer %i\n", i+1);
            n_ignored++;
        }
        if (!data.k_a[i] && data.ob_data->a[i] >= 0) {
            printf("Ignore: 'desent_angle'   for observer %i\n", i+1);
            n_ignored++;
        }
    }
    printf("Total ingored: %i\n\n", n_ignored);

    // Exit
    return 0;
}
