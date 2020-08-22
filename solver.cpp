#include <cstdio>
#include "data.h"
#include "hyperparams.h"


/*
 * Run binary-tree-like search in 3D value-space.
 */
vec3d_t btree_flash_search(data_t &data) {

    // For the search, latitude and longitude
    // are defined between -10 and 10 degrees.
    // This 10 degrees are calculated with:
    //      a = acos(R/(R+h)),
    // where R - radius of Earth
    //       h - maximum height of flash
    //
    // An offset (mean of observer's values)
    // is added to get real value.
    vec3d_t min_val {-PI/18, -PI/18, z0_min};
    vec3d_t max_val { PI/18,  PI/18, z0_max};
    vec3d_t flash_pos = {0, 0, (min_val.z + max_val.z) / 2};
    vec3d_t flash_pos_offset {data.mean_lat, data.mean_lon, 0};

    for (int i = 0; i < FLASH_SEARCH_N; i++) {

        // X
        for (int j = 0; j < FLASH_SEARCH_DEPTH; j++) {
            flash_pos.x = (min_val.x + max_val.x) / 2;

            vec3d_t correction = {(max_val.x - min_val.x) / 100, 0, 0};
            double e1 = data.rate_flash_pos(flash_pos + flash_pos_offset - correction);
            double e2 = data.rate_flash_pos(flash_pos + flash_pos_offset + correction);

            if (e1 < e2)
                max_val.x = (max_val.x + flash_pos.x) / 2;
            else
                min_val.x = (min_val.x + flash_pos.x) / 2;
        }

        // Y
        for (int j = 0; j < FLASH_SEARCH_DEPTH; j++) {
            flash_pos.y = (min_val.y + max_val.y) / 2;

            vec3d_t correction {0, (max_val.y - min_val.y) / 100, 0};
            double e1 = data.rate_flash_pos(flash_pos + flash_pos_offset - correction);
            double e2 = data.rate_flash_pos(flash_pos + flash_pos_offset + correction);

            if (e1 < e2)
                max_val.y = (max_val.y + flash_pos.y) / 2;
            else
                min_val.y = (min_val.y + flash_pos.y) / 2;
        }

        // Z
        for (int j = 0; j < FLASH_SEARCH_DEPTH; j++) {
            flash_pos.z = (min_val.z + max_val.z) / 2;

            vec3d_t correction = {0, 0, (max_val.z - min_val.z) / 100};
            double e1 = data.rate_flash_pos(flash_pos + flash_pos_offset - correction);
            double e2 = data.rate_flash_pos(flash_pos + flash_pos_offset + correction);

            if (e1 < e2)
                max_val.z = (max_val.z + flash_pos.z) / 2;
            else
                min_val.z = (min_val.z + flash_pos.z) / 2;
        }

        // Reset bounds
        min_val = {-PI/18, -PI/18, z0_min};
        max_val = { PI/18,  PI/18, z0_max};

    }
    return flash_pos + flash_pos_offset;
}
vec3d_t btree_traj_search(data_t &data, const vec3d_t flash_pos) {
    vec3d_t min_val {-1, -1, -1};
    vec3d_t max_val {1, 1, 1};
    vec3d_t flash_traj;

    for (int i = 0; i < TRAJ_SEARCH_DEPTH; i++) {
        flash_traj.x = (min_val.x + max_val.x) / 2;
        flash_traj.y = (min_val.y + max_val.y) / 2;
        flash_traj.z = (min_val.z + max_val.z) / 2;

        // Found combination of parameter's
        // changes that gives least error.
        double min_err = INFINITY;
        int best_index_x=0, best_index_y=0, best_index_z=0;
        vec3d_t correction;
        for (int xi = -1; xi <= 1; xi+=2) {
            correction.x = (max_val.x - min_val.x) / 10 * xi;
            for (int yi = -1; yi <= 1; yi+=2) {
                correction.y = (max_val.y - min_val.y) / 10 * yi;
                for (int zi = -1; zi <= 1; zi+=2) {
                    correction.z = (max_val.z - min_val.z) / 10 * zi;
                    double error = data.rate_flash_traj(flash_pos, (flash_traj + correction).normalized());
                    if (error < min_err) {
                        min_err = error;
                        best_index_x = xi;
                        best_index_y = yi;
                        best_index_z = zi;
                    }
                }
            }
        }

        if (best_index_x == -1)
            max_val.x = (max_val.x +flash_traj.x) / 2;
        else
            min_val.x = (min_val.x + flash_traj.x) / 2;

        if (best_index_y == -1)
            max_val.y = (max_val.y +flash_traj.y) / 2;
        else
            min_val.y = (min_val.y + flash_traj.y) / 2;

        if (best_index_z == -1)
            max_val.z = (max_val.z +flash_traj.z) / 2;
        else
            min_val.z = (min_val.z + flash_traj.z) / 2;
    }
    return flash_traj.normalized();
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

    vec3d_t flash_pos = btree_flash_search(data);
    data.eliminate_inconsistent_flash_data(flash_pos);
    flash_pos = btree_flash_search(data);
    double flash_error = data.rate_flash_pos(flash_pos);

    printf("\nSummary on finding flash position:\n");
    printf("    Total square-error         : %#9.6f rad\n", flash_error);
    double n = data.data_Ne * 2 - data.k_count_flash;
    printf("    Standard deviation of meaan: %#9.6f°\n", sqrt(flash_error / n / (n-1))/PI*180);


    ///////////////////////////////
    ///// Find flash trajectory ///
    ///////////////////////////////

    vec3d_t flash_traj = btree_traj_search(data, flash_pos);
    data.eliminate_inconsistent_traj_data(flash_pos, flash_traj);
    flash_traj = btree_traj_search(data, flash_pos);
    double traj_error = data.rate_flash_traj(flash_pos, flash_traj);
    vec3d_t flash_vel = data.get_flash_vel(flash_pos, flash_traj);
    data.normalize_t(flash_vel);

    printf("\nSummary on finding flash trajectory:\n");
    printf("    Total square-error:          %#9.6f rad\n", traj_error);
    n = data.data_Ne * 3 - data.k_count_traj;
    printf("    Standard deviation of meaan: %#9.6f°\n", sqrt(traj_error / n / (n-1))/PI*180);


    // Print answer
    printf("\nAnswer:\n");
    printf("    lat =  %8.5f°\n", flash_pos.x*180/PI);
    printf("    lon =  %8.5f°\n", flash_pos.y*180/PI);
    printf("    z   =  %8.5f(km)\n", flash_pos.z/1000);
    printf("    vx  =  %8.5f(km/s)\n", flash_vel.x/1000);
    printf("    vy  =  %8.5f(km/s)\n", flash_vel.y/1000);
    printf("    vz  =  %8.5f(km/s)\n", flash_vel.z/1000);

    // Print processed answer for each observer
    printf("\n                  i     z0        h0        zb        hb         A         t   \n");
    for (int i = 0; i < data.data_N; i++) {
        printf("Processed Answer (%i): %#7.3f | %7.3f | %7.3f | %7.3f | %7.3f | %f\n", i+1,
                data.ex_data->z0[i]*180/PI, data.ex_data->h0[i]*180/PI,
                data.ex_data->zb[i]*180/PI, data.ex_data->hb[i]*180/PI,
                data.ex_data->a[i]*180/PI, data.ex_data->t[i]);
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
