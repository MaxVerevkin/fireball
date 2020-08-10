#include <cstdio>
#include "data.h"

#include "hyperparams.h"


/////////////////
/// Functions ///
/////////////////

/*
 * Run binary-tree search in 3D value-space.
 */
vec3d_t btree_flash_search(data_t &data) {
    vec3d_t flash_pos;
    vec3d_t min_val {lat0_min/180*PI, lon0_min/180*PI, z0_min};
    vec3d_t max_val {lat0_max/180*PI, lon0_max/180*PI, z0_max};
    flash_pos.x = (min_val.x + max_val.x) / 2;
    flash_pos.y = (min_val.y + max_val.y) / 2;
    flash_pos.z = (min_val.z + max_val.z) / 2;

    for (int i = 0; i < FLASH_SEARCH_DEPTH; i++) {

        // Found combination of parameter's
        // changes that gives least error.
        double min_err = INFINITY;
        int best_index_x, best_index_y, best_index_z;
        vec3d_t correction;
        for (int xi = 0; xi < 3; xi++) {
            correction.x = (max_val.x - min_val.x) / 100 * (xi-1);
            for (int yi = 0; yi < 3; yi++) {
                correction.y = (max_val.y - min_val.y) / 100 * (yi-1);
                for (int zi = 0; zi < 3; zi++) {
                    correction.z = (max_val.z - min_val.z) / 1000 * (zi-1);
                    double error = data.rate_flash_pos(flash_pos + correction, data.ex_data);
                    if (error < min_err) {
                        min_err = error;
                        best_index_x = xi;
                        best_index_y = yi;
                        best_index_z = zi;
                    }
                }
            }
        }

        //printf("Answer: %8.5f %8.5f %8.5f\n", flash_pos.x*180/PI, flash_pos.y*180/PI, flash_pos.z);
        //printf("Min: %8.5f %8.5f %8.5f\n", min_val.x*180/PI, min_val.y*180/PI, min_val.z);
        //printf("Max: %8.5f %8.5f %8.5f\n", max_val.x*180/PI, max_val.y*180/PI, max_val.z);

        if (best_index_x == 0) {
            max_val.x = (max_val.x +flash_pos.x) / 2;
            flash_pos.x = (flash_pos.x + min_val.x) / 2;
        } else if (best_index_x == 2) {
            min_val.x = (min_val.x + flash_pos.x) / 2;
            flash_pos.x = (flash_pos.x + max_val.x) / 2;
        }

        if (best_index_y == 0) {
            max_val.y = (max_val.y +flash_pos.y) / 2;
            flash_pos.y = (flash_pos.y + min_val.y) / 2;
        } else if (best_index_y == 2) {
            min_val.y = (min_val.y + flash_pos.y) / 2;
            flash_pos.y = (flash_pos.y + max_val.y) / 2;
        }

        if (best_index_z == 0) {
            max_val.z = (max_val.z +flash_pos.z) / 2;
            flash_pos.z = (flash_pos.z + min_val.z) / 2;
        } else if (best_index_z == 2){
            min_val.z = (min_val.z + flash_pos.z) / 2;
            flash_pos.z = (flash_pos.z + max_val.z) / 2;
        }
    }
    return flash_pos;
}
//vec3d_t btree_traj_search(data_t &data, const vec3d_t &flash) {
    //vec3d_t flash_traj;
    //vec3d_t min_val {kx_min, ky_min, kz_min};
    //vec3d_t max_val {kx_max, ky_max, kz_max};
    //flash_traj.x = (min_val.x + max_val.x) / 2;
    //flash_traj.y = (min_val.y + max_val.y) / 2;
    //flash_traj.z = (min_val.z + max_val.z) / 2;

    //for (int i = 0; i < TRAJ_SEARCH_DEPTH; i++) {
        //double e1 = data.rate_flash_traj(flash, {flash_traj.x+10,flash_traj.y+10,flash_traj.z+10}, data.ex_data);
        //double e2 = data.rate_flash_traj(flash, {flash_traj.x+10,flash_traj.y+10,flash_traj.z-10}, data.ex_data);
        //double e3 = data.rate_flash_traj(flash, {flash_traj.x+10,flash_traj.y-10,flash_traj.z+10}, data.ex_data);
        //double e4 = data.rate_flash_traj(flash, {flash_traj.x+10,flash_traj.y-10,flash_traj.z-10}, data.ex_data);
        //double e5 = data.rate_flash_traj(flash, {flash_traj.x-10,flash_traj.y+10,flash_traj.z+10}, data.ex_data);
        //double e6 = data.rate_flash_traj(flash, {flash_traj.x-10,flash_traj.y+10,flash_traj.z-10}, data.ex_data);
        //double e7 = data.rate_flash_traj(flash, {flash_traj.x-10,flash_traj.y-10,flash_traj.z+10}, data.ex_data);
        //double e8 = data.rate_flash_traj(flash, {flash_traj.x-10,flash_traj.y-10,flash_traj.z-10}, data.ex_data);
        //double min = min_8d(e1, e2, e3, e4, e5, e6, e7, e8);

        //if (min == e1 || min == e2 || min == e3 || min == e4) {
            //min_val.x = flash_traj.x;
            //flash_traj.x = (flash_traj.x + max_val.x) / 2;
        //} else {
            //max_val.x = flash_traj.x;
            //flash_traj.x = (flash_traj.x + min_val.x) / 2;
        //}

        //if (min == e1 || min == e2 || min == e5 || min == e6) {
            //min_val.y = flash_traj.y;
            //flash_traj.y = (flash_traj.y + max_val.y) / 2;
        //} else {
            //max_val.y = flash_traj.y;
            //flash_traj.y = (flash_traj.y + min_val.y) / 2;
        //}

        //if (min == e1 || min == e3 || min == e5 || min == e7) {
            //min_val.z = flash_traj.z;
            //flash_traj.z = (flash_traj.z + max_val.z) / 2;
        //} else {
            //max_val.z = flash_traj.z;
            //flash_traj.z = (flash_traj.z + min_val.z) / 2;
        //}
    //}
    //return flash_traj;
//}


/*
 * main
 */
int main() {
    data_t data;
    printf("Data is initialized\n");

    
    ///////////////////////////
    /// Find flash position ///
    ///////////////////////////

    vec3d_t flash_pos = btree_flash_search(data);
    data.eliminate_inconsistent_flash_data(flash_pos);
    flash_pos = btree_flash_search(data);
    double flash_error = data.rate_flash_pos(flash_pos, data.ex_data);

    printf("\nSummary on finding flash position:\n");
    printf("    Total square-error (rad): %#9.6f\n", flash_error);
    printf("    Mean square-error  (rad): %#9.6f\n", flash_error / (data.data_Ne * 2 - data.k_count));
    printf("    Standard error     (deg): %#9.6f\n", sqrt(flash_error / (data.data_Ne * 2 - data.k_count))/PI*180);


    ///////////////////////////////
    ///// Find flash trajectory ///
    ///////////////////////////////

    //vec3d_t flash_traj = btree_traj_search(data, flash_pos);
    //data.eliminate_inconsistent_traj_data(flash_pos, flash_traj);
    //flash_traj = btree_traj_search(data, flash_pos);
    //double traj_error = data.rate_flash_traj(flash_pos, flash_traj, data.ex_data);

    //printf("\nSummary on finding flash trajectory:\n");
    //printf("    Total square-error (rad): %#9.6f\n", traj_error);
    //printf("    Mean square-error  (rad): %#9.6f\n", traj_error / (data.data_Ne * 3 - data.k_count));
    //printf("    Standard error     (deg): %#9.6f\n", sqrt(traj_error / (data.data_Ne * 3 - data.k_count))/PI*180);


    // Print answer
    printf("\nAnswer: %8.5f %8.5f %8.5f\n",
            flash_pos.x*180/PI, flash_pos.y*180/PI, flash_pos.z);

    // Print processed answer for each observer
    printf("\n                  i     z0        h0        zb        hb         a         T   \n");
    for (int i = 0; i < data_N; i++) {
        printf("Processed Answer (%i): %#7.3f | %7.3f\n", i+1,
                data.ex_data.z0[i]*180/PI, data.ex_data.h0[i]*180/PI);
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
            printf("Ignore: 'desent_angle'   for observer %i\n", i+1);
    }

    // Exit
    printf("\n");
    return 0;
}
