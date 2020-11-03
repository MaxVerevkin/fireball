#include <cstdio>
#include "data.h"
#include "hyperparams.h"

#include <random>
std::default_random_engine e(time(0));



/*
 * Run binary-tree-like search.
 */
vec2d2_t btree_2d_search(data_t &data, data_t::rate_2d_fn err_fn, data_t::sigma_2d_fn sig_fn) {
    double range = PI/36.;
    vec2d_t min_val = { data.mean_lat - range, data.mean_lon - range };
    vec2d_t max_val = { data.mean_lat + range, data.mean_lon + range };
    vec2d_t pos = (min_val + max_val) * .5;

    for (int i = 0; i < FLASH_2D_SEARCH_N; i++) {

        // Lat
        for (int j = 0; j < FLASH_2D_SEARCH_DEPTH; j++) {
            vec2d_t correction = {(max_val.x - min_val.x) * .01, 0};
            double e1 = (data.*err_fn)(pos - correction);
            double e2 = (data.*err_fn)(pos + correction);

            if (e1 < e2)
                max_val.x = pos.x;
            else
                min_val.x = pos.x;

            pos.x = (min_val.x + max_val.x) * .5;
        }
        // Lon
        for (int j = 0; j < FLASH_2D_SEARCH_DEPTH; j++) {
            vec2d_t correction = {0, (max_val.y - min_val.y) * .01};
            double e1 = (data.*err_fn)(pos - correction);
            double e2 = (data.*err_fn)(pos + correction);

            if (e1 < e2)
                max_val.y = pos.y;
            else
                min_val.y = pos.y;
            
            pos.y = (min_val.y + max_val.y) * .5;
        }
        
        // Reset bounds
        min_val = { data.mean_lat - PI/18, data.mean_lon - PI/18};
        max_val = { data.mean_lat + PI/18, data.mean_lon + PI/18};
    }

    return {pos, (data.*sig_fn)(pos)};
}


/*
 * Returns a pair vectors:
 *  1) The middle point of a cloud
 *  2) The standard deviation.
 */
inline vec2d2_t btree_end_cloud(data_t &data) {
    return btree_2d_search(data, &data_t::rate_z0, &data_t::sigma_z0);
}
inline vec2d2_t btree_start_cloud(data_t &data) {
    return btree_2d_search(data, &data_t::rate_zb, &data_t::sigma_zb);
}


/*
 * Takes in two clouds: of the start and the end of the trajectory.
 * Returns trajectory that passes through those two clouds and
 * fits the observer's data best.
 */
line3d_t btree_traj_search(data_t &data, const vec2d2_t &start_cloud, const vec2d2_t &end_cloud) {

    std::normal_distribution<double> start_lat(start_cloud.v1.x, start_cloud.v2.x);
    std::normal_distribution<double> start_lon(start_cloud.v1.y, start_cloud.v2.y);

    std::normal_distribution<double> end_lat(end_cloud.v1.x, end_cloud.v2.x);
    std::normal_distribution<double> end_lon(end_cloud.v1.y, end_cloud.v2.y);

    std::normal_distribution<double> start_height(100000., 30000);
    std::normal_distribution<double> end_height(40000., 15000);

    double best_error = INFINITY;
    line3d_t best_traj;
    line3d_t traj;

    for (int i = 0; i < 50000; i++) {
        traj.start = {start_lat(e), start_lon(e), start_height(e)};
        traj.end = {end_lat(e), end_lon(e), end_height(e)};

        double error = data.rate_traj(traj);

        if (error < best_error) {
            best_error = error;
            best_traj = traj;
        }
    }

    return best_traj;
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

    
    ///////////////////////////////
    /// Find the flash position ///
    ///////////////////////////////
    
    // (lat, lon) coordinates
    vec2d2_t end_cloud = btree_end_cloud(data);
    vec2d2_t start_cloud = btree_start_cloud(data);
    vec2d2_t flash_2d = end_cloud;
    for (int i = 0; i < 7; i++) {
        data.eliminate_inconsistent_z0(flash_2d.v1);
        flash_2d = btree_end_cloud(data);
    }

    // height
    vec2d_t flash_height = data.calc_flash_height(flash_2d.v1); // (x => value, y => sigma)


    /////////////////////////////////
    /// Find the flash trajectory ///
    /////////////////////////////////

    line3d_t traj;
    //traj = {
        //{0.586361, 0.596406, 136777.302050},
        //{0.586076, 0.583937, 41169.047112}
    //};
    for (int i = 0; i < 10; i++) {
        traj = btree_traj_search(data, start_cloud, end_cloud);
        data.eliminate_inconsistent_traj_data(traj);
    }
    //printf("start: %f, %f, %f\nend: %f, %f, %f\n", traj.start.x, traj.start.y, traj.start.z, traj.end.x, traj.end.y, traj.end.z);
    data.process_traj(traj);
    vec3d_t flash_vel = global_to_local(traj.vec(), flash_2d.v1);


    // Print answer
    printf("\nAnswer:\n");
    printf("  Location:\n");
    printf("    lat =  %8.4f ± %6.4f(°)\n", flash_2d.v1.x*180/PI, flash_2d.v2.x*180/PI);
    printf("    lon =  %8.4f ± %6.4f(°)\n", flash_2d.v1.y*180/PI, flash_2d.v2.y*180/PI);
    printf("    z   =  %8.4f ± %6.4f(km)\n", flash_height.x/1000, flash_height.y/1000);
    printf("  Velocity (normalized):\n");
    printf("    Local:\n");
    printf("      v_East  =  %7.4f\n", flash_vel.x);
    printf("      v_North =  %7.4f\n", flash_vel.y);
    printf("      v_z     =  %7.4f\n", flash_vel.z);
    printf("    Global:\n");
    printf("      v_x = %7.4f\n", traj.vec().x);
    printf("      v_y = %7.4f\n", traj.vec().y);
    printf("      v_z = %7.4f\n", traj.vec().z);

    // Print processed answer for each observer
    printf("\nProcessed answer:\n");
    printf("       |            Start           |  |             End            |  |    Desent angle   |\n");
    printf(" i     |   z°         h°     E°   U |  |   z°         h°     E°   U |  |   A°       E°   U |\n");
    for (int i = 0; i < data.data_N; i++) {

        char start_u = ' ';
        if (data.traj_accept_start[i] == 0 || data.ob_data->hb[i] < 0 || data.ob_data->zb[i] < 0)
            start_u = 'x';
        else if (data.k_traj_start[i] == 1)
            start_u = '*';

        char end_u = ' ';
        if (data.traj_accept_end[i] == 0 || data.ob_data->h0[i] < 0 || data.ob_data->z0[i] < 0)
            end_u = 'x';
        else if (data.k_traj_end[i] == 1)
            end_u = '*';

        char a_u = ' ';
        if (data.ob_data->a[i] < 0)
            a_u = 'x';
        else if (data.k_traj_a[i] == 1)
            a_u = '*';

        printf("%3i |  | %3.0f/%3.0f | %3.0f/%3.0f | %2.0f | %c |  | %3.0f/%3.0f | %3.0f/%3.0f | %2.0f | %c |  | %3.0f/%3.0f | %3.0f | %c |\n", i+1,

                DEG(data.ex_data->zb[i]),
                DEG(data.ob_data->zb[i]),
                DEG(data.ex_data->hb[i]),
                DEG(data.ob_data->hb[i]),
                DEG(sqrt(data.traj_error_start[i])),
                start_u,

                DEG(data.ex_data->z0[i]),
                DEG(data.ob_data->z0[i]),
                DEG(data.ex_data->h0[i]),
                DEG(data.ob_data->h0[i]),
                DEG(sqrt(data.traj_error_end[i])),
                end_u,

                DEG(data.ex_data->a[i]),
                DEG(data.ob_data->a[i]),
                abs(DEG(angle_delta(data.ob_data->a[i], data.ex_data->a[i]))),
                a_u
                );
    }

    // Print stats on data filtering
    printf("\nStatistics on data filtering:\n");
    int z0_total = 0, z0_used = 0;
    int h0_total = 0, h0_used = 0;
    int start_total = 0, start_used = 0;
    int end_total = 0, end_used = 0;
    int a_total = 0, a_used = 0;
    for (int i = 0; i < data.data_N; i++) {
        if (data.ob_data->z0[i] >= 0) {
            z0_total++;
            if (data.k_z0[i])
                z0_used++;
        }
        if (data.ob_data->h0[i] >= 0) {
            h0_total++;
            if (data.k_h0[i])
                h0_used++;
        }
        if (data.ob_data->zb[i] >= 0 && data.ob_data->hb[i] >= 0) {
            start_total++;
            if (data.k_traj_start[i])
                start_used++;
        }
        if (data.ob_data->z0[i] >= 0 && data.ob_data->h0[i] >= 0) {
            end_total++;
            if (data.k_traj_end[i])
                end_used++;
        }
        if (data.ob_data->a[i] >= 0) {
            a_total++;
            if (data.k_traj_a[i])
                a_used++;
        }
    }
    printf("Azimuth End  - %i used from %i\n", z0_used, z0_total);
    printf("Altitude End - %i used from %i\n", h0_used, h0_total);
    printf("Start        - %i used from %i\n", start_used, start_total);
    printf("End          - %i used from %i\n", end_used, end_total);
    printf("Desent angle - %i used from %i\n", a_used, a_total);
    printf("Total - %i used from %i\n\n",
            z0_used + h0_used + start_used + end_used + a_used,
            z0_total + h0_total + start_total + end_total + a_total);

    // Exit
    return 0;
}
