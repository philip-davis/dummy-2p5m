#ifndef __DUMMY_2P5_H__
#define __DUMMY_2P5_H__

// TODO: make sure axes are right, i.e. x:y::longitude:latitute
//    Because they are NOT right somewhere resulting in transposition.
#define LAT_IDX 0
#define LONG_IDX 1

struct sim_args {
    double min_coord[2]; // lower bounds of sim
    double max_coord[2]; // upper bounds of sim
    double grid_deltas[2]; // dx, dy
    double plume_source[2]; // source location
    double wind[2]; // wind (convection) vector
    /*
        baseline has three effects:
            1. sets conditions at f(x,0)
            2. sets boundary conditons for tail end of convection
            3. source injects 10 x baseline every t
    */
    double baseline;
    int steps; // how many steps to simulate
    int out_steps; // how often to create output file
    char *input_file;
    char *sim_out_dir;
};

struct sim_grid {
    int nx, ny; // number of columns and rows
    int plumex, plumey; // plume/source matrix indices
    double **data; // the data matrix
    struct sim_args *args; // inputed simulation parameters
};

/*
        Print the grid paramters (not the grid data)
*/
void printf_grid(struct sim_grid *grid);

/*
    Do a convection-diffusion step with explicit advance (Euler's method)
    The diffusion and convection coefficients are trivial, and numerical
    instability is VERY LIKELY for many values that result in a snappy 
    run time on a laptop.

    Return the largest value in the new timestep (artifact of troubleshooting)
*/
double convect_diffuse(struct sim_grid *grid, double dt);

int parse_conf(const char *filename, struct sim_args *args);

int create_dir(const char *dir);

#endif // __DUMMY_2P5_H__