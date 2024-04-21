#ifndef __DUMMY_2P5_H__
#define __DUMMY_2P5_H__

// TODO: make sure axes are right, i.e. x:y::longitude:latitute
//    Because they are NOT right somewhere resulting in transposition.
#define LAT_IDX 0
#define LONG_IDX 1

struct sim_args {
    double plume_source[2]; // source location
    double wind[2];         // wind (convection) vector
    /*
        baseline has three effects:
            1. sets conditions at f(x,0)
            2. sets boundary conditons for tail end of convection
            3. source injects 10 x baseline every t
    */
    double baseline;
    int steps; // how many steps to simulate
    double dt;

    double min_coord[2];   // lower bounds of sim
    double max_coord[2];   // upper bounds of sim
    double grid_deltas[2]; // dx, dy

    int out_steps;     // how often to create output file
    char *sim_out_dir; // simulation output folder

    char *input_file; // shared toml config file
};

struct sensor_args {
    int count;    // number of sensors
    int interval; // atomic interval for sensor checking
    int rate;     // percent chance a given sensor reports in an interval
    int mobile;   // how many sensors are mobile
    int speed; // how many steps for a mobile sensor to cross a pixel - must be
               // >1!
    double noise; // how much noise in a reading

    // environmental arguments that fit with sensors for now
    double wind_shift[2]; // magnitude by which wind can vary
    double variation;     // magnitude of particulate variation
};

struct sim_grid {
    int nx, ny;                // number of columns and rows
    int plumex, plumey;        // plume/source matrix indices
    double **data;             // the data matrix
    struct sim_args *args;     // inputed simulation parameters
    struct sensor_args *sargs; // inputed sensor parameters
};

/*
    Initialize a grid based on passed arguments. sargs can be NULL, but
    args must be present. Allocate and populate a 2D grid data structure, etc.
*/
struct sim_grid *init_grid(struct sim_args *args, struct sensor_args *sargs);

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
double convect_diffuse(struct sim_grid *grid, struct sensor_args *sargs);

/*
    Parse a .toml file to get startup arguments.
*/
int parse_conf(const char *filename, struct sim_args *args,
               struct sensor_args *sargs);

/*
    A helper function to create a directory. Succeeds if the directory already
   exists.
*/
int create_dir(const char *dir);

/*
    A helper function to uniformly sample in an interval
*/
double uniform_sample_d_interval(double x0, double x1);

#endif // __DUMMY_2P5_H__