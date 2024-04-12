#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
};

struct sim_grid {
    int nx, ny; // number of columns and rows
    int plumex, plumey; // plume/source matrix indices
    double **data; // the data matrix
    struct sim_args *args; // inputed simulation parameters
};

/*
    Write the simulation grid to a file
        - grid: the grid
        - filename: the filname to write to
        - step: which step is this (stored in the header)
*/
void fwrite_grid_data(struct sim_grid *grid, const char *filename, int step) {
    FILE *file = fopen(filename, "wb");
    
    if (file != NULL) {
        fwrite(&step, sizeof(step), 1, file);
        fwrite(&grid->nx, sizeof(grid->nx), 1, file);
        fwrite(&grid->ny, sizeof(grid->ny), 1, file);
        fwrite(grid->args->min_coord, sizeof(*grid->args->min_coord), 2, file);
        fwrite(grid->args->max_coord, sizeof(*grid->args->min_coord), 2, file);
        fwrite(grid->args->plume_source, sizeof(*grid->args->plume_source), 2, file);
        fwrite(grid->args->wind, sizeof(*grid->args->wind), 2, file);

        // Write the array data to the file
        fwrite(grid->data[0], sizeof(*grid->data[0]), grid->nx * grid->ny, file);

        fclose(file); // Close the file
    } else {
        fprintf(stderr, "Unable to open file '%s'.\n", filename);
    }
}

/*
        Print the grid paramters (not the grid data)
*/
void printf_grid(struct sim_grid *grid)
{
    struct sim_args *passed_args = grid->args;

    printf("-=-Model Setup-=-\n\n");

    printf("--Passed Arguments--\n");
    printf("Min Value: %lf, %lf\n", passed_args->min_coord[LAT_IDX], passed_args->min_coord[LONG_IDX]);
    printf("Max Value: %lf, %lf\n", passed_args->max_coord[LAT_IDX], passed_args->max_coord[LONG_IDX]);
    printf("Grid Deltas: %lf, %lf\n", passed_args->grid_deltas[LAT_IDX], passed_args->grid_deltas[LONG_IDX]);
    if(passed_args->plume_source[LAT_IDX] != 0 && passed_args->plume_source[LONG_IDX] != 0) {
        printf("Plume at %lf, %lf\n", passed_args->plume_source[LAT_IDX], passed_args->plume_source[LONG_IDX]);
    } else {
        printf("No Plume coordinates set.\n");
    }
    printf("Wind vector: %lf, %lf\n", passed_args->wind[LAT_IDX], passed_args->wind[LONG_IDX]);
    printf("Baseline particulate: %lf\n", passed_args->baseline);

    printf("--Grid Parameters--\n");
    printf("Dimensions: (%i, %i)\n", grid->nx, grid->ny);
    if(grid->plumex != -1 && grid->plumey != -1) {
        printf("Plume grid point: (%i, %i)\n", grid->plumex, grid->plumey);
    } else {
        printf("No plume\n");
    }


}

/*
    Do a convection-diffusion step with explicit advance (Euler's method)
    The diffusion and convection coefficients are trivial, and numerical
    instability is VERY LIKELY for many values that result in a snappy 
    run time on a laptop.

    Return the largest value in the new timestep (artifact of troubleshooting)
*/
double convect_diffuse(struct sim_grid *grid, double dt)
{
    static double **new_data = NULL;
    double i_bias, j_bias;
    double w_lat = grid->args->wind[LAT_IDX];
    double w_long = grid->args->wind[LONG_IDX];
    double dx = grid->args->grid_deltas[LONG_IDX];
    double dy = grid->args->grid_deltas[LAT_IDX];
    double **data = grid->data;
    double conv_du, diff_du;
    double max = -1;
    int i, j;

    // keep reusing the same static buffer
    if(!new_data) {
        new_data = malloc(sizeof(*new_data) * grid->ny);
        new_data[0] = malloc(sizeof(*new_data[0]) * grid->nx * grid->ny);
        for(i = 1; i < grid->ny; i++) {
            new_data[i] = &(new_data[0][i * grid->nx]);
        }
    }
    
    // TODO: fiddle with coefficients? Better put in some stability checks at least.
    for(i = 0; i < grid->ny; i++) {
        /* the point of i_bias and j_bias is that we can be rational about convection
            at the boundaries towards which the wind is blowing, but we can't in the
            direction the wind is coming from (thar be monsters). We bias the indices
            by the corresponding components of the wind vector, and don't try to calculate
            the convection component if that puts us out of bounds.

            Very small values in the wind vector might cause a round-off error bug.
        */ 
        i_bias = (double)i - w_long;
        for(j = 0; j < grid->nx; j++) {
            // Convection component
            j_bias = (double)j - w_lat;
            conv_du = 0;
            if(i_bias > 0 && i_bias < grid->ny-1 && j_bias > 0 && j_bias < grid->nx-1) {
                // choose forward or backwards difference to align with the wind direction
                if(w_long > 0) {
                    conv_du += (w_long * (data[i][j] - data[i-1][j])) / (2 * dy);
                } else {
                    conv_du += (w_long * (data[i+1][j] - data[i][j])) / (2 * dy);
                }
                if(w_lat > 0) {
                    conv_du += (w_lat * (data[i][j] - data[i][j-1])) / (2 * dx);
                } else {
                    conv_du += (w_lat * (data[i][j+1] - data[i][j])) / (2 * dx);
                }
            }

            // Diffusion Component
            diff_du = 0;
            if(i > 0 && i < grid->ny-1 && j > 0 && j < grid->nx-1) {
                diff_du += (data[i+1][j] - 2 * data[i][j] + data[i-1][j]) / (dy * dy);
                diff_du += (data[i][j+1] - 2 * data[i][j] + data[i][j-1]) / (dx * dx);
            }

            new_data[i][j] = data[i][j] + dt * (diff_du - conv_du);
            if(new_data[i][j] > max) {
                max = new_data[i][j];
            }
        }
    }

    // plume source
    if(grid->plumex != -1 && grid->plumey != -1) {
        new_data[grid->plumex][grid->plumey] += 10 * grid->args->baseline;
    }

    memcpy(data[0], new_data[0], sizeof(*new_data[0]) * grid->nx * grid->ny);

    return(max);
}

int parse_arguments(int argc, char *argv[], struct sim_args *args)
{
     // Define long and short options
    struct option long_options[] = {
        {"min", required_argument, 0, 'm'},
        {"max", required_argument, 0, 'M'},
        {"grid_deltas", required_argument, 0, 'd'},
        {"plume", optional_argument, 0, 'p'},
        {"wind", optional_argument, 0, 'w'},
        {"baseline", optional_argument, 0, 'b'},
        {"steps", required_argument, 0, 't'},
        {"out_steps", optional_argument, 0, 'o'},
        {0, 0, 0, 0}
    };

    // Parse command line arguments
    int option;
    int option_index = 0;
    while ((option = getopt_long(argc, argv, "m:M:p:w:d:b:t:o:", long_options, &option_index)) != -1) {
        switch (option) {
            case 'd':
                sscanf(optarg, "%lf,%lf", &args->grid_deltas[LAT_IDX], &args->grid_deltas[LONG_IDX]);
                break;
            case 'm':
                sscanf(optarg, "%lf,%lf", &args->min_coord[LAT_IDX], &args->min_coord[LONG_IDX]);
                break;
            case 'M':
                sscanf(optarg, "%lf,%lf", &args->max_coord[LAT_IDX], &args->max_coord[LONG_IDX]);
                break;
            case 'p':
                sscanf(optarg, "%lf,%lf", &args->plume_source[LAT_IDX], &args->plume_source[LONG_IDX]);
                break;
            case 'w':
                sscanf(optarg, "%lf,%lf", &args->wind[LAT_IDX], &args->wind[LONG_IDX]);
                break;
            case 'b':
                sscanf(optarg, "%lf", &args->baseline);
                break;
            case 't':
                sscanf(optarg, "%i", &args->steps);
                break;
            case 'o':
                sscanf(optarg, "%i", &args->out_steps);
                break;
            default:
                printf("Usage: %s --min <float,float> --max <float,float> --grid_deltas <float, float> --steps <int> [-- out_steps <int> --plume <float, float>] [--wind <float, float>] [--baseline float]\n", argv[0]);
                return EXIT_FAILURE;
        }
    }

    if(!args->steps) {
        fprintf(stderr, "No steps argument found!\n\n");
        printf("Usage: %s --min <float,float> --max <float,float> --grid_deltas <float, float> --steps <int> [-- out_steps <int> --plume <float, float>] [--wind <float, float>] [--baseline float]\n", argv[0]);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

struct sim_grid *init_grid(struct sim_args *args)
{
    struct sim_grid *grid;
    int i, j;

    // Need a few more than this!
    if (fabs(args->min_coord[LAT_IDX]) + fabs(args->min_coord[LONG_IDX]) == 0.0
            || fabs(args->max_coord[LAT_IDX]) + fabs(args->max_coord[LONG_IDX]) == 0.0
            || fabs(args->grid_deltas[LAT_IDX]) + fabs(args->grid_deltas[LONG_IDX]) == 0.0) {
        printf("Error: --min, --max, and --grid_deltas arguments are required.\n");
        return NULL;
    }

    // Try to pick up some obvious errors
    if(args->min_coord[LAT_IDX] > args->max_coord[LAT_IDX]
        || args->min_coord[LONG_IDX] > args->max_coord[LONG_IDX]) {
        printf("Minimum coordinates should be smaller than maximum coordinates.\n");
        return NULL;
    }

    if(args->grid_deltas[LAT_IDX] <= 0 || args->grid_deltas[LONG_IDX] <= 0) {
        printf("Grid deltas must be positive.\n");
        return NULL;
    }

    grid = malloc(sizeof(*grid));
    grid->nx = 1 + (int)((args->max_coord[LONG_IDX] - args->min_coord[LONG_IDX]) / args->grid_deltas[LONG_IDX]);
    grid->ny = 1 + (int)((args->max_coord[LAT_IDX] - args->min_coord[LAT_IDX]) / args->grid_deltas[LAT_IDX]);

    grid->args = args;

    grid->plumex = -1;
    grid->plumey = -1;
    // This is a bad check, since (0,0) is an obvious coordinate to choose
    if(args->plume_source[LAT_IDX] != 0 || args->plume_source[LONG_IDX] != 0) {
        if(args->plume_source[LAT_IDX] <= args->min_coord[LAT_IDX] ||
            args->plume_source[LAT_IDX] >= args->max_coord[LAT_IDX] ||
            args->plume_source[LONG_IDX] <= args->min_coord[LONG_IDX]) {
                printf("Plume must be between min and max coordinates.\n");
                return NULL;
        } else {
            grid->plumex = (int)((args->plume_source[LONG_IDX] - args->min_coord[LONG_IDX])/args->grid_deltas[LONG_IDX]);
            grid->plumey = (int)((args->plume_source[LAT_IDX] - args->min_coord[LAT_IDX])/args->grid_deltas[LAT_IDX]);
        }
    }

    grid->data = malloc(sizeof(*grid->data) * grid->ny);
    grid->data[0] = malloc(sizeof(*grid->data[0]) * grid->nx * grid->ny);
    for(i = 0; i < grid->ny; i++) {
        if(i) {
            grid->data[i] = &(grid->data[0][i * grid->nx]);
        }
        for(j = 0; j < grid->nx; j++) {
            grid->data[i][j] = grid->args->baseline;
        }
    }

    return(grid);
}

int main(int argc, char *argv[]) {
    struct sim_args args = {0};
    struct sim_grid *grid;
    char outfile[100];
    double max;
    double too_big = 5000;
    int t;

   
    if(parse_arguments(argc, argv, &args)) {
        return EXIT_FAILURE;
    }

    grid = init_grid(&args);
    if(!grid) {
        return EXIT_FAILURE;
    }

    printf_grid(grid);

    for(t = 0; t < args.steps; t++) {
        max = convect_diffuse(grid, .00001);
        if((args.out_steps && t % args.out_steps == 0) || max > too_big) {
            printf("Step %i\n", t);
            sprintf(outfile, "out.%i.dat", t);
            fwrite_grid_data(grid, outfile, t);
        }
        if(max > too_big) {
            break;
        }
    }

    return EXIT_SUCCESS;
}
