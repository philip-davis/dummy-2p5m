#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"

struct sensor {
    int loc[2];
    int fine_loc[2];
    int mobile;
};

struct reading {
    int id;
    double loc[2];
    double value;
};

/*
    Do a convection-diffusion step with explicit advance (Euler's method)
    The diffusion and convection coefficients are trivial, and numerical
    instability is VERY LIKELY for many values that result in a snappy
    run time on a laptop.

    Return the largest value in the new timestep (artifact of troubleshooting)
*/
double convect_diffuse(struct sim_grid *grid, struct sensor_args *sargs)
{
    static double **new_data = NULL;
    struct sim_args *args = grid->args;
    double i_bias, j_bias;
    double w_lat = args->wind[LAT_IDX];
    double w_long = args->wind[LONG_IDX];
    double dx = args->grid_deltas[LONG_IDX];
    double dy = args->grid_deltas[LAT_IDX];
    double **data = grid->data;
    double conv_du, diff_du;
    double max = -1;
    double dt = args->dt;
    int nx = grid->nx;
    int ny = grid->ny;
    int i, j;

    // keep reusing the same static buffer
    if(!new_data) {
        new_data = malloc(sizeof(*new_data) * ny);
        new_data[0] = malloc(sizeof(*new_data[0]) * nx * ny);
        for(i = 1; i < ny; i++) {
            new_data[i] = &(new_data[0][i * nx]);
        }
    }

    // TODO: fiddle with coefficients? Better put in some stability checks at
    // least.
    for(i = 0; i < ny; i++) {
        /* the point of i_bias and j_bias is that we can be rational about
           convection at the boundaries towards which the wind is blowing, but
           we can't in the direction the wind is coming from (thar be monsters).
           We bias the indices by the corresponding components of the wind
           vector, and don't try to calculate the convection component if that
           puts us out of bounds.

            Very small values in the wind vector might cause a round-off error
           bug.
        */
        i_bias = (double)i - w_long;
        for(j = 0; j < nx; j++) {
            if(sargs) {
                // wind varies stochastically if wind shift is set
                w_lat = uniform_sample_d_interval(
                    args->wind[LAT_IDX] - sargs->wind_shift[LAT_IDX],
                    args->wind[LAT_IDX] + sargs->wind_shift[LAT_IDX]);
                w_long = uniform_sample_d_interval(
                    args->wind[LONG_IDX] - sargs->wind_shift[LONG_IDX],
                    args->wind[LONG_IDX] + sargs->wind_shift[LONG_IDX]);
                i_bias = (double)i - w_long;
            }
            // Convection component
            j_bias = (double)j - w_lat;
            conv_du = 0;
            if(i_bias > 0 && i_bias < ny - 1 && j_bias > 0 && j_bias < nx - 1) {
                // choose forward or backwards difference to align with the wind
                // direction
                if(w_long > 0) {
                    conv_du +=
                        (w_long * (data[i][j] - data[i - 1][j])) / (2 * dy);
                } else {
                    conv_du +=
                        (w_long * (data[i + 1][j] - data[i][j])) / (2 * dy);
                }
                if(w_lat > 0) {
                    conv_du +=
                        (w_lat * (data[i][j] - data[i][j - 1])) / (2 * dx);
                } else {
                    conv_du +=
                        (w_lat * (data[i][j + 1] - data[i][j])) / (2 * dx);
                }
            }

            // Diffusion Component
            diff_du = 0;
            if(i > 0 && i < ny - 1 && j > 0 && j < nx - 1) {
                diff_du += (data[i + 1][j] - 2 * data[i][j] + data[i - 1][j]) /
                           (dy * dy);
                diff_du += (data[i][j + 1] - 2 * data[i][j] + data[i][j - 1]) /
                           (dx * dx);
            }

            new_data[i][j] =
                data[i][j] + dt * (args->diffusivity * diff_du - conv_du);
            if(sargs && sargs->variation) {
                new_data[i][j] = uniform_sample_d_interval(
                    new_data[i][j] * (1 - sargs->variation),
                    new_data[i][j] * (1 + sargs->variation));
            }
            if(new_data[i][j] > max) {
                max = new_data[i][j];
            }
        }
    }

    // plume source
    if(grid->plumex != -1 && grid->plumey != -1) {
        new_data[grid->plumex][grid->plumey] += args->source;
    }

    memcpy(data[0], new_data[0], sizeof(*new_data[0]) * nx * ny);

    return (max);
}


/*
    allocate a sensor array. Each sensor has a locaiton (nearby
    gridpoint) and a "fine location" (a positive offest from that location.)
    Some are mobile, and wander the simulation grid.
*/
void init_sensors(struct sim_grid *grid, struct sensor **sarray)
{
    struct sensor_args *sargs = grid->sargs;
    int count = sargs->count;
    struct sensor *s;
    int i;

    *sarray = malloc(sizeof(**sarray) * count);
    for(i = 0; i < count; i++) {
        s = &(*sarray)[i];
        s->loc[LAT_IDX] = rand() % (grid->ny - 1);
        s->loc[LONG_IDX] = rand() % (grid->nx - 1);
        s->fine_loc[LAT_IDX] = rand() % sargs->speed;
        s->fine_loc[LONG_IDX] = rand() % sargs->speed;
        if(i < sargs->mobile) {
            s->mobile = 1;
        }
    }
}

int parse_arguments(int argc, char *argv[], struct sim_args *args,
                    struct sensor_args *sargs)
{
    int nchar;
    // Define long and short options
    struct option long_options[] = {{"min", optional_argument, 0, 'm'},
                                    {"max", optional_argument, 0, 'M'},
                                    {"grid_deltas", optional_argument, 0, 'd'},
                                    {"plume", optional_argument, 0, 'p'},
                                    {"wind", optional_argument, 0, 'w'},
                                    {"baseline", optional_argument, 0, 'b'},
                                    {"steps", optional_argument, 0, 't'},
                                    {"out_steps", optional_argument, 0, 'o'},
                                    {"input_file", optional_argument, 0, 'i'},
                                    {0, 0, 0, 0}};

    // Parse command line arguments
    int option;
    int option_index = 0;
    while((option = getopt_long(argc, argv, "m:M:p:w:d:b:t:o:i:", long_options,
                                &option_index)) != -1) {
        switch(option) {
        case 'd':
            sscanf(optarg, "%lf,%lf", &args->grid_deltas[LAT_IDX],
                   &args->grid_deltas[LONG_IDX]);
            break;
        case 'm':
            sscanf(optarg, "%lf,%lf", &args->min_coord[LAT_IDX],
                   &args->min_coord[LONG_IDX]);
            break;
        case 'M':
            sscanf(optarg, "%lf,%lf", &args->max_coord[LAT_IDX],
                   &args->max_coord[LONG_IDX]);
            break;
        case 'p':
            sscanf(optarg, "%lf,%lf", &args->plume_source[LAT_IDX],
                   &args->plume_source[LONG_IDX]);
            args->plume = 1;
            break;
        case 'w':
            sscanf(optarg, "%lf,%lf", &args->wind[LAT_IDX],
                   &args->wind[LONG_IDX]);
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
        case 'i':
            sscanf(optarg, "%*s%n", &nchar);
            args->input_file = malloc(nchar + 1);
            sscanf(optarg, "%s", args->input_file);
            break;
        default:
            fprintf(stderr,
                    "Usage: %s (--input-file <filename> | --min <float,float> "
                    "--max <float,float> --grid_deltas <float, float> --steps "
                    "<int> [--out_steps <int> --plume <float, float>] [--wind "
                    "<float, float>] [--baseline float])\n",
                    argv[0]);
            return EXIT_FAILURE;
        }
    }

    if((!args->steps ||
        (args->min_coord[LAT_IDX] == args->max_coord[LAT_IDX]) ||
        (args->min_coord[LONG_IDX] == args->max_coord[LONG_IDX]) ||
        !args->grid_deltas[LAT_IDX] || !args->grid_deltas[LONG_IDX]) &&
       !args->input_file) {
        fprintf(stderr, "Missing required arguments!\n\n");
        fprintf(stderr,
                "Usage: %s (--input-file <filename> | --min <float,float> "
                "--max <float,float> --grid_deltas <float, float> --steps "
                "<int> [--out_steps <int> --plume <float, float>] [--wind "
                "<float, float>] [--baseline float])\n",
                argv[0]);
        return EXIT_FAILURE;
    }

    args->diffusivity = 1.0;

    if(args->input_file) {
        if(parse_conf(args->input_file, args, sargs) < 0) {
            free(args->input_file);
            args->input_file = NULL;
        }
    }

    if(!args->dt) {
        args->dt = .00001;
    }

    return EXIT_SUCCESS;
}

#define SENSOR_MOVE_N 0
#define SENSOR_MOVE_S 1
#define SENSOR_MOVE_E 2
#define SENSOR_MOVE_W 3

/*
    Move all the mobile sensors. There's a 50% chance of a given mobile
    sensor moving, 12.5% in each cardinal direction. If a move would take
    it out of bounds, it does not move. Movement is by 1 "tick". Speed sets
    the inverse size of a tick; the grid points are spaced speed ticks apart.
*/
void move_sensors(struct sim_grid *grid, struct sensor *sarray)
{
    struct sensor_args *sargs = grid->sargs;
    struct sim_args *args = grid->args;

    int i;

    for(i = 0; i < sargs->mobile; i++) {
        switch(rand() % 8) {
        case SENSOR_MOVE_N:
            sarray[i].fine_loc[LAT_IDX] =
                (sarray[i].fine_loc[LAT_IDX] + 1) % sargs->speed;
            if(sarray[i].fine_loc[LAT_IDX] == 0) {
                sarray[i].loc[LAT_IDX]++;
            } else if(sarray[i].loc[LAT_IDX] == (grid->ny - 1)) {
                // beyond northern extereme
                sarray[i].fine_loc[LAT_IDX] = 0;
            }
            break;
        case SENSOR_MOVE_S:
            sarray[i].fine_loc[LAT_IDX] =
                (sargs->speed + sarray[i].fine_loc[LAT_IDX] - 1) % sargs->speed;
            if(sarray[i].fine_loc[LAT_IDX] == sargs->speed - 1) {
                sarray[i].loc[LAT_IDX]--;
                if(sarray[i].loc < 0) {
                    // beyond southern extreme
                    sarray[i].loc[LAT_IDX] = 0;
                    sarray[i].fine_loc[LAT_IDX] = 0;
                }
            }
            break;
        case SENSOR_MOVE_E:
            sarray[i].fine_loc[LONG_IDX] =
                (sarray[i].fine_loc[LONG_IDX] + 1) % sargs->speed;
            if(sarray[i].fine_loc[LONG_IDX] == 0) {
                sarray[i].loc[LONG_IDX]++;
            } else if(sarray[i].loc[LONG_IDX] == (grid->ny - 1)) {
                // beyond eastern extereme
                sarray[i].fine_loc[LONG_IDX] = 0;
            }
            break;
        case SENSOR_MOVE_W:
            sarray[i].fine_loc[LONG_IDX] =
                (sargs->speed + sarray[i].fine_loc[LONG_IDX] - 1) %
                sargs->speed;
            if(sarray[i].fine_loc[LONG_IDX] == sargs->speed - 1) {
                sarray[i].loc[LONG_IDX]--;
                if(sarray[i].loc < 0) {
                    // beyond western extreme
                    sarray[i].loc[LONG_IDX] = 0;
                    sarray[i].fine_loc[LONG_IDX] = 0;
                }
            }
            break;
        default:
            // 50% chance of no attempted movement
            break;
        }
    }
}

/*
    Sample the sensor, with noise. The sample is a bilinear
    interpolation of nearby gridpoints.
*/
double sample_sensor(struct sim_grid *grid, struct sensor *s)
{
    struct sensor_args *sargs = grid->sargs;
    struct sim_args *args = grid->args;
    double true_val;
    int x, y;
    int speed = sargs->speed;
    double noise = sargs->noise;
    double frac;
    double q[6];

    x = s->loc[LONG_IDX];
    y = s->loc[LAT_IDX];
    q[0] = grid->data[y][x];

    // I handle this by cases because I'm afraid of out of bounds memory access
    // if we're at the far N or E boundaries. Not very elegant - I'm sure
    // there's a nicer way to do this.
    if(s->fine_loc[LONG_IDX] == 0) {
        if(s->fine_loc[LAT_IDX] == 0) {
            // right on a grid point
            true_val = q[0];
        } else {
            // on a N/S grid line; linear interpolation
            q[1] = grid->data[y + 1][x];
            frac = (double)(s->fine_loc[LAT_IDX]) / speed;
            true_val = q[1] * frac + q[0] * (1. - frac);
        }
    } else {
        if(s->fine_loc[LAT_IDX] == 0) {
            // on an E/W grid line; liner interpolation
            q[1] = grid->data[y][x + 1];
            frac = (double)(s->fine_loc[LONG_IDX]) / speed;
            true_val = q[1] * frac + q[0] * (1. - frac);
        } else {
            // bilinear interpolation
            /*
                q[0] SW corner
                q[1] SE corner
                q[2] NW corner
                q[3] NE corner
                q[4] S end of final interp
                q[5] N end of final interp
            */
            q[1] = grid->data[y][x + 1];
            q[2] = grid->data[y + 1][x];
            q[3] = grid->data[y + 1][x + 1];
            frac = (double)(s->fine_loc[LONG_IDX]) / speed;
            q[4] = q[1] * frac + q[0] * (1 - frac);
            q[5] = q[3] * frac + q[2] * (1 - frac);
            frac = (double)(s->fine_loc[LAT_IDX]) / speed;
            true_val = q[4] * frac + q[5] * (1 - frac);
        }
    }

    return (uniform_sample_d_interval(true_val * (1 - noise),
                                      true_val + (1 + noise)));
}

int gather_readings(struct sim_grid *grid, struct sensor *sarray,
                    struct reading **readings)
{
    static struct reading *sensor_out = NULL;
    struct sim_args *args = grid->args;
    struct sensor_args *sargs = grid->sargs;
    int i, n = 0;

    // keep reusing the same buffer for readings
    if(!sensor_out) {
        sensor_out = calloc(sizeof(*sensor_out), sargs->count);
    }
    *readings = sensor_out;

    for(i = 0; i < sargs->count; i++) {
        if(random() % 100 < sargs->rate) {
            sensor_out[n].id = i;
            sensor_out[n].value = sample_sensor(grid, &sarray[i]);
            sensor_out[n].loc[LAT_IDX] =
                args->min_coord[LAT_IDX] +
                args->grid_deltas[LAT_IDX] * sarray[i].loc[LAT_IDX] +
                (double)(sarray[i].fine_loc[LAT_IDX]) *
                    (args->grid_deltas[LAT_IDX] / sargs->speed);
            sensor_out[n].loc[LONG_IDX] =
                args->min_coord[LONG_IDX] +
                args->grid_deltas[LONG_IDX] * sarray[i].loc[LONG_IDX] +
                (double)(sarray[i].fine_loc[LONG_IDX]) *
                    (args->grid_deltas[LONG_IDX] / sargs->speed);
            n++;
        }
    }

    return (n);
}

int main(int argc, char *argv[])
{
    struct sim_args args = {0};
    struct sim_grid *grid;
    struct sensor_args sargs;
    struct sensor *sarray;
    float max;
    char outbase[100], outfile[100];
    struct reading *readings;
    int num_reads;
    int i, t;

    if(parse_arguments(argc, argv, &args, &sargs)) {
        return EXIT_FAILURE;
    }

    grid = init_grid(&args, &sargs);
    if(!grid) {
        return EXIT_FAILURE;
    }

    init_sensors(grid, &sarray);
    if(!sarray) {
        return EXIT_FAILURE;
    }

    for(t = 0; t < args.steps; t++) {
        max = convect_diffuse(grid, &sargs);
        if(sargs.interval && t % sargs.interval == 0) {
            move_sensors(grid, sarray);
            num_reads = gather_readings(grid, sarray, &readings);
            for(i = 0; i < num_reads; i++) {
                printf("{ t:%i, loc: (%lf, %lf), value: %lf }\n", t,
                       readings[i].loc[0], readings[i].loc[1],
                       readings[i].value);
            }
        }
    }
}