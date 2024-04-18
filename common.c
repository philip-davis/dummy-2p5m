#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<sys/stat.h>

#include "common.h"
#include "toml.h"

void printf_grid(struct sim_grid *grid)
{
    struct sim_args *passed_args = grid->args;

    if(passed_args->input_file) {
        printf("-=-Input Source-=-\n\n");
        printf("%s (overrides CLI arguments)\n\n", passed_args->input_file);
    }

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

    printf("\n--Grid Parameters--\n");
    printf("Dimensions: (%i, %i)\n", grid->nx, grid->ny);
    if(grid->plumex != -1 && grid->plumey != -1) {
        printf("Plume grid point: (%i, %i)\n", grid->plumex, grid->plumey);
    } else {
        printf("No plume\n");
    }
    printf("\n\n");
}

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

static void get_toml_pair(toml_table_t *t, const char *key, double d[2])
{
    toml_array_t *arr;
    toml_datum_t dat;

    arr = toml_array_in(t, key);
    if(arr) {
        dat = toml_double_at(arr, 0);
        d[LAT_IDX] = dat.u.d;
        dat = toml_double_at(arr, 1);
        d[LONG_IDX] = dat.u.d;
    }
}

static void get_toml_int(toml_table_t *t, const char *key, int *i)
{
    toml_datum_t dat;

    if(toml_key_exists(t, key)) {
        dat = toml_int_in(t, key);
        *i = dat.u.i;
    }
}

static void get_toml_double(toml_table_t *t, const char *key, double *d)
{
    toml_datum_t dat;

    if(toml_key_exists(t, key)) {
        dat = toml_double_in(t, key);
        *d = dat.u.d;
    }
}

static void get_toml_str(toml_table_t *t, const char *key, char **str)
{
    toml_datum_t dat;

    if(toml_key_exists(t, key)) {
        dat = toml_string_in(t, key);
        *str = strdup(dat.u.s);
    }
}

int parse_conf(const char *filename, struct sim_args *args)
{
     FILE *f;
     toml_table_t *conf, *model, *sim, *sensors, *grid;
     toml_datum_t dat;
     toml_array_t *arr;
     char errbuf[200];

    f = fopen(filename, "r");
    if(!f) {
        fprintf(stderr, "WARNING: could not open configuration file '%s'.\n",
                filename);
        return(-1);
    }

    conf = toml_parse_file(f, errbuf, sizeof(errbuf));
    if(!conf) {
        fprintf(stderr, "could not parse %s, %s.\n", filename, errbuf);
        return -1;
    }
    fclose(f);

    model = toml_table_in(conf, "model");
    get_toml_pair(model, "plume", args->plume_source);
    get_toml_pair(model, "wind", args->wind);
    get_toml_double(model, "baseline", &args->baseline);
    get_toml_int(model, "steps", &args->steps);

    grid = toml_table_in(model, "grid");
    get_toml_pair(grid, "min", args->min_coord);
    get_toml_pair(grid, "max", args->max_coord);
    get_toml_pair(grid, "delta", args->grid_deltas);

    sim = toml_table_in(conf, "sim");
    get_toml_int(sim, "out_steps", &args->out_steps);
    get_toml_str(sim, "out_dir", &args->sim_out_dir);

    toml_free(conf);

    return(0);
}

int create_dir(const char *dir)
{
    struct stat st = {0};

    // Check if the directory already exists
    if(stat(dir, &st) == -1) {
        if(mkdir(dir, 0700) == -1) {
            perror("mkdir");
            return 0;
        } else {
            return 1;
        }
    } else {
        // Directory already exists
        return 1;
    }
}