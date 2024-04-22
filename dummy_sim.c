#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "toml.h"

/*
    Write the simulation grid to a file
        - grid: the grid
        - filename: the filname to write to
        - step: which step is this (stored in the header)
*/
void fwrite_grid_data(struct sim_grid *grid, const char *filename, int step)
{
    FILE *file = fopen(filename, "wb");

    if(file != NULL) {
        fwrite(&step, sizeof(step), 1, file);
        fwrite(&grid->nx, sizeof(grid->nx), 1, file);
        fwrite(&grid->ny, sizeof(grid->ny), 1, file);
        fwrite(grid->args->min_coord, sizeof(*grid->args->min_coord), 2, file);
        fwrite(grid->args->max_coord, sizeof(*grid->args->min_coord), 2, file);
        fwrite(grid->args->plume_source, sizeof(*grid->args->plume_source), 2,
               file);
        fwrite(grid->args->wind, sizeof(*grid->args->wind), 2, file);

        // Write the array data to the file
        fwrite(grid->data[0], sizeof(*grid->data[0]), grid->nx * grid->ny,
               file);

        fclose(file); // Close the file
    } else {
        fprintf(stderr, "Unable to open file '%s'.\n", filename);
    }
}

int parse_arguments(int argc, char *argv[], struct sim_args *args)
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

    if(args->input_file) {
        if(parse_conf(args->input_file, args, NULL) < 0) {
            free(args->input_file);
            args->input_file = NULL;
        }
    }

    args->diffusivity = 1.0;

    if(!args->dt) {
        args->dt = .00001;
    }

    return EXIT_SUCCESS;
}

int main(int argc, char *argv[])
{
    struct sim_args args = {0};
    struct sim_grid *grid;
    char outbase[100], outfile[100];
    double max;
    double too_big = 5000;
    int t;

    if(parse_arguments(argc, argv, &args)) {
        return EXIT_FAILURE;
    }

    grid = init_grid(&args, NULL);
    if(!grid) {
        return EXIT_FAILURE;
    }

    printf_grid(grid);

    if(args.sim_out_dir && create_dir(args.sim_out_dir)) {
        sprintf(outbase, "%s/out.", args.sim_out_dir);
    } else {
        sprintf(outbase, "out.");
    }

    printf("Simulate %i steps.", args.steps);
    if(args.out_steps) {
        printf("Output every %i steps.\n", args.out_steps);
        if(args.sim_out_dir) {
            printf("Output will be placed in %s/\n", args.sim_out_dir);
        }
    }
    printf("\n");

    for(t = 0; t < args.steps; t++) {
        max = convect_diffuse(grid, NULL);
        if((args.out_steps && t % args.out_steps == 0) || max > too_big) {
            printf("Step %i\n", t);
            sprintf(outfile, "%s%i.dat", outbase, t);
            fwrite_grid_data(grid, outfile, t);
        }
        if(max > too_big) {
            break;
        }
    }

    return EXIT_SUCCESS;
}
