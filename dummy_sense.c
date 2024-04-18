#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "common.h"

int parse_arguments(int argc, char *argv[], struct sim_args *args)
{
    // Define long and short options
    struct option long_options[] = {{"min", required_argument, 0, 'm'},
                                    {"max", required_argument, 0, 'M'},
                                    {"grid_deltas", required_argument, 0, 'd'},
                                    {"plume", optional_argument, 0, 'p'},
                                    {"wind", optional_argument, 0, 'w'},
                                    {"baseline", optional_argument, 0, 'b'},
                                    {"steps", required_argument, 0, 't'},
                                    {"out_steps", optional_argument, 0, 'o'},
                                    {"sensors", required_argument, 0, 's'},
                                    {0, 0, 0, 0}};

    // Parse command line arguments
    int option;
    int option_index = 0;
    while((option = getopt_long(argc, argv, "m:M:p:w:d:b:t:o:s:", long_options,
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
        case 's':
            break;
        default:
            printf("Usage: %s --min <float,float> --max <float,float> "
                   "--grid_deltas <float, float> --steps <int> [-- out_steps "
                   "<int> --plume <float, float>] [--wind <float, float>] "
                   "[--baseline float]\n",
                   argv[0]);
            return EXIT_FAILURE;
        }
    }

    if(!args->steps) {
        fprintf(stderr, "No steps argument found!\n\n");
        printf(
            "Usage: %s --min <float,float> --max <float,float> --grid_deltas "
            "<float, float> --steps <int> [-- out_steps <int> --plume <float, "
            "float>] [--wind <float, float>] [--baseline float]\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) { return (0); }