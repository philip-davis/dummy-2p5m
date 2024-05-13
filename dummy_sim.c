#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "common.h"
#include "toml.h"

struct sim_pgrid {
    int nx, ny;                // number of columns and rows
    int ox, oy;                // offset within the global grid
    int gx, gy;                // global grid size
    int plumex, plumey;        // plume/source matrix indices
    double *ghost[4];        // ghost regions 
    double **data;             // the data matrix
    double **fuel;             // available fuel
    double **fire;             // fire intensity
    struct sim_args *args;     // inputed simulation parameters
    struct sensor_args *sargs; // inputed sensor parameters
};

#define WEST_RANK 0
#define EAST_RANK 1
#define NORTH_RANK 2
#define SOUTH_RANK 3

struct sim_app {
    int (*bounds)[4];
    struct sim_args args;
    struct sim_pgrid pgrid;
    MPI_Comm comm;
    int nranks[4];             // neighbor ranks
    int rank;
    int size;
    int px, py;                // process dimensions
    int x, y;                  // this rank's process coordinates
};

/*
    Write the simulation grid to a file
        - grid: the grid
        - filename: the filname to write to
        - step: which step is this (stored in the header)
*/
void fwrite_grid_data(struct sim_pgrid *pgrid, double *data, const char *filename, int step)
{
    FILE *file = fopen(filename, "wb");

    if(file != NULL) {
        fwrite(&step, sizeof(step), 1, file);
        fwrite(&pgrid->gx, sizeof(pgrid->gx), 1, file);
        fwrite(&pgrid->gy, sizeof(pgrid->gy), 1, file);
        fwrite(pgrid->args->min_coord, sizeof(*pgrid->args->min_coord), 2, file);
        fwrite(pgrid->args->max_coord, sizeof(*pgrid->args->min_coord), 2, file);
        fwrite(pgrid->args->plume_source, sizeof(*pgrid->args->plume_source), 2,
               file);
        fwrite(pgrid->args->wind, sizeof(*pgrid->args->wind), 2, file);

        // Write the array data to the file
        fwrite(data, sizeof(*data), pgrid->gx * pgrid->gy,
               file);

        fclose(file); // Close the file
    } else {
        fprintf(stderr, "Unable to open file '%s'.\n", filename);
    }
}

int parse_arguments(int argc, char *argv[], struct sim_app *app)
{
    struct option long_options[] = {{"ranks", required_argument, 0, 'r'},
                                    {"input_file", required_argument, 0, 'i'},
                                    {0, 0, 0, 0}};
    int option;
    int option_index = 0;
    int nchar;

    while((option = getopt_long(argc, argv, "i:r:", long_options,
                                &option_index)) != -1) {
        switch(option) {
        case 'i':
            sscanf(optarg, "%*s%n", &nchar);
            app->args.input_file = malloc(nchar + 1);
            sscanf(optarg, "%s", app->args.input_file);
            break;
        case 'r':
            sscanf(optarg, "%i,%i", &app->px, &app->py);
            break;
        default:
            fprintf(stderr, "Usage: %s --input_file <filename> --ranks <integer, integer>\n", argv[0]);
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

int validate_args(struct sim_app *app)
{
    if(app->size != app->px * app->py) {
        if(!app->rank) {
            fprintf(stderr, "ERROR: %i rank(s), but %i expected.\n", app->size, (app->px * app->py));
        }
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

int validate_config(struct sim_app *app) 
{
    struct sim_args *args = &app->args;
    int config_bad = 0;

    if(args->steps <= 0) {
        if(!app->rank) {
            fprintf(stderr, "INVALID CONFIG: model.steps should be present and positive.\n");
        }
        config_bad = 1;
    }
    if((args->min_coord[LAT_IDX] == args->max_coord[LAT_IDX]) ||
        (args->min_coord[LONG_IDX] == args->max_coord[LONG_IDX])) {
        if(!app->rank) {
            fprintf(stderr, "INVALID CONFIG: model.grid.min should be less than modle.grid.max.\n");
        }
        config_bad = 1;
    }
    if((args->grid_deltas[LAT_IDX] * args->grid_deltas[LONG_IDX]) == 0) {
        if(!app->rank) {
            fprintf(stderr, "INVALID CONFIG: model.grid.delta must be populated and positive in both dimensions.\n");
        }
        config_bad = 1;
    }
    if(args->dt <= 0) {
        if(!app->rank) {
            fprintf(stderr, "INVALID CONFIG: model.dt must be populated and positive.\n");
        }
        config_bad = 1;
    }
    if(args->diffusivity == 0) {
         if(!app->rank) {
            fprintf(stderr, "INVALID CONFIG: diffusivity zero is invalid...changing to 1.0\n");
        }
        args->diffusivity = 1.0;
    }
    if(args->plume && (args->plume_source[LAT_IDX] <= args->min_coord[LAT_IDX] ||
        args->plume_source[LAT_IDX] >= args->max_coord[LAT_IDX] ||
        args->plume_source[LONG_IDX] <= args->min_coord[LONG_IDX] ||
        args->plume_source[LONG_IDX] >= args->max_coord[LONG_IDX])) {
        if(!app->rank) {
            fprintf(stderr, "INVALID CONFIG: plume must be between min and max coordinates.\n");
        }
        config_bad = 1;
    }
    if(config_bad) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

int x_y_within_coords(int x, int y, int lbx, int lby, int lx, int ly)
{
    if(x >= lbx && x < lbx + lx && y >= lby && y < lby + ly) {
        return(1);
    }
    return(0);
}

int init_sim_grid(struct sim_app *app)
{

    struct sim_pgrid *pgrid = &app->pgrid;
    struct sim_args *args = &app->args;
    int gplumex, gplumey;
    int bounds[4];
    int i, j;

    pgrid->args = &app->args;

    pgrid->gx =
        1 + (int)((args->max_coord[LONG_IDX] - args->min_coord[LONG_IDX]) /
                  args->grid_deltas[LONG_IDX]);
    pgrid->gy = 1 + (int)((args->max_coord[LAT_IDX] - args->min_coord[LAT_IDX]) /
                         args->grid_deltas[LAT_IDX]);

    pgrid->nx = pgrid->gx / app->px;
    pgrid->ny = pgrid->gy / app->py;
    pgrid->ox = pgrid->nx * app->x;
    pgrid->oy = pgrid->ny * app->y;
    if((pgrid->ox + pgrid->nx) >= pgrid->gx) {
        pgrid->nx = (pgrid->gx - pgrid->ox);
    }
    if((pgrid->oy + pgrid->ny) >= pgrid->gy) {
        pgrid->ny = (pgrid->gy - pgrid->oy);
    }

    pgrid->plumex = -1;
    pgrid->plumey = -1;
    
    if(args->plume) {
        gplumex = (int)((args->plume_source[LONG_IDX] -
                                  args->min_coord[LONG_IDX]) /
                                 args->grid_deltas[LONG_IDX]);
        gplumey =
            (int)((args->plume_source[LAT_IDX] - args->min_coord[LAT_IDX]) /
                    args->grid_deltas[LAT_IDX]);
        if(x_y_within_coords(gplumex, gplumey, pgrid->ox, pgrid->oy, pgrid->nx, pgrid->ny)) {
            pgrid->plumex = gplumex - pgrid->ox;
            pgrid->plumey = gplumey - pgrid->oy;
        } else {
            args->plume = 0;
        }
    }

    // expand for ghosts
    pgrid->data = malloc(sizeof(*pgrid->data) * pgrid->ny+2);
    pgrid->data[0] = malloc(sizeof(*pgrid->data[0]) * (pgrid->nx+2) * (pgrid->ny+2));
    for(i = 0; i < pgrid->ny+2; i++) {
        if(i) {
            pgrid->data[i] = &(pgrid->data[0][i * pgrid->nx]);
        }
        for(j = 0; j < pgrid->nx+2; j++) {
            //pgrid->data[i][j] = args->baseline;
            pgrid->data[i][j] = i * pgrid->nx * i + j;
        }
    }

    for(i = 0; i < 4; i++) {
        if(app->nranks[i] > -1) {
            pgrid->ghost[i] = malloc(sizeof(*pgrid->ghost[i]) *
                                ((i==WEST_RANK||i==EAST_RANK)?pgrid->ny:pgrid->nx));
        } else {
            pgrid->ghost[i] = NULL;
        }
    }

    if(!app->rank) {
        app->bounds = malloc(sizeof(*app->bounds) * app->size);
    }
    bounds[0] = pgrid->ox;
    bounds[1] = pgrid->oy;
    bounds[2] = pgrid->nx;
    bounds[3] = pgrid->ny;
    MPI_Gather(bounds, 4, MPI_INT, app->bounds, 4, MPI_INT, 0, app->comm);

    return EXIT_SUCCESS;
}

void coords_from_rank(int px, int rank, int *x, int *y)
{
    *x = rank % px;
    *y = rank / px;
}

int rank_from_coords(int px, int py, int x, int y)
{
    if(x < 0 || x >= px || y < 0 || y >= py) {
        return(-1);
    }
    return(y * px + x);
}

void init_neighbors(struct sim_app *app)
{
    int nanks[4];

    coords_from_rank(app->px, app->rank, &app->x, &app->y);
    app->nranks[WEST_RANK] = rank_from_coords(app->px, app->y, app->x-1, app->y);
    app->nranks[EAST_RANK] = rank_from_coords(app->px, app->y, app->x+1, app->y);
    app->nranks[SOUTH_RANK] = rank_from_coords(app->px, app->y, app->x, app->y-1);
    app->nranks[NORTH_RANK] = rank_from_coords(app->px, app->y, app->x, app->y+1);
}

int init_sim_app(int argc, char **argv, struct sim_app *app)
{
    if(!app) {
        return EXIT_FAILURE;
    }
    app->comm = MPI_COMM_WORLD;
    MPI_Comm_rank(app->comm, &app->rank);
    MPI_Comm_size(app->comm, &app->size);
    parse_arguments(argc, argv, app);
    if(validate_args(app)) {
        return EXIT_FAILURE;
    }

    init_neighbors(app);

    if(!app->rank) {
        if(parse_conf(app->args.input_file, &app->args, NULL)) {
            app->args.steps = -1; // signal bad conf
        }
    }

    MPI_Bcast(&app->args, sizeof(app->args), MPI_BYTE, 0, app->comm);
    if(app->args.steps == -1) {
        //Config was bad
        return EXIT_FAILURE;
    }

    if(validate_config(app)) {
        return EXIT_FAILURE;
    }

    if(init_sim_grid(app)) {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

void exchange_ghosts(struct sim_app *app)
{
    static double *ghost_buf[4] = {NULL, NULL};
    MPI_Request send_req[4];
    double *out_buf;
    struct sim_pgrid *pgrid = &app->pgrid;
    int xstart, ystart;
    int xmult, ymult;
    int ghost_size;

    int i, j;

    for(i = 0; i < 4; i++) {
        if(app->nranks[i] > -1) {
            switch(i) {
                case WEST_RANK:
                    ghost_size = pgrid->ny;
                    if(ghost_buf[0] == NULL) {
                        ghost_buf[0] = malloc(sizeof(**ghost_buf) * ghost_size);
                    }
                    out_buf = ghost_buf[0];
                    for(j = 0; j < ghost_size; j++) {
                        out_buf[j] = pgrid->data[j+1][1];
                    }
                    break;
                case EAST_RANK:
                    ghost_size = pgrid->ny;
                    if(ghost_buf[1] == NULL) {
                        ghost_buf[1] = malloc(sizeof(**ghost_buf) * ghost_size);
                    }
                    out_buf = ghost_buf[1];
                    for(j = 0; j < ghost_size; j++) {
                        out_buf[j] = pgrid->data[j+1][(pgrid->nx)-2];
                    }
                    break;
                case NORTH_RANK:
                    ghost_size = pgrid->nx;
                    out_buf = &pgrid->data[(pgrid->ny)-2][1];
                    break;
                case SOUTH_RANK:
                    ghost_size = pgrid->nx;
                    out_buf = &pgrid->data[1][1];
                    break;
            }
            MPI_Isend(out_buf, ghost_size, MPI_DOUBLE, app->nranks[i], 0, app->comm, &send_req[i]);
        }
    }

    for(i = 0; i < 4; i++) {
        if(app->nranks[i] > -1) {
            if(app->nranks[i] > -1) {
                ghost_size = (i==EAST_RANK||i==WEST_RANK)?pgrid->ny:pgrid->nx;
                MPI_Recv(pgrid->ghost[i], ghost_size, MPI_DOUBLE, app->nranks[i], 0, app->comm, MPI_STATUS_IGNORE);
            }
            switch(i) {
                case WEST_RANK:
                    ystart = 1;
                    ymult = 1;
                    xstart = 0;
                    xmult = 0;
                    break;
                case EAST_RANK:
                    ystart = 1;
                    ymult = 1;
                    xstart = pgrid->nx-1;
                    xmult = 0;
                    break;
                case NORTH_RANK:
                    ystart = pgrid->ny-1;
                    ymult = 0;
                    xstart = 1;
                    xmult = 1;
                    break;
                case SOUTH_RANK:
                    ystart = 0;
                    ymult = 0;
                    xstart = 1;
                    xmult = 1;
                    break;
            }
            for(j = 0; j < ghost_size; j++) {
                pgrid->data[ystart + j * ymult][xstart + j * xmult] = pgrid->ghost[i][j];
            }
        }
    }
}

/*
    Do a convection-diffusion step with explicit advance (Euler's method)
    The diffusion and convection coefficients are trivial, and numerical
    instability is VERY LIKELY for many values that result in a snappy
    run time on a laptop.

    Return the largest value in the new timestep (artifact of troubleshooting)
*/
double convect_diffuse(struct sim_app *app)
{
    static double **new_data = NULL;
    struct sim_args *args = &app->args;
    struct sim_pgrid *pgrid = &app->pgrid;
    double i_bias, j_bias;
    double w_lat = args->wind[LAT_IDX];
    double w_long = args->wind[LONG_IDX];
    double dx = args->grid_deltas[LONG_IDX];
    double dy = args->grid_deltas[LAT_IDX];
    double **data = pgrid->data;
    double conv_du, diff_du;
    double max = -1;
    double dt = args->dt;
    int nx = pgrid->nx;
    int ny = pgrid->ny;
    int gx = pgrid->gx;
    int gy = pgrid->gy;
    int ox = pgrid->ox;
    int oy = pgrid->oy;
    int i, j;

    // keep reusing the same static buffer
    if(!new_data) {
        new_data = malloc(sizeof(*new_data) * (ny+2));
        new_data[0] = malloc(sizeof(*new_data[0]) * (nx+2) * (ny+2));
        for(i = 1; i < ny+2; i++) {
            new_data[i] = &(new_data[0][i * (nx+2)]);
        }
    }

    // TODO: fiddle with coefficients? Better put in some stability checks at
    // least.
    for(i = 1; i < ny-1; i++) {
        /* the point of i_bias and j_bias is that we can be rational about
           convection at the boundaries towards which the wind is blowing, but
           we can't in the direction the wind is coming from (thar be monsters).
           We bias the indices by the corresponding components of the wind
           vector, and don't try to calculate the convection component if that
           puts us out of bounds.

            Very small values in the wind vector might cause a round-off error
           bug.
        */
        i_bias = (double)(oy + i) - w_long;
        for(j = 1; j < nx-1; j++) {
            // Convection component
            j_bias = (double)(ox + j) - w_lat;
            conv_du = 0;
            if(i_bias > 0 && i_bias < gy - 1 && j_bias > 0 && j_bias < gx - 1) {
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
            if(i > 0 && i < gy - 1 && j > 0 && j < gx - 1) {
                diff_du += (data[i + 1][j] - 2 * data[i][j] + data[i - 1][j]) /
                           (dy * dy);
                diff_du += (data[i][j + 1] - 2 * data[i][j] + data[i][j - 1]) /
                           (dx * dx);
            }

            new_data[i][j] =
                data[i][j] + dt * (args->diffusivity * diff_du - conv_du);
            if(new_data[i][j] > max) {
                max = new_data[i][j];
            }
        }
    }

    // plume source
    if(pgrid->plumex != -1 && pgrid->plumey != -1) {
        new_data[pgrid->plumex+1][pgrid->plumey+1] += args->source;
    }

    //memcpy(data[0], new_data[0], sizeof(*new_data[0]) * nx * ny);
    for(i = 1; i < pgrid->ny-1; i++) {
        //memcpy(&data[i][1], &new_data[i][1], sizeof(*new_data[0]) * pgrid->nx-2);
    }

    return (max);
}

void fwrite_pgrid_data(struct sim_app *app, const char *filename, int step)
{
    static double *gdata = NULL;
    static double *recv_buf = NULL;
    static int *sizes = NULL;
    static int *dspls = NULL;
    static MPI_Datatype send_type;
    struct sim_pgrid *pgrid = &app->pgrid;
    int arrsizes[2], subsizes[2], starts[2];
    int ox, oy, nx, ny;
    double *rdata;
    int gidx, lidx;
    MPI_Request req;
    int i, j, k;

    if(!app->rank && !gdata) {
        gdata = malloc(sizeof(*gdata) * pgrid->gx * pgrid->gy);
        recv_buf = malloc(sizeof(*recv_buf) * pgrid->gx * pgrid->gy);
        dspls = malloc(sizeof(*dspls) * app->size);
        sizes = malloc(sizeof(*sizes) * app->size);
        dspls[0] = 0;
        for(i = 0; i < app->size - 1; i++) {
            sizes[i] = app->bounds[i][2] * app->bounds[i][3];
            dspls[i+1] = dspls[i] + sizes[i];
        }
        sizes[i] = app->bounds[i][2] * app->bounds[i][3];
        arrsizes[1] = pgrid->nx+2;
        arrsizes[0] = pgrid->ny+2;
        subsizes[1] = pgrid->nx;
        subsizes[0] = pgrid->ny;
        starts[0] = 1;
        starts[1] = 1;
        MPI_Type_create_subarray(2, arrsizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &send_type);
        MPI_Type_commit(&send_type);
    } 
    MPI_Isend(pgrid->data[0], 1, send_type, 0, 0, app->comm, &req);
    if(!app->rank) {
        for(i = 0; i < app->size; i++) {
            rdata = &recv_buf[dspls[i]];
            MPI_Recv(rdata, sizes[i], MPI_DOUBLE, i, 0, app->comm, MPI_STATUS_IGNORE);
            ox = app->bounds[i][0];
            oy = app->bounds[i][1];
            nx = app->bounds[i][2];
            ny = app->bounds[i][3];
            for(j = 0; j < ny; j++) {
                for(k = 0; k < nx; k++) {
                    lidx = j * nx + k;
                    gidx = ((oy + j) * pgrid->nx) + (ox + k);
                    gdata[gidx] = rdata[lidx];
                }
            }
        }
    }
    if(!app->rank) {
        fwrite_grid_data(pgrid, gdata, filename, step);
    }
    MPI_Wait(&req, MPI_STATUS_IGNORE);
}

int main(int argc, char *argv[])
{
    struct sim_app app = {0};
    struct sim_args *args = &app.args;
    struct sim_pgrid *pgrid = &app.pgrid;
    char outbase[100], outfile[100];
    double max;
    int t;

    MPI_Init(&argc, &argv);

    if(init_sim_app(argc, argv, &app)) {
        goto err_out;
    }

    //TODO: update for parallel grid
    //printf_grid(grid);

    if(!app.rank) {
        if(args->sim_out_dir && create_dir(args->sim_out_dir)) {
            sprintf(outbase, "%s/out.", args->sim_out_dir);
        } else {
            sprintf(outbase, "out.");
        }
        printf("Simulate %i steps.", args->steps);
        if(args->out_steps) {
            printf("Output every %i steps.\n", args->out_steps);
            if(args->sim_out_dir) {
                printf("Output will be placed in %s/\n", args->sim_out_dir);
            }
        }
        printf("\n");
    }

    for(t = 0; t < args->steps; t++) {
        exchange_ghosts(&app);
        convect_diffuse(&app);
        if(args->out_steps && t % args->out_steps == 0) {
            if(!app.rank) {
                printf("Step %i\n", t);
                sprintf(outfile, "%s%i.dat", outbase, t);
            }
            fwrite_pgrid_data(&app, outfile, t);
        }
    }

    MPI_Finalize();
    return EXIT_SUCCESS;

err_out:
    MPI_Finalize();
    return EXIT_FAILURE;

}
