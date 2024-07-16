#include "dummy_benesh.h"

#include<stdio.h>
#include<stdlib.h>

int benesh_init(const char *name, const char *conf, MPI_Comm gcomm, int dummy, int wait,
                struct benesh_handle **bnh)
{
    return(0);
}

struct bnh_storage_def *benesh_new_grid_def(const int ndim, int *dims, bnh_order_t order)
{
    struct bnh_storage_def *gdef = malloc(sizeof(*gdef));
    int i;

    gdef->type = BNH_GRID;
    gdef->order = order;
    gdef->ndim = ndim;
    gdef->dims = malloc(sizeof(*gdef->dims) * ndim);
    gdef->ghosts = malloc(sizeof(*gdef->ghosts) * ndim * 2);
    for(i = 0; i < ndim; i++) {
        gdef->dims[i] = dims[i];
    }

    return(gdef);
}

void benesh_grid_ghosts(struct bnh_storage_def *gdef, int *depths)
{
    int i;

    if(gdef->type != BNH_GRID) {
        fprintf(stderr, "WARNING: %s: trying to set ghosts for a non-grid definition.\n", __func__);
        return;
    }

    for(i = 0; i < gdef->ndim * 2; i++) {
        gdef->ghosts[i] = depths[i];
    }

}

void benesh_grid_ghosts_uniform(struct bnh_storage_def *gdef, int depth)
{
    int depths[gdef->ndim * 2];
    int i;

    for(i = 0; i < gdef->ndim * 2; i++) {
        depths[i] = depth;
    }

    benesh_grid_ghosts(gdef, depths);
}