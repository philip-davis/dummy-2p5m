#include "dummy_benesh.h"

#include<stdio.h>
#include<stdlib.h>

int benesh_init(const char *name, const char *conf, MPI_Comm gcomm, int dummy, int wait,
                struct benesh_handle **bnh)
{
    return(0);
}

int benesh_fini(benesh_app_id bnh)
{
    return(0);
}

void benesh_touchpoint(benesh_app_id bnh, const char *tpname)
{
    return;
}

struct bnh_storage_def *benesh_grid_def(const int ndim, int *dims, bnh_order_t order)
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

struct benesh_domain_fragment *benesh_domain_geotile_decompose(struct benesh_handle *bnh, const char *domain, double lb[2], double ub[2])
{
    struct benesh_domain_fragment *dfrag;
    // Validate against configuration
    // Precompute offsets, etc

    dfrag = malloc(sizeof(*dfrag));
    dfrag->type = BNH_GEO_TILE; 
    dfrag->gtile.lb[0] = lb[0];
    dfrag->gtile.lb[1] = lb[1];
    dfrag->gtile.ub[0] = ub[0];
    dfrag->gtile.ub[1] = ub[1];

    return(dfrag);
}

int benesh_bind_var_with_type(struct benesh_handle *bnh, const char *var_name, struct bnh_storage_def *sdef, struct benesh_domain_fragment *dfrag, benesh_datatype dtype)
{
    return 0;
}

int benesh_bind_var_with_size(struct benesh_handle *bnh, const char *var_name, struct bnh_storage_def *sdef, struct benesh_domain_fragment *dfrag, size_t dsize)
{
    //call to benesh_bind_var_with_size()
    return(0);
}