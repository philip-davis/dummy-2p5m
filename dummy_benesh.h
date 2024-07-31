#ifndef __DUMMY_BENESH__
#define __DUMMY_BENESH__

#include<mpi.h>
#include<stdint.h>
#include<stdlib.h>

typedef struct benesh_handle *benesh_app_id;
typedef void *benesh_arg;

typedef int (*benesh_method)(benesh_app_id, benesh_arg);

int benesh_init(const char *name, const char *conf, MPI_Comm gcomm, int dummy, int wait,
                struct benesh_handle **bnh);

int benesh_bind_method(struct benesh_handle *bnh, const char *name,
                       benesh_method method, void *user_arg);

int benesh_bind_var(struct benesh_handle *bnh, const char *var_name, void *buf);

void *benesh_bind_var_mesh(struct benesh_handle *bnh, const char *var_name, int *idx, unsigned int idx_len);

void *benesh_bind_field_mpient(struct benesh_handle *bnh, const char *var_name, int idx, const char *rcn_file, MPI_Comm comm, void *buffer, int length, int participates);

void *benesh_bind_field_dummy(struct benesh_handle *bnh, const char *var_name, int idx, int participates);

typedef enum bnh_order {
    BNH_ROW_MAJOR,
    BNH_COL_MAJOR
} bnh_order_t;

typedef enum bnh_storage_type {
    BNH_GRID,
    BNH_MESH,
    BNH_SET
} bnh_storage_type;

typedef struct bnh_storage_def {
    bnh_storage_type type;
    union {
        struct {
            enum bnh_order order;
            int ndim;
            int *dims;
            int *ghosts;
        };
    };
} *bnh_storage_def_t;

bnh_storage_def_t benesh_grid_def(const int ndim, int *dims, bnh_order_t order);

void benesh_grid_ghosts(bnh_storage_def_t gdef, int *depths);

void benesh_grid_ghosts_uniform(bnh_storage_def_t gdef, int depth);

typedef enum benesh_datatype {
    BNH_INT,
    BNH_LONG,
    BNH_FLOAT,
    BNH_DOUBLE,
    BNH_BOOL,
    BNH_BYTES
} benesh_datatype;

enum benesh_fragment_type {
    BNH_TILE,
    BNH_DISC_TILE,
    BNH_GEO_TILE,
    BNH_SPANS
};

struct benesh_tile_fragment {
    int ndim;
    double *lb;
    double *ub;
};

struct benesh_discrete_fragment {
    int ndim;
    uint64_t *lb;
    uint64_t *ub;
};

struct benesh_geotile_fragment {
    double lb[2];
    double ub[2];
};

struct benensh_spans_fragment {
    int count;
    uint64_t *starts;
    uint64_t *spans;
};

typedef struct benesh_domain_fragment {
     enum benesh_fragment_type type;
     union {
       struct benesh_tile_fragment tile;
       struct benesh_discrete_fragment dtile;
       struct benesh_geotile_fragment gtile;
       struct benensh_spans_fragment spans;
    };
} *benesh_domain_fragment_t;

benesh_domain_fragment_t benesh_domain_geotile_decompose(benesh_app_id bnh, const char *domain, double lb[2], double ub[2]);

int benesh_bind_var_with_size(benesh_app_id bnh, const char *var_name, bnh_storage_def_t sdef, benesh_domain_fragment_t dfrag, size_t dsize);

int benesh_bind_var_with_type(benesh_app_id bnh, const char *var_name, bnh_storage_def_t sdef, benesh_domain_fragment_t dfrag, benesh_datatype dtype);

int benesh_bind_grid_domain(benesh_app_id bnh, const char *dom_name,
                            double *grid_offset, double *grid_dims,
                            uint64_t *grid_points, int alloc);

int benesh_bind_mesh_domain(benesh_app_id bnh, const char *dom_name,
                            const char *grid_file, const char *cpn_file, int alloc);

void benesh_touchpoint(benesh_app_id bnh, const char *tpname);

int benesh_fini(benesh_app_id bnh);

int benesh_get_var_domain(benesh_app_id bnh, const char *var_name,
                          char **dom_name, int *ndim, double **lb, double **ub);

void *benesh_get_var_buf(benesh_app_id bnh, const char *var_name, uint64_t *size);

double benesh_get_var_val(benesh_app_id bnh, const char *var_name);

typedef long benesh_int_t;
typedef double benesh_real_t;

#endif //__DUMMY_BENESH__