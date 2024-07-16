#ifndef __DUMMY_BENESH__
#define __DUMMY_BENESH__

#include<mpi.h>
#include<stdint.h>

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

bnh_storage_def_t benesh_new_grid_def(const int ndim, int *dims, bnh_order_t order);

void benesh_grid_ghosts(struct bnh_storage_def *gdef, int *depths);

void benesh_grid_ghosts_uniform(struct bnh_storage_def *gdef, int depth);

typedef enum benesh_datatype {
    BNH_INT,
    BNH_LONG,
    BNH_FLOAT,
    BNH_DOUBLE,
    BNH_BOOL,
    BNH_BYTES
} benesh_datatype;

int benesh_bind_var_with_type(struct benesh_handle *bnh, const char *var_name, 

int benesh_bind_grid_domain(struct benesh_handle *bnh, const char *dom_name,
                            double *grid_offset, double *grid_dims,
                            uint64_t *grid_points, int alloc);

int benesh_bind_mesh_domain(struct benesh_handle *bnh, const char *dom_name,
                            const char *grid_file, const char *cpn_file, int alloc);

void benesh_tpoint(struct benesh_handle *bnh, const char *tpname);

int benesh_fini(struct benesh_handle *bnh);

int benesh_get_var_domain(struct benesh_handle *bnh, const char *var_name,
                          char **dom_name, int *ndim, double **lb, double **ub);

void *benesh_get_var_buf(struct benesh_handle *bnh, const char *var_name, uint64_t *size);

double benesh_get_var_val(struct benesh_handle *bnh, const char *var_name);

typedef long benesh_int_t;
typedef double benesh_real_t;

#endif //__DUMMY_BENESH__