#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h> 
#include <stdbool.h>
typedef unsigned long ul;
typedef struct stack
    {
        ul *sites;
        ul length;
        ul top;
        ul NH;
    }
stack;

/* Stack functions */
int push(stack *s, ul site);
ul pop(stack *s, int *error);
stack *setup_stack(ul NH);
void reset_visited(bool visited[], ul length);

/* Hypercube functions */
ul DFS_hypercube(stack *s, bool visited[], float p, ul N, ul start_state, gsl_rng *RNG, int *error);
ul *clusters_hypercube(ul N, ul NR, float p, gsl_rng *RNG, int *error);

/* PXP functions */
void populate_sites_PXP(ul *sites, ul N);
ul *clusters_PXP(ul N, ul NR, float p, gsl_rng *RNG, int *error);
bool PXP_flip_allowed(ul u, ul i, ul N);
ul DFS_PXP(stack *s, ul *sites, bool visited[], float p, ul N, ul start_state, gsl_rng *RNG, int *error);

/* Maths functions */
ul intpower(ul base, ul exponent);
ul binomialCoeff(ul n, ul r);
ul index_site(ul *sites, ul site, ul left, ul right, int *idx_flag);
ul fibonacci(ul n);

/* Misc functions */
bool check_args(ul N, ul NR, float p);



