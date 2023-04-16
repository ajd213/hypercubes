typedef unsigned long ul;
typedef struct stack
    {
        ul *sites;
        ul length;
        ul top;
        ul NH;
    }
stack;

// an instance of the GNU scientific library's random number generator
extern gsl_rng *RNG;

int push(stack *s, ul site);
ul pop(stack *s, int *error);
ul intpower(ul base, ul exponent);
ul binomialCoeff(ul n, ul r);
ul index_site(ul *sites, ul site, ul left, ul right, int *idx_flag);
void reset_visited(bool visited[], ul length);
void populate_sites_PXP(ul *sites, ul N);
ul DFS_hypercube(stack *s, bool visited[], float p, ul N, ul start_state, gsl_rng *RNG, int *error);
ul *clusters_hypercube(ul N, ul NR, float p, int *error);
bool check_args(ul N, ul NR, float p);
ul fibonacci(ul n);
ul *clusters_PXP(ul N, ul NR, float p, int *error);
bool PXP_flip_allowed(ul u, ul i, ul N);
ul DFS_PXP(stack *s, ul *sites, bool visited[], float p, ul N, ul start_state, gsl_rng *RNG, int *error);
stack *setup_stack(ul NH);
