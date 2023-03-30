typedef unsigned long ul;
typedef struct stack
    {
        ul *sites;
        ul length;
        ul top;
        ul NH;
    }
stack;

int push(stack *s, ul site);
ul pop(stack *s, int *error);
ul intpower(ul base, ul exponent);
ul binomialCoeff(ul n, ul k);
ul index_site(ul *sites, ul site, ul left, ul right, int *idx_flag);
void reset_visited(bool visited[], ul length);
void test_stack(stack *s);
void populate_sites_XXZ(ul *sites, ul N, int UP);
