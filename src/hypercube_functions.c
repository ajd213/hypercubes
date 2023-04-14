/* Core functions for running a depth-first search algorithm on a hypercube of arbitrary
dinemsion for the percolation problem. */

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <gsl/gsl_rng.h> 
#include "hypercube_functions.h"



/*
 * Function:  DFS_PXP
 * --------------------
 * grows a percolation cluster on a Fibonacci cube (the HS of the PXP model) using a depth-
 * first search.
 *
 *  s: a pointer to a stack, defined in functions.h
 *  sites: a pointer to an array of sites, which contains the Fibonacci cube nodes
 *  visited: an array of bools of length NH, to record whether each site has been visited
 *  p: the percolation strength. 0 <= p <= 1
 *  N: the dimension of the hypercube. (E.g. N=3 is a regular cube.)
 *  start_state: which site on the hypercube to grow the cluster from
 *  RNG: a random number generator from the gsl library
 *  error: a pointer to an error flag, in case something goes wrong
 *
 *  returns: the size of the cluster. I.e., the number of sites visited by the DFS algorithm.
 */
ul DFS_PXP(stack *s, ul *sites, bool visited[], float p, ul N, ul start_state, gsl_rng *RNG, int *error)
{
    ul size = 0; // cluster size
    ul NH = fibonacci(N+2);
    ul u, v, idx_u, idx_v;
    int idx_flag = 0;

    if (s->top != 0)
    {
        printf("Error! Stack not empty!\n");
        *error = -1; 
        return 0;
    }
    push(s, start_state);

    while (s->top > 0)
    {

        u = pop(s, error);
        if (*error == -1)
        {
            // error flag has been set, so return.
            return 0;
        }

        idx_u = index_site(sites, u, 0, NH-1, &idx_flag);
        if (idx_flag == -1)
        {
            printf("Error! Site u not found!\n");
            *error = -1;
            return 0;
        }

        if (visited[idx_u])
        {
            continue;
        }
        visited[idx_u] = true;
        size++;

        // loop from right to left over the spins
        for (ul i = 0; i < N; i++)
        {

            if (PXP_flip_allowed(u, i, N))
            {
                // do the flip
                v = u ^ (1UL << i);

                // find the index
                idx_v = index_site(sites, v, 0, NH-1, &idx_flag);
                if (idx_flag == -1)
                {
                    printf("Error! Site v not found!\n");
                    *error = -1;
                    return 0;
                }

                // and push to the stack
                if (!visited[idx_v] && (gsl_rng_uniform(RNG) < p))
                {
                    if (push(s, v) == 1) 
                    { 
                        // stack error!
                        *error = -1; 
                        return 0;
                    }
                }
            }

        }
    }
    return size;
}

/*
* Function: PXP_flip_allowed
* ---------------------------
* For a node in the PXP model (Fibonacci cube), is flipping the spin (bit) at index i
* a valid flip? I.e., do the PXP constraints allow it? I.e., does it lead to another
* node within the Fibonacci cube?
*
* u : the state
* i : the index of the bit in u to be flipped
* N : the cube dimension == the length of the chain
*/
bool PXP_flip_allowed(ul u, ul i, ul N)
{
    ul left_bit, right_bit;
    bool flip_allowed = false;

    // if the bit is the first site
    if (i == 0)
    {
        // the bit to the left of bit i must be zero
        left_bit = u & (1UL << 1);
        if (left_bit == 0UL)
        {
            flip_allowed = true;
        }
    }
    // if the bit is the last site
    else if (i == N - 1)
    {
        // the bit to the right of bit i must be zero
        right_bit = u & (1UL << (N - 2));
        if (right_bit == 0UL)
        {
            flip_allowed = true;
        }
    }
    // if the bit is in the middle 
    else
    {
        // both adjacent bits must be zero
        left_bit = u & (1UL << (i + 1));
        right_bit = u & (1UL << (i - 1));
        if ((left_bit == 0UL) && (right_bit == 0UL))
        {
            flip_allowed = true;
        }
    }
    return flip_allowed;
}


/*
 * Function:  DFS_hypercube
 * --------------------
 * grows a percolation cluster on a hypercube with a depth-first search:
 *
 *  s: a pointer to a stack, defined in functions.h
 *  visited: an array of bools of length NH, to record whether each site has been visited
 *  p: the percolation strength. 0 <= p <= 1
 *  N: the dimension of the hypercube. (E.g. N=3 is a regular cube.)
 *  start_state: which site on the hypercube to grow the cluster from
 *  RNG: a random number generator from the gsl library
 *  error: a pointer to an error flag, in case something goes wrong
 *
 *  returns: the size of the cluster. I.e., the number of sites visited by the DFS algorithm.
 */
ul DFS_hypercube(stack *s, bool visited[], float p, ul N, ul start_state, gsl_rng *RNG, int *error)
{
    ul size = 0; // cluster size
    ul u, v;

    if (s->top != 0)
    {
        printf("Error! Stack not empty!\n");
        *error = -1; 
        return 0;
    }
    push(s, start_state);

    while (s->top > 0)
    {

        u = pop(s, error);
        if (*error == -1)
        {
            // error flag has been set, so return.
            return 0;
        }

        if (visited[u])
        {
            continue;
        }
        visited[u] = true;
        size++;

        for (ul i = 0; i < N; i++)
        {
            // flip the ith bit
            v = u ^ (1UL << i);

            if (!visited[v] && (gsl_rng_uniform(RNG) < p))
            {
                if (push(s, v) == 1) 
                { 
                    // stack error!
                    *error = -1; 
                    return 0;
                }
            }
        }
    }
    return size;
}

/*
 * Function:  check_args
 * --------------------
 *  check user arguments N, NR and p.
 *
 *  N: the dimension of the hypercube: 1 <= N <= 32
 *  NR: the Number of Realisations: number of clusters to grow. NR >= 1
 *  p: the percolation concentration. 0 <= p <= 1
 *
 *  returns: true if arguments are OK, else false.
 */
bool check_args(ul N, ul NR, float p)
{
    if (N < 1 || N > 32)
    {
        printf("N must be between 1 and 32.\n");
        return false;
    }

    if (p < 0 || p > 1)
    {
        printf("p is a probability, and must satisfy 0 <= p <= 1.\n");
        return false;
    }

    if (NR < 1)
    {
        printf("Please use one or more realisations: NR > 1.\n");
        return false;
    }
    return true;
}

stack *setup_stack(NH)
{
    stack *s = malloc(sizeof(stack));
    if (!s)
    {
        return NULL;
    }
    s->top = 0;
    s->length = NH;
    s->NH = NH;
    s->sites = malloc(s->length*sizeof(ul));
    if (!s->sites)
    {
        return NULL;
    }

    return s;
}

ul *clusters_PXP(ul N, ul NR, float p, int *error)
{

    if (!check_args(N, NR, p)) {*error = 3; return NULL;}

    // set and seed the random number generator
    gsl_rng *RNG = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RNG, time(NULL));

    // the size of the graph
    ul NH = fibonacci(N+2);

    stack *s = setup_stack(NH);
    if (!s)
    {
        printf("Error using malloc for s.sites!\n");
        *error = 1;
        return NULL;
    }

    // a list of the sites
    ul *sitelist = malloc(sizeof(ul)*NH);
    if (sitelist == NULL)
    {
        printf("Error using malloc for sitelist!\n");
        free(s->sites);
        free(s);
        *error = 1;
        return NULL;
    }
    populate_sites_PXP(sitelist, N);

    // keep track of visited nodes
    bool *visited = malloc(sizeof(bool)*NH);
    if (!visited)
    {
        printf("Error using malloc for visited!\n");
        free(s->sites);
        free(s);
        *error = 1;
        return NULL;
    }

    int error_flag = 0;

    ul *cluster_sizes = malloc(NR*sizeof(ul));
    if (!cluster_sizes)
    {
        printf("Error using malloc for cluster_sizes!\n");
        *error = 1;
        free(visited);
        free(s->sites);
        free(s);
        return NULL;
    }

    // where we will grow the cluster from
    ul start_site;

    // run the DFS algorithm over NR realisations
    for (ul i = 0; i < NR; i++)
    {
        // set all nodes to not visited
        reset_visited(visited, NH); 
        start_site = sitelist[gsl_rng_uniform_int(RNG, NH)];
        // run DFS algorithm, get a cluster size
        cluster_sizes[i] = DFS_PXP(s, sitelist, visited, p, N, start_site, RNG, &error_flag);
        if (error_flag == -1)
        {
            // error!
            printf("Error! Exiting DFS algorithm and closing file...\n");
            free(s->sites);
            free(s);
            free(visited);
            free(cluster_sizes);
            *error = 2;
            return NULL;
        }
    }

    // free heap memory, except cluster sizes
    free(s->sites);
    free(s);
    free(visited);
    return cluster_sizes;
}

/*
 * Function:  clusters_hypercube
 * --------------------
 *  driver code for running DFS_hypercube many times and returning a pointer to the cluster sizes.
 *  Checks the user input.
 *
 *  N: the dimension of the hypercube
 *  NR: the Number of Realisations: number of clusters to grow
 *  p: the percolation concentration. 0 <= p <= 1
 *  error: a pointer to an error flag in case of probems.
 *
 *  returns: a pointer to an array of NR cluster sizes, of type ul (unsigned long)
 */
ul *clusters_hypercube(ul N, ul NR, float p, int *error)
{

    if (!check_args(N, NR, p)) {*error = 3; return NULL;}

    // set and seed the random number generator
    gsl_rng *RNG = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RNG, time(NULL));

    // the size of the graph
    ul NH = intpower(2, N); 

    stack *s = setup_stack(NH);
    if (!s)
    {
        printf("Error using malloc for s.sites!\n");
        *error = 1;
        return NULL;
    }

    // keep track of visited nodes
    bool *visited = malloc(sizeof(bool)*NH);
    if (!visited)
    {
        printf("Error using malloc for visited!\n");
        free(s->sites);
        free(s);
        *error = 1;
        return NULL;
    }

    // where to grow the cluster
    ul start_site = 0; 
    int error_flag = 0;

    ul *cluster_sizes = malloc(NR*sizeof(ul));
    if (!cluster_sizes)
    {
        printf("Error using malloc for cluster_sizes!\n");
        *error = 1;
        free(visited);
        free(s->sites);
        free(s);
        return NULL;
    }

    // run the DFS algorithm over NR realisations
    for (ul i = 0; i < NR; i++)
    {
        // set all nodes to not visited
        reset_visited(visited, NH); 

        // run DFS algorithm, get a cluster size
        cluster_sizes[i] = DFS_hypercube(s, visited, p, N, start_site, RNG, &error_flag);
        if (error_flag == -1)
        {
            // error!
            printf("Error! Exiting DFS algorithm and closing file...\n");
            free(s->sites);
            free(s);
            free(visited);
            free(cluster_sizes);
            *error = 2;
            return NULL;
        }
    }

    // free heap memory, except cluster sizes
    free(s->sites);
    free(s);
    free(visited);
    return cluster_sizes;
}

/*
 * Function:  push
 * --------------------
 * push a site to the stack, and re-size the stack with realloc() if necessary.
 * Stack struct defined in header file.
 *
 *  s: a pointer to a stack, defined in functions.h
 *  site: the graph site to be pushed to the stack
 *
 *  returns: 0 if push successful, 1 otherwise.
 */
int push(stack *s, ul site)
{
    if (s->top == s->length - 1)
    {
        // Overflow! Increase the length of the stack.
        ul *new_stackmem = realloc(s->sites, (s->length + s->NH)*sizeof(ul));
        if (new_stackmem == NULL)
        {
            printf("Error with realloc in stack!\n");
            return 1;
        }
        else
        {
            s->sites = new_stackmem;
            s->length += s->NH;
        }
    }
    s->sites[s->top] = site;
    s->top++;

    // pushed item
    return 0;
}

/*
 * Function:  pop
 * --------------------
 * pop a site from the stack. Throw an error if the stack is empty.
 *
 *  s: a pointer to a stack, defined in functions.h
 *  error: a pointer to an error flag
 *
 *  returns: the site popped from the stack, if successful. Otherwise 0 (and error raised).
 */
ul pop(stack *s, int *error)
{
    if (s->top == 0)
    {
        // cannot pop if the stack is empty
        *error = -1;
        return 0;
    }

    s->top--;
    return s->sites[s->top];
}

/*
 * Function:  reset_visited
 * --------------------
 * loop through visited[], which tracks whether each site has been visited by the DFS algorithm.
 *
 *  visited: a list of bools equal to NH, the number of nodes in the graph.
 *  length: the length of visited (= NH)
 */
void reset_visited(bool visited[], ul length)
{
    for (ul i = 0; i < length; i++)
    {
        visited[i] = false;
    }
}

/*
 * Function:  index_site
 * --------------------
 * search an ordered list of sites for a particular site using a binary search.
 *
 *  sites: a pointer to an ordered list of sites
 *  site: the site to search for
 *  left: the left index of the list (use 0 when calling)
 *  right: the right index of the list (use len(sites) - 1 when calling)
 *  idx_flag: a pointer to a flag which is set to -1 if the site cannot be found
 *
 *  returns: the index of the site, if found. If not found, return 0 and set idx_flag to -1.
 */
ul index_site(ul *sites, ul site, ul left, ul right, int *idx_flag)
{
    if (right >= left) {
        int mid = left + (right - left) / 2;
 
        if (sites[mid] == site)
            return mid;
 
        if (sites[mid] > site)
            return index_site(sites, site, left, mid - 1, idx_flag);
 
        return index_site(sites, site, mid + 1, right, idx_flag);
    }
    // not found
    *idx_flag = -1;
    return 0;
}

/*
 * Function:  populate_sites_PXP
 * --------------------
 * For the PXP model, populate an ordered list of the graph nodes. The largest sector of the HS is called the Fibonacci cube, and consists
 * of all nodes which have no adjacent set bits.
 *
 *  sites: a pointer to an array of sites (unsigned long integers)
 *  N: the dimension of the graph

 */
void populate_sites_PXP(ul *sites, ul N)
{
    ul full_HS = intpower(2, N);
    for (ul i = 0, counter = 0; i < full_HS; i++)
    {
        if ((i & (i >> 1)) == 0UL)
        {
            sites[counter] = i;
            counter++;
        }
    }
}

/*
 * Function:  intpower
 * --------------------
 * recursive integer exponentiation for ul (unsigned long) type. Used seldomly.
 *
 *  base: the number to be exponentiated
 *  exponent: the power which the base is raised to
 *
 *  returns: base ** exponent (base to the power of exponent)
 */
ul intpower(ul base, ul exponent)
{
    if (exponent == 1)
    {
        return base;
    }
    else
    {
        return base * intpower(base, exponent - 1);
    }
}

/*
 * Function:  binomialCoeff
 * --------------------
 * Efficiently compute binomial Coefficient C(n, r). Note that 
 * we must evaluate the k! term `backwards' to retain exact divisibiity.
 * O(r) time complexity and O(1) space complexity.
 *
 *  n: integer in C(n, r)
 *  r: integer in C(n, r)
 *
 *  returns: C(n, r) = n!/(r!(n-r)!)
 */
ul binomialCoeff(ul n, ul r)
{
    ul result = 1;
 
    // Shorten the loop: C(n, r) = C(n, n-r)
    if (r > n - r)
        r = n - r;

    for (ul i = 0; i < r; i++) {
        result *= (n - i); // n, (n-1), ..., (n-r+1)
        result /= (i + 1); // 1, 2, ..., r
    }
 
    return result;
}

/*
* Function: fibonacci

*/
ul fibonacci(ul n)
{
    assert(n >= 0);

    if (n == 0)
    {
        return 0;
    }

    ul a = 0;
    ul b = 1;
    while (n-- > 1) {
        int t = a;
        a = b;
        b += t;
    }
    return b;
}
