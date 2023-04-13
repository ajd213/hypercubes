/* Core functions for running a depth-first search algorithm on a hypercube of arbitrary
dinemsion for the percolation problem. */

#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>
#include <gsl/gsl_rng.h> 
#include "hypercube_functions.h"


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
 * Function:  clusters_hypercube
 * --------------------
 *  driver code for running DFS_hypercube many times and returning a pointer to the cluster sizes.
 *  Checks the user input.
 *
 *  N: the dimension of the hypercube
 *  NR: the Number of Realisations: number of clusters to grow
 *  p: the percolation concentration. 0 <= 0 <= 1
 *  error: a pointer to an error flag in case of probems.
 *
 *  returns: a pointer to an array of NR cluster sizes, of type ul (unsigned long)
 */
ul *clusters_hypercube(ul N, ul NR, float p, int *error)
{

    if (N < 1 || N > 32)
    {
        printf("N must be between 1 and 32.\n");
        *error = 3;
        return NULL;
    }

    if (p < 0 || p > 1)
    {
        printf("p is a probability, and must satisfy 0 <= p <= 1.\n");
        *error = 3;
        return NULL;
    }

    if (NR < 1)
    {
        printf("Please use one or more realisations: NR > 1.\n");
        *error = 3;
        return NULL;
    }


    // set and seed the random number generator
    gsl_rng *RNG = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RNG, time(NULL));

    // the size of the graph. Equal to the Hilbert-space dimension.
    ul NH = intpower(2, N); 

    // setup stack. Size will be increased dynamically if needed.
    stack s;
    s.top = 0;
    s.length = NH;
    s.NH = NH;
    s.sites = malloc(s.length*sizeof(ul));
    if (!s.sites)
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
        free(s.sites);
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
        free(s.sites);
        return NULL;
    }

    // run the DFS algorithm over NR realisations
    for (ul i = 0; i < NR; i++)
    {
        // set all nodes to not visited
        reset_visited(visited, NH); 

        // run DFS algorithm, get a cluster size
        cluster_sizes[i] = DFS_hypercube(&s, visited, p, N, start_site, RNG, &error_flag);
        if (error_flag == -1)
        {
            // error!
            printf("Error! Exiting DFS algorithm and closing file...\n");
            free(s.sites);
            free(visited);
            free(cluster_sizes);
            *error = 2;
            return NULL;
        }
    }

    // free heap memory, except cluster sizes
    free(s.sites);
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
 * Function:  populate_sites_XXZ
 * --------------------
 * For the XXZ graph, populate a list with the ordered nodes of the graph. That is, the numbers from 0 to 2**N - 1 which have a specified number of bits set.
 *
 *  sites: a pointer to an array of sites (unsigned long integers)
 *  N: the dimension of the graph
 *  UP: the number of bits we require each site to have set. Alternatively, the (conserved) number of UP spins in the spin chain which has the the XXZ graph as its Hilbert space.

 */
void populate_sites_XXZ(ul *sites, ul N, int UP)
{
    ul full_HS = intpower(2, N);
    for (ul i = 0, counter = 0; i < full_HS; i++)
    {
        if (__builtin_popcount(i) == UP)
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
 * binomial Coefficient C(n, k). Inefficient but simple algorithm: used seldomly.
 *
 *  n: integer in C(n, k)
 *  k: integer in C(n, k)
 *
 *  returns: C(n, k) = n!/(k!(n-k)!)
 */
ul binomialCoeff(ul n, ul k)
{
    // Base Cases
    if (k > n)
        return 0;

    if (k == 0 || k == n)
        return 1;
 
    // Pascal's triangle
    return binomialCoeff(n - 1, k - 1) + binomialCoeff(n - 1, k);
}

