/* Core functions for running a depth-first search algorithm on a hypercube of arbitrary
dimension for the percolation problem. */

#include "functions.h"
#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>


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
 *  p: the percolation concentration. 0 <= p <= 1
 *  error: a pointer to an error flag in case of probems.
 *
 *  returns: a pointer to an array of NR cluster sizes, of type ul (unsigned long)
 */
ul *clusters_hypercube(ul N, ul NR, float p, gsl_rng *RNG, int *error)
{

    if (!check_args(N, NR, p)) {*error = 3; return NULL;}

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
