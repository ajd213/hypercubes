/* Core algorithms for growing percolation clusters on the PXP model graph, also called
the Fibonacci cube. */
#define NO_IMPORT_ARRAY
#include "functions.h"
#include <stdlib.h>
#include <stdio.h>



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
ul DFS_PXP(stack *s, ul *sites, bool visited[], const float p, const ul N, const ul start_state, gsl_rng *RNG, int *error)
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
