/* Use a depth-first search algorithm to explore the percolation problem on a hypercube.
The user specifies the dimension of the hypercube, the percolation concentration, and the
number of times to run the DFS algorithm (I.E., the number of clusters to grow).
The sizes of the clusters are then saved to a file, which can be analysed to extract
the percolation transition. */

#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include <gsl/gsl_rng.h> 
#include "functions.h"

gsl_rng *RNG; // random number generator

/*
 * Function:  main
 * --------------------
 *  Grows many clusters on the hypercube graph and saves their sizes to a file.
 * 
 * ARGS:
 * 
 *  N: (long int) the dimension of the hypercube. (E.g. N=3 is a regular cube.)
 *  NR: (long int) the number of disorder realisations to use. = the number of clusters to grow.
 *  p: (float) the percolation strength. 0 <= p <= 1
 *  filename: where to save the cluster sizes
 *
 *  returns: 0 if successful, otherwise a nonzero integer.
 */
int main(int argc, char **argv)
{

    if (argc != 5)
    {
        printf("Usage: ./main N NR p <filename>\n");
        return 1;
    }

    // set and seed the RNG
    RNG = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RNG, time(NULL));

    
    ul N = (ul) atoi(argv[1]); // dimension of hypercube
    int NR = atoi(argv[2]); // number of realisations
    float p = atof(argv[3]); // percolation concentration

    char *filename = argv[4];

    int error_flag = 0;

    ul *cluster_sizes = clusters_hypercube(N, NR, p, &error_flag);
    if (error_flag != 0 || !cluster_sizes)
    {
        // error! Can add handling of different errors via error_flag
        return 2;
    }

    FILE *to_save = fopen(filename, "w");
    if (!to_save)
    {
        printf("Error opening file: %s!\n", filename);
        free(cluster_sizes);
        return 3;
    }

    for (int i = 0; i < NR; i++)
    {
        fprintf(to_save, "%lu\n", cluster_sizes[i]);
    }

    free(cluster_sizes);
    fclose(to_save);
    return 0;
}
