#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include <gsl/gsl_rng.h> 
#include "functions.h"

gsl_rng *RNG; // random number generator

int main(int argc, char **argv)
{

    if (argc != 2)
    {
        printf("Usage: ./main <filename>\n");
        return 1;
    }

    // set and seed the RNG
    RNG = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RNG, time(NULL));

    
    ul N = 16; // dimension of hypercube
    int NR = 100; // number of realisations
    float p = 1; // percolation concentration

    char *filename = argv[1];

    int error_flag = 0;

    ul *cluster_sizes = clusters_hypercube(N, NR, p, &error_flag);
    if (error_flag != 0)
    {
        //error!
        return 2;
    }

    FILE *to_save = fopen(filename, "w");
    if (!to_save)
    {
        printf("Error opening file: %s!\n", filename);
        return 3;
    }

    for (int i = 0; i < NR; i++)
    {
        fprintf(to_save, "%lu\n", cluster_sizes[i]);
    }

    free(cluster_sizes);
    fclose(to_save);
}
