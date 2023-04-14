#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include <gsl/gsl_rng.h> 
#include "hypercube_functions.h"

int main(int argc, char **argv)
{   

    if (argc != 5)
    {
        printf("Usage: ./main N NR p filename\n");
        return 4;
    }

    ul N = (ul) atoi(argv[1]); // chain length
    if (N > 32)
    {
        printf("Error! Max chain length = 32.\n");
        return 4;
    }

    int NR = (ul) atoi(argv[2]); // number of realisations
    float p = (float) atof(argv[3]); // percolation concentration
    char *filename = argv[4]; // where to save data

    FILE *outfile = fopen(filename, "w");
    if (!outfile)
    {
        printf("Error opening file: %s!\n", filename);
        return 3;
    }


    int errorflag = 0;
    // get the clusters
    ul *cs = clusters_PXP(N, NR, p, &errorflag);
    if (errorflag != 0 || !cs)
    {
        return 4;
    }

    for (int i = 0; i < NR; i++)
    {
        fprintf(outfile, "%lu\n", cs[i]);
    }

    free(cs);

}
