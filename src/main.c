#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include <gsl/gsl_rng.h> 
#include "functions.h"

int main(int argc, char **argv)
{   

    if (argc != 6)
    {
        printf("Usage: ./main h/p N NR p filename\n");
        return 1;
    }

    char arg = argv[1][0];
    ul N = (ul) atoi(argv[2]); // chain length
    int NR = atoi(argv[3]); // number of realisations
    float p = (float) atof(argv[4]); // percolation concentration
    char *filename = argv[5]; // where to save data

    
    // execute code for either hypercubes or PXP

    int errorflag = 0;
    ul *cs;
    switch (arg) 
    {
        case 'h':
            cs = clusters_hypercube(N, NR, p, &errorflag);
            break;
        case 'p':
            cs = clusters_PXP(N, NR, p, &errorflag);
            break;
        default:
            printf("Error: First argument must be 'h' or 'p'\n");
            printf("h: hypercubes, p: PXP/Fibonacci cube\n");
            return 1;
    }
    if (errorflag != 0 || !cs)
    {
        return 2;
    }

    FILE *outfile = fopen(filename, "w");
    if (!outfile)
    {
        printf("Error opening file: %s!\n", filename);
        return 3;
    }
    
    for (int i = 0; i < NR; i++)
    {
        fprintf(outfile, "%lu\n", cs[i]);
    }
    free(cs);
    fclose(outfile);
    return 0;
}
