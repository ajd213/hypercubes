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

    // dimension of hypercube
    ul N = 10; 
    // HS dim
    ul NH = intpower(2, N); 
    // number of realisations
    int NR = 100;
    // percolation concentration
    float p = 0.5;
    // where to save data
    char *filename = argv[1];

    // setup stack. Size will be increased dynamically if needed.
    stack s;
    s.top = 0;
    s.length = NH;
    s.NH = NH;
    s.sites = malloc(s.length*sizeof(ul));
    if (!s.sites)
    {
        printf("Error using malloc for s.sites!\n");
        return 2;
    }

    // keep track of visited nodes
    bool *visited = malloc(sizeof(bool)*NH);
    if (!visited)
    {
        printf("Error using malloc for visited!\n");
        return 2;
    }

    // open the output file
    FILE *sizes = fopen(filename, "w");
    if (!sizes)
    {
        printf("Error opening file: %s!\n", filename);
        return 3;
    }

    // where to grow the cluster
    ul start_site = 0; 
    ul size;
    int error_flag = 0;


    // run the DFS algorithm over NR realisations
    for (int i = 0; i < NR; i++)
    {
        // set all nodes to not visited
        reset_visited(visited, NH); 
        // run DFS algorithm, get a cluster size
        size = DFS_hypercube(&s, visited, p, N, start_site, RNG, &error_flag);
        if (error_flag == -1)
        {
            // error!
            printf("Error! Exiting DFS algorithm and closing file...\n");
            free(s.sites);
            free(visited);
            fclose(sizes);
            return 1;
        }
        fprintf(sizes, "%lu\n", size);
    }

    // free heap memory, close file
    free(s.sites);
    free(visited);
    fclose(sizes);
    return 0;
}
