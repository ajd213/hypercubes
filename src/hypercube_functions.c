/* Core functions for running a depth-first search algorithm on a hypercube of arbitrary
dimension for the percolation problem. */
#define NO_IMPORT_ARRAY
#include "functions.h"
#include <stdlib.h>
#include <stdio.h>


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
ul DFS_hypercube(stack *s, bool visited[], const float p, const ul N, const ul start_state, gsl_rng *RNG, int *error)
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
 * Function:  hypercube_clusters
 * --------------------
 *  driver code for running DFS_hypercube many times and returning a NumPy array of cluster sizes.
 *  Checks the user input.
 *
 *  N: the dimension of the hypercube
 *  NR: the Number of Realisations: number of clusters to grow
 *  p: the percolation concentration. 0 <= p <= 1
 *  error: a pointer to an error flag in case of probems.
 *
 *  returns: a PyObject * (NumPy array) of NR cluster sizes, of type NPY_ULONG (unsigned long)
 */
PyObject *hypercube_clusters(PyObject *self, PyObject *args)
{
    // set and seed the RNG
    gsl_rng *RNG2 = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RNG2, time(NULL));

    ul N; // hypercube dimension
    int NR; // Number of Realisations
    float p; // percolation concentration

    if (!PyArg_ParseTuple(args, "kif", &N, &NR, &p)) goto error;

    ul NH = intpower(2, N); // the size of the graph

    if (!check_args(N, NR, p)) 
    {
        PyErr_SetString(PyExc_ValueError, "Invalid input arguments");
        goto error;
    }

    stack *s = setup_stack(NH);
    if (!s)
    {
        PyErr_SetString(PyExc_RuntimeError, "Error setting up stack");
        goto error;
    }

    bool *visited = malloc(sizeof(bool)*NH);
    if (!visited)
    {
        PyErr_SetString(PyExc_RuntimeError, "Error setting up visited");
        goto error;
    }


    ul start_site = 0; 

    // Create a NumPy array of uls
    npy_intp dimensions[1] = {NR};
    PyArrayObject *numpy_array = (PyArrayObject *) PyArray_SimpleNew(1, dimensions, NPY_ULONG);
    if (!numpy_array)
    {
        PyErr_SetString(PyExc_RuntimeError, "Unable to create NumPy array in hypercube_clusters");
        goto error;
    }

    int error_flag = 0;
    // run the DFS algorithm over NR realisations
    npy_intp index[1];
    for (npy_intp i = 0; i < NR; i++)
    {
        // set all nodes to not visited
        reset_visited(visited, NH); 

        // run DFS algorithm, get a cluster size
        index[0] = i;
        ul *array_ptr = (ul *) PyArray_GetPtr(numpy_array, index);
        *array_ptr = DFS_hypercube(s, visited, p, N, start_site, RNG2, &error_flag);

        if (error_flag == -1)
        {
            // error!
            PyErr_SetString(PyExc_RuntimeError, "Error in DFS algorithm!");
            goto error;
        }
    }

    // free heap memory, except cluster sizes
    free(s->sites);
    free(s);
    free(visited);
    gsl_rng_free(RNG2);

    return numpy_array;

    error:
        if (s)
        {
            free(s->sites);
            free(s);
        }
        if (visited) free(visited);
        if (numpy_array) Py_DECREF(numpy_array);
        if (RNG2) gsl_rng_free(RNG2);
        return NULL;

}
