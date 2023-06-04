/* Core functions for running a depth-first search algorithm on a hypercube of arbitrary
dimension for the percolation problem. */
#define NO_IMPORT_ARRAY
#include "functions.h"
#include <stdlib.h>
#include <stdio.h>

extern gsl_rng *RNG;

/*
 * Function:  H_hypercube
 * --------------------
 * build the Hamiltonian (adjacency matrix) for the hypercube, as a NumPy array.
 *
 *  N: the dimension of the hypercube
 *  p: the percolation concentration
 *
 *  returns: a pointer to the Ndarray (Hamiltonian matrix).
 */
PyObject* H_hypercube(PyObject *self, PyObject *args)
{
    
    PyObject *py_N = NULL; // N as a Python object
    ul N; // hypercube dimension
    float p; // percolation concentration

    if (!PyArg_ParseTuple(args, "Of", &py_N, &p)) goto error;

    N = pyobject_to_ul(py_N);
    // Check for overflow
    if (PyErr_Occurred()) goto error;

    if (!check_args(N, 1, p)) 
    {
        PyErr_SetString(PyExc_ValueError, "Invalid input arguments");
        goto error;
    }

    // the size of the graph
    ul NH = intpower(2, N); 

    // Create a new NumPy array of integers with the same dimensions
    npy_intp dimensions[2] = {NH, NH};
    PyArrayObject *numpy_array = (PyArrayObject *) PyArray_ZEROS(2, dimensions, NPY_INT, 0);
    if (!numpy_array)
    {
        PyErr_SetString(PyExc_RuntimeError, "Unable to create NumPy array in H_hypercube");
        goto error;
    }

    int connected = 1;
    int disconnected = 0;
    // Loop over the matrix elements
    for (ul row = 0; row < NH; row++) 
    {
        for (ul i = 0; i < N; i++)
        {
            // flip the ith bit
            ul col = row ^ (1UL << i);

            int *array_ptr = (int *) PyArray_GETPTR2(numpy_array, row, col);
            // ptr to the transpose element
            int *array_ptr_T = (int *) PyArray_GETPTR2(numpy_array, col, row);

            // with probability p, create a link
            if (gsl_rng_uniform(RNG) < p)
            {
                
                *array_ptr = connected;
                *array_ptr_T = connected;
            }
            else
            {
                *array_ptr = disconnected;
                *array_ptr_T = disconnected;
            }
        }
    }

    return numpy_array;

    error:
        if (!PyErr_Occurred()) PyErr_SetString(PyExc_RuntimeError, "Fatal error occurred");
        if (numpy_array) Py_DECREF(numpy_array);
        return NULL;
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
ul DFS_hypercube(stack *s, bool visited[], const float p, const ul N, const ul start_state, int *error)
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
 *
 *  returns: a PyObject * (NumPy array) of NR cluster sizes, of type NPY_ULONG (unsigned long)
 */
PyObject *hypercube_clusters(PyObject *self, PyObject *args)
{

    PyObject *py_N = NULL; // N as a Python object
    ul N; // hypercube dimension
    int NR; // Number of Realisations
    float p; // percolation concentration

    if (!PyArg_ParseTuple(args, "Oif", &py_N, &NR, &p)) goto error;

    N = pyobject_to_ul(py_N);
    // Check for overflow
    if (PyErr_Occurred()) goto error;

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
        *array_ptr = DFS_hypercube(s, visited, p, N, start_site, &error_flag);

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

    return numpy_array;

    error:
        if (s)
        {
            free(s->sites);
            free(s);
        }
        if (visited) free(visited);
        if (numpy_array) Py_DECREF(numpy_array);
        return NULL;

}
