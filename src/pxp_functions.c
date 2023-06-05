/* Core algorithms for growing percolation clusters on the PXP model graph, also called
the Fibonacci cube. */
#define NO_IMPORT_ARRAY
#include "functions.h"
#include <stdlib.h>
#include <stdio.h>

extern gsl_rng *RNG;

/*
 * Function:  PXP_sites
 * --------------------
 * construct the (ordered) list of nodes/sites for the PXP graph
 *
 *  N: the dimension of the Fibonacci cube
 *
 *  returns: a NumPy array of sites of type NPY_ULONGLONG
 */
PyObject *PXP_sites(PyObject *self, PyObject *args)
{
    PyObject *py_N = NULL; // N as a Python object

    // parse and check arguments
    if (!PyArg_ParseTuple(args, "O", &py_N)) goto error;

    // Fibonacci cube dimension
    ul N = pyobject_to_ul(py_N);
    // Check for overflow and invalid arguments
    if (PyErr_Occurred() || !check_args(N, 1, 1)) goto error;

    // the size of the graph
    ul NH = fibonacci(N+2);

    ul *sitelist = construct_PXP_sitelist(N);
    if (!sitelist) goto error;
    
    return CArrayToNumPyArray(sitelist, NH);

error:
    if (!PyErr_Occurred()) PyErr_SetString(PyExc_RuntimeError, "Fatal error occurred");
    if (sitelist) free(sitelist);
    return NULL;
}


/*
 * Function:  H_PXP
 * --------------------
 * build the Hamiltonian (adjacency matrix) for the PXP model, as a NumPy array.
 *
 *  N: the dimension of the Fibonacci cube
 *  p: the percolation concentration
 *
 *  returns: a pointer to the Ndarray (Hamiltonian matrix).
 */
PyObject *PXP_H(PyObject *self, PyObject *args)
{

    PyObject *py_N = NULL; // N as a Python object
    ul N; // Fibonacci cube dimension
    float p; // percolation concentration

    // parse and check arguments
    if (!PyArg_ParseTuple(args, "Of", &py_N, &p)) 
    {
        PyErr_SetString(PyExc_ValueError, "Invalid input arguments in H_PXP(N, p)");
        goto error;
    }
    N = pyobject_to_ul(py_N);
    // Check for overflow/arguements
    if (PyErr_Occurred() || !check_args(N, 1, p)) goto error;

    // the size of the graph
    ul NH = fibonacci(N+2);

    // a list of the PXP sites
    ul *sitelist = construct_PXP_sitelist(N);
    if (!sitelist) goto error;

    // Create a zeroed NH x NH Ndarray of ints
    npy_intp dimensions[2] = {NH, NH};
    PyArrayObject *numpy_array = (PyArrayObject *) PyArray_ZEROS(2, dimensions, NPY_INT, 0);
    if (!numpy_array)
    {
        PyErr_SetString(PyExc_RuntimeError, "Unable to create NumPy array in H_PXP");
        return NULL;
    }

    int error = 0;
    int connected = 1;
    int disconnected = 0;
    ul row, col_index, flipped;

    // Loop over the PXP nodes
    for (ul j = 0; j < NH; j++) 
    {
        row = sitelist[j];

        for (ul i = 0; i < N; i++)
        {
            // flip the ith bit
            if (PXP_flip_allowed(row, i, N))
            {
                flipped = row ^ (1UL << i);
                col_index = index_site(sitelist, flipped, 0, NH-1, &error);
                if (error != 0) { goto error; }

                
                // ptr to matrix[row, col]
                int *array_ptr = (int *) PyArray_GETPTR2(numpy_array, j, col_index);

                // ptr to the transpose element
                int *array_ptr_T = (int *) PyArray_GETPTR2(numpy_array, col_index, j);

                // with probability p, create a link
                // of course, we double the number of operations: but this is needed 
                // to preserve Hermiticity!
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
    }
    

    return numpy_array;

    error:
        if (sitelist) free(sitelist);
        if (numpy_array) Py_DECREF(numpy_array);
        return NULL;
}

/*
 * Function:  PXP_clusters
 * --------------------
 *  driver code for running DFS_PXP many times and returning a pointer to the cluster sizes.
 *  Checks the user input. For each realisation, starts from a random site
 *
 *  N: the dimension of the hypercube
 *  NR: the Number of Realisations: number of clusters to grow
 *  p: the percolation concentration. 0 <= p <= 1
 *  error: a pointer to an error flag in case of probems.
 *
 *  returns: a pointer to an NumPy array of NR cluster sizes, of type ul (unsigned long)
 */
PyObject *PXP_clusters(PyObject *self, PyObject *args)
{

    PyObject *py_N = NULL; // N as a Python object
    ul N; // hypercube dimension
    int NR; // Number of Realisations
    float p; // percolation concentration

    if (!PyArg_ParseTuple(args, "Oif", &py_N, &NR, &p)) goto error;

    N = pyobject_to_ul(py_N);
    // Check for overflow & arguments
    if (PyErr_Occurred() || !check_args(N, NR, p)) goto error;

    // the size of the graph
    ul NH = fibonacci(N+2);

    stack *s = setup_stack(NH);
    if (!s) goto error;

    // a list of the sites
    ul *sitelist = construct_PXP_sitelist(N);
    if (!sitelist) goto error;

    // keep track of visited nodes
    bool *visited = malloc(sizeof(bool)*NH);
    if (!visited)
    {
        PyErr_SetString(PyExc_RuntimeError, "Error setting up visited");
        goto error;
    }

    // where we will grow the cluster from
    ul start_site;

    // Create a NumPy array of uls
    npy_intp dimensions[1] = {NR};
    PyArrayObject *numpy_array = (PyArrayObject *) PyArray_SimpleNew(1, dimensions, NPY_ULONG);
    if (!numpy_array)
    {
        PyErr_SetString(PyExc_RuntimeError, "Unable to create NumPy array in PXP_clusters");
        goto error;
    }

    int error_flag = 0;
    // run the DFS algorithm over NR realisations
    npy_intp index[1];
    for (npy_intp i = 0; i < NR; i++)
    {
        // set all nodes to not visited
        reset_visited(visited, NH); 
        // choose a random start site
        start_site = sitelist[gsl_rng_uniform_int(RNG, NH)];

        // run DFS algorithm, get a cluster size
        index[0] = i;
        ul *array_ptr = (ul *) PyArray_GetPtr(numpy_array, index);
        *array_ptr =  DFS_PXP(s, sitelist, visited, p, N, start_site, &error_flag);

        if (error_flag == -1) goto error;
    }

    // free heap memory, except cluster sizes
    free(s->sites);
    free(s);
    free(visited);
    free(sitelist);

    return numpy_array;

    error:
        if (s)
        {
            free(s->sites);
            free(s);
        }
        if (visited) free(visited);
        if (numpy_array) Py_DECREF(numpy_array);
        if (sitelist) free(sitelist);
        return NULL;

}

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
 *  error: a pointer to an error flag, in case something goes wrong
 *
 *  returns: the size of the cluster. I.e., the number of sites visited by the DFS algorithm.
 */
ul DFS_PXP(stack *s, ul *sites, bool visited[], const float p, const ul N, const ul start_state, int *error)
{
    ul size = 0; // cluster size
    ul NH = fibonacci(N+2);
    ul u, v, idx_u, idx_v;
    int idx_flag = 0;

    if (s->top != 0)
    {
        PyErr_SetString(PyExc_RuntimeError, "Error! Stack not empty in DFS_PXP");
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

/*
 * Function:  construct_PXP_sitelist
 * --------------------
 * Allocate memory for and construct a list of the sites (nodes) of the PXP graph (Fibonacci cube)
 *
 *  N: the dimension of the graph
 * 
 * returns: a pointer to an array of uls containing the sites

 */
ul *construct_PXP_sitelist(ul N)
{
    // the size of the graph
    ul NH = fibonacci(N+2);

    // a list of the PXP sites
    ul *sitelist = malloc(sizeof(ul)*NH);
    if (sitelist == NULL)
    {
        PyErr_SetString(PyExc_RuntimeError, "Error setting up sitelist");
        return NULL;
    }

    // fill sitelist with PXP nodes
    populate_sites_PXP(sitelist, N);

    return sitelist;
}
