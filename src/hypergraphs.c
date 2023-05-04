#include "functions.h"
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

gsl_rng *RNG; // random number generator
PyObject *CArrayToNumPyArray(ul *cs, ul NR);


static PyObject* H_hypercube(PyObject *self, PyObject *args)
{
    // set and seed the RNG
    RNG = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RNG, time(NULL));

    ul N; // hypercube dimension
    float p; // percolation concentration

    if (!PyArg_ParseTuple(args, "kf", &N, &p))
    {
        return NULL;
    }

    int error = 0;
    if (!check_args(N, 1, p)) 
    {
        PyErr_SetString(PyExc_ValueError, "Invalid input arguments");
        return NULL;
    }

    // the size of the graph
    ul NH = intpower(2, N); 

    // Create a new NumPy array of integers with the same dimensions
    npy_intp dimensions[2] = {NH, NH};
    PyArrayObject *numpy_array = (PyArrayObject *) PyArray_ZEROS(2, dimensions, NPY_INT, 0);
    if (!numpy_array)
    {
        PyErr_SetString(PyExc_RuntimeError, "Unable to create NumPy array in MatrixToNumPyArray");
        return NULL;
    }

    int connected = 1;
    // Loop over the matrix elements
    for (ul row = 0; row < NH; row++) 
    {
        for (ul i = 0; i < N; i++)
        {
            // flip the ith bit
            ul col = row ^ (1UL << i);

            // with probability p, create a link
            if (gsl_rng_uniform(RNG) < p)
            {
                int *array_ptr = (int *) PyArray_GETPTR2(numpy_array, row, col);
                *array_ptr = connected;
            }
        }
    }

    gsl_rng_free(RNG);
    return numpy_array;
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
static PyObject *hypercube_clusters(PyObject *self, PyObject *args)
{
    // set and seed the RNG
    RNG = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RNG, time(NULL));

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
    ul *cluster_sizes = malloc(NR*sizeof(ul));
    if (!cluster_sizes)
    {
        PyErr_SetString(PyExc_RuntimeError, "Error setting up cluster_sizes");
        goto error;
    }

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
        *array_ptr = DFS_hypercube(s, visited, p, N, start_site, RNG, &error_flag);

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
    gsl_rng_free(RNG);

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

static PyObject *PXP_clusters(PyObject *self, PyObject *args)
{
    // set and seed the RNG
    RNG = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RNG, time(NULL));

    ul N; // pxp model dimension
    ul NR; // Number of Realisations
    float p; // percolation concentration

    if (!PyArg_ParseTuple(args, "kkf", &N, &NR, &p))
    {
        return NULL;
    }


    // get the clusters
    int errorflag = 0;
    ul *cs = clusters_PXP(N, NR, p, RNG, &errorflag);
    if (errorflag != 0 || !cs)
    {
        
        return NULL;
    }


    PyObject* cs_python = CArrayToNumPyArray(cs, NR);

    gsl_rng_free(RNG);
    return cs_python;

}

PyObject *CArrayToNumPyArray(ul *cs, ul NR)
{
    npy_intp dims[] = {NR};
    PyObject *numpy_array = PyArray_SimpleNewFromData(1, dims, NPY_ULONGLONG, (void *)cs);

    if (!numpy_array)
    {
        PyErr_SetString(PyExc_RuntimeError, "Error: Unable to create NumPy array in CArrayToNumPyArray.");
        return NULL;
    }

    return numpy_array;
}

static PyObject *version(PyObject *self)
{
    return Py_BuildValue("s", "Version 1.1");
}

static PyMethodDef myMethods[] = {
    {"hypercube_clusters", hypercube_clusters, METH_VARARGS, "Computes the sizes of NR clusters on a hypercube of dimension N with concentration p."},
    {"PXP_clusters", PXP_clusters, METH_VARARGS, "Computes the sizes of NR clusters on a PXP graph (Fibonacci cube) of dimension N with concentration p."},
    {"H_hypercube", H_hypercube, METH_VARARGS, "Compute the Hamiltonian for the Hypercube with concentration p."},
    {"version", (PyCFunction) version, METH_NOARGS, "returns the version."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef hypergraphs = {
    PyModuleDef_HEAD_INIT, 
    "hypergraphs",
    "Module for percolation on hypergraphs",
    -1,
    myMethods
};

PyMODINIT_FUNC PyInit_hypergraphs(void)
{
    PyObject *module = PyModule_Create(&hypergraphs);
    // Initialize the NumPy C API
    import_array();

    return module;
}
