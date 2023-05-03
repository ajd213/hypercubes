#define PY_SSIZE_T_CLEAN
#include "functions.h"
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdio.h>

gsl_rng *RNG; // random number generator
PyObject *ClustersToPyList(ul *cs, ul NR);


static PyObject* H_hypercube(PyObject *self, PyObject *args)
{
    // set and seed the RNG
    RNG = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RNG, time(NULL));

    ul N; // hypercube dimension
    float p; // percolation concentration

    // per https://docs.python.org/3/c-api/arg.html#numbers 
    // parse two long unsigneds and a float

    if (!PyArg_ParseTuple(args, "kf", &N, &p))
    {
        return NULL;
    }

    int error = 0;
    if (!check_args(N, 1, p)) {error = 3; return NULL;}

    // the size of the graph
    ul NH = intpower(2, N); 

    // Create a new NumPy array of integers with the same dimensions
    npy_intp dimensions[2] = {NH, NH};
    PyObject *numpy_array = PyArray_ZEROS(2, dimensions, NPY_INT, 0);
    if (!numpy_array)
    {
        printf("Error: Unable to create NumPy array in MatrixToNumPyArray.\n");
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
                *((int *)PyArray_GETPTR2((PyArrayObject *) numpy_array, row, col)) = connected;
            }
        }
    }

    gsl_rng_free(RNG);
    return numpy_array;
}

static PyObject *hypercube_clusters(PyObject *self, PyObject *args)
{
    // set and seed the RNG
    RNG = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RNG, time(NULL));

    ul N; // hypercube dimension
    ul NR; // Number of Realisations
    float p; // percolation concentration

    // per https://docs.python.org/3/c-api/arg.html#numbers 
    // parse two long unsigneds and a float
    if (!PyArg_ParseTuple(args, "kkf", &N, &NR, &p))
    {
        return NULL;
    }


    // get the clusters
    int errorflag = 0;
    ul *cs = clusters_hypercube(N, NR, p, RNG, &errorflag);
    if (errorflag != 0 || !cs)
    {
        return NULL;
    }


    PyObject* cs_python = ClustersToPyList(cs, NR);

    free(cs);
    gsl_rng_free(RNG);
    return cs_python;

}

static PyObject *PXP_clusters(PyObject *self, PyObject *args)
{
    // set and seed the RNG
    RNG = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RNG, time(NULL));

    ul N; // pxp model dimension
    ul NR; // Number of Realisations
    float p; // percolation concentration

    // per https://docs.python.org/3/c-api/arg.html#numbers 
    // parse two long unsigneds and a float
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


    PyObject* cs_python = ClustersToPyList(cs, NR);
    free(cs);
    gsl_rng_free(RNG);
    return cs_python;

}


PyObject *ClustersToPyList(ul *cs, ul NR)
{
    PyObject* cs_python = PyList_New(NR);
    if (!cs_python)
    {
        return NULL;
    }

    for (ul i = 0; i < NR; i++)
    {
        PyObject* python_int = Py_BuildValue("k", cs[i]);
        if (!python_int)
        {
            // error! Perhaps OOB.
            // Release the items that were already added to the list
            for (ul j = 0; j < i; j++)
            {
                Py_DECREF(PyList_GetItem(cs_python, j));
            }
            // Release the list itself
            Py_DECREF(cs_python);
            return NULL;
        }

        PyList_SetItem(cs_python, i, python_int);
    }

    return cs_python;
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
    // Initialize the NumPy C API
    import_array();

    return PyModule_Create(&hypergraphs);
}
