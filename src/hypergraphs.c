#include "functions.h"
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

gsl_rng *RNG;

static PyObject *version(PyObject *self)
{
    return Py_BuildValue("s", "Version 1.1");
}

static PyMethodDef myMethods[] = {
    {"hypercube_clusters", hypercube_clusters, METH_VARARGS, "Computes the sizes of NR clusters on a hypercube of dimension N with concentration p."},
    {"PXP_clusters", PXP_clusters, METH_VARARGS, "Computes the sizes of NR clusters on a PXP graph (Fibonacci cube) of dimension N with concentration p."},
    {"hypercube_H", hypercube_H, METH_VARARGS, "Compute the Hamiltonian for the Hypercube with concentration p."},
    {"hypercube_H_SC", hypercube_H_SC, METH_VARARGS, "Compute the Hamiltonian for a single cluster of the Hypercube with concentration p."},
    {"hypercube_H_LC", hypercube_H_LC, METH_VARARGS, "Compute the Hamiltonian for the LARGEST cluster of the Hypercube with concentration p."},
    {"hypercube_dijkstra", hypercube_dijkstra, METH_VARARGS, "Compute the shortest paths between site 0 and all other sites in the same cluster"},
    {"hypercube_dijkstra_LC", hypercube_dijkstra_LC, METH_VARARGS, "Compute the shortest paths within the largest cluster"},
    {"PXP_H", PXP_H, METH_VARARGS, "Compute the Hamiltonian for the PXP model with concentration p."},
    {"PXP_sites", PXP_sites, METH_VARARGS, "Return an array containing the basis sites/nodes of the PXP graph."},
    {"Hamming_distance", Hamming_distance, METH_VARARGS, "Return the Hamming distance between two integers."},
    {"RNG_test", RNG_test, METH_NOARGS, "Return a random int between 0 and INT_MAX"},
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
    
    import_array(); // Initialize the NumPy C API

    RNG = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(RNG, time(NULL));

    return module;
}
