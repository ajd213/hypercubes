#include "functions.h"
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

static PyObject *version(PyObject *self)
{
    return Py_BuildValue("s", "Version 1.1");
}

static PyMethodDef myMethods[] = {
    {"hypercube_clusters", hypercube_clusters, METH_VARARGS, "Computes the sizes of NR clusters on a hypercube of dimension N with concentration p."},
    {"PXP_clusters", PXP_clusters, METH_VARARGS, "Computes the sizes of NR clusters on a PXP graph (Fibonacci cube) of dimension N with concentration p."},
    {"H_hypercube", H_hypercube, METH_VARARGS, "Compute the Hamiltonian for the Hypercube with concentration p."},
    {"H_PXP", H_PXP, METH_VARARGS, "Compute the Hamiltonian for the PXP model with concentration p."},
    {"PXP_sites", PXP_sites, METH_VARARGS, "Return an array containing the basis sites/nodes of the PXP graph."},
    {"Hamming_distance", Hamming_distance, METH_VARARGS, "Return the Hamming distance between two integers."},
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
