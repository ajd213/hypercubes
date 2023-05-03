#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <stdio.h>
#include <gsl/gsl_rng.h> 
#include "functions.h"

gsl_rng *RNG; // random number generator

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
    ul *cs = clusters_hypercube(N, NR, p, &errorflag);
    if (errorflag != 0 || !cs)
    {
        return NULL;
    }


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
            return NULL;
        }

        PyList_SetItem(cs_python, i, python_int);
    }

    free(cs);
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
    ul *cs = clusters_PXP(N, NR, p, &errorflag);
    if (errorflag != 0 || !cs)
    {
        
        return NULL;
    }


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
            return NULL;
        }

        PyList_SetItem(cs_python, i, python_int);
    }

    free(cs);
    return cs_python;

}


static PyObject *version(PyObject *self)
{
    return Py_BuildValue("s", "Version 1.1");
}

static PyMethodDef myMethods[] = {
    {"hypercube_clusters", hypercube_clusters, METH_VARARGS, "Computes the sizes of NR clusters on a hypercube of dimension N with concentration p."},
    {"PXP_clusters", PXP_clusters, METH_VARARGS, "Computes the sizes of NR clusters on a PXP graph (Fibonacci cube) of dimension N with concentration p."},
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
    return PyModule_Create(&hypergraphs);
}
