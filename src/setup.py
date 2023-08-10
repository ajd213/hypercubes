from distutils.core import setup, Extension
import numpy

module = Extension("hypergraphs", sources = ["hypergraphs.c", "hypercube_functions.c", "pxp_functions.c", "helper_functions.c"],
libraries = ["m", "gsl", "gslcblas"])
setup(name = "hypergraphs", 
version = "1.0", 
description = "Percolation on hypergraphs", 
include_dirs=[numpy.get_include()],
ext_modules = [module])
