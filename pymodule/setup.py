from distutils.core import setup, Extension

module = Extension("hypercubes", sources = ["hypercubes.c", "hypercube_functions.c", "pxp_functions.c", "helper_functions.c"],
libraries = ["m", "gsl", "gslcblas"])
setup(name = "hypercubes", 
version = "1.0", 
description = "Percolation on hypercubes", 
ext_modules = [module])
