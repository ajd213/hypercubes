from distutils.core import setup, Extension

module = Extension("Hypercubes", sources = ["module.c", "functions.c"],
libraries = ["m", "gsl", "gslcblas"])
setup(name = "Hypercubes", 
version = "1.0", 
description = "Percolation on hypercubes.", 
ext_modules = [module])
