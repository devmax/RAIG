from distutils.core import setup
from Cython.Build import cythonize
import numpy


setup(
    name = 'viterbi app',
    ext_modules = cythonize("Viterbi_cython.pyx"),

    include_dirs=[numpy.get_include()]
)
