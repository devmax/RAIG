from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'viterbi app',
  ext_modules = cythonize("Viterbi_cython.pyx"),
)
