from distutils.core import setup
from distutils.extension import Extension

from Cython.Distutils import build_ext

include_dirs = []
ext_modules = [Extension("g4hunter", ["g4hunter.pyx"], include_dirs)]

setup(name="optimisation cython", cmdclass={"build_ext": build_ext}, ext_modules=ext_modules)
