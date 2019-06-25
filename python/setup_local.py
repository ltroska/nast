from setuptools import setup, Extension
from Cython.Build import cythonize

ext = [Extension("nast", sources=["nast.pyx"], include_dirs=["../include"])]

setup(name = "nast",
      version = "0.1",
      ext_modules = cythonize(ext)
      )
