from setuptools import setup#, find_namespace_packages
from Cython.Build import cythonize

setup(setup_requires=["pbr"], 
      ext_modules = cythonize(["admixfrog/*.pyx"], annotate=True),
      pbr=True)

