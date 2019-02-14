from setuptools import setup#, find_namespace_packages
from distutils.extension import Extension


USE_CYTHON = False

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [Extension("admixfrog",
                        ["admixfrog/distributions"+ext,
                         "admixfrog/hmm_updates"+ext
                         ])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

    setup(setup_requires["pbr"],
        ext_modules = extensions,
          pbr=True)

setup(setup_requires=["pbr"], 
      ext_modules = extensions,#cythonize(["admixfrog/*.pyx"], annotate=True),
      pbr=True)

