from setuptools import setup#, find_namespace_packages
from distutils.extension import Extension


USE_CYTHON = False

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [Extension("admixfrog.distributions", ["admixfrog/distributions"+ext]),
              Extension("admixfrog.hmm_updates", ["admixfrog/hmm_updates"+ext])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

    setup(
        setup_requires =["pbr>=1.9", 'setuptools>=17.1'],
        ext_modules = extensions,
          pbr=True)

setup(
    setup_requires=["pbr>=1.9", 'setuptools>=17.1'],
    python_requires=">=3.6",
    ext_modules = extensions,#cythonize(["admixfrog/*.pyx"], annotate=True),
      pbr=True)

