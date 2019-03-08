from setuptools import setup#, find_namespace_packages
from distutils.extension import Extension


USE_CYTHON = True

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [Extension("admixfrog.distributions", ["admixfrog/cydistributions"+ext]),
              Extension("admixfrog.read_emissions", ["admixfrog/cyread_emissions"+ext])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

    setup(
        setup_requires =["pbr>=1.9", 'setuptools>=17.1'],
        ext_modules = extensions,
          pbr=True)
else:

    setup(
        setup_requires=["pbr>=1.9", 'setuptools>=17.1'],
        python_requires=">=3.5",
        ext_modules = extensions,#cythonize(["admixfrog/*.pyx"], annotate=True),
          pbr=True)

