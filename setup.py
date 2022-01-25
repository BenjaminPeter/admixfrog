from setuptools import setup#, find_namespace_packages
from distutils.extension import Extension


USE_CYTHON = False

ext = '.pyx' if USE_CYTHON else '.c'

extensions = [Extension("admixfrog.utils.distributions", ["admixfrog/utils/distributions"+ext]),
              Extension("admixfrog.gll.read_emissions", ["admixfrog/gll/read_emissions"+ext])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

    setup(
        setup_requires =["pbr>=1.9", 'setuptools>=17.1', 
                         'pytest-runner'
                         ],
        tests_require = ['pytest',
                         'pytest-console-scripts',
                         'pytest-cov'],
        ext_modules = extensions,
          pbr=True)
else:

    setup(
        setup_requires=["pbr>=1.9", 'setuptools>=17.1', 'pytest-runner'],
        python_requires=">=3.6",
        tests_require = ['pytest',
                         'pytest-console-scripts',
                         'pytest-cov'],
        ext_modules = extensions,#cythonize(["admixfrog/*.pyx"], annotate=True),
          pbr=True)

