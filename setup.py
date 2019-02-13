from setuptools import setup, find_namespace_packages
from Cython.Build import cythonize
import Cython.Compiler.Options
from Cython.Compiler.Options import get_directive_defaults
directive_defaults = get_directive_defaults()
#directive_defaults['linetrace'] = True
#directive_defaults['binding'] = True



setup(
    name='admixfrog',
    version='0.2',
    description='HMM to call fragments from contaminated genomes',
    author='Ben Peter',
    author_email='benjamin_peter@eva.mpg.de',
    ext_modules = cythonize(["admixfrog/*.pyx"], annotate=True),
#    ext_modules = cythonize(["admixfrog/*.pyx"], annotate=False),
    packages=find_namespace_packages(),
    setup_requires=['scipy', 'cython'],
    install_requires=[
                  'numba',
                  'numpy',
                  'scipy',
                  'pandas>=0.24',
                  'cython'
      ],
    entry_points={
        'console_scripts': [
            'admixfrog=admixfrog:run'
        ]
    }
)

