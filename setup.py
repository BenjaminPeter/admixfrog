from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import Cython.Compiler.Options
from Cython.Compiler.Options import get_directive_defaults
directive_defaults = get_directive_defaults()
directive_defaults['linetrace'] = True
directive_defaults['binding'] = True



extensions = [
        Extension(
                    "admixfrog.hmm",
                    ["admixfrog/hmm.pyx", "admixfrog/hmm.pyx"],
                    extra_compile_args=['-fopenmp'],
                    extra_link_args=['-fopenmp'],
                ),
    
]


setup(
    name='admixfrog',
    version='0.1',
    description='HMM to call fragments from contaminated genomes',
    author='Ben Peter',
    author_email='benjamin_peter@eva.mpg.de',
    ext_modules = cythonize(extensions, annotate=True),
    packages = ["admixfrog"]
)

