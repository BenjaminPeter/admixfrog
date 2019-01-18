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
                    "hmm",
                    ["hmm2.pyx"],
                    extra_compile_args=['-fopenmp'],
                    extra_link_args=['-fopenmp'],
                )
]


setup(
    ext_modules = cythonize(extensions, annotate=True)
#    ext_modules = cythonize(extensions, annotate=True)
)

