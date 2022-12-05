#from setuptools import setup#, find_namespace_packages
from setuptools import Extension
from Cython.Build import cythonize


extensions = [Extension("admixfrog.utils.distributions",
                        ["admixfrog/utils/distributions.pyx"]),
              Extension("admixfrog.gll.read_emissions",
                        ["admixfrog/gll/read_emissions.pyx"])]

extensions = cythonize(extensions)

def build(setup_kwargs):
    setup_kwargs.update({
        'ext_modules' : cythonize(extensions, language_level=3),
        'zip_safe' : False})
