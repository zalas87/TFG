from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
ext_modules1 = Extension("tfgmain", ["tfgmain.pyx"])
ext_modules2 = Extension("TrafficModels", ["TrafficModels.pyx"])

setup(
	name = '1',
	cmdclass = {'build_ext': build_ext},
	ext_modules = [ext_modules1, ext_modules2]
)
