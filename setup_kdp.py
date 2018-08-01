# run this script via
# 		$ python3 setup_kdp.py build_ext --inplace
# (or via command python2 if you want a python2 library).
# This creates a shared library kdp.so, whose objects (Vec2, Vec3, Vec4) can be imported via
# 		>>> from kdp import *
#
# I have configured my Linux + GNU system to have the following directory structure
#    ~/local/include/package/xyz.hpp      <= softlinks to my library headers, segregated (so I #include "package/xyz.hpp")
#    ~/local/lib/xyz.so                   <= softlinks to my libraries (not segregrated, all in the same folder)
#    ~/local/pyLib/xyz.so                 <= softlinks to my Python libraries (also not segregated)
# Then I have to set a few environmental variables (in my .bash_profile)
#    LD_LIBRARY_PATH = $HOME/local/lib/
#    PYTHONPATH = $HOMES/local/pyLib
# This allows me to easily build against and link to the C++ libraries, 
# and to "import module" from any running Python shell

from sys import version_info
from distutils.core import setup
from distutils.extension import Extension

flags = ['-std=c++11', '-mfpmath=sse', '-mieee-fp', '-march=native', '-ftree-vectorize'] # -O2 is default

from subprocess import check_output
# Add CPU specific flags to accelerate vector math
extraFlags = check_output(["sh", "getSSE_AVX.sh"]).split()

# convert 'bytes' to string in python3, and 
if (version_info > (3, 0)):
	for i in range(len(extraFlags)):	
		extraFlags[i] = extraFlags[i].decode()		
		
flags += extraFlags

try:
	from Cython.Distutils import build_ext

	# This build script leverages the pre-compiled kdp library (libkdp.so), 
	# a library which must be accessible from anywhere you intend to use pYqRand.
	# This keeps the C-compilation of pYqRand.pyx to a minimal (so the Python module is smaller),
	# and allows different compile flags to be used for the Cython build.

	setup(
	  name = "kdp",
	  ext_modules=[
		 Extension('kdp',
					  sources=['source/kdp.pyx'],
					  include_dirs = ['./include/'],
					  libraries = ['kdp'],
					  library_dirs = ['./lib'],
					  extra_compile_args=flags,
					  language='c++')
			 ], cmdclass = {'build_ext': build_ext})
except ImportError:
	setup(
	  name = "kdp",
	  ext_modules=[
		 Extension('kdp',
					  sources=['source/kdp.cpp'],
					  include_dirs = ['./include/'],
					  libraries = ['kdp'],
					  library_dirs = ['./lib'],
					  extra_compile_args = flags,
					  language='c++')
			 ], cmdclass = {'build_ext': build_ext})
