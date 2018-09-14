#!/usr/bin/env python

"""
setup.py file for epg
"""

 #extra_compile_args = [ "-std=gnu99", "-O3", "-march=native", "-ffast-math", "-fstrict-aliasing",],

from distutils.core import setup, Extension
import os
import sys
import numpy

if __name__ == "__main__":

    [setup_py_directory, filename] = os.path.split(__file__)

    eigen_include = os.path.join( setup_py_directory, 'include', 'Eigen334')
    numpy_include1 = numpy.get_include()
    numpy_include2 = os.path.join( numpy_include1, 'numpy')


    if "win32" == sys.platform:


    	epg_module = Extension('_epg',
    				   sources=['epg.i',  'epg_cpmg.cpp',],
    				   swig_opts= ['-c++',],
    				   include_dirs=[numpy_include1, numpy_include2, ],
    				    extra_compile_args = [ "/GS-",  "/Zc:inline", "/fp:fast", ],
    				   )
    elif "linux" == sys.platform:

    	epg_module = Extension('_epg',
    				   sources=['epg.i',  'epg_cpmg.cpp',],
    				   swig_opts= ['-c++',],
    				   include_dirs=[numpy_include1, numpy_include2, ],
    				    extra_compile_args = [  "-std=gnu++11", "-O3", "-march=native", "-ffast-math", "-fstrict-aliasing",],

    				   )


    setup (name = 'epg',
           version = '0.1',
           author      = "Eric Hughes",
           author_email = 'erichughesabc@gmail.com',
           description = """EPG Extended Phase Graph CPMG module""",
	   keywords="EPG MRI python Extended Phase Graph T2 CPMG",
	   license="MIT",
           ext_modules = [epg_module],
           py_modules = ["epg"],
           )
