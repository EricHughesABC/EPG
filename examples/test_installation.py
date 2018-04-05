# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 12:21:51 2018

@author: neh69
"""
import numpy
import os
import sys

if __name__ == "__main__":
    root = __file__
    print(numpy.get_include())
    print(root)
    print(os.path.split(root))
    print(os.path.join(numpy.get_include(), 'numpy'))

    include_dirs=[numpy.get_include(), os.path.join(numpy.get_include(), 'numpy'), os.path.join( root, 'include', 'Eigen334')]
    
    print(include_dirs)
    
    [setup_py_directory, filename] = os.path.split(root)
    
    print(setup_py_directory)