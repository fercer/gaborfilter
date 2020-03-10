#!/usr/bin/env python
import os
import numpy as np
from distutils.core import setup, Extension

setup(name='gabor',
      version='1.0',
      description='Patch extraction tools for coronary angiograms',
      author='Fernando Cervantes',
      author_email='iie.fercer@gmail.com',
      
      ext_modules=[Extension('gabor', ['include/gabor.c'],
                             define_macros=[('BUILDING_PYTHON_MODULE','BUILDING_GABOR_DLL'), ('NDEBUG',)],
                             include_dirs=[np.get_include(), os.environ['FFTW_INCLUDE_PATH']],
                             library_dirs=[os.environ['FFTW_LIBS_PATH']],
                             libraries=['fftw3'],
                             )],
      )
