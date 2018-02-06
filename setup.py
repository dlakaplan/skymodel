from __future__ import print_function
#from setuptools import setup, Extension
import numpy, os, sys, os.path
import glob
# need this to use f2py from within numpy
from numpy.distutils.core import setup, Extension

if '--use-cython' in sys.argv:
    USE_CYTHON = True
    sys.argv.remove('--use-cython')
else:
    USE_CYTHON = False
ext = '.pyx' if USE_CYTHON else '.c'

# standard YMW16 .c files
ymw16_sources=['dmdtau.c',
               'dora.c',
               'fermibubble.c',
               'frb_d.c',
               'galcen.c',
               'gum.c',
               'lmc.c',
               'localbubble.c',
               'ne_crd.c','nps.c',
               'smc.c',
               'spiral.c',
               'thick.c',
               'thin.c',
               'ymw16par.c']

# check for required modules
try:
    import astropy
except ImportError as e:
    print("ERROR! You are missing a dependency!")
    print(e)
    raise

try:
    import pint
except ImportError:
    print('You are missing PINT: certain functionality may not work')

try:
    import pygsm
    import healpy
except ImportError:
    print('You are missing PyGSM: certain functionality may not work')
    
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
packages = ['skymodel','skymodel.ne2001','skymodel.ymw16','skymodel.data']


ymw16_sources = ['skymodel/ymw16/src/' + x for x in ["ymw16_dmdtau" + ext] + ymw16_sources]


extensions=[
    Extension("skymodel.ymw16.ymw16_dmdtau",
              ymw16_sources,
              language='c',
              libraries=['m'],
              extra_compile_args=['-fPIC','-fcommon'],
              extra_link_args=['-fPIC']),
    Extension(name='skymodel.ne2001.ne2001',
              sources=['skymodel/ne2001/src/dmdsm.NE2001.f',
                       'skymodel/ne2001/src/density.NE2001.f',
                       'skymodel/ne2001/src/NE2001.f',
                       'skymodel/ne2001/src/neclumpN.NE2001.f',
                       'skymodel/ne2001/src/neLISM.NE2001.f',
                       'skymodel/ne2001/src/nevoidN.NE2001.f',
                       'skymodel/ne2001/src/scattering98.f'],                       
              f2py_options=['']
              )
]


if USE_CYTHON:
    from Cython.Build import cythonize, build_ext
    extensions = cythonize(extensions)

setup(
    name="skymodel",
    author = "D. Kaplan/NANOGrav",
    author_email = "kaplan@uwm.edu",
    description = ("Set of tools for sky models"),
    license = "BSD",
    #cmdclass={ "version": Version},
    #keywords = "MWA radio",
    #url = "http://mwa-lfd.haystack.mit.edu",
    packages=packages,
    #package_dir={'mwapy':'mwapy','':'configs'},
    scripts=glob.glob('scripts/*.py'),
    ext_modules=extensions,
    long_description=read('README.md'),
    # put the data files in same place as the python module
    data_files=[('skymodel/data/', glob.glob('data/*.dat') + glob.glob('data/*.inp') + glob.glob('data/*.txt'))]
    )
    
