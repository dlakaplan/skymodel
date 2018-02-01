#from setuptools import setup, Extension
import numpy, os, sys, os.path
import glob
# need this to use f2py from within numpy
from numpy.distutils.core import setup, Extension

# check for required modules
try:
    import pygsm
    import astropy
    import healpy
except ImportError,e:
    print "ERROR! You are missing a dependency!"
    print e
    raise

try:
    import pint
except ImportError:
    print 'You are missing PINT: certain functionality may not work'
    
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
packages = ['skymodel','skymodel.ne2001','skymodel.data']
#pythonlib=get_python_lib()
#pythonlib=EXEC_PREFIX

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
    scripts=glob.glob('scripts/*'),
    ext_modules=[
    Extension(name='skymodel.ne2001',
              sources=glob.glob('skymodel/ne2001/*NE2001.f') + ['skymodel/ne2001/scattering98.f'],
              f2py_options=['']
              )
    ],
    long_description=read('README.md'),
    # put the data files in same place as the python module
    data_files=[('skymodel/data/', glob.glob('data/*.dat') + glob.glob('data/*.inp'))]
    )
    
