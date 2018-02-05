from __future__ import print_function
#import sys;sys.path.append("./lib/python2.7/site-packages/")
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units as u, constants as c
from skymodel import skymodel

source = SkyCoord("05h35m17.3s -05d23m28s", frame='icrs')

d=skymodel.SkyModel()
print(d.distance(source, 10)[0])
print(d.DM(source, 659*u.pc)[0])
print(d.Tsky(source))

print(d.distance('0038-2501.par', 10)[0])
print(d.DM('0038-2501.par', 449*u.pc)[0])
print(d.Tsky('0038-2501.par'))

ra=np.linspace(0,240,15).reshape((5,3))
dec=np.linspace(-30,30,15).reshape((5,3))

source=SkyCoord(ra,dec,unit='deg')
try:
    print(d.Tsky(source))
except TypeError:
    print('Error converting to Galactic coords: probably astropy/numpy bug')
try:
    print(d.DM(source, 1*u.kpc)[0])
except TypeError:
    print('Error converting to Galactic coords: probably astropy/numpy bug')
try:
    print(d.distance(source, 10)[0])
except TypeError:
    print('Error converting to Galactic coords: probably astropy/numpy bug')




print(skymodel.Tsky(source))
print(skymodel.DM(source, 10*u.pc)[0])
print(skymodel.distance(source, 10)[0])

