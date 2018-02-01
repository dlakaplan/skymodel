# SkyModel
generic/pulsar-specific interface to the sky model for computing:
 * sky temperatures (Global Sky Model)
 * dispersion measure (from distance)
 * distance (from dispersion measure)

Requires:
 * astropy
 * global sky model python interface [pygsm](https://github.com/telegraphic/PyGSM)
 * numpy/f2py
Optional (to work with .par files):
 * pulsar timing interface [PINT](https://github.com/nanograv/PINT)
 
Usage:
```python
from astropy.coordinates import SkyCoord
from astropy import units as u, constants as c
from skymodel import skymodel

source = SkyCoord("05h35m17.3s -05d23m28s", frame='icrs')

d=skymodel.SkyModel()
print d.distance(source, 10)
print d.DM(source, 659*u.pc)
print d.Tsky(source)

print d.distance('0038-2501.par', 10)
print d.DM('0038-2501.par', 449*u.pc)
print d.Tsky('0038-2501.par')

ra=np.linspace(0,240,15).reshape((5,3))
dec=np.linspace(-30,30,15).reshape((5,3))

source=SkyCoord(ra,dec,unit='deg')
print d.Tsky(source)
print d.DM(source, 1*u.kpc)
print d.distance(source, 10)

```

 * When model is initialized the frequency (for Tsky) can be passed as a number (MHz assumed) or an astropy quantity.
 * The Tsky model can be 'GSM2008' or 'GSM2016'
 * parfiles are parsed by PINT
 * Currently only NE2001 is supported for DM models
 * Source positions can be vectorized
 * There are convenience functions as well:
 ```python
 print skymodel.Tsky(source)
 print skymodel.DM(source, 10*u.pc)
 print skymodel.distance(source, 10)
```

To Do:
 * scripts for command-line access
 * vectorize distance/DM too?
 * add YMW16 model
 


