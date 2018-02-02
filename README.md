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

And there are command-line scripts:
```
plock[test]% bin/NE2001_distance.py --par=1949+3426.par 10
For source (RA,Dec)=19h49m13.6713s +34d26m33.8648s [1949+3426.par]
(l,b)=69.7201 4.28843
and DM=10.000 pc/cm**3
0.890833377838 kpc

plock[test]% bin/NE2001_DM.py --par=1949+3426.par 10
For source (RA,Dec)=19h49m13.6713s +34d26m33.8648s [1949+3426.par]
(l,b)=69.7201 4.28843
and d=10.0 kpc
230.689941406 pc / cm3

plock[test]% bin/NE2001_DM.py --par=1949+3426.par --unit=pc 10
For source (RA,Dec)=19h49m13.6713s +34d26m33.8648s [1949+3426.par]
(l,b)=69.7201 4.28843
and d=10.0 pc
0.0499999858439 pc / cm3

plock[test]% bin/GSM2008_Tsky.py --par=1949+3426.par -f 350
For source (RA,Dec)=19h49m13.6713s +34d26m33.8648s [1949+3426.par]
(l,b)=69.7201 4.28843
and freq=350.0 MHz
76.676389301 K
```

To Do:
 * vectorize distance/DM too?
 * add YMW16 model
 * access scattering quantities


