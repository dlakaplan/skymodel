# pulsar_gsm
pulsar-specific interface to the global sky model for computing sky temperatures

Requires:
 * astropy
 * global sky model python interface [pygsm](https://github.com/telegraphic/PyGSM)
 * pulsar timing interface [PINT](https://github.com/nanograv/PINT)
 
Usage:
```python
import pulsar_Tsky
from astropy import units as u
Tsky=pulsar_Tsky.pulsar_Tsky(<parfile>, freq=1400*u.MHz, model='2008')
```

 * The frequency can be passed as a number (MHz assumed) or an astropy quantity.
 * The model can be '2008' or '2016'
 * parfiles are parsed by PINT
 
Or on the command line:
```
plock[gbncc]% python ./pulsar_Tsky.py -f 350 FYRE/0038-2501.par
WARNING: Unrecognized parfile line 'EPHVER         5' [pint.models.timing_model]
WARNING: Unrecognized parfile line 'T2CMETHOD      IAU2000B' [pint.models.timing_model]
WARNING: Unrecognized parfile line 'NE_SW          4.000' [pint.models.timing_model]
WARNING: Unrecognized parfile line 'CHI2R          0.0000 215' [pint.models.timing_model]
WARNING: Unrecognized parfile line 'EPHVER         5' [pint.models.timing_model]
WARNING: Unrecognized parfile line 'T2CMETHOD      IAU2000B' [pint.models.timing_model]
WARNING: Unrecognized parfile line 'NE_SW          4.000' [pint.models.timing_model]
WARNING: Unrecognized parfile line 'CHI2R          0.0000 215' [pint.models.timing_model]
PSR J0038-25: 27.1 K (at 350 MHz)
```
 

