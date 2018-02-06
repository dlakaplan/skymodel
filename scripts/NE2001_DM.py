#!/usr/bin/env python

from __future__ import print_function
from astropy.coordinates import SkyCoord
from astropy import units as u, constants as c
import os,sys
import numpy as np
from optparse import OptionParser,OptionGroup

from skymodel import skymodel



def main():
    

    
    usage="Usage: %prog [options] DM\n"
    usage+='\tTranslates a distance into a DM using NE2001\n'
    
    parser = OptionParser(usage=usage)
    parser.add_option('-p','--par',dest="parfile",default=None,
                      help="Get positions from a parfile (PINT required)")
    parser.add_option('-r','--ra',default=None,
                      help="Right Ascension")
    parser.add_option('-d','--dec',default=None,
                      help="Declanation")
    parser.add_option('-l','--l','--gall',dest='l',default=None,
                      help="Galactic longitude")
    parser.add_option('-b','--b','--galb',dest='b',default=None,
                      help="Galactic latitude")
    parser.add_option('-u','--unit',default='kpc',
                      help='Units for distance [default=%default]')

    (options, args) = parser.parse_args()

    if options.ra is not None:
        if options.dec is None:
            raise ValueError, "Must supply Dec with RA"
        try:
            ra=float(options.ra)
            dec=float(options.dec)
            source=SkyCoord(ra,dec,unit='deg')
        except:
            if ':' in options.ra:
                try:
                    source=SkyCoord(options.ra,options.dec,unit=('hour','deg'))
                except:
                    raise ValueError,'Cannot parse (RA,Dec)=(%s,%s)' % (options.ra,options.dec)
            else:        
                source=SkyCoord(options.ra,options.dec)
                                
    elif options.l is not None:
        if options.b is None:
            raise ValueError, "Must supply b with l"
        try:
            l=float(options.l)
            b=float(options.b)
            source=SkyCoord(l,b,unit='deg',frame='galactic').icrs
        except:
            source=SkyCoord(options.l,options.b,frame='galactic').icrs
    elif options.parfile is not None:
        source=skymodel.parfile2SkyCoord(options.parfile)
        
    m=skymodel.SkyModel(dmmodel='NE2001')
    for distance in args:
        d=eval(distance + '*u.' + options.unit)
        s='For source (RA,Dec)=%s' % (source.to_string('hmsdms'))
        if options.parfile is not None:
            s+=' [%s]' % options.parfile
        print(s)
        print('(l,b)=%s' % source.galactic.to_string('decimal'))
        print('and d=%s' % d)
        print(m.DM(source, d)[0])
        

 






######################################################################

if __name__=="__main__":
    main()
