from pygsm import GlobalSkyModel,GlobalSkyModel2016
import numpy as np
from astropy.coordinates import SkyCoord
import healpy
from astropy import units as u
import astropy
from optparse import OptionParser
import os,sys
try:    
    import pint.models as models
    _usePINT=True
except ImportError:
    _usePINT=False

import ne2001

import data


def parfile2SkyCoord(parfile):
    """
    psr=parfile2SkyCoord(parfile)
    uses PINT to interpret a parfile
    and return a SkyCoord object
    """
    
    if not _usePINT:
        raise 'PINT is not available: cannot use parfiles'
    m=models.get_model(parfile)
    try:
        if m.RAJ.unitsuffix=='h':
            psr=SkyCoord(m.RAJ.value,m.DECJ.value,unit=('hour','deg'))
        else:
            psr=SkyCoord(m.RAJ.value,m.DECJ.value,unit=('deg','deg'))            
    except:
        try:
            psr=SkyCoord(m.ELONG.value*u.deg,m.ELAT.value*u.deg,frame='pulsarecliptic')
        except:
            raise KeyError,'Cannot find RAJ,DECJ or ELONG,ELAT in:\n%s' % m.as_parfile()

    return psr
##################################################
class SkyModel:
    """
    class interface to Sky models for Tsky and DM/distance
    
    so that structures/data stay in memory

    m=SkyModel(freq=350*u.MHz, tskymodel='2008', dmmodel='NE2001', datadir=data._datadir)
    

    can return distance, DM, or Tsky
    input can be SkyCoord or parfile

    m.DM(source, distance)
    m.distance(source, DM)
    m.Tsky(source)
    
    """

    def __init__(self, freq=350*u.MHz, tskymodel='2008', dmmodel='NE2001', datadir=data._datadir):
        """
        m=SkyModel(freq=350*u.MHz, tskymodel='2008', dmmodel='NE2001', datadir=data._datadir)
        """

        assert dmmodel.lower() in ['ne2001','ymw16']
        assert str(tskymodel) in ['2008','2016']


        self.datadir=datadir
        self.dmmodel=dmmodel
        if not isinstance(freq, astropy.units.quantity.Quantity):
            # assume MHz
            freq=freq*u.MHz

        if str(tskymodel)=='2008':
            self.gsm = GlobalSkyModel()
        elif str(tskymodel)=='2016':
            self.gsm=GlobalSkyModel2016()
        self.map=self.gsm.generate(freq.to(u.MHz).value)
        self.freq=freq
        self.tskymodel=tskymodel

    ##################################################
    def Tsky(self, source):
        """
        T=m.Tsky(source)

        returns sky temperature in K for the given model
        given a source, either a SkyCoord object or a parfile
        """
        
        if not isinstance(source, astropy.coordinates.sky_coordinate.SkyCoord):
            if isinstance(source,str):
                # assume .par file
                source=parfile2SkyCoord(source)
            else:
                raise TypeError, 'Do not know how to interpret an object of type %s' % source.__class__

        T=healpy.pixelfunc.get_interp_val(self.map,
                                          source.galactic.l.value,
                                          source.galactic.b.value,
                                          lonlat=True)
        return T*u.K


    ##############################
    def DM(self, source, distance, smweight='uniform'):
        """
        DM,SM=m.DM(source, distance, smweight='uniform')

        returns DM in standard units and scattering measure for the given model
        given a source, either a SkyCoord object or a parfile
        and a distance

        scattering measure can be weighted by:
        'uniform': uniform weighting [default]
        'tau': weighted for pulse broadening
        'theta': weighted for angular broadening of galactic sources
        'iso':  appropriate for calculating the isoplanatic angle at the source's location


        if no units supplied kpc assumed for distance
        """

        assert smweight.lower() in ['uniform','tau','theta','iso']

        if not isinstance(distance, astropy.units.quantity.Quantity):
            # assume kpc
            distance*=u.kpc            
        if distance <= 0:
            raise ValueError,'distance must be > 0'

        if not isinstance(source, astropy.coordinates.sky_coordinate.SkyCoord):
            if isinstance(source,str):
                # assume .par file
                source=parfile2SkyCoord(source)
            else:
                raise TypeError, 'Do not know how to interpret an object of type %s' % source.__class__
        source=source.icrs


        if len(source.ra.shape)==0:
            results=ne2001.dmdsm(self.datadir,
                                 np.radians(source.galactic.l.value),
                                 np.radians(source.galactic.b.value),
                                 -1,
                                 0,
                                 distance.to(u.kpc).value)
            if results[2]=='>':
                raise ValueError,'DM returned a lower limit'
            if smweight.lower() == 'uniform':
                SM=results[3]*u.kpc/u.m**(20./3)
            elif smweight.lower() == 'tau':
                SM=results[4]*u.kpc/u.m**(20./3)
            elif smweight.lower() == 'theta':
                SM=results[5]*u.kpc/u.m**(20./3)
            elif smweight.lower() == 'iso':
                SM=results[6]*u.kpc/u.m**(20./3)

            return results[0]*u.pc/u.cm**3,SM
        else:
            dm=np.zeros_like(source.ra.value)
            SM=np.zeros_like(source.ra.value)
            it = np.nditer(source.ra, flags=['multi_index'])
            while not it.finished:
                results=ne2001.dmdsm(self.datadir,
                                     np.radians(source[it.multi_index].galactic.l.value),
                                     np.radians(source[it.multi_index].galactic.b.value),
                                     -1,
                                     0,
                                     distance.to(u.kpc).value)
                if results[2]=='>':
                    raise ValueError,'DM returned a lower limit'
                dm[it.multi_index]=results[0]
                if smweight.lower() == 'uniform':
                    SM[it.multi_index]=results[3]
                elif smweight.lower() == 'tau':
                    SM[it.multi_index]=results[4]
                elif smweight.lower() == 'theta':
                    SM[it.multi_index]=results[5]
                elif smweight.lower() == 'iso':
                    SM[it.multi_index]=results[6]

                it.iternext()
            return dm*u.pc/u.cm**3,SM*u.kpc/u.m**(20./3)
            
    ##################################################
    def distance(self, source, DM, smweight='uniform'):
        """
        d,SM=m.distance(source, DM, smweight='uniform')

        returns distance in kpc and scattering measure for the given model
        given a source, either a SkyCoord object or a parfile
        and a DM

        scattering measure can be weighted by:
        'uniform': uniform weighting [default]
        'tau': weighted for pulse broadening
        'theta': weighted for angular broadening of galactic sources
        'iso':  appropriate for calculating the isoplanatic angle at the source's location

        if no units supplied standard DM units assumed for DM
        """

        assert smweight.lower() in ['uniform','tau','theta','iso']

        if not isinstance(DM, astropy.units.quantity.Quantity):
            # assume DM unit
            DM*=u.pc/u.cm**3
        if DM <= 0:
            raise ValueError, 'DM must be > 0'
        if not isinstance(source, astropy.coordinates.sky_coordinate.SkyCoord):
            if isinstance(source,str):
                # assume .par file
                source=parfile2SkyCoord(source)
            else:
                raise TypeError, 'Do not know how to interpret an object of type %s' % source.__class__
        source=source.icrs

        if len(source.ra.shape)==0:
            
            results=ne2001.dmdsm(self.datadir,
                                 np.radians(source.galactic.l.value),
                                 np.radians(source.galactic.b.value),
                                 1,
                                 DM.to(u.pc/u.cm**3).value,
                                 0)
            if results[2]=='>':
                raise ValueError,'distance returned a lower limit'
            distance=results[1]*u.kpc
            if smweight.lower() == 'uniform':
                SM=results[3]*u.kpc/u.m**(20./3)
            elif smweight.lower() == 'tau':
                SM=results[4]*u.kpc/u.m**(20./3)
            elif smweight.lower() == 'theta':
                SM=results[5]*u.kpc/u.m**(20./3)
            elif smweight.lower() == 'iso':
                SM=results[6]*u.kpc/u.m**(20./3)
            return distance,SM
        else:
            distance=np.zeros_like(source.ra.value)
            SM=np.zeros_like(source.ra.value)
            it = np.nditer(source.ra, flags=['multi_index'])
            dm=DM.to(u.pc/u.cm**3).value
            while not it.finished:
                results=ne2001.dmdsm(self.datadir,
                                     np.radians(source[it.multi_index].galactic.l.value),
                                     np.radians(source[it.multi_index].galactic.b.value),
                                     1,
                                     dm,
                                     0)
                if results[2]=='>':
                    raise ValueError,'distance returned a lower limit'
                distance[it.multi_index]=results[1]
                if smweight.lower() == 'uniform':
                    SM[it.multi_index]=results[3]
                elif smweight.lower() == 'tau':
                    SM[it.multi_index]=results[4]
                elif smweight.lower() == 'theta':
                    SM[it.multi_index]=results[5]
                elif smweight.lower() == 'iso':
                    SM[it.multi_index]=results[6]
                it.iternext()
            return distance*u.kpc,SM*u.kpc/u.m**(20./3)



##############################
# convenience functions
# provide a single interface
##############################
def Tsky(source, freq=350*u.MHz, model='2008'):
    """
    T=Tsky(source, freq=350*u.MHz)
    takes an astropy SkyCoord object or parfile and computes the sky temperature for that position
    using pyGSM
    scaled to the given frequency

    If frequency has no units, MHz is assumed

    model='2008' or '2016' for GSM2008 (Oliveira-Costa et al., 2008) or GSM2016 (Zheng et al., 2016)

    returns sky temperature in K
    """

    if not isinstance(source, astropy.coordinates.sky_coordinate.SkyCoord):
        if isinstance(source,str):
            # assume .par file
            source=parfile2SkyCoord(source)
        else:
            raise TypeError, 'Do not know how to interpret an object of type %s' % source.__class__


    m=SkyModel(freq=freq, tskymodel=model)
    return m.Tsky(source)


##############################
def DM(source, distance, model='NE2001'):
    """
    DM=DM(source, distance, model='NE2001')
    takes an astropy SkyCoord object or parfile
    return the DM associated with source at a given distance, for the specified DM model

    if distance has no units, kpc assumed

    return DM in DM units
    """
    
    d=SkyModel(dmmodel=model)
    return d.DM(source, distance)

##############################
def distance(source, DM, model='NE2001'):
    """
    distance=distance(source, DM, model='NE2001')
    takes an astropy SkyCoord object or parfile
    return the distance associated with source at a given DM, for the specified DM model

    if DM has no units, standard units assumed (pc/cm**3)

    return distance in kpc
    """

    d=SkyModel(dmmodel=model)
    return d.distance(source, DM)

