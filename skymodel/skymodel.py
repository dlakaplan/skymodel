import numpy as np
from astropy.coordinates import SkyCoord

from astropy import units as u
import astropy
from optparse import OptionParser
import os,sys
try:    
    import pint.models as models
    _usePINT=True
except ImportError:
    _usePINT=False

try:
    from pygsm import GlobalSkyModel,GlobalSkyModel2016
    import healpy
    _usePyGSM=True
except ImportError:
    _usePyGSM=False
    
from . import ne2001
from . import ymw16

from . import data

def parfile2SkyCoord(parfile):
    """
    psr=parfile2SkyCoord(parfile)
    uses PINT to interpret a parfile
    and return a SkyCoord object
    """
    
    if not _usePINT:
        raise ImportError('PINT is not available: cannot use parfiles')
    m=models.get_model(parfile)
    return SkyCoord(m.get_psr_coords())


"""
dmdtau_c usage:
result=ymw16.dmdtau_c(l, b, dordm, ndir, dirname)

l,b: Galactic coordinates in degrees
dordm is input coordinate (distance in pc)
ndir=1 for DM -> Distance
ndir=2 for Distance -> DM
dirname is location of data files

returns dmord

"""                      


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
    
    tskymodel can be '2008' or '2016'
    dmmodel can be 'NE2001' or 'YMW16'

    if DM/d is a lower limit, returned value will be < 0 (NE2001-only)

    """

    def __init__(self, freq=350*u.MHz, tskymodel='2008', dmmodel='NE2001', datadir=data._datadir):
        """
        m=SkyModel(freq=350*u.MHz, tskymodel='2008', dmmodel='NE2001', datadir=data._datadir)
        """

        assert dmmodel.lower() in ['ne2001','ymw16']
        assert str(tskymodel) in ['2008','2016']


        self.datadir=datadir
        self.dmmodel=dmmodel.lower()
        if not isinstance(freq, astropy.units.quantity.Quantity):
            # assume MHz
            freq=freq*u.MHz

        if _usePyGSM:
            if str(tskymodel)=='2008':
                self.gsm = GlobalSkyModel()
            elif str(tskymodel)=='2016':
                self.gsm=GlobalSkyModel2016()
            self.map=self.gsm.generate(freq.to(u.MHz).value)
            self.tskymodel=tskymodel
        else:
            self.tskymodel=None
        self.freq=freq

    ##################################################
    def Tsky(self, source):
        """
        T=m.Tsky(source)

        returns sky temperature in K for the given model
        given a source, either a SkyCoord object or a parfile
        """

        if not _usePyGSM:
            raise ImportError('PyGSM is not available: cannot access sky temperatures')
        if not isinstance(source, astropy.coordinates.sky_coordinate.SkyCoord):
            if isinstance(source,str):
                # assume .par file
                source=parfile2SkyCoord(source)
            else:
                raise TypeError('Do not know how to interpret an object of type %s' % source.__class__)

        source=source.galactic
        T=healpy.pixelfunc.get_interp_val(self.map,
                                          source.l.value,
                                          source.b.value,
                                          lonlat=True)
        return T*u.K
    ##############################
    def DM(self, source, distance, smweight='uniform'):
        """
        DM,SM=m.DM(source, distance, smweight='uniform')

        returns DM in standard units and scattering measure (NE2001 only, else None) for the specified model
        given a source, either a SkyCoord object or a parfile
        and a distance

        scattering measure can be weighted by:
        'uniform': uniform weighting [default]
        'tau': weighted for pulse broadening
        'theta': weighted for angular broadening of galactic sources
        'iso':  appropriate for calculating the isoplanatic angle at the source's location


        if no units supplied kpc assumed for distance
        """

        if self.dmmodel == 'ne2001':
            return self.DM_NE2001(source, distance, smweight=smweight)
        elif self.dmmodel == 'ymw16':
            return self.DM_YMW16(source, distance)


    ##############################
    def DM_NE2001(self, source, distance, smweight='uniform'):
        """
        DM,SM=m.DM_NE2001(source, distance, smweight='uniform')

        returns DM in standard units and scattering measure for NE2001
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
            distance=distance*u.kpc            
        if (len(distance.shape)>0 and distance.value.any() <= 0) or (len(distance.shape)==0 and distance.value < 0):
            raise ValueError('distance must be > 0')

        if not isinstance(source, astropy.coordinates.sky_coordinate.SkyCoord):
            if isinstance(source,str):
                # assume .par file
                source=parfile2SkyCoord(source)
            else:
                raise TypeError('Do not know how to interpret an object of type %s' % source.__class__)
        source=source.galactic


        if len(source.l.shape)==0:
            results=ne2001.dmdsm(self.datadir,
                                 np.radians(source.l.value),
                                 np.radians(source.b.value),
                                 -1,
                                 0,
                                 distance.to(u.kpc).value)
            sign=1
            if results[2]=='>':
                #raise ValueError('DM returned a lower limit')
                sign=-1
            if smweight.lower() == 'uniform':
                SM=results[3]*u.kpc/u.m**(20./3)
            elif smweight.lower() == 'tau':
                SM=results[4]*u.kpc/u.m**(20./3)
            elif smweight.lower() == 'theta':
                SM=results[5]*u.kpc/u.m**(20./3)
            elif smweight.lower() == 'iso':
                SM=results[6]*u.kpc/u.m**(20./3)

            return sign*results[0]*u.pc/u.cm**3,SM
        else:
            dm=np.zeros_like(source.l.value)
            SM=np.zeros_like(source.l.value)
            it = np.nditer(source.l, flags=['multi_index'])
            if len(dm.shape)==0:
                dm_touse=dm
            else:
                dm_touse=dm[it.multi_index]
            while not it.finished:
                if len(distance.shape)==0:
                    d_touse=distance
                else:
                    d_touse=distance[it.multi_index]
                results=ne2001.dmdsm(self.datadir,
                                     np.radians(source[it.multi_index].l.value),
                                     np.radians(source[it.multi_index].b.value),
                                     -1,
                                     0,
                                     d_touse.to(u.kpc).value)
                sign=1
                if results[2]=='>':
                    #raise ValueError('DM returned a lower limit')
                    sign=-1
                dm[it.multi_index]=results[0]*sign
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
   ##############################
    def DM_YMW16(self, source, distance):
        """
        DM,None=m.DM_YMW16(source, distance)

        returns DM in standard units for YMW16
        given a source, either a SkyCoord object or a parfile
        and a distance

        if no units supplied kpc assumed for distance
        """


        if not isinstance(distance, astropy.units.quantity.Quantity):
            # assume kpc
            distance=distance*u.kpc            
        if (len(distance.shape)>0 and distance.value.any() <= 0) or (len(distance.shape)==0 and distance.value < 0):
            raise ValueError('distance must be > 0')

        if not isinstance(source, astropy.coordinates.sky_coordinate.SkyCoord):
            if isinstance(source,str):
                # assume .par file
                source=parfile2SkyCoord(source)
            else:
                raise TypeError('Do not know how to interpret an object of type %s' % source.__class__)
        source=source.galactic


        if len(source.l.shape)==0:
            results=ymw16.dmdtau_c(source.l.value,
                                   source.b.value,
                                   distance.to(u.pc).value,
                                   2,
                                   self.datadir)

            return results*u.pc/u.cm**3,None
        else:
            dm=np.zeros_like(source.l.value)
            it = np.nditer(source.l, flags=['multi_index'])
            if not (len(distance.shape)==0 or distance.shape==source.l.shape):
                raise IndexError('Shape of distance must be scalar or the same as shape of coordinates')
            d=distance.to(u.pc).value
            while not it.finished:
                if len(d.shape)==0:
                    d_touse=d
                else:
                    d_touse=d[it.multi_index]
                results=ymw16.dmdtau_c(source[it.multi_index].l.value,
                                       source[it.multi_index].b.value,
                                       d_touse,
                                       2,
                                       self.datadir)
                    
                dm[it.multi_index]=results
                it.iternext()
                
            return dm*u.pc/u.cm**3,None
            
    ##################################################
    def distance(self, source, DM, smweight='uniform'):
        """
        d,SM=m.distance(source, DM, smweight='uniform')

        returns distance in kpc and scattering measure (NE2001 only, else None) for the specified model
        given a source, either a SkyCoord object or a parfile
        and a DM

        scattering measure can be weighted by:
        'uniform': uniform weighting [default]
        'tau': weighted for pulse broadening
        'theta': weighted for angular broadening of galactic sources
        'iso':  appropriate for calculating the isoplanatic angle at the source's location

        if no units supplied standard DM units assumed for DM
        """
        
        if self.dmmodel == 'ne2001':
            return self.distance_NE2001(source, DM, smweight=smweight)
        elif self.dmmodel == 'ymw16':
            return self.distance_YMW16(source, DM)


    ##################################################
    def distance_NE2001(self, source, DM, smweight='uniform'):
        """
        d,SM=m.distance_NE2001(source, DM, smweight='uniform')

        returns distance in kpc and scattering measure for NE2001
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
            DM=DM*u.pc/u.cm**3
        if (len(DM.shape)>0 and DM.value.any() <= 0) or (len(DM.shape)==0 and DM.value < 0):
            raise ValueError('DM must be > 0')
        if not isinstance(source, astropy.coordinates.sky_coordinate.SkyCoord):
            if isinstance(source,str):
                # assume .par file
                source=parfile2SkyCoord(source)
            else:
                raise TypeError('Do not know how to interpret an object of type %s' % source.__class__)
        source=source.galactic

        if len(source.l.shape)==0:
            
            results=ne2001.dmdsm(self.datadir,
                                 np.radians(source.l.value),
                                 np.radians(source.b.value),
                                 1,
                                 DM.to(u.pc/u.cm**3).value,
                                 0)
            sign=1
            if results[2]=='>':
                #raise ValueError('distance returned a lower limit')
                sign=-1
            distance=results[1]*u.kpc*sign
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
            distance=np.zeros_like(source.l.value)
            SM=np.zeros_like(source.l.value)
            it = np.nditer(source.l, flags=['multi_index'])
            dm=DM.to(u.pc/u.cm**3).value
            if not (len(dm.shape)==0 or dm.shape==source.l.shape):
                raise IndexError('Shape of DM must be scalar or the same as shape of coordinates')
            while not it.finished:
                if len(dm.shape)==0:
                    dm_touse=dm
                else:
                    dm_touse=dm[it.multi_index]
                results=ne2001.dmdsm(self.datadir,
                                     np.radians(source[it.multi_index].l.value),
                                     np.radians(source[it.multi_index].b.value),
                                     1,
                                     dm_touse,
                                     0)
                sign=1
                if results[2]=='>':
                    #raise ValueError('distance returned a lower limit')
                    sign=-1
                distance[it.multi_index]=results[1]*sign
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
    ##################################################
    def distance_YMW16(self, source, DM):
        """
        d=m.distance_YMW16(source, DM)

        returns distance in kpc for YMW16
        given a source, either a SkyCoord object or a parfile
        and a DM

        if no units supplied standard DM units assumed for DM
        """

        if not isinstance(DM, astropy.units.quantity.Quantity):
            # assume DM unit
            DM=DM*u.pc/u.cm**3
        if (len(DM.shape)>0 and DM.value.any() <= 0) or (len(DM.shape)==0 and DM.value < 0):
            raise ValueError('DM must be > 0')
        if not isinstance(source, astropy.coordinates.sky_coordinate.SkyCoord):
            if isinstance(source,str):
                # assume .par file
                source=parfile2SkyCoord(source)
            else:
                raise TypeError('Do not know how to interpret an object of type %s' % source.__class__)
        source=source.galactic

        if len(source.l.shape)==0:
            
            results=ymw16.dmdtau_c(source.l.value,
                                   source.b.value,                                   
                                   DM.to(u.pc/u.cm**3).value,
                                   1,
                                   self.datadir)
            distance=results*u.pc
            return distance,None
        else:
            distance=np.zeros_like(source.l.value)
            it = np.nditer(source.l, flags=['multi_index'])
            dm=DM.to(u.pc/u.cm**3).value
            if not (len(dm.shape)==0 or dm.shape==source.l.shape):
                raise IndexError('Shape of DM must be scalar or the same as shape of coordinates')
            while not it.finished:
                if len(dm.shape)==0:
                    dm_touse=dm
                else:
                    dm_touse=dm[it.multi_index]
                results=ymw16.dmdtau_c(source[it.multi_index].l.value,
                                       source[it.multi_index].b.value,                                   
                                       dm_touse,
                                       1,
                                       self.datadir)
                distance[it.multi_index]=results
                it.iternext()
            return distance*u.pc,None



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
            raise TypeError('Do not know how to interpret an object of type %s' % source.__class__)


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

