from constants import COSMO, SNAP_MAP, REFERENCE_Z
from astropy import constants as astroconst
from astropy import units as u
from astropy.io import fits
import numpy as np

class SkirtCube():

    def __init__(self, filename, base_path, z=None, z0=0.002265, reference_size=190,
                 FOV=60, units=None, out_units=u.Jy * u.Hz / u.micron):

        self.filename = filename
        self.snap, self.subfind, self.orientation = self.filename.split('_')
        self.base_path = base_path
        self.cube_path = f'{base_path}/DATACUBE/{filename}_total.fits'
        self.sed_path = f'{base_path}/SED/{filename}_sed.dat'
        self.z0 = z0
        self.reference_size = reference_size
        self.FOV = FOV
        self.out_units = out_units
        self.pixel_size = self.FOV * u.kpc * COSMO.arcsec_per_kpc_proper(self.z0)/self.reference_size

        if units is None:
            self.units = u.W / u.m**2 / u.arcsecond**2 * self.pixel_size**2
        else:
            self.units = units

        if z is None:
            self.z = SNAP_MAP[1][np.where(SNAP_MAP[0] == int(self.snap))][0]
        else:
            self.z = z

        if self.z < REFERENCE_Z:
            self.target_z = REFERENCE_Z
        else:
            self.target_z = self.z

        #self.dimming = 

        self._load()
            

    def _load(self):

        self.wav, self.flux = np.loadtxt(self.sed_path).T
        self.wav *= u.micron
        self.wav_z = self.wav * (1+self.target_z)
        self.flux *= u.W / u.m**2

        with fits.open(self.cube_path, mode='readonly') as hdul:
            self.header = hdul[0].header
            self.data = hdul[0].data.T * self.units
            self.data = (self.data / self.wav).to(self.out_units).T #flux to flux density
    
    
    def integrate_filter(self, band):
        """
            Based on code from Connor Bottrell on RealSim for SDSS.
            https://github.com/cbottrell/RealSim
        """
        
        filter_wav, filter_res = np.loadtxt(band).T
        filter_wav = (filter_wav * u.AA).to(u.micron)

        dwl = np.median(np.diff(self.wav_z))

        filter_res = np.interp(xp=filter_wav, x=self.wav_z, fp=filter_res, left=0,right=0)
        wl_pivot2 = np.sum(filter_res*self.wav_z*dwl)/np.sum(filter_res*dwl/self.wav_z)
        
        f_wl = np.sum(self.wav_z*filter_res*self.data.T*dwl,axis=2)/np.sum(self.wav_z*filter_res*dwl)
        f_jy = ((f_wl * wl_pivot2 / astroconst.c)).to(u.Jy)

        return f_jy