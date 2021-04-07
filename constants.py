from astropy.io import fits
from astropy import units as u
from astropy.cosmology import LambdaCDM
import numpy as np 

_h = 0.6774
SKIRT_Z0 = 0.002265
SKIRT_OUTSIZE = 319
FOV = 60

COSMO = LambdaCDM(H0=100*_h, Om0=0.3089, Ob0=0.0486, Ode0=0.6911, Tcmb0=2.73)

ps = FOV * u.kpc * COSMO.arcsec_per_kpc_comoving(SKIRT_Z0)/SKIRT_OUTSIZE


""" UNITS """
SIM_UNITS = u.W / u.m**2 / u.arcsecond**2 * ps**2 # units of the datacube produced by SKIRT
OUT_UNITS = u.Jy * u.Hz / u.micron




SNAP_MAP = np.load('snapTNG.npy')

bg606w = fits.getdata('/home/ppxlf2/CANDELS/bgf606w_Jy.fits', memmap=True)
bg814w = fits.getdata('/home/ppxlf2/CANDELS/bgf814w_Jy.fits', memmap=True)
bg125w = fits.getdata('/home/ppxlf2/CANDELS/bgf125w_Jy.fits', memmap=True)
bg160w = fits.getdata('/home/ppxlf2/CANDELS/bgf160w_Jy.fits', memmap=True)

f606w_header = fits.getheader('/home/ppxlf2/CANDELS/hlsp_candels_hst_acs_cos-tot_f606w_v1.0_drz.fits', memmap=True)
f814w_header = fits.getheader('/home/ppxlf2/CANDELS/hlsp_candels_hst_acs_gs-tot_f814w_v1.0_drz.fits', memmap=True)
f125w_header = fits.getheader('/home/ppxlf2/CANDELS/hlsp_candels_hst_wfc3_cos-tot_f125w_v1.0_drz.fits', memmap=True)
f160w_header = fits.getheader('/home/ppxlf2/CANDELS/hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_drz.fits', memmap=True)



psf606w = fits.getdata('/home/ppxlf2/psfs/PSFs/f606w.fits')
psf814w = fits.getdata('/home/ppxlf2/psfs/PSFs/f814w.fits')
psf125w = fits.getdata('/home/ppxlf2/psfs/f125wpsf.fits')
psf160w = fits.getdata('/home/ppxlf2/psfs/f160wpsf.fits')


HEADERS = [f606w_header, f814w_header, f125w_header, f160w_header]
BGS = [bg606w, bg814w, bg125w, bg160w]
PSFS = [psf606w, psf814w, psf125w, psf160w]


REFERENCE_Z = 0.5