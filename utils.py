from astropy.cosmology import LambdaCDM

import  h5py
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from PIL import Image
import matplotlib
matplotlib.use('Agg')

h = 0.6774

from astropy.io import fits
from astropy import constants as const
from scipy.interpolate import interp1d
from scipy.integrate import simps
import pyphot
from astropy.convolution import convolve
from astropy import constants as const

from scipy import ndimage as ndi
from scipy.stats import median_abs_deviation
from pyphot import Sun, unit
lib = pyphot.get_library()

import sys  
import glob

skirt_z0 = 0.002265
cosmo = LambdaCDM(H0=100*h, Om0=0.3089, Ob0=0.0486, Ode0=0.6911, Tcmb0=2.73)

ps = 60*u.kpc * cosmo.arcsec_per_kpc_comoving(0.002265)/319 /u.pix

def to_Jansky(l, cube):
    return (cube.T * (ps.value**2 * u.W / u.m**2 / (const.c / (l*u.micron).to(u.m) ) ).to(u.Jy)).T

def to_ergs(l, cube):
    return (cube.T * (ps.value**2 * u.W / u.m**2 /  (u.AA)).to(u.erg / u.cm**2 / u.AA / u.s)).T
    
def to_eps(band, header):
    band = band.to(u.erg / u.cm**2 / u.AA / u.s, equivalencies=u.spectral_density(u.AA, header['PHOTPLAM']))
    return band / (header['PHOTFLAM'] * u.erg / u.cm**2 / u.AA / u.ct) 


def Jy_to_eps(cube, header):
    cube_in_ergs = cube * u.Jy.to(u.erg/u.cm**2 / u.AA / u.s, equivalencies=u.spectral_density(u.AA, header['PHOTPLAM']))
    return cube_in_ergs / header['PHOTFLAM'] # [e/s]

def z_rescale(data, z1, z2, d1, d2, ps1, ps2):
    '''                                                      
        Rescale data from z1 to z2, and from pixel scale ps1 
        to pixel scale ps2. This simulates the change of resolution 
        if a galaxy at z1 would be observed at z2.                 
                                                                   
        http://esoads.eso.org/abs/2008ApJS..175..105B              
    '''                                                            
    scale_factor = (d1 * (1+z2) * ps1) / (d2 * (1+z1) * ps2)                                
   
    flux = data.sum()
    rescaled = ndi.zoom(data, scale_factor, prefilter=True)
    rescaled /= rescaled.sum()
    return rescaled * flux 

    def to_Jansky(l, cube):
        return (cube.T * (ps.value**2 * u.W / u.m**2 / (const.c / (l*u.micron).to(u.m) ) ).to(u.Jy)).T.value

def to_ergs(l, cube):
    return (cube.T * (ps.value**2 * u.W / u.m**2 /  (l*u.micron)).to(u.erg / u.cm**2 / u.AA / u.s)).T

def fix_size(img, max_size):   
    if img.shape[0] < max_size:
        imgflux = img.sum()
        img = ndi.zoom(img, max_size/img.shape[0])
        img /= img.sum()
        img *= imgflux
    
    return img



def plot_RGB(img):

    import numpy as np
    from astropy.io import fits
    from astropy.visualization import SqrtStretch, LogStretch
    from astropy.visualization import LinearStretch, AsinhStretch
    from astropy.visualization import make_lupton_rgb

    forCasting = np.float_()
    
    g = img[0]
    r = img[1]
    i = img[2]

    max_size = np.max([g.shape[0], r.shape[0], i.shape[0]])

    image_b = fix_size(g, max_size)
    image_g = fix_size(r, max_size)
    image_r = fix_size(i, max_size)


    norm = np.nanpercentile(image_b,99)
    image_r/=norm
    image_g/=norm
    image_b/=norm

    image = make_lupton_rgb(image_r, image_g, image_b, Q=8, stretch=0.5,minimum=image_b.min())
    # Cut the top rows - contains black pixels

    return image





def integrate_cube(cube, wl, band, header):
    
    # useful constants / quantities
    speed_of_light = 2.998e8 # [m/s]
    speed_of_light = speed_of_light*1e10 # [Angstrom/s]

    filter_wl, filter_res = np.loadtxt(band).T
    filter_wl = (filter_wl * u.AA).to(u.micron)

    dwl = np.median(np.diff(wl))

    filter_res = np.interp(xp=filter_wl, x=wl, fp=filter_res, left=0,right=0)
    wl_pivot2 = np.sum(filter_res*wl*dwl)/np.sum(filter_res*dwl/wl)
    
    f_wl = np.sum(wl*filter_res*cube.T*dwl,axis=2)/np.sum(wl*filter_res*dwl)
    f_jy = (f_wl * wl_pivot2 / const.c).to(u.Jy)

    return f_jy

def skybg_iterative(data, max_iter=100, old_med=0, old_mad=0 ):
    '''
    1. Defines initial median and mad based on mask (initially the whole image). 
    2. Select those pixels below the limit;
    3. calculates new median and mad on 2.
    4. repeat until median and mad converges or max_iter is reached
    
    Based on idea and code from Leonardo Ferreira, 2015. 
    '''

    skymed = np.median(data)
    skymad = median_abs_deviation(data, axis=None)
    if max_iter==0 or (skymed==old_med and skymad==old_mad):
        return (skymed, skymad, data.mean(), data.std())
    else:
        mask = (data < skymed + 3.*skymad)
        return skybg_iterative(data[mask], max_iter-1, skymed, skymad)

def dict_to_hdu(primary, dict, hdrs):

    hdul = [primary]
    for band in list(dict.keys()):
        
        stacked = np.dstack([dict[band]['final'],
                                dict[band]['clean'], 
                                dict[band]['bg'], 
                                dict[band]['segmap'],
                                dict[band]['clean_segmap'],
                                dict[band]['overlap_segmap']])

        hdul.append(fits.ImageHDU(stacked, header=hdrs[band]))

    return fits.HDUList(hdul)

import config as cfg
def config_to_header(hdr):

    hdr['F_RFRAME'] = str(cfg.restframe)
    hdr['F_SHOTN'] = str(cfg.shotnoise)
    hdr['F_BG'] = str(cfg.bg)
    hdr['F_CONV'] = str(cfg.convolve)
    hdr['F_ZSCALE'] = str(cfg.zscale)
    hdr['F_FSCALE'] = str(cfg.fix_scaling)
    hdr['F_BGMODE'] = cfg.bg_mode

    hdr['FA_h'] = float(cfg.areia_config.h)
    hdr['FA_DIM'] = str(cfg.areia_config.dimming)
    hdr['FA_SCORR'] = str(cfg.areia_config.size_correction)
    hdr['FA_ZEVO'] = str(cfg.areia_config.size_correction)
    hdr['FA_OSIZE'] = int(cfg.areia_config.output_size)

    return hdr
    