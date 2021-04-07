from areia import areia

from astropy.cosmology import LambdaCDM
from constants import *
from background import find_background_sky_patch
from SkirtCube import SkirtCube

import  h5py
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import h5py
from PIL import Image
import matplotlib
matplotlib.use('Agg')


from astropy.io import fits
from astropy import constants as const
from scipy.interpolate import interp1d
from scipy.integrate import simps
import pyphot
from astropy.convolution import convolve

from scipy import ndimage as ndi
from scipy.stats import median_abs_deviation
from pyphot import Sun, unit

from utils import *

import sys  
import glob
import os


import config as cfg

if __name__ == '__main__':

    """
        TNG Data is stored as {snapshot}_{subfindID}. Example: 99_234510
    """
    rootname = sys.argv[1]

    for orientation in cfg.ALLOWED_ORIENTATIONS:
        filename = f'{rootname}_{orientation}'

        base_folder = cfg.SOURCE_FOLDER
        if not os.path.exists(f'{base_folder}/DATACUBE/{filename}_total.fits'):
            base_folder = cfg.SOURCE_FOLDER_FB

        cube = SkirtCube(filename, base_folder)
    
        dimming = ((cosmo.comoving_distance(SKIRT_Z0) / cosmo.comoving_distance(cube.target_z))**2).value # Dimming calculated from the ratio of the luminosities of target and source frames

        broadbands = [cube.integrate_filter(band) * dimming  for band in cfg.FILTERS]  
        broadbands = [to_eps(band, header).value for band, header in zip(broadbands, HEADERS)]  

        bg_sections, coords = find_background_sky_patch(broadbands)
   
        restframe = areia.ObservationFrame(REFERENCE_Z, pixelscale=0.03, exptime=1)

        if cfg.restframe:
            obs_frames = [restframe for tps in cfg.TARGET_PS]
        else:
            obs_frames = [areia.ObservationFrame(cube.target_z, pixelscale=tps, exptime=exptime) for tps, exptime in zip(cfg.TARGET_PS, cfg.EXPTIME)]

        final = [areia.ArtificialRedshift(band, psf, bg, restframe, obs_frame, MAG=None, config=cfg.areia_config).final for band, psf, bg, obs_frame in zip(broadbands, PSFS, bg_sections, obs_frames)]

        for i, band in enumerate(cfg.BANDS):
            output_path = f'{cfg.OUTPUT_FOLDER}/{band}/{rootname}_{orientation}.fits'
            dirname = os.path.dirname(output_path)
            if not os.path.exists(dirname):
                os.makedirs(dirname)
        
            fits.writeto(output_path, final[i], overwrite=True)
        


"""        if cfg.zscale:
            d1 = cosmo.luminosity_distance(REFERENCE_Z).value
            d2 = cosmo.luminosity_distance(z).value
            broadbands = [z_rescale(b, REFERENCE_Z, z, d1, d2, 0.03, tps) for b, tps in zip(broadbands, cfg.TARGET_PS)]

              final = []
        for i, (facei, bg, psf, header, exptime) in enumerate(zip(broadbands, bg_sections, PSFS, HEADERS, cfg.EXPTIME)):


            if cfg.convolve:
                facei = convolve(facei.value, psf) * u.Jy
            
            #convert from Jy to e/s
            ergs = facei.to(u.erg / u.cm**2 / u.AA / u.s, equivalencies=u.spectral_density(u.AA, header['PHOTPLAM']))
            cts = ergs / (header['PHOTFLAM'] * u.erg / u.cm**2 / u.AA / u.ct) 
            
            #add error and back to Jy
            if cfg.shotnoise:
                #img = (((cts * u.s / u.ct) + np.sqrt(abs(cts*exptime * u.s / u.ct))*np.random.randn(cts.shape[0], cts.shape[1])/exptime) * (header['PHOTFLAM'] * u.erg / u.cm**2 / u.AA / u.s)).to(u.Jy, equivalencies=u.spectral_density(u.AA, header['PHOTPLAM']))
                img = cts * u.s / u.ct + np.sqrt(abs(cts*exptime * u.s / u.ct))*np.random.randn(cts.shape[0], cts.shape[1])/exptime
                #print(img)
            else:
                img = ergs

            if cfg.bg:
                size = facei.shape[0]
                candels_face = bg[0:int(size)*2, 0:int(size)*2].copy()
                bg_size = candels_face.shape[0]
                xi = int(bg_size/2 - size/2)
                xf = int(bg_size/2 + size/2)
                
                
                img += candels_face[xi:xf, xi:xf] 

                
                #img = img.to(u.erg / u.cm**2 / u.AA / u.s, equivalencies=u.spectral_density(u.AA, header['PHOTPLAM']))
                #img = img / (header['PHOTFLAM'] * u.erg / u.cm**2 / u.AA / u.ct) 
            
        
            final.append(img)
"""        
        