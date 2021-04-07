import sys  
import glob
import os
import matplotlib
matplotlib.use('Agg')

import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import config as cfg

from astropy.io import fits
from areia import areia
from constants import *
from background import find_background_sky_patch
from SkirtCube import SkirtCube
from utils import *


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