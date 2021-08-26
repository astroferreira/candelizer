import sys  
import glob
import os
import matplotlib
matplotlib.use('Agg')

import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import config as cfg

from sklearn.model_selection import ParameterGrid

from astropy.io import fits
from areia import areia
from constants import *
from background import find_background_sky_patch, _find_bg_section, find_CANDELS_sky_patch
from SkirtCube import SkirtCube
from scipy.ndimage import zoom
from utils import *

from pixelstatistics import segmentate

import logging
logging.getLogger("astropy").setLevel(logging.WARNING)


if __name__ == '__main__':

    """
        TNG Data is stored as {snapshot}_{subfindID}. Example: 99_234510
    """
    rootname = sys.argv[1]

    seed = int(''.join(rootname.split('_')))
    
    
    """
        Create a ParameterGrid to avoid 3 nested loops
        The runs are independent on FIELD/REDSHIFT/ORIENTATION combination.
    """
    parameter_space = ParameterGrid({
                        'FIELDS': list(FIELDS.keys()),
                        'orientations' : cfg.ALLOWED_ORIENTATIONS,
                        'redshifts' : cfg.TARGET_REDSHIFTS
                      })

    for i, params in enumerate(parameter_space):

        seed += i*787
        np.random.seed(seed)
        
        field = params['FIELDS']
        orientation = params['orientations']
        target_redshift = params['redshifts']

        print(f'Running {rootname} in {field}/{orientation}')
        
        primary_header = fits.Header()
    
        filename = f'{rootname}_{orientation}'

        primary_header['ROOTNAME'] = rootname
        primary_header['FILENAME'] = filename
        primary_header['NP_SEED'] = seed
        primary_header['SKIRTORI'] = orientation
        primary_header['CANDELSF'] = field
        primary_header['TARGET_Z'] = target_redshift
        primary_header['BANDS'] = ','.join(cfg.BANDS)

        output_headers = {}
        for band in cfg.BANDS:
            hdr = fits.Header()
            hdr['BAND'] = band
            output_headers[band] = hdr

        cube = SkirtCube(filename, cfg.SOURCE_FOLDER)

        if target_redshift != 'full':
            cube.target_z = target_redshift

        dimming = (1/(1+cube.target_z))*((cosmo.luminosity_distance(SKIRT_Z0) / cosmo.luminosity_distance(cube.target_z))**2).value # Dimming calculated from the ratio of the luminosities of target and source frames
        cube.data *= dimming 

        primary_header['DIMMFACT'] = dimming

        broadbands = [cube.integrate_filter(band)  for band in cfg.FILTERS]  
        broadbands = [to_eps(band, header).value    for band, header in zip(broadbands, cfg.HEADERS)]  

        if (cfg.fix_scaling) & (cube.z < REFERENCE_Z):
            correct_scale = COSMO.scale_factor(0.5)/COSMO.scale_factor(cube.z)

            for i, band in enumerate(broadbands):
                flux = band.sum()
                band = zoom(band, correct_scale)
                band = (band/band.sum())*flux
                broadbands[i] = band


        if cfg.bg:
            bg_sections, _, coord = find_CANDELS_sky_patch(broadbands, field)
            primary_header['RA_BG'] = coord.ra.value
            primary_header['DEC_BG'] = coord.dec.value
        else:
            bg_sections = [None] * len(cfg.FILTERS)

        restframe = areia.ObservationFrame(REFERENCE_Z, pixelscale=0.05, exptime=1)

        if cfg.restframe:
            obs_frames = [restframe for tps in cfg.TARGET_PS]
        else:
            obs_frames = [areia.ObservationFrame(cube.target_z, pixelscale=tps, exptime=exptime) for tps, exptime in zip(cfg.TARGET_PS, cfg.EXPTIME)]

        result = {}
        for band, psf, bg, obs_frame, bandname in zip(broadbands, cfg.PSFS, bg_sections, obs_frames, cfg.BANDS):
            AR = areia.ArtificialRedshift(band, psf, bg, restframe, obs_frame, MAG=None, bg_position=None, config=cfg.areia_config)
            result[bandname] = {}
            result[bandname]['final'] = AR.final_crop
            result[bandname]['bg'] = AR.background
            result[bandname]['clean'] = AR.convolved

        for i, band in enumerate(cfg.BANDS):
            keyword = f'BG{band}'
            
            
            segmap, flux_per_pixel, flux, num_pix = segmentate(result[band]['bg'])
            galsegmap, _, _, _ = segmentate(result[band]['clean'])
            
            overlap_index = np.where((segmap >= 1) & (galsegmap >= 1))
            overlapping = overlap_index[0].shape[0]

            overlap_map = np.zeros_like(galsegmap)
            overlap_map[overlap_index] = 1

            coverage = np.round(100 * overlapping / len(galsegmap[galsegmap >= 1]), 2)

            output_headers[band]['FLUX_PP'] = flux_per_pixel
            output_headers[band]['FLUX_INT'] = flux
            output_headers[band]['NUM_PIX'] = num_pix
            output_headers[band]['OVERLAP'] = overlapping
            output_headers[band]['COVERAGE'] = coverage

            result[band]['segmap'] = segmap
            result[band]['clean_segmap'] = galsegmap
            result[band]['overlap_segmap'] = overlap_map
        
        if cfg.shotnoise:
            for filter in FIELDS[field]:
                result[filter]['final'] = result[filter]['final'] + np.random.normal(0, FIELDS[field][filter]['skystd'], size=result[filter]['clean'].shape)
                result[filter]['clean'] = result[filter]['clean'] + np.random.normal(0, FIELDS[field][filter]['skystd'], size=result[filter]['clean'].shape)

        primary_header = config_to_header(primary_header)

        hdul = dict_to_hdu(fits.PrimaryHDU(header=primary_header), result, output_headers)

        output_path = f'{cfg.OUTPUT_FOLDER}/{rootname}_{orientation}_{target_redshift}_{field}.fits'

        dirname = os.path.dirname(output_path)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        
        hdul.writeto(output_path, overwrite=True)

        #for i, (band, hdr) in enumerate(zip(cfg.BANDS, output_headers)):
        #    output_path = f'{cfg.OUTPUT_FOLDER}/{band}/{rootname}_{orientation}_{target_redshift}_{field}.fits'
        #    dirname = os.path.dirname(output_path)
        #    if not os.path.exists(dirname):
        #        os.makedirs(dirname)
        #
        #    fits.writeto(output_path, final[i], overwrite=True)







