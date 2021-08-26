import config as cfg

import pandas as pd
import numpy as np
import pickle

from matplotlib import pyplot as plt
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.io import fits
from astropy import units as u
from photutils import CircularAperture
from photutils import aperture_photometry
from photutils.segmentation import deblend_sources
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils.segmentation import detect_sources, detect_threshold
 

from constants import FIELDS, WCS_MEMORY

gds = pd.read_pickle(f'{cfg.CANDELS_FOLDER}/gs_complete.pk')

cats = ['/home/ppxlf2/CANDELS/PANDAS/cos_complete.pk',
        '/home/ppxlf2/CANDELS/PANDAS/egs_complete.pk',
        '/home/ppxlf2/CANDELS/PANDAS/uds_complete.pk',
        '/home/ppxlf2/CANDELS/PANDAS/gn_complete.pk',
        '/home/ppxlf2/CANDELS/PANDAS/gs_complete.pk']

sources_dfs = {}
for i, field in enumerate(FIELDS):
    sources_dfs[field] = pd.read_pickle(cats[i])

"""
field606 = fits.open(f'{cfg.CANDELS_FOLDER}/hlsp_candels_hst_acs_gs-tot_f606w_v1.0_drz.fits', memmap=True)
field814 = fits.open(f'{cfg.CANDELS_FOLDER}/hlsp_candels_hst_acs_gs-tot_f814w_v1.0_drz.fits.1', memmap=True)
field125 = fits.open(f'{cfg.CANDELS_FOLDER}/hlsp_candels_hst_wfc3_gs-tot_f125w_v1.0_drz.fits', memmap=True)
field160 = fits.open(f'{cfg.CANDELS_FOLDER}/hlsp_candels_hst_wfc3_gs-tot_f160w_v1.0_drz.fits', memmap=True)  

field_data = [field606, field814, field125, field160]

wcs160 = WCS(field160[0].header)
wcs125 = WCS(field125[0].header)
wcs814 = WCS(field814[0].header)
wcs606 = WCS(field606[0].header)

wcs_fields = [wcs606, wcs814, wcs125, wcs160]

center_coord = SkyCoord(ra=53.12272274*u.degree, dec=-27.805064*u.degree)
ra_min = (center_coord.ra - 120*u.arcsec).to(u.degree).value
ra_max = (center_coord.ra + 120*u.arcsec).to(u.degree).value
dec_min = (center_coord.dec - 120*u.arcsec).to(u.degree).value
dec_max = (center_coord.dec + 120*u.arcsec).to(u.degree).value
subset = np.where((gds.DEC >= dec_min) & (gds.DEC <= dec_max) & (gds.RA >= ra_min) & (gds.RA <= ra_max) & (gds.FLAGS == 0))
RAs = gds.RA.values*u.degree
DECs = gds.DEC.values*u.degree

source_coords = SkyCoord(ra=RAs, dec=DECs)
"""

def _find_bg_section(bgs, output_size, position=None):
    
    a, b = bgs.shape
    c, d = output_size
    c = int(c)
    d = int(d)
    if position is None:
        if((c <= a) & (d <= b)):
            x = np.random.randint(0, a-c)
            y = np.random.randint(0, b-d)
            #print(x, y)
        else:
            raise(IndexError('Galaxy Image larger than BG'))
    else:
        x, y = position

    bgs_section = bgs[x:(x+c), y:(y+d)]

    return bgs_section, (x, y)

#@profile
def find_bg_section_from_CANDELS(size, bg_coord=None, field='COSMOS'):
    
    
    offset = size[-1]
    bg_sections = []
    if bg_coord is None:

        reiterate = True
        num_iters = 0
        while(reiterate):    
            bg_sections = []
            x_c = np.random.randint(offset, field160[0].data.shape[0]-offset)
            y_c = np.random.randint(offset, field160[0].data.shape[0]-offset)

            x_i = int(x_c-offset/2)
            x_f = int(x_c+offset/2)
            y_i = int(y_c-offset/2)
            y_f = int(y_c+offset/2)
            
            #bg_coord_idx = np.random.choice(np.arange(pre_selected_bg_coords.shape[0]), 1)[0]
            #bg_coord = pre_selected_bg_coords[bg_coord_idx]
            bg_coord = SkyCoord(ra=wcs160.wcs_pix2world(x_c, y_c, 0)[0]*u.degree, dec=wcs160.wcs_pix2world(x_c, y_c, 0)[1]*u.degree)

            BREAK = False
            for wcs_i in wcs_fields:
                if not bg_coord.contained_by(wcs_i):
                    BREAK = True
                    continue

            if BREAK:
                continue

            try:
                bg160 = Cutout2D(field160[0].data, position=bg_coord, size=offset*u.pixel, wcs=wcs160)
            except:
                continue

            distances = bg_coord.separation(source_coords[subset]).to(u.arcsec)/0.06
            within_reach = np.where(distances.value < offset/2)
            too_close = np.where(distances[within_reach].value <= offset/2)[0]

            #if(distances.min().value > offset/2):
            #    continue
    
            aperture = CircularAperture((offset/2, offset/2), r=offset/3)
            LT = aperture_photometry(bg160.data, aperture)['aperture_sum'][0]
            light_outside = bg160.data.sum() - LT
            #print(LT, light_outside)
            if(LT == 0 or LT > light_outside):
                continue

            FALTY_BG = False
            for field_i, size_i, wcs_i in zip(field_data, size, wcs_fields):
                cutout = Cutout2D(field_i[0].data, position=bg_coord, size=size_i*u.pixel, wcs=wcs_i) 
                bg_sections.append(cutout.data)
                if cutout.data.sum() == 0:
                    FALTY_BG = True
                    
            if FALTY_BG:
                continue

            if(len(too_close)==0):
                reiterate = False

            num_iters += 1
    else:
        for field_i, size_i, wcs_i in zip(field_data, size, wcs_fields):
            cutout = Cutout2D(field_i[0].data, position=bg_coord, size=size_i*u.pixel, wcs=wcs_i) 
            bg_sections.append(cutout.data)

        
    return bg_sections, bg_coord


from astropy.stats import gaussian_fwhm_to_sigma


def search_input_position(cutout):

    sigma = 3.0 * gaussian_fwhm_to_sigma 
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    threshold = detect_threshold(cutout.data, nsigma=2.)
    segm = detect_sources(cutout.data, threshold, npixels=5, filter_kernel=kernel)

    segm_deblend = deblend_sources(cutout.data, segm, npixels=5,
                                filter_kernel=kernel, nlevels=5,
                                contrast=0.001)

    idx = np.random.choice(np.arange(len((np.where(segm_deblend.data==0)[0]))), 1)
    x = np.where(segm_deblend.data==0)[0][idx][0]
    y = np.where(segm_deblend.data==0)[0][idx][0]

    return cutout.wcs.pixel_to_world(x, y)

def find_CANDELS_sky_patch(bands, field):

    looking_for_input = True
    while looking_for_input:

        skypatches = []
        input_coord = None
        try:
            sources = sources_dfs[field]
            source = sources.iloc[np.random.randint(0, sources.shape[0])]
            coord = SkyCoord(ra=source.RA * u.deg, dec=source.DEC* u.deg)
            
            for band, filter in zip(bands, FIELDS[field]):

                #gal_img = fits.getdata(f'/data/captain/MERGERS/SKIRT/TNG50-1/PMxSF/TESTS/{filter}/30_100453_oct_1.fits')

                wcs  = WCS_MEMORY[field][filter]#header = fits.getheader(f'/home/ppxlf2/CANDELS/{field}/{FIELDS[field][filter]["fits_filename"]}.fits')
                data = FIELDS[field][filter]['data']#fits.getdata(f'/home/ppxlf2/CANDELS/{field}/{FIELDS[field][filter]["fits_filename"]}.fits', memmap=True)
                
                #wcs = WCS(header)
                #wcs.sip = None

                if input_coord is None:

                    cutout = Cutout2D(data, position=coord, size=512, wcs=wcs)

                    if np.all(cutout.data == 0):
                        raise ValueError('Empty Patch of the Sky')

                    input_coord = search_input_position(cutout)


                cutout = Cutout2D(data, position=input_coord, size=256, wcs=wcs)

                skypatches.append(cutout.data)
        except ValueError as e:
            print(e)
            continue
        
        break

    return skypatches, field, input_coord
    
    
def find_background_sky_patch(broadbands):
    if cfg.bg_mode == 'CLEAN':
        bg_sections = _find_bg_section(BGS, broadbands)
        coord = None
    elif cfg.bg_mode == 'FIELD':
        bg_coord = None
        if cfg.pre_loaded_coords:
            coords = np.load('good_coords.npy', allow_pickle=True)
            bg_coord = coords[np.random.randint(len(coords))]
    
        sizes = []
        for b in broadbands:
            sizes.append(b.shape[0])

        bg_sections, coord = find_bg_section_from_CANDELS(sizes, bg_coord=bg_coord, field=cfg.bg_field)
    else:
        bg_sections = np.zeros_like(broadbands)
        coord = None

    return bg_sections, coord