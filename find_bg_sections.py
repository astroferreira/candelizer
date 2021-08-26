"""
from background import find_bg_section_from_CANDELS
import numpy as np

sizes = [int(319*1.5), int(319*1.5), 319, 319]
coords = []
bg_coord=None
for i in range(1000):
    bg_sections, coord = find_bg_section_from_CANDELS(sizes, bg_coord=bg_coord, field='GOODS')
    coords.append(coord)
    print(coord, i)

np.save('coords.npy', np.array(coords), allow_pickle=True)
"""
import matplotlib
matplotlib.use('Agg')

class Config(): 

    SOURCE_FOLDER = '/home/ppxlf2/MERGERS/SKIRT/TNG50-1/HALO/'
    SOURCE_FOLDER_FB = '/home/ppxlf2/MERGERS/SKIRT/TNG50-1/SUBHALO/'
    OUTPUT_FOLDER = '/home/ppxlf2/MERGERS/SKIRT/TNG50-1/AMY_LARGE/'
    FILTER_FOLDER = '/home/ppxlf2/sources/skirt_to_candels/filters/'
    CANDELS_FOLDER = '/home/ppxlf2/BGS'

    FILTERS = [f'{FILTER_FOLDER}/HST_ACS_WFC.F606W.dat', 
            f'{FILTER_FOLDER}/HST_ACS_WFC.F814W.dat', 
            f'{FILTER_FOLDER}/HST_WFC3_IR.F125W.dat', 
            f'{FILTER_FOLDER}/HST_WFC3_IR.F160W.dat']

    BANDS = ['F606W', 'F814W', 'F125W', 'F160W']
    EXPTIME = [3049,5761, 1963, 3480]

    TARGET_PS = [0.03, 0.03, 0.06, 0.06]

    ALLOWED_ORIENTATIONS = ['xy2', 'xz', 'yz', 'oct']



    from areia.areia import Config as AreiaConfig

    #FLAGS
    restframe=False
    shotnoise=True
    bg=True
    convolve=True
    pre_loaded_coords = False
    zscale = True
    bg_mode = 'FIELD' # CLEAN, FIELD
    bg_field = 'COSMOS'
    
cfg = Config()

import pandas as pd
import numpy as np

from matplotlib import pyplot as plt
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.io import fits
from astropy import units as u
from photutils import CircularAperture
from photutils import aperture_photometry
 

gds = pd.read_pickle(f'{cfg.CANDELS_FOLDER}/gs_complete.pk')

field606 = fits.open(f'{cfg.CANDELS_FOLDER}/hlsp_candels_hst_acs_gs-tot_f606w_v1.0_drz.fits', memmap=True)
field814 = fits.open(f'{cfg.CANDELS_FOLDER}/hlsp_candels_hst_acs_gs-tot_f814w_v1.0_drz.fits.1', memmap=True)
field125 = fits.open(f'{cfg.CANDELS_FOLDER}/hlsp_candels_hst_wfc3_gs-tot_f125w_v1.0_drz.fits', memmap=True)
field160 = fits.open(f'{cfg.CANDELS_FOLDER}/hlsp_candels_hst_wfc3_gs-tot_f160w_v1.0_drz.fits', memmap=True)  

field_data = [field814, field125, field160]

wcs160 = WCS(field160[0].header)
wcs125 = WCS(field125[0].header)
wcs814 = WCS(field814[0].header)
wcs606 = WCS(field606[0].header)

wcs_fields = [wcs814, wcs125, wcs160]

center_coord = SkyCoord(ra=53.12272274*u.degree, dec=-27.805064*u.degree)
ra_min = (center_coord.ra - 120*u.arcsec).to(u.degree).value
ra_max = (center_coord.ra + 120*u.arcsec).to(u.degree).value
dec_min = (center_coord.dec - 120*u.arcsec).to(u.degree).value
dec_max = (center_coord.dec + 120*u.arcsec).to(u.degree).value
subset = np.where((gds.DEC >= dec_min) & (gds.DEC <= dec_max) & (gds.RA >= ra_min) & (gds.RA <= ra_max) & (gds.FLAGS == 0))
RAs = gds.RA.values*u.degree
DECs = gds.DEC.values*u.degree

source_coords = SkyCoord(ra=RAs, dec=DECs)

def _find_bg_section(bgs, output_size, position=None):
    
    a, b = bgs.shape
    c, d = output_size
    c = int(c*1.5)
    d = int(d*1.5)
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


    
def find_background_sky_patch(broadbands):
    if cfg.bg_mode == 'CLEAN':
        bg_sections = _find_bg_section(BGS, broadbands)
        coord = None
    elif cfg.bg_mode == 'FIELD':
        bg_coord = None
        if cfg.pre_loaded_coords:
            coords = np.load('coords.npy', allow_pickle=True)
            bg_coord = coords[np.random.randint(1000)]
    
        sizes = []
        for b in broadbands:
            sizes.append(broadbands[0].shape[0]*2)

        bg_sections, coord = find_bg_section_from_CANDELS(sizes, bg_coord=bg_coord, field=cfg.bg_field)
    else:
        bg_sections = np.zeros_like(broadbands)
        coord = None

    return bg_sections, coord


selected_coords = np.load('coords.npy', allow_pickle=True)

keeps = []
for index, coord in enumerate(selected_coords):
    bg_sections, _ = find_bg_section_from_CANDELS([(1000, 1000),(500, 500),(500, 500)], bg_coord=coord)
    f, axs = plt.subplots(1, 3, figsize=(600/80, 200/80), dpi=80)
    for bg, axi in zip(bg_sections, axs):
        axi.imshow(np.log10(bg))
    
    print(index)
    plt.savefig(f'BGS/{index}.png')
    plt.close()
    del f, axs, bg_sections