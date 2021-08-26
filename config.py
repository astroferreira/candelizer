from astropy import units as u
from astropy.io import fits
from constants import COSMO, _h

SOURCE_FOLDER = '/home/ppxlf2/MERGERS/SKIRT/TNG100-1'
SOURCE_FOLDER_FB = '/home/ppxlf2/MERGERS/SKIRT/TNG100-1'
OUTPUT_FOLDER = '/home/ppxlf2/MERGERS/SKIRT/TNG100-1/PAPER_v2'
FILTER_FOLDER = 'resources/filters'
CANDELS_FOLDER = '/home/ppxlf2/BGS'

FILTERS = [f'{FILTER_FOLDER}/HST_ACS_WFC.F814W.dat', 
           f'{FILTER_FOLDER}/HST_WFC3_IR.F125W.dat', 
           f'{FILTER_FOLDER}/HST_WFC3_IR.F160W.dat']

BANDS = ['F814W', 'F125W', 'F160W']
EXPTIME = [7500, 2100, 2100]

TARGET_PS = [0.05, 0.10, 0.10]

ALLOWED_ORIENTATIONS = ['xy2', 'xz', 'yz', 'oct']
TARGET_REDSHIFTS = ['full']


#bg606w = fits.getdata('/home/ppxlf2/CANDELS/bgf606w.fits', memmap=True)
bg814w = fits.getdata('/home/ppxlf2/CANDELS/bgf814w.fits', memmap=True)
bg125w = fits.getdata('/home/ppxlf2/CANDELS/bgf125w.fits', memmap=True)
bg160w = fits.getdata('/home/ppxlf2/CANDELS/bgf160w.fits', memmap=True)

#f606w_header = fits.getheader('/home/ppxlf2/CANDELS/hlsp_candels_hst_acs_cos-tot_f606w_v1.0_drz.fits', memmap=True)
f814w_header = fits.getheader('/home/ppxlf2/CANDELS/hlsp_candels_hst_acs_gs-tot_f814w_v1.0_drz.fits', memmap=True)
f125w_header = fits.getheader('/home/ppxlf2/CANDELS/hlsp_candels_hst_wfc3_cos-tot_f125w_v1.0_drz.fits', memmap=True)
f160w_header = fits.getheader('/home/ppxlf2/CANDELS/hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_drz.fits', memmap=True)



#psf606w = fits.getdata('/home/ppxlf2/psfs/PSFs/f606w.fits')
psf814w = fits.getdata('/home/ppxlf2/psfs/TT814.fits')
psf125w = fits.getdata('/home/ppxlf2/psfs/TT125.fits')
psf160w = fits.getdata('/home/ppxlf2/psfs/TT160.fits')


HEADERS = [f814w_header, f125w_header, f160w_header]
BGS = [bg814w, bg125w, bg160w]
PSFS = [psf814w, psf125w, psf160w]


from areia.areia import Config as AreiaConfig

#FLAGS
restframe=False
shotnoise=True
bg=True
convolve=True
pre_loaded_coords = False
zscale = True
fix_scaling = True
bg_mode = 'FIELD' # CLEAN, FIELD
bg_field = 'COSMOS'


areia_config = AreiaConfig()
areia_config.h = _h
areia_config.make_cutout = False 
areia_config.dimming = False
areia_config.size_correction = False
areia_config.evo = False
areia_config.rebinning = zscale
areia_config.shot_noise = False
areia_config.convolve_with_psf = convolve
areia_config.add_background = bg
areia_config.output_size = 256
