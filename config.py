from astropy import units as u
from constants import COSMO, _h

SOURCE_FOLDER = '/home/ppxlf2/MERGERS/SKIRT/TNG50-1/HALO/'
SOURCE_FOLDER_FB = '/home/ppxlf2/MERGERS/SKIRT/TNG50-1/SUBHALO/'
OUTPUT_FOLDER = '/home/ppxlf2/MERGERS/SKIRT/TNG50-1/CLEAN/'
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
pre_loaded_coords = True
zscale = True
bg_mode = 'FIELD' # CLEAN, FIELD
bg_field = 'COSMOS'

areia_config = AreiaConfig()
areia_config.h = _h
areia_config.make_cutout = False
areia_config.dimming = False
areia_config.size_correction = False
areia_config.evo = False
areia_config.rebinning = zscale
areia_config.shot_noise = shotnoise
areia_config.convolve_with_psf = convolve
areia_config.add_background = bg