import pickle

from astropy.io import fits
from astropy import units as u
from astropy.cosmology import LambdaCDM
import numpy as np 

FIELDS = {
            'COSMOS' : {
                        'F814W' : {
                                    'exptime' : 6900,
                                    'zpt' : 25.947,
                                    'lim_mag' : 28.4,
                                    'skystd' : 0.000721,
                                    'fits_filename' : 'hlsp_candels_hst_acs_cos-tot_f814w_v1.0_drz',
                                },
                        'F125W' : {
                                    'exptime' : 1900,
                                    'zpt'     : 26.2303,
                                    'lim_mag' : 27,
                                    'skystd' : 0.003130,
                                    'fits_filename' : 'hlsp_candels_hst_wfc3_cos-tot_f125w_v1.0_drz'
                                },
                        'F160W' : {
                                    'exptime' : 3200,
                                    'zpt'     : 25.9463,
                                    'lim_mag' : 26.9,
                                    'skystd'  : 0.002652,
                                    'fits_filename' : 'hlsp_candels_hst_wfc3_cos-tot_f160w_v1.0_drz'
                                }
                    },
            'EGS' : {
                        'F814W' : {
                                    'exptime' : 6900,
                                    'zpt' : 25.94333,
                                    'lim_mag' : 28.4,
                                    'skystd' : 0.000691,
                                    'fits_filename' : 'hlsp_candels_hst_acs_egs-tot-30mas_f814w_v1.0_drz'
                                },
                        'F125W' : {
                                    'exptime' : 1900,
                                    'zpt'     : 26.25,
                                    'lim_mag' : 27,
                                    'skystd' : 0.003176,
                                    'fits_filename' : 'hlsp_candels_hst_wfc3_egs-tot-60mas_f125w_v1.0_drz'
                                },
                        'F160W' : {
                                    'exptime' : 3200,
                                    'zpt'     : 25.96,
                                    'lim_mag' : 26.9,
                                    'skystd'  : 0.002624,
                                    'fits_filename' : 'hlsp_candels_hst_wfc3_egs-tot-60mas_f160w_v1.0_drz'
                                }
                    },
            'UDS' : {
                        'F814W' : {
                                    'exptime' : 5700,
                                    'zpt' : 25.94333,
                                    'lim_mag' : 28.4,
                                    'skystd' : 0.000629,
                                    'fits_filename' : 'hlsp_candels_hst_acs_uds-tot_f814w_v1.0_drz'
                                },
                        'F125W' : {
                                    'exptime' : 1900,
                                    'zpt'     : 26.25,
                                    'lim_mag' : 27,
                                    'skystd' : 0.003180,
                                    'fits_filename' : 'hlsp_candels_hst_wfc3_uds-tot_f125w_v1.0_drz'
                                },
                        'F160W' : {
                                    'exptime' : 3300,
                                    'zpt'     : 25.96,
                                    'lim_mag' : 27.1,
                                    'skystd'  : 0.002155,
                                    'fits_filename' : 'hlsp_candels_hst_wfc3_uds-tot_f160w_v1.0_drz'
                                }
                    },
                'GOODS_N' : {
                        'F814W' : {
                                    'exptime' : 1200,
                                    'zpt' : 25.947,
                                    'lim_mag' : 28.8,
                                    'skystd' : 0.000487,
                                    'fits_filename' : 'hlsp_candels_hst_acs_gn-tot-30mas_f814w_v1.0_drz'
                                },
                        'F125W' : {
                                    'exptime' : 1900,
                                    'zpt'     : 26.2303,
                                    'lim_mag' : 27.2,
                                    'skystd' : 0.002788,
                                    'fits_filename' : 'hlsp_candels_hst_wfc3_gn-tot-60mas_f125w_v1.0_drz'
                                },
                        'F160W' : {
                                    'exptime' : 2800,
                                    'zpt'     : 25.9463,
                                    'lim_mag' : 27,
                                    'skystd'  : 0.002371,
                                    'fits_filename' : 'hlsp_candels_hst_wfc3_gn-tot-60mas_f160w_v1.0_drz'
                                }
                    },

                'GOODS_S' : {
                        'F814W' : {
                                    'exptime' : 7500,
                                    'zpt' : 25.947,
                                    'lim_mag' : 28.6,
                                    'skystd' : 0.000549,
                                    'fits_filename' : 'hlsp_candels_hst_acs_gs-tot_f814w_v1.0_drz'
                                },
                        'F125W' : {
                                    'exptime' : 2700,
                                    'zpt'     : 26.2687,
                                    'lim_mag' : 27.2,
                                    'skystd' : 0.002564,
                                    'fits_filename' : 'hlsp_candels_hst_wfc3_gs-tot_f125w_v1.0_drz'
                                },
                        'F160W' : {
                                    'exptime' : 2100,
                                    'zpt'     : 25.9463,
                                    'lim_mag' : 26.7,
                                    'skystd'  : 0.003260,
                                    'fits_filename' : 'hlsp_candels_hst_wfc3_gs-tot_f160w_v1.0_drz'
                                }
                    }
        }


for field in list(FIELDS.keys()):
    for filter in FIELDS[field]:
        FIELDS[field][filter]['data'] = fits.getdata(f'/home/ppxlf2/CANDELS/{field}/{FIELDS[field][filter]["fits_filename"]}.fits', memmap=True)

with open('resources/wcs_memory.pickle', 'rb') as handle:
    WCS_MEMORY = pickle.load(handle)


_h = 0.6774
SKIRT_Z0 = 0.002265
SKIRT_OUTSIZE = 190
FOV = 60

COSMO = LambdaCDM(H0=100*_h, Om0=0.3089, Ob0=0.0486, Ode0=0.6911, Tcmb0=2.73)
ps = FOV * u.kpc * COSMO.arcsec_per_kpc_comoving(SKIRT_Z0)/SKIRT_OUTSIZE


""" UNITS """
SIM_UNITS = u.W / u.m**2 / u.arcsecond**2 * ps**2 # units of the datacube produced by SKIRT
OUT_UNITS = u.Jy * u.Hz / u.micron



SNAP_MAP = np.load('resources/snapTNG.npy')
REFERENCE_Z = 0.5