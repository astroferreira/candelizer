import numpy as np

from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from photutils.segmentation import detect_sources
from photutils.segmentation import detect_threshold 

from scipy.ndimage import binary_dilation, zoom                  

sigma = 2.0 * gaussian_fwhm_to_sigma  # FWHM = 3.

def segmentate(bg):                   
    
    threshold = detect_threshold(bg, nsigma=2.)
    
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm_data = detect_sources(bg, threshold, npixels=10, filter_kernel=kernel)


    flux = 0
    num_pix = 1
    if segm_data is not None:
        circ_kernel = generate_circular_kernel(7)
        segm_data.data = binary_dilation(segm_data.data, circ_kernel)

        binarymap = segm_data.data.copy()
        segm_data = segm_data.data
        binarymap[segm_data > 0] = 1
        num_pix = np.where(binarymap == 1)[0].shape[0]
        flux = bg[binarymap==1].sum()
    else:
        segm_data = np.zeros_like(bg)
        
    flux_per_pixel =  flux / num_pix
    #print(flux_per_pixel, flux, num_pix)
    return segm_data, flux_per_pixel, flux, num_pix


def generate_circular_kernel(d):
    '''
        Generates a circular kernel to be used with
        dilation transforms used in GalClean. This
        is done with an 2D matrix of 1s, flipping to
        0 when x^2 + y^2 > r^2.

        Parameters
        ----------
        d : int
            The diameter of the kernel. This should be an odd
            number, but the function handles it in either case.

        Returns
        -------
        circular_kernel : 2D numpy.array bool
            A 2D boolean array with circular shape.
    '''

    d = int(d)

    if (d % 2) == 0:
        d = d + 1

    mask = np.ones((d, d))
    r = np.round(d/2)
    x0 = r
    y0 = r
    for i in range(0, d, 1):
        for j in range(0, d, 1):
            xx = (i-x0)**2
            yy = (j-y0)**2
            rr = r**2
            if xx + yy > rr:
                mask[i][j] = 0

    return mask