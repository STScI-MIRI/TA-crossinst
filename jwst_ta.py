import numpy as np
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import pdb



def rclip(x, xmin, xmax):
    dum = x
    if dum<xmin:
        dum = xmin
    elif dum>xmax:
        dum = xmax
    
    return dum


#=====================================================

def centroid(infile=None, ext=0, checkbox=5, cwin=6, incoord=(0., 0.), roi=None, bgcorr=-1, out=None):
    
    '''
    Implementation of the JWST GENTALOCATE algorithm. Parameters key:
    
    - infile:       FITS filename
    - ext:          extension number of the FITS file containing the science data (default = 0)
    - checkbox:     the FULL size of the checkbox, in pixels, for coarse centroiding (default = 5)
    - cwin:         the FULL size of the centroid window, in pixels, for fine centroiding (default = 6)
    - incoord:      (x,y) input coordinates of the source position (optional)
    - roi:          size of a region of interest to be used for the centroiding (optional). If not set, full image will be used for coarse centroiding
    - bgcorr:       background correction parameter. set to:
                        * negative value for NO background subtraction (default)
                        * 0 < bgcorr < 1 for fractional background subtraction
                        * bgcorr > 1 for constant background subtraction number (this number will be subtracted from the entire image)
    - out:          enter a filename for output of the fit results to a file (default = None)
    '''
    
    hdu = fits.open(infile)
    im = hdu[ext].data
    h = hdu[ext].header
    
    ndim = np.ndim(im)
    #pdb.set_trace()
    
    n = [np.size(im, axis=i) for i in range(ndim)]
    
    xin = h['X000']
    yin = h['Y000']
    print 'Input coordinates = ({0}, {1})'.format(xin, yin)
    
    # consider first the simple case where the image is 2-D and we perform just 1 iteration, and no background subtraction
    
    
    
    return 0
    
#=====================================================


    
    
    
    