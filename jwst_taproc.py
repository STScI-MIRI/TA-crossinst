import numpy as np
import astropy.io.fits as fits
import astropy.io.ascii as ascii



def rclip(x, xmin, xmax):
	dum = x
    if (dum < xmin):
        dum = xmin
    elif (dum > xmax):
        dum = xmax
    return dum


#=====================================================

def jwst_centroid(infile=None, ext=0, checkbox=3, cwin=3, bgcorr=0.0, out=None):
    
    hdu = fits.open(infile)
    im = hdu[ext].data
    h = hdu[ext].header
    
    ndim = np.ndim(im)
    nx, ny, nz = [np.size(im, axis=i) for i in range(3))]
    
    xin = h['X000']
    yin = h['Y000']
    
    
    
    
    
    