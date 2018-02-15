import numpy as np
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import pdb
import matplotlib.pyplot as plt


from jwst_ta import centroid

# TEST SCRIPT FOR THE TA ALGORITHM
plt.close('all')

f = 'det_image_1_MIRIMAGE_F560Wexp1_TAImage_sim1_def3.fits'

h = fits.getheader(f)
xin = h['X000']
yin = h['Y000']

x = centroid(infile=f, cbox=1, incoord=(xin, yin), roi=None, bgcorr=-1, out=None)

    