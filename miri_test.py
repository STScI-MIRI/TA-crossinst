import numpy as np
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import pdb
import matplotlib.pyplot as plt
import os
import glob


from jwst_ta import centroid

# TEST SCRIPT FOR THE TA ALGORITHM
plt.close('all')

f = 'det_image_1_MIRIMAGE_F560Wexp1_TAImage_sim1_def3.fits'

cdir= os.getcwd()
os.chdir('../TA_sims/TA_ims/')
files = glob.glob('*.fits')

for file in files:
    print file
    h = fits.getheader(file)
    xin = h['X000']
    yin = h['Y000']
    x,y = centroid(infile=file, cbox=5, incoord=(xin, yin), roi=64., bgcorr=-1, out=None)

os.chdir(cdir)
  