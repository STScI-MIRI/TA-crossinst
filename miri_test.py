import numpy as np
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import pdb
import matplotlib.pyplot as plt
import os
import glob


from jwst_ta import centroid, make_ta_image

# TEST SCRIPT FOR THE TA ALGORITHM
plt.close('all')

f = 'MIRI_5489_55_S_20170308-045954_SCE3.fits'
#f = 'TA_img_for_MIRI_5489_49_S_20170308-042414_SCE3.fits'

cdir= os.getcwd()
#os.chdir('../TA_sims/TA_ims/')
#files = glob.glob('*.fits')
files = [f]

for file in files:
    print(file)
    # h = fits.getheader(file)
    # these coordinates were obtained with DS9
    xin = 389.313
    yin = 792.000
    x,y = centroid(infile=file, input_type='ramp', cbox=5, incoord=(xin, yin), roi=64, bgcorr=-1, thresh=0.5)
    
   


  