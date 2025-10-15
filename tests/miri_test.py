import numpy as np
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import pdb
import matplotlib.pyplot as plt
import os
import glob
from importlib import reload


from jwst_ta import jwst_ta

# TEST SCRIPT FOR THE TA ALGORITHM
plt.close('all')

f = 'jw01053001001_02101_00001_mirimage_uncal.fits'
#f = 'TA_img_for_MIRI_5489_49_S_20170308-042414_SCE3.fits'

cdir= os.getcwd()
#os.chdir('../TA_sims/TA_ims/')
#files = glob.glob('*.fits')
files = [f]

tarefx = 412.5
tarefy = 288.5

for file in files:
    print(file)
    # h = fits.getheader(file)
    # these coordinates were obtained with DS9
    x,y = jwst_ta.centroid(infile=file, input_type='ramp', cbox=5, incoord=(tarefx, tarefy), roi=48, bgcorr=0.3, thresh=0.5)
    
   


  