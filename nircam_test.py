#! /usr/bin/env python

"""
Run the centroiding algorithm on some simulated NIRCam data
"""

import os
import glob
import numpy as np
import astropy.io.fits as fits
from jwst_ta import centroid

#files = glob.glob('*.fits')
#files = ['nrca5_TA_timeseries_NRCFLATA5GRTS_mag09.4_center_15.5_15.5_noiserealiz_1_uncal.fits']
files = ['nrca5_TA_timeseries_NRCFLATA5GRTS_mag09.4_center_15.25_15.25_noiserealiz_1_uncal.fits', 'nrca5_TA_timeseries_NRCFLATA5GRTS_mag09.4_center_15.75_15.75_noiserealiz_1_uncal.fits']
#flatfile = 'NRCFLATA5GRTS_flat.fits'
flatfile = 'NRCFLATA5GRTS_flat_withbadpix.fits'

for file in files:
    print(file)
    with fits.open(file) as h:
        yd, xd = h[1].data.shape[-2:]
    xin = 15.
    yin = 15.
    # NOTE: The NIRCam checkbox is 3x3 and the fine centroiding box
    # is 9x9. Making the checkbox larger (e.g. 5x5) leads to very
    # poor centroiding results.
    x,y = centroid(infile=file, input_type='ramp', ext=1,
                   cbox=3, cwin=9, incoord=(xin, yin), bgcorr=-1, thresh=0.5)
