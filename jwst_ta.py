#! /usr/bin/env python

import os
from copy import copy
import numpy as np
import astropy.io.fits as fits
import astropy.io.ascii as ascii
import matplotlib.pyplot as plt
import pdb

def rclip(x, xmin, xmax):
    dum = x
    if dum<xmin:
        dum = xmin
    elif dum>xmax:
        dum = xmax
    
    return dum


#=====================================================
def checkbox(data, box, bgcorr):
    
    ''' 
    This function performs the coarse centroiding on the data array provided.
    
    - data:         a 2D array
    - box:          the size over which each element of the checkbox image will be computed. has to be smaller than the size of data
    - bgcorr:       inherit the background correction parameter so can apply here.
    
    '''
    
    print('performing coarse centroiding on an array of size {0}'.format(np.shape(data)))
    
    
    hbox = np.int(np.floor(box/2.))
    #print hbox
    psum = 0.
    ix = 0.
    iy = 0.
    
    x = hbox + np.arange(np.size(data, axis=1)-box)
    y = hbox + np.arange(np.size(data, axis=0)-box)
    if (box == 1.):
        chbox = np.copy(data)
        maxpx = np.argmax(chbox)
        ix, iy = np.unravel_index(maxpx, np.shape(data))    
    else:
        chbox = np.zeros_like(data)
        for i in x:
            for j in y:
                # remember to add 1 so python sums all 5 pixels
                psumu = np.sum(data[int(j)-hbox:int(j)+hbox+1, int(i)-hbox:int(i)+hbox+1])
                chbox[j,i] = psumu
                if (psumu > psum):
                    ix = i
                    iy = j
                    psum = psumu
    
    #plt.figure()
    #plt.imshow(chbox, origin='lower', aspect='equal', cmap='gist_heat', interpolation='None')
    #plt.plot(iy, ix, marker='x', mew=2., ms = 10.)
    #plt.title('Checkbox image with coarse centroid')
    #plt.show()
    
    
    return ix, iy
    
    

#=====================================================
def fine_centroid(data, cwin, xc, yc):
    
    sumx = 0.
    sumy = 0.
    sump = 0.
    
    # define rwin, half the cwin setting
    rwin = np.floor(cwin/2.)
    
    # remember to add 1 so we get all cwin pixels
    x = (xc - np.floor(cwin/2.)) + np.arange(cwin)
    y = (yc - np.floor(cwin/2.)) + np.arange(cwin) 
    #pdb.set_trace()
    
    # write weights out to an array - for testing
    #wtarr = np.zeros_like(data)
    
    for i in x:
        for j in y:
            wx = rclip(rwin-abs(i-xc)+0.5, 0.0, 1.0)
            wy = rclip(rwin-abs(j-yc)+0.5, 0.0, 1.0)
            ww = wx*wy
            # for testing - delete once tested
            #wtarr[i,j] = ww
            
            sumx += ww * data[int(j),int(i)] * i
            sumy += ww * data[int(j),int(i)] * j
            sump += ww * data[int(j),int(i)]
    
    
    # plot of pixel weights, for testing
    #plt.figure()
    #plt.imshow(wtarr[xc-5:xc+5, yc-5:yc+5], origin='lower', aspect='equal', interpolation='None')
    #plt.colorbar()
    #plt.show()
    # end test plot
    
    xc_old = xc
    yc_old = yc
    
    xc = sumx/sump
    yc = sumy/sump
    
    

    return xc, yc
#=====================================================
def bgrsub(data, val):
    
    '''
    Does the background correction step, if needed. Method is determined from the val parameter:
    
    - val <= 0:     no background subtraction (shouldn't come to this function, but let's check anyway)
    - 0 < val < 1:  fractional background. sorts all pixels and subtracts the value of this percentile
    - val > 1:      this value is subtracted uniformly from all pixels
    
    '''
    
    if (val <= 0.):
        print('No background subtracted')
        outdata = data
    elif (val > 0.) & (val < 1.):
        
        # test plotting code - cumulative distribution function of the pixels in data
        
        bgrval = np.percentile(data, val*100.)

        # Histogram of signal values
        plt.figure()
        bins = np.arange(-100,1000,50)
        bins2 = np.linspace(1000,np.max(data),5)
        bins = np.concatenate([bins,bins2])
        plt.hist(np.ravel(data), bins=bins, lw = 1.5, histtype='step', cumulative=True, normed=True)
        plt.axvline(bgrval, ymin=0., ymax=1., color='r',
                    lw = 1.5, label='Background Level')
        plt.xlim(np.min(data), np.min(data) + np.abs(bgrval)*4)
        plt.xlabel("Signal (ADU)")
        plt.ylabel("Cumulative Fraction of Pixels")
        plt.legend(loc='lower right')
        plt.show()

        # Subtract background level
        outdata = data - bgrval
        print('subtracting {0} from image'.format(bgrval))

    else:
        outdata = data - val
    
    return outdata
    
    
#=====================================================

def centroid(infile=None, input_type='image', ext=0, cbox=5, cwin=5, incoord=(0., 0.), roi=None, bgcorr=-1, flat=None, flatext=0, out=None):
    
    '''
    Implementation of the JWST GENTALOCATE algorithm. Parameters key:
    
    - infile:       FITS filename
    - input_type:   description of input data: 'image' or 'ramp'. If 'ramp'
                    then make_ta_image functin is run. If 'image' (default)
                    centroiding is performed directly on the data in the
                    input file
    - ext:          extension number of the FITS file containing the science data (default = 0)
    - checkbox:     the FULL size of the checkbox, in pixels, for coarse centroiding (default = 5)
    - cwin:         the FULL size of the centroid window, in pixels, for fine centroiding (default = 5)
    - incoord:      (x,y) input coordinates of the source position
    - roi:          size of a region of interest to be used for the centroiding (optional). If not set, full image will be used for coarse centroiding
                        * setting an ROI also requires input coordinates
    - bgcorr:       background correction parameter. set to:
                        * negative value for NO background subtraction (default)
                        * 0 < bgcorr < 1 for fractional background subtraction
                        * bgcorr > 1 for constant background subtraction number (this number will be subtracted from the entire image)
    - out:          enter a filename for output of the fit results to a file (default = None)
    '''

    # Read in data. Create the TA image if requested
    if input_type.lower() == 'image':
        hdu = fits.open(infile)
        im = hdu[ext].data
        h = hdu[ext].header
    elif input_type.lower() == 'ramp':
        im = make_ta_image(infile, ext=ext, useframes=3)
        # Save TA image for code testing
        h0 = fits.PrimaryHDU(im)
        hl = fits.HDUList([h0])
        indir, inf = os.path.split(infile)
        tafile = os.path.join(indir, 'TA_img_for_'+inf)
        hl.writeto(tafile, overwrite=True)

    # Do background correction first
    if bgcorr > 0.:
        im = bgrsub(im, bgcorr)

    # Apply flat field
    if flat is not None:
        # Read in flat
        with fits.open(flat) as ff:
            flatfield = ff[flatext].data
            
        # Flat must be the same size as the data
        ffshape = flatfield.shape
        dshape = im.shape
        if dshape != ffshape:
            raise RunTimeError(("WARNING: flat field shape ({}) does "
                                "not match data shape ({})!"
                                .format(ffshape,dshape)))
        # Apply flat
        im = apply_flat_field(im, flatfield)
        
    ndim = np.ndim(im)
    #pdb.set_trace()
    
    n = [np.size(im, axis=i) for i in range(ndim)]
    
    # NOTE: in python the x-coord is axis 1, y-coord is axis 0
    xin = incoord[1]
    yin = incoord[0]
    print('Input coordinates = ({0}, {1})'.format(xin, yin))
    
    # Extract the ROI
    if (roi is not None):
        # first check that the ROI is a sensible number.
        # if it's bigger than the size of the array, use the
        # full array instead
        if (roi >= n[0])|(roi >= n[1]):
            print('ROI size is bigger than the image; using full image instead')
            roi_im = im
            xoffset = 0
            yoffset = 0
            #xc, yc = checkbox(im, cbox, bgcorr)
        else:
            roi_im = im[np.round(yin-(roi/2.)):np.round(yin+(roi/2.)+1),
                        np.round(xin-(roi/2.)):np.round(xin+(roi/2.))]
            #xc, yc = checkbox(im[np.round(yin-(roi/2.)):np.round(yin+(roi/2.)+1), np.round(xin-(roi/2.)):np.round(xin+(roi/2.))], cbox, bgcorr)
            #xc += np.round(xin-(roi/2.))
            #yc += np.round(yin-(roi/2.))
            xoffset = np.round(xin-(roi/2.))
            yoffset = np.round(yin-(roi/2.))
    else:
        #xc, yc = checkbox(im, cbox, bgcorr)
        roi_im = im
        xoffset = 0
        yoffset = 0

    # Perform coarse centroiding. Pay attention to coordinate
    # offsets
    xc, yc = checkbox(roi_im, cbox, bgcorr)
    xc += xoffset
    yc += yoffset
    print('Coarse centroid found at ({0}, {1})'.format(xc, yc))

    # Iterate fine centroiding
    iter_thresh = 0.1
    nconv = 0
    while nconv == 0:
        xf, yf = fine_centroid(im, cwin, xc, yc)
        err = np.sqrt((xf-xin)**2 + (yf-yin)**2)
        print(("Fine centroid found at (x, y) = ({0:.4f}, {1:.4f}). "
               "Rms error = {2:.4f}".format(xf, yf, err)))
        if (abs(xf-xc) <= iter_thresh) & (abs(yf-yc) <= iter_thresh):
            nconv = 1
        xc = xf
        yc = yf
    

    return xf, yf
    
#=====================================================

def make_ta_image(infile, ext=0, useframes=3):
    """
    Create the image on which to perform the centroiding
    given a fits file containing an exposure with 
    multiple groups. In the case of multiple integrations
    in the exposure, use only the first, and produce a 
    single TA image.

    Parameters:
    -----------
    infile -- Name of FITS file containing exposure. Must be
              an exposure containing at least one integration
              and the integration must have an odd number of groups
    ext -- Extension number containing the data.
           (Change to use 'SCI' rather than number?)
    useframes -- Number of frames to use for the calculation.
                 Most of the time 3 frames are used, these 
                 being the 1st, (N+1)/2, and Nth groups
                 of the integration. 

                 There are, apparently, ways to do the calculation
                 with larger numbers of groups, but I haven't 
                 seen the math for those yet.

    Returns:
    --------
    2D numpy ndarray of the TA image
    """

    # Read in data. Convert to floats
    with fits.open(infile) as h:
        data = h[ext].data
    data = data * 1.
        
    shape = data.shape
    if len(shape) <= 2:
        raise RunTimeError(("Warning: Input target acq exposure must "
                            "have multiple groups!"))
    elif len(shape) == 4:
        # If there are multiple integrations, use only the first
        data = data[0, :, :, :]
        
    ngroups = shape[-3]
    if ngroups % 2 == 0:
        raise RunTimeError(("Warning: Input target acq exposure "
                            "must have an odd number of groups!"))

    # Group numbers to use. Adjust the values to account for
    # python being 0-indexed
    if useframes == 3:
        frames = [0, np.int((ngroups-1)/2), ngroups-1]
        diff21 = data[frames[1], :, :] - data[frames[0], :, :]
        diff32 = data[frames[2], :, :] - data[frames[1], :, :]
        ta_img = np.minimum(diff21, diff32)
    elif useframes == 5:
        something_else

    return ta_img


#=====================================================
def apply_flat_field(image, flat):
    """
    Apply flat field to TA image. Assume the flat 
    has the format matching those to be used on 
    board by GENTALOCATE. Pixel values are multiplied
    by 1000 relative to traditional flat field files.
    (i.e. flat is normalized to a value of 1000).
    Bad pixels have a value of 65535. Bad pixels
    receive a value that is interpolated from 
    nearest neighbors.

    Parameters:
    -----------
    image -- TA image. 2D ndarray
    flat -- Flat field image. 2D ndarray.
    
    Returns:
    --------
    Flat fielded image -- 2D ndarray
    """
    # Make sure flat field values are floats
    flat = flat * 1.
    
    # Find bad pixels and set to NaN
    bad = flat == 65535
    print("Found {} bad pixels in the flat.".format(np.sum(bad)))
    flat[bad] = np.nan
    
    # Apply flat
    image /= (flat/1000.)

    # Use surrounding pixels to set bad pixel values
    # NOT SURE IF THIS IS IMPLEMENTED IN THE REAL
    # GENTALOCATE OR NOT...
    if np.any(bad):
        image = fixbadpix(image)

    return image



#======================================================
def fixbadpix(data, maxstampwidth=3, method='median'):
    """
    Replace the values of bad pixels in the TA image
    with interpolated values from neighboring good
    pixels. Bad pixels are identified as those with a
    value of 65535 in the flat field file

    Parameters:
    -----------
    data -- 2D array containing the TA image
    bpix -- Tuple of lists of bad pixel coordinates. 
            (output from np.where on the 2D flat field image)
    maxstampwidth -- Maximum width of area centered on a bad pixel
                     to use when calculating median to fix the bad 
                     pixel. Must be an odd integer.
    method -- The bad pixel is fixed by taking the median or the 
              mean of the surrounding pixels. This is a string
              that can be either 'median' or 'mean'.

    Returns:
    --------
    2D image with fixed bad pixels
    """
    ny, nx = data.shape

    # Set up the requested calculation method
    if method == "median":
        mmethod = np.nanmedian
    elif method == "mean":
        mmethod = np.nanmean
    else:
        print("Invalid method. Must be either 'median' or 'mean'.")
        sys.exit()
    
    if (maxstampwidth % 2) == 0:
        print("maxstampwidth must be odd. Adding one to input value.")
        maxstampwidth += 1 

    half = np.int((maxstampwidth - 1)/2)

    bpix = np.isnan(data)
    bad = np.where(bpix)
    # Loop over the bad pixels and correct
    for bady, badx in zip(bad[0], bad[1]):

        print('Bad pixel:',bady,badx)
        
        substamp = np.zeros((maxstampwidth, maxstampwidth))
        substamp[:,:] = np.nan
        minx = badx - half
        maxx = badx + half + 1
        miny = bady - half
        maxy = bady + half + 1

        # Offset between data coordinates and stamp
        # coordinates
        dx = copy(minx)
        dy = copy(miny)

        # Check for stamps that fall off the edges
        # of the data array
        sminx = 0
        sminy = 0
        smaxx = maxstampwidth
        smaxy = maxstampwidth
        if minx < 0:
            sminx = 0 - minx
            minx = 0
        if miny < 0:
            sminy = 0 - miny
            miny = 0
        if maxx > nx:
            smaxx = maxstampwidth - (maxx - nx)
            maxx = nx
        if maxy > ny:
            smaxy = maxstampwidth - (maxy - ny)
            maxy = ny
            
        substamp[sminy:smaxy, sminx:smaxx] = data[miny:maxy, minx:maxx]

        # First try the mean of only the 4 adjacent pixels
        neighborsx = [half, half+1, half, half-1]
        neighborsy = [half+1, half, half-1, half]
        if np.sum(np.isnan(substamp[neighborsx, neighborsy])) < 4:
            data[bady, badx] = mmethod(substamp[neighborsx, neighborsy])
            print(("Good pixels within nearest 4 neighbors. Mean: {}"
                   .format(mmethod(substamp[neighborsx, neighborsy]))))
            continue

        # If the adjacent pixels are all NaN, expand to include corners
        else:
            neighborsx.append([half-1, half+1, half+1, half-1])
            neighborsy.append([half+1, half+1, half-1, half-1])
            if np.sum(np.isnan(substamp[neighborsx, neighborsy])) < 8:
                data[bady, badx] = mmethod(substamp[neighborsx, neighborsy])
                print(("Good pixels within 8 nearest neighbors. Mean: {}"
                       .format(mmethod(substamp[neighborsx, neighborsy]))))
                continue

        # If all pixels are still NaN, iteratviely expand to include
        # more rings of pixels until the entire stamp image is used
        # (This step not included in Goudfrooij's original bad pixel
        # correction script).
        delta = 2
        while delta <= half:
            newy = np.arange(half-(delta-1), half+delta)
            newx = np.repeat(half - delta, len(newy))
            neighborsx.extend(newx)
            neighborsy.extend(newy)
            newx = np.repeat(half + delta, len(newy))
            neighborsx.extend(newx)
            neighborsy.extend(newy)
            newx = np.arange(half-delta, half+delta+1)
            newy = np.repeat(half - delta, len(newx))
            neighborsx.extend(newx)
            neighborsy.extend(newy)
            newy = np.repeat(half + delta, len(newx))
            neighborsx.extend(newx)
            neighborsy.extend(newy)
            if np.sum(np.isnan(substamp[neighborsx, neighborsy])) < (len(neighbosrsx)):
                data[bady, badx] = mmethod(substamp[neighborsx, neighborsy])
                print("Expanding to {} rows".format(delta))
                continue
            else:
                delta += 1
        print(("Warning: all pixels within {} rows/cols of the bad pixel at ({},{}) "
               "are also bad. Cannot correct this bad pixel with this stamp image"
               "size.".format(delta, badx, bady)))

    return data

        
#====================================================== 
    
    
    
    
