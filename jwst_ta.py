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
def checkbox(data, box):
    
    ''' 
    This function performs the coarse centroiding on the data array provided.
    
    - data:         a 2D array
    - box:          the size over which each element of the checkbox image will be computed. has to be smaller than the size of data
    
    '''
    
    print 'performing coarse centroiding on an array of size {0}'.format(np.shape(data))
    
    
    hbox = np.int(np.floor(box/2.))
    #print hbox
    psum = 0.
    ix = 0.
    iy = 0.
    
    x = hbox + np.arange(np.size(data, axis=0)-box)
    y = hbox + np.arange(np.size(data, axis=1)-box)
    #print np.min(x), np.max(x), np.min(y), np.max(y)
    
    
    if (box == 1.):
        chbox = np.copy(data)
        maxpx = np.argmax(chbox)
        ix, iy = np.unravel_index(maxpx, np.shape(data))
        
    
    else:
        chbox = np.zeros_like(data)
        for i in x:
            for j in y:
                # remember to add 1 so python sums all 5 pixels
                psumu = np.sum(data[int(i)-hbox:int(i)+hbox+1, int(j)-hbox:int(j)+hbox+1])
                chbox[i,j] = psumu
                if (psumu > psum):
                    ix = i
                    iy = j
                    psum = psumu
            
    
    
    plt.figure()
    plt.imshow(chbox, origin='lower', aspect='equal', cmap='gist_heat')
    plt.plot(iy, ix, marker='x', mew=2., ms = 10.)
    plt.title('Checkbox image with coarse centroid')
    plt.show()
    
    print 'Coarse centroid found at ({0}, {1})'.format(ix, iy)
    
    
   
   
   
    return ix, iy
    
    

#=====================================================

def centroid(infile=None, ext=0, cbox=5, cwin=6, incoord=(0., 0.), roi=None, bgcorr=-1, out=None):
    
    '''
    Implementation of the JWST GENTALOCATE algorithm. Parameters key:
    
    - infile:       FITS filename
    - ext:          extension number of the FITS file containing the science data (default = 0)
    - checkbox:     the FULL size of the checkbox, in pixels, for coarse centroiding (default = 5)
    - cwin:         the FULL size of the centroid window, in pixels, for fine centroiding (default = 6)
    - incoord:      (x,y) input coordinates of the source position
    - roi:          size of a region of interest to be used for the centroiding (optional). If not set, full image will be used for coarse centroiding
                        * setting an ROI also requires input coordinates
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
    
    # NOTE: in python the x-coord is axis 1, y-coord is axis 0
    xin = incoord[1]
    yin = incoord[0]
    print 'Input coordinates = ({0}, {1})'.format(xin, yin)
    
    # consider first the simple case where the image is 2-D and we perform just 1 iteration, and no background subtraction
    if (roi != None):
        # first check that the ROI is a sensible number. if it's bigger than the size of the array, use the full array instead
        if (roi >= n[0])|(roi >= n[1]):
            print 'ROI size is bigger than the image; using full image instead'
            xc, yc = checkbox(im, cbox)
            print xc, yc
        else:
            print np.round(xin-(roi/2.))
            xc, yc = checkbox(im[np.round(xin-(roi/2.)):np.round(xin+(roi/2.)), np.round(yin-(roi/2.)):np.round(yin+(roi/2.))], cbox)
            xc += np.round(xin-(roi/2.))
            yc += np.round(yin-(roi/2.))
            print xc, yc
    else:
        xc, yc = checkbox(im, cbox)
        print xc, yc
            
            
    
    
    return 0
    
#=====================================================


    
    
    
    