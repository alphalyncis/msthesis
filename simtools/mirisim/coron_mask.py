# Mask central star to 1%
import numpy as np
import glob                 # glob is used to find the output directory
import os                   # for listing directory contents
import sys
from astropy.io import fits # for reading FITS file contents


def mask(fname):
    '''
    e.g.
    fname = '10jup50au/cube5npix23lam1020inc0.fits'
    lams = ['1020','1500','2100']
    '''
    lam = fname[fname.find('lam')+3:fname.find('inc')]
    datadir = '/home/ipa/ipaszulagyi/users/xuechen/thesis/RADMC3D/'

    with fits.open(datadir+fname) as hdulist:
        img = hdulist[0].data
        hdr = hdulist[0].header

    a = img[0].shape[0]
    if (a % 2) == 0:
        img[:,int(a/2-1):int(a/2+1), int(a/2-1):int(a/2+1)] *= 0.01
    else:
        img[:,int((a-1)/2),int((a-1)/2)] *= 0.01

    hdu = fits.PrimaryHDU(img)
    hdulist = fits.HDUList([hdu])
    hdulist[0].header = hdr
    print('writing to '+fname[:-5]+'_coron.fits ...')
    hdu.writeto(datadir+fname[:-5]+'_coron.fits', overwrite=True)

if __name__ == "__main__":
    fname = sys.argv[1]
    mask(fname)
