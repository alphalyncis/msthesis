from __future__ import absolute_import
from __future__ import print_function
import os
import sys
import astropy.units as u

import numpy as np
from astropy.io import fits as pf
import matplotlib.pylab as plt

# Constants 
cc  = 2.9979245800000e10      # Light speed             [cm/s]
pc  = 3.08572e18              # Parsec                  [cm]
au  = 1.49598e13              # Astronomical Unit       [cm]

class radmc3dImage(object):
    """
    RADMC-3D image class

    Attributes
    ----------

    image       : ndarray
                  The image as calculated by radmc3d (the values are intensities in erg/s/cm^2/Hz/ster)

    imageJyppix : ndarray
                  The image with pixel units of Jy/pixel

    x           : ndarray
                  x coordinate of the image [cm]

    y           : ndarray
                  y coordinate of the image [cm]

    nx          : int
                  Number of pixels in the horizontal direction

    ny          : int
                  Number of pixels in the vertical direction

    sizepix_x   : float
                  Pixel size in the horizontal direction [cm]

    sizepix_y   : float
                  Pixel size in the vertical direction [cm]

    nfreq       : int
                  Number of frequencies in the image cube

    freq        : ndarray
                  Frequency grid in the image cube

    nwav        : int
                  Number of wavelengths in the image cube (same as nfreq)

    wav         : ndarray
                  Wavelength grid in the image cube

    """

    def __init__(self):
        self.image = None
        self.x = None
        self.y = None
        self.nx = 0
        self.ny = 0
        self.sizepix_x = 0
        self.sizepix_y = 0
        self.nfreq = 0
        self.freq = None
        self.nwav = 0
        self.wav = None
        self.stokes = False
        self.psf = {}
        self.fwhm = []
        self.pa = 0
        self.dpc = 0     
        self.filename = 'image.out'

    # --------------------------------------------------------------------------------------------------
   
    def readImage(self, fname=None):
        """Reads an image calculated by RADMC-3D

        Parameters
        ----------

        fname   : str, optional
                 File name of the radmc3d output image (if omitted 'image.out' is used)

        """

        # Look for the image file
        print('Reading '+ fname)

        self.filename = fname
        with open(fname, 'r') as rfile:

            dum = ''

            # Format number
            iformat = int(rfile.readline())

            # Nr of pixels
            dum = rfile.readline()
            dum = dum.split()
            self.nx = int(dum[0])
            self.ny = int(dum[1])
            # Nr of frequencies
            self.nfreq = int(rfile.readline())
            self.nwav = self.nfreq
            # Pixel sizes
            dum = rfile.readline()
            dum = dum.split()
            self.sizepix_x = float(dum[0])
            self.sizepix_y = float(dum[1])
            # Wavelength of the image
            self.wav = np.zeros(self.nwav, dtype=np.float64)
            for iwav in range(self.nwav):
                self.wav[iwav] = float(rfile.readline())
            self.wav = np.array(self.wav)
            self.freq = cc / self.wav * 1e4

            # If we have a normal total intensity image
            self.stokes = False

            self.image = np.zeros([self.nx, self.ny, self.nwav], dtype=np.float64)
            for iwav in range(self.nwav):
                # Blank line
                dum = rfile.readline()
                for iy in range(self.ny):
                    for ix in range(self.nx):
                        self.image[ix, iy, iwav] = float(rfile.readline())
        '''
        # Conversion from erg/s/cm/cm/Hz/ster to Jy/pixel
        conv1 = self.sizepix_x * self.sizepix_y / (dpc * pc)**2. * 1e23
        conv2 = 1e23 / 4.25e10
        self.imageJyppix = self.image * conv1
        self.imageJyarcsec2 = self.image * conv2
        '''

        self.x = ((np.arange(self.nx, dtype=np.float64) + 0.5) - self.nx / 2) * self.sizepix_x
        self.y = ((np.arange(self.ny, dtype=np.float64) + 0.5) - self.ny / 2) * self.sizepix_y


def readImage(fname=None):
    """Reads an image calculated by RADMC-3D.
       This function is an interface to radmc3dImage.readImage().

    Parameters
    ----------
        fname   : str, optional
                 File name of the radmc3d output image (if omitted 'image.out' is used)
    """

    dum = radmc3dImage()
    dum.readImage(fname=fname)
    return dum


def write_fits_mirisim(fname):
    '''
    fname = 'npix23lam1020inc0,out'
    '''
    dpc = 100
    coord='03h10m05s -10d05m30s'

    im = readImage(fname)
    data = np.ones([im.nfreq+1, im.nx, im.ny], dtype=float)

    # Calculate pixel scale
    dx = im.sizepix_x / au / dpc   # from cm to arcsec
    dy = im.sizepix_y / au / dpc

    # Conversion from erg/s/cm/cm/Hz/sr to Jy/arcsec^2
    conv = 1e23 / 4.25e10

    for inu in range(im.nfreq):
        data[inu, :, :] = im.image[:, :, inu] * conv
    
    # remove last frame of ones
    data = data[:-1,:,:]

    naxis = len(data.shape)
    hdu = pf.PrimaryHDU(data.swapaxes(naxis - 1, naxis - 2))

    ### Create fits header ###
    dum = coord
    ra = []
    delim = ['h', 'm', 's']
    for i in delim:
        ind = dum.find(i)
        ra.append(float(dum[:ind]))
        dum = dum[ind + 1:]

    dec = []
    delim = ['d', 'm', 's']
    for i in delim:
        ind = dum.find(i)
        dec.append(float(dum[:ind]))
        dum = dum[ind + 1:]

    target_ra = (ra[0] + ra[1] / 60. + ra[2] / 3600.) * 15.
    if dec[0] >= 0:
        target_dec = (dec[0] + dec[1] / 60. + dec[2] / 3600.)
    else:
        target_dec = (dec[0] - dec[1] / 60. - dec[2] / 3600.)

    lam = cc/im.freq * 10000 # Hz to cm to um

    hdu.header['CRVAL1'] = target_ra
    hdu.header['CRPIX1'] = (im.nx + 1.) / 2.
    hdu.header['CDELT1'] = dx / 3600
    hdu.header['CTYPE1'] = "RA---TAN"
    hdu.header['CUNIT1'] = "deg"

    hdu.header['CRVAL2'] = target_dec
    hdu.header['CRPIX2'] = (im.ny + 1.) / 2.
    hdu.header['CDELT2'] = dy / 3600
    hdu.header['CTYPE2'] = "DEC--TAN"
    hdu.header['CUNIT2'] = "deg"

    hdu.header['CRVAL3'] = lam[2]
    hdu.header['CRPIX3'] = 5
    hdu.header['CDELT3'] = lam[1]-lam[0]
    hdu.header['CTYPE3'] = 'WAVE'
    hdu.header['CUNIT3'] = 'um'

    hdu.header['UNITS'] = "Jy / arcsec^2"

    outFileName = fname[:-4]+'.fits'

    if os.path.isfile(outFileName):
        os.remove(outFileName)
    hdu.writeto(outFileName)


if __name__ == "__main__":
    fname = sys.argv[1]
    write_fits_mirisim(fname)
