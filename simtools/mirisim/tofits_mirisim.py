"""This module contains classes/functions to create and read images with RADMC-3D and to calculate
interferometric visibilities and write fits files
For help on the syntax or functionality of each function see the help of the individual functions
"""
from __future__ import absolute_import
from __future__ import print_function
import traceback
import copy
import subprocess as sp
import os
import astropy.units as u

try:
    import numpy as np
except ImportError:
    np = None
    print(' Numpy cannot be imported ')
    print(' To use the python module of RADMC-3D you need to install Numpy')
    print(traceback.format_exc())

try:
    import scipy.special as spc
except ImportError:
    spc = None
    print('scipy.special cannot be imported ')
    print('This module is required to be able to calculate Airy-PSFs. Now PSF calculation is limited to Gaussian.')
    print(traceback.format_exc())

try:
    from astropy.io import fits as pf
except ImportError:
    print('astropy.io.fits cannot be imported trying pyfits')
    try:
        import pyfits as pf
    except ImportError:
        pf = None
        print('pyfits cannot be imported. Either of these modules is needed to write RADMC-3D images '
              + 'to FITS format. The rest of radmc3dPy can be used but fits output is now disabled.')
        print(traceback.format_exc())

try:
    import matplotlib.pylab as plt
except ImportError:
    plt = None
    print('Warning')
    print('matplotlib.pyplot cannot be imported')
    print('Without matplotlib you can use the python module to set up a model but you will not be able to plot things')
    print('or display images')

import natconst as nc

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
        self.imageJyppix = None
        self.imageJyarcsec2 = None
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
    def writeFits(self, fname='', dpc=1., coord='03h10m05s -10d05m30s', bandwidthmhz=2000.0,
                  casa=False, nu0=0., stokes='I', fitsheadkeys=[], ifreq=None):
        """Writes out a RADMC-3D image data in fits format.

        Parameters
        ----------

        fname        : str
                        File name of the radmc3d output image (if omitted 'image.fits' is used)

        dpc          : float
                        Distance of the source in pc

        coord        : str
                        Image center coordinates

        bandwidthmhz : float
                        Bandwidth of the image in MHz (equivalent of the CDELT keyword in the fits header)

        casa         : bool
                        If set to True a CASA compatible four dimensional image cube will be written

        nu0          : float
                        Rest frequency of the line (for channel maps)

        stokes       : {'I', 'Q', 'U', 'V', 'PI'}
                       Stokes parameter to be written if the image contains Stokes IQUV (possible
                       choices: 'I', 'Q', 'U', 'V', 'PI' -Latter being the polarized intensity)

        fitsheadkeys : dictionary
                        Dictionary containing all (extra) keywords to be added to the fits header. If
                        the keyword is already in the fits header (e.g. CDELT1) it will be updated/changed
                        to the value in fitsheadkeys, if the keyword is not present the keyword is added to
                        the fits header.

        ifreq        : int
                       Frequency index of the image array to write. If set only this frequency of a multi-frequency
                       array will be written to file.
        """
        # --------------------------------------------------------------------------------------------------
        istokes = 0


        if self.stokes:
            if fname == '':
                fname = 'image_stokes_' + stokes.strip().upper() + '.fits'

            if stokes.strip().upper() == 'I':
                istokes = 0
            if stokes.strip().upper() == 'Q':
                istokes = 1
            if stokes.strip().upper() == 'U':
                istokes = 2
            if stokes.strip().upper() == 'V':
                istokes = 3
        else:
            if fname == '':
                fname = 'image.fits'

        # Decode the image center cooridnates
        # Check first whether the format is OK
        dum = coord

        ra = []
        delim = ['h', 'm', 's']
        for i in delim:
            ind = dum.find(i)
            if ind <= 0:
                msg = 'coord keyword has a wrong format. The format should be coord="0h10m05s -10d05m30s"'
                raise ValueError(msg)
            ra.append(float(dum[:ind]))
            dum = dum[ind + 1:]

        dec = []
        delim = ['d', 'm', 's']
        for i in delim:
            ind = dum.find(i)
            if ind <= 0:
                msg = 'coord keyword has a wrong format. The format should be coord="0h10m05s -10d05m30s"'
                raise ValueError(msg)
            dec.append(float(dum[:ind]))
            dum = dum[ind + 1:]

        target_ra = (ra[0] + ra[1] / 60. + ra[2] / 3600.) * 15.
        if dec[0] >= 0:
            target_dec = (dec[0] + dec[1] / 60. + dec[2] / 3600.)
        else:
            target_dec = (dec[0] - dec[1] / 60. - dec[2] / 3600.)

        if len(self.fwhm) == 0:
            # Conversion from erg/s/cm/cm/ster to Jy/pixel
            conv = self.sizepix_x * self.sizepix_y / (dpc * nc.pc)**2. * 1e23

            # Conversion from Jy/pixel to Jy/arcsec2
            fluxconv = 1 / (dx ** 2) 


        else:
            # If the image has already been convolved with a gaussian psf then self.image has
            # already the unit of erg/s/cm/cm/beam, so we need to multiply it by 10^23 to get
            # to Jy/beam
            conv = 1e23

        # Create the data to be written
        if casa:
            # Put the stokes axis to the 4th dimension
            # data = np.zeros([1, self.nfreq, self.ny, self.nx], dtype=float)
            data = np.zeros([1, self.nfreq, self.ny, self.nx], dtype=float)
            if self.nfreq == 1:
                data[0, 0, :, :] = self.image[:, :] * conv

            else:
                for inu in range(self.nfreq):
                    data[inu, 0, :, :] = self.image[:, :, inu] * conv
        else:
            data = np.zeros([self.nfreq, self.nx, self.ny], dtype=float)
            if self.stokes:
                if stokes.strip().upper() != 'PI':
                    if self.nfreq == 1:
                        data[0, :, :] = self.image[:, :, istokes, 0] * conv

                    else:
                        for inu in range(self.nfreq):
                            data[inu, :, :] = self.image[:, :, istokes, inu] * conv
                else:
                    if self.nfreq == 1:
                        data[0, :, :] = np.sqrt(self.image[:, :, 1, 0]**2 + self.image[:, :, 2, 0]**2) * conv

                    else:
                        for inu in range(self.nfreq):
                            data[inu, :, :] = np.sqrt(
                                self.image[:, :, 1, inu]**2 + self.image[:, :, 2, inu]**2) * conv

            else:
                if self.nfreq == 1:
                    data[0, :, :] = self.image[:, :, 0] * conv * fluxconv

                else:
                    for inu in range(self.nfreq):
                        data[inu, :, :] = self.image[:, :, inu] * conv * fluxconv

        if ifreq is not None:
            if len(data.shape) == 3:
                data = data[ifreq, :, :]

        naxis = len(data.shape)
        hdu = pf.PrimaryHDU(data.swapaxes(naxis - 1, naxis - 2))
        hdulist = pf.HDUList([hdu])

        hdulist[0].header.set('CRPIX1', (self.nx + 1.) / 2., ' ')
        hdulist[0].header.set('CDELT1', -self.sizepix_x / nc.au / dpc / 3600., '')
        # hdulist[0].header.set('CRVAL1', self.sizepix_x/1.496e13/dpc/3600.*0.5+target_ra, '')
        hdulist[0].header.set('CRVAL1', target_ra, '')
        hdulist[0].header.set('CUNIT1', '     DEG', '')
        hdulist[0].header.set('CTYPE1', 'RA---TAN', '')

        hdulist[0].header.set('CRPIX2', (self.ny + 1.) / 2., '')
        hdulist[0].header.set('CDELT2', self.sizepix_y / nc.au / dpc / 3600., '')
        # hdulist[0].header.set('CRVAL2', self.sizepix_y/1.496e13/dpc/3600.*0.5+target_dec, '')
        hdulist[0].header.set('CRVAL2', target_dec, '')
        hdulist[0].header.set('CUNIT2', '     DEG', '')
        hdulist[0].header.set('CTYPE2', 'DEC--TAN', '')

        # For ARTIST compatibility put the stokes axis to the 4th dimension
        if casa:
            hdulist[0].header.set('CRPIX4', 1., '')
            hdulist[0].header.set('CDELT4', 1., '')
            hdulist[0].header.set('CRVAL4', 1., '')
            hdulist[0].header.set('CUNIT4', '        ', '')
            hdulist[0].header.set('CTYPE4', 'STOKES  ', '')

            if self.nwav == 1:
                hdulist[0].header.set('CRPIX3', 1.0, '')
                hdulist[0].header.set('CDELT3', bandwidthmhz * 1e6, '')
                hdulist[0].header.set('CRVAL3', self.freq[0], '')
                hdulist[0].header.set('CUNIT3', '      HZ', '')
                hdulist[0].header.set('CTYPE3', 'FREQ-LSR', '')

            else:
                if ifreq is None:
                    hdulist[0].header.set('CRPIX3', 1.0, '')
                    hdulist[0].header.set('CDELT3', (self.freq[1] - self.freq[0]), '')
                    hdulist[0].header.set('CRVAL3', self.freq[0], '')
                    hdulist[0].header.set('CUNIT3', '      HZ', '')
                    hdulist[0].header.set('CTYPE3', 'FREQ-LSR', '')
                    hdulist[0].header.set('RESTFRQ', self.freq[0])
                else:
                    hdulist[0].header.set('CRPIX3', 1.0, '')
                    hdulist[0].header.set('CDELT3', bandwidthmhz * 1e6, '')
                    hdulist[0].header.set('CRVAL3', self.freq[ifreq], '')
                    hdulist[0].header.set('CUNIT3', '      HZ', '')
                    hdulist[0].header.set('CTYPE3', 'FREQ-LSR', '')

        else:
            if self.nwav == 1:
                hdulist[0].header.set('CRPIX3', 1.0, '')
                hdulist[0].header.set('CDELT3', bandwidthmhz * 1e6, '')
                hdulist[0].header.set('CRVAL3', self.freq[0], '')
                hdulist[0].header.set('CUNIT3', '      HZ', '')
                hdulist[0].header.set('CTYPE3', 'FREQ-LSR', '')
            else:
                if ifreq is None:
                    hdulist[0].header.set('CRPIX3', 1.0, '')
                    hdulist[0].header.set('CDELT3', (self.freq[1] - self.freq[0]), '')
                    hdulist[0].header.set('CRVAL3', self.freq[0], '')
                    hdulist[0].header.set('CUNIT3', '      HZ', '')
                    hdulist[0].header.set('CTYPE3', 'FREQ-LSR', '')
                else:
                    hdulist[0].header.set('CRPIX3', 1.0, '')
                    hdulist[0].header.set('CDELT3', bandwidthmhz * 1e6, '')
                    hdulist[0].header.set('CRVAL3', self.freq[ifreq], '')
                    hdulist[0].header.set('CUNIT3', '      HZ', '')
                    hdulist[0].header.set('CTYPE3', 'FREQ-LSR', '')

        if nu0 > 0:
            hdulist[0].header.set('RESTFRQ', nu0, '')
        else:
            if self.nwav == 1:
                hdulist[0].header.set('RESTFRQ', self.freq[0], '')

        if len(self.fwhm) == 0:
            hdulist[0].header.set('BUNIT', 'Jy / pixel', '')
        else:
            hdulist[0].header.set('BUNIT', 'JY/BEAM', '')
            hdulist[0].header.set('BMAJ', self.fwhm[0] / 3600., '')
            hdulist[0].header.set('BMIN', self.fwhm[1] / 3600., '')
            hdulist[0].header.set('BPA', -self.pa, '')

        hdulist[0].header.set('BTYPE', 'INTENSITY', '')
        hdulist[0].header.set('BZERO', 0.0, '')
        hdulist[0].header.set('BSCALE', 1.0, '')

        hdulist[0].header.set('EPOCH', 2000.0, '')
        hdulist[0].header.set('LONPOLE', 180.0, '')

        if fitsheadkeys:
            if len(fitsheadkeys.keys()) > 0:
                for ikey in fitsheadkeys.keys():
                    # hdulist[0].header.update(ikey, fitsheadkeys[ikey], '')
                    hdulist[0].header.set(ikey, fitsheadkeys[ikey], '')

        if os.path.exists(fname):
            print(fname + ' already exists')
            dum = input('Do you want to overwrite it (yes/no)?')
            if (dum.strip()[0] == 'y') | (dum.strip()[0] == 'Y'):
                os.remove(fname)
                hdu.writeto(fname)
            else:
                print('No image has been written')
        else:
            hdu.writeto(fname)
            # --------------------------------------------------------------------------------------------------
    
    def readImage(self, fname=None, binary=False, old=False):
        """Reads an image calculated by RADMC-3D

        Parameters
        ----------

        fname   : str, optional
                 File name of the radmc3d output image (if omitted 'image.out' is used)

        old     : bool
                 If set to True it reads old radmc-2d style image

        binary  : bool, optional
                 False - the image format is formatted ASCII if True - C-compliant binary (omitted if old=True)
        """
        if old:
            if fname is None:
                fname = 'image.dat'

            self.filename = fname
            print('Reading ' + fname)


            with open(fname, 'r') as rfile:

                dum = rfile.readline().split()
                self.nx = int(dum[0])
                self.ny = int(dum[1])
                self.nfreq = int(dum[2])
                self.nwav = self.nfreq

                dum = rfile.readline().split()
                self.sizepix_x = float(dum[0])
                self.sizepix_y = float(dum[1])
                self.wav = np.zeros(self.nwav, dtype=float) - 1.
                self.freq = np.zeros(self.nwav, dtype=float) - 1.

                self.stokes = False
                self.image = np.zeros([self.nx, self.ny, self.nwav], dtype=np.float64)
                for iwav in range(self.nwav):
                    dum = rfile.readline()
                    for iy in range(self.ny):
                        for ix in range(self.nx):
                            self.image[ix, iy, iwav] = float(rfile.readline())

        else:
            if binary:
                if fname is None:
                    fname = 'image.bout'

                self.filename = fname

                dum = np.fromfile(fname, count=4, dtype=int)
                iformat = dum[0]
                self.nx = dum[1]
                self.ny = dum[2]
                self.nfreq = dum[3]
                self.nwav = self.nfreq
                dum = np.fromfile(fname, count=-1, dtype=np.float64)

                self.sizepix_x = dum[4]
                self.sizepix_y = dum[5]
                self.wav = dum[6:6 + self.nfreq]
                self.freq = nc.cc / self.wav * 1e4

                if iformat == 1:
                    self.stokes = False
                    self.image = np.reshape(dum[6 + self.nfreq:], [self.nfreq, self.ny, self.nx])
                    self.image = np.swapaxes(self.image, 0, 2)
                elif iformat == 3:
                    self.stokes = True
                    self.image = np.reshape(dum[6 + self.nfreq:], [self.nfreq, 4, self.ny, self.nx])
                    self.image = np.swapaxes(self.image, 0, 3)
                    self.image = np.swapaxes(self.image, 1, 2)

            else:

                # Look for the image file

                if fname is None:
                    fname = 'image.out'

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
                    self.freq = nc.cc / self.wav * 1e4

                    # If we have a normal total intensity image
                    if iformat == 1:
                        self.stokes = False

                        self.image = np.zeros([self.nx, self.ny, self.nwav], dtype=np.float64)
                        for iwav in range(self.nwav):
                            # Blank line
                            dum = rfile.readline()
                            for iy in range(self.ny):
                                for ix in range(self.nx):
                                    self.image[ix, iy, iwav] = float(rfile.readline())

                    # If we have the full stokes image
                    elif iformat == 3:
                        self.stokes = True
                        self.image = np.zeros([self.nx, self.ny, 4, self.nwav], dtype=np.float64)
                        for iwav in range(self.nwav):
                            # Blank line
                            dum = rfile.readline()
                            for iy in range(self.ny):
                                for ix in range(self.nx):
                                    dum = rfile.readline().split()
                                    imstokes = [float(i) for i in dum]
                                    self.image[ix, iy, 0, iwav] = float(dum[0])
                                    self.image[ix, iy, 1, iwav] = float(dum[1])
                                    self.image[ix, iy, 2, iwav] = float(dum[2])
                                    self.image[ix, iy, 3, iwav] = float(dum[3])

        # Conversion from erg/s/cm/cm/Hz/ster to Jy/pixel
        conv1 = self.sizepix_x * self.sizepix_y / (dpc * nc.pc)**2. * 1e23
        conv2 = 1e23 / 4.25e10
        self.imageJyppix = self.image * conv1
        self.imageJyarcsec2 = self.image * conv2

        self.x = ((np.arange(self.nx, dtype=np.float64) + 0.5) - self.nx / 2) * self.sizepix_x
        self.y = ((np.arange(self.ny, dtype=np.float64) + 0.5) - self.ny / 2) * self.sizepix_y


def readImage(fname=None, binary=False, old=False):
    """Reads an image calculated by RADMC-3D.
       This function is an interface to radmc3dImage.readImage().

    Parameters
    ----------
        fname   : str, optional
                 File name of the radmc3d output image (if omitted 'image.out' is used)

        old     : bool
                 If set to True it reads old radmc-2d style image

        binary  : bool, optional
                 False - the image format is formatted ASCII if True - C-compliant binary (omitted if old=True)
    """

    dum = radmc3dImage()
    dum.readImage(fname=fname, binary=binary, old=old)
    return dum




dpc = 100
coord='03h10m05s -10d05m30s'
ifreq = None


im = readImage("image50au10_25lams.out")
data = np.ones([im.nfreq+1, im.nx, im.ny], dtype=float)

#spec = readSpectrum("spectrum.out")



# Calculate pixel scale
dx = im.sizepix_x / nc.au / dpc   # from cm to arcsec
dy = im.sizepix_y / nc.au / dpc

# Conversion from Jy/pixel to Jy/arcsec^2
fluxconv = 1 / (dx * dy) 

for inu in range(im.nfreq):
    data[inu, :, :] = im.imageJyppix[:, :, inu] * fluxconv 


naxis = len(data.shape)

hdu = pf.PrimaryHDU(data.swapaxes(naxis - 1, naxis - 2))


###Create fits header###

dum = coord

ra = []
delim = ['h', 'm', 's']
for i in delim:
    ind = dum.find(i)
    if ind <= 0:
        msg = 'coord keyword has a wrong format. The format should be coord="0h10m05s -10d05m30s"'
        raise ValueError(msg)
    ra.append(float(dum[:ind]))
    dum = dum[ind + 1:]

dec = []
delim = ['d', 'm', 's']
for i in delim:
    ind = dum.find(i)
    if ind <= 0:
        msg = 'coord keyword has a wrong format. The format should be coord="0h10m05s -10d05m30s"'
        raise ValueError(msg)
    dec.append(float(dum[:ind]))
    dum = dum[ind + 1:]

target_ra = (ra[0] + ra[1] / 60. + ra[2] / 3600.) * 15.
if dec[0] >= 0:
    target_dec = (dec[0] + dec[1] / 60. + dec[2] / 3600.)
else:
    target_dec = (dec[0] - dec[1] / 60. - dec[2] / 3600.)

lam = nc.cc/im.freq * 10000 # Hz to cm to um

hdu.header['CRVAL1'] = target_ra
hdu.header['CRPIX1'] = (im.nx + 1.) / 2.
hdu.header['CDELT1'] = dx / 3600
hdu.header['CTYPE1'] = "RA---TAN"
hdu.header['CUNIT1'] = "DEG"

hdu.header['CRVAL2'] = target_dec
hdu.header['CRPIX2'] = (im.ny + 1.) / 2.
hdu.header['CDELT2'] = dy / 3600
hdu.header['CTYPE2'] = "DEC--TAN"
hdu.header['CUNIT2'] = "DEG"

hdu.header['CRVAL3'] = lam[12]
hdu.header['CRPIX3'] = 12
hdu.header['CDELT3'] = lam[1]-lam[0]
hdu.header['CTYPE3'] = 'WAVE'
hdu.header['CUNIT3'] = 'um'

hdu.header['UNITS'] = "mJy / arcsec2"

outFileName = 'testfile_25lams.fits'
if os.path.isfile(outFileName):
 os.remove(outFileName)
hdu.writeto(outFileName)


