# GMTIFS imager simulation wrapper
import sys
import numpy as np

import poppy
from astropy import units as u
from astropy.io import fits
from astropy import constants as const
from scipy.signal import fftconvolve
from scipy.ndimage import gaussian_filter

def sim_wrapper(fname, starflag):
    '''
    e.g.
    fname = '10jup50au/npix520lam214inc0.fits'
    lams = ['125','164','214']
    '''
    lam = fname[fname.find('lam')+3:fname.find('inc')]
    lam_ref = int(lam)/100
    inc = fname[fname.find('inc')+3:fname.find('.fits')]
    npix = int(fname[fname.find('npix')+4:fname.find('lam')])
    dir = fname.split('/')[0]
    datadir = '/home/ipa/ipaszulagyi/users/xuechen/thesis/RADMC3D/'

    ########## Set sim parameters ##########

    D = 25.4 * u.m  # telescope effective primary aperture
    QE = 0.9        # quantum efficiency
    gain = 2        # sensor gain
    d = 100         # distance in pc
    t_exp = 3600    # exposure time in total
    coron = False

    # FOV
    if '50au' in dir:
        fov_arcsec = 2.5
    elif '30au' in dir:
        fov_arcsec = 1.5

    # filter to use
    if lam_ref == 1.25:
        filter = 'hJ'
        dlam = (1.35 - 1.10) * u.um
    elif lam_ref == 1.64:
        filter = 'hH'
        dlam = (1.80 - 1.47) * u.um
    elif lam_ref == 2.14:
        filter = 'hK'
        dlam = (2.51-2.04) * u.um
    
    print('Center wavelength is {} um.'.format(lam_ref))
    print('Using filter {}...'.format(filter))

    ########## Build GMT aperture ##########

    primary = poppy.optics.MultiCircularAperture(rings=1, segment_radius=8.365*u.m/2, gap=0.293*u.m)
    secondary = poppy.SecondaryObscuration(secondary_radius=3.2*u.m/2)
    gmt_aperture = poppy.CompoundAnalyticOptic(opticslist = [primary, secondary])

    ########## Compute GMT PSF ##########

    osys = poppy.OpticalSystem()
    osys.add_pupil(gmt_aperture)
    osys.add_detector(pixelscale=0.005, fov_arcsec=fov_arcsec)  # in arcsec
    gmt_psf = osys.calc_psf(lam_ref * 1e-6)

    ########## Create source image ##########

    # Import fits image
    with fits.open(datadir+fname) as inhdul:
        img = inhdul[0].data[0]
        hdr = inhdul[0].header

    # if coron, set star to 1%
    if starflag == 'coron':
        a = len(img[0])
        if (a % 2) == 0:
            img[int(a/2-1):int(a/2+1), int(a/2-1):int(a/2+1)] *= 0.01
        else:
            img[int((a-1)/2),int((a-1)/2)] *= 0.01

    try:
        print("Image unit is:", hdr['UNITS'])
    except:
        print("Image unit is:", hdr['BUNIT'])    

    # Check image pixel scale
    try:
        pixscale = hdr['CDELT1'] * -3600
    except:
        pixscale = hdr['CD1_1'] * -3600

    print("pixel scale is {} arcsec/pix".format(pixscale))

    # calc conv factor from Jy to photon count
    lam_ref = lam_ref * u.um
    E_ph = const.h * const.c / lam_ref.to(u.m)
    Jy2SI = (1 * u.Jy).to(u.J * u.s**-1 * u.m**-2 * u.um**-1, equivalencies=u.spectral_density(2* u.um))
    Jy2DN = (Jy2SI / E_ph) * (D/2)**2 * dlam * QE * (1 / gain)
    print("1 Jy at {0:5.2f} corresponds to {1:5.2f}".format(lam_ref, Jy2DN))

    # convert img to DN/s
    imgDN = img * Jy2DN.value * t_exp

    '''
    # coronagraph - dim the star by 99%
    if coron:
        coron_mask = np.ones((npix,npix))
        cen = int(npix/2)
        rmask = 1
        coron_mask[cen-rmask:cen+rmask, cen-rmask:cen+rmask] = 0.01

        imgDN = np.multiply(imgDN, coron_mask)
        print('Coronagraphy simulation applied.') 
    '''

    ########## Convolve image with PSF ##########

    print('Convolving image with PSF...')
    conv_img = fftconvolve(imgDN, gmt_psf[0].data, mode='same')

    # Add seeing halo
    conv_img = gaussian_filter(conv_img, sigma=5)


    ########## Add noise to concolved image ##########
    print('Adding noise... type: photon, background, det')
    # Add photon shot noise
    noise_phot = np.random.poisson(lam=conv_img)

    # Add background and detector noise
    noise_det = np.random.poisson(lam=np.ones((npix,npix)) * conv_img.max()/5)

    noise = noise_det + noise_phot
    det_img = noise + conv_img
    

    ########## Save final image ##########
    print("Saving image...")

    hdu = fits.PrimaryHDU([det_img, noise_phot, noise_det])

    hdu.header = hdr
    hdu.header['BUNIT'] = 'DN'
    hdu.header['EXPTIME'] = t_exp
    hdu.header['OTHER'] = '[det_img, noise_phot, noise_det]'
    if starflag == 'coron':
        hdu.writeto('./image{}_sim/{}/{}_lam{}inc{}_coron.fits'.format(lam, dir, filter, lam, inc), overwrite=True)
    elif starflag == 'star':
        hdu.writeto('./image{}_sim/{}/{}_lam{}inc{}.fits'.format(lam, dir, filter, lam, inc), overwrite=True)


if __name__ == "__main__":
    fname = sys.argv[1]
    try:
        starflag = sys.argv[2]  # 'star', 'nostar', 'coron'
    except:
        starflag = 'star'
    sim_wrapper(fname, starflag)
