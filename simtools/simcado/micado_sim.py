# MICADO imager simulation wrapper
import os
import sys
import numpy as np

from astropy.io import fits

import simcado
from simcado.source import source_from_image
from simcado import __data_dir__
print("SimCADO data directory is", __data_dir__)


def readSpectrum(specfname):
    with open(specfname, 'r') as rfile:
        # Read the format number
        dum = rfile.readline()
        # Read the number of wavelengths
        nwav = int(rfile.readline())
        # Read a blank line
        dum = rfile.readline()

        res = np.zeros([nwav, 2], dtype=np.float64)
        for iwav in range(nwav):
            dum = rfile.readline().split()
            res[iwav, 0] = float(dum[0])
            res[iwav, 1] = float(dum[1])
    return res

def sim_wrapper(fname, starflag):
    '''
    e.g.
    fname = '10jup50au/npix1733lam214inc0.fits'
    lams = ['125','164','214']
    '''
    lam = fname[fname.find('lam')+3:fname.find('inc')]
    lam_ref = int(lam)/100
    inc = fname[fname.find('inc')+3:fname.find('.fits')]
    npix = fname[fname.find('npix')+4:fname.find('lam')]
    dir = fname.split('/npix')[0]
    datadir = '/home/ipa/ipaszulagyi/users/xuechen/thesis/RADMC3D/'

    ########## Set sim parameters ##########

    d = 100 # pc
    layout = 'centre'
    mode, micado_pixscale = 'zoom', 0.0015
    t_long = 3600  # second
    lam_min = 0.78
    lam_max = 2.5

    # filter to use
    if lam_ref == 1.25:
        filter = 'J'

    elif lam_ref == 1.64:
        filter = 'H'

    elif lam_ref == 2.14:
        filter = 'Ks'

    
    ########## Create source image ##########

    # Import fits image
    with fits.open(datadir+fname) as inhdul:
        sourceimg = inhdul[0].data[0]
        sourcehdr = inhdul[0].header

    # if coron, set star to 1%
    if starflag == 'coron':
        a = len(sourceimg[0])
        if (a % 2) == 0:
            sourceimg[int(a/2-1):int(a/2+1), int(a/2-1):int(a/2+1)] *= 0.01
        else:
            sourceimg[int((a-1)/2),int((a-1)/2)] *= 0.01

    try:
        print("Pixel unit:", sourcehdr['UNITS'])
    except:
        print("Pixel unit:", sourcehdr['BUNIT'])    

    # Check image pixel scale
    try:
        pixscale = sourcehdr['CDELT1'] * -3600
    except:
        pixscale = sourcehdr['CD1_1'] * -3600

    # If image pixscale not the same as the detector pixscale, need to specify oversample to interpolate between the two scales. 
    oversample = pixscale / micado_pixscale

    print("Image pixel scale:", pixscale, "arcsec/pixel, oversample:", oversample)

    # Import spectrum
    spectrum = readSpectrum(datadir+dir+'/sed_star.out')
    wav = spectrum[:,0]
    spec = spectrum[:,1] * 1e23 / d ** 2      # Conversion from erg s-1 cm-2 Hz-1 to Jy at d pc

    # Convert flux to photon counts
    from astropy import units as u
    lam_ref = lam_ref * u.um
    flam_ref = (1 * u.Jy).to(u.ph * u.s**-1 * u.m**-2 * u.um**-1, 
                            equivalencies=u.spectral_density(lam_ref))
    print("1 Jy at {0:5.2f} corresponds to {1:5.2f}".format(lam_ref, flam_ref))

    dlam = np.zeros(len(wav)-1)  # [um], spectral bin width
    for i in range(len(wav)-1):
        dlam[i] = wav[i+1] - wav[i]
    flux = spec[:-1] * flam_ref.value * dlam
    lams = wav[:-1]

    # Creating the source object
    src_disk = source_from_image(sourceimg, lams, flux, units="ph/s/m2", # spectra in [ph/s/m2/bin]
                                        plate_scale=pixscale, # The plate scale of the simulated image
                                        oversample=oversample, conserve_flux=True) 


    ########## Run full simcado simulation ##########
    sim_image = simcado.run(src_disk, mode=mode, detector_layout=layout, filter_name=filter, OBS_NDIT=1, OBS_DIT=t_long, SIM_LAM_MIN=lam_min, SIM_LAM_MAX=lam_max)
    #sim_image_2s = simcado.run(src_disk, mode=mode, detector_layout=layout, filter_name=filter, OBS_NDIT=1, OBS_DIT=t_short)

    
    ########## Crop final image ##########
    size = sim_image[0].data.shape[0]
    cen = int(size/2)
    r = int(int(npix)/2)
    if cen-r > 0:
        sim_image[0].data = sim_image[0].data[cen-r:cen+r,cen-r:cen+r]

    print("Saving image...")
    if starflag == 'coron':
        sim_image.writeto('./image{}_sim/{}/{}_lam{}inc{}_coron.fits'.format(lam, dir, filter, lam, inc), overwrite=True)
    elif starflag == 'star':
        sim_image.writeto('./image{}_sim/{}/{}_lam{}inc{}.fits'.format(lam, dir, filter, lam, inc), overwrite=True)
    print("ALL DONE.")


if __name__ == "__main__":
    fname = sys.argv[1]
    try:
        starflag = sys.argv[2]  # 'star', 'nostar', 'coron'
    except:
        starflag = 'star'
    sim_wrapper(fname, starflag)
