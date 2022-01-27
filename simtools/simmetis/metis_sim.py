# METIS imager simulation wrapper
import os
import sys
import numpy as np

from astropy.io import fits

import simmetis as sim


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
    fname = '10jup50au/cube25npix23lam1020inc0.fits'
    lams = ['1020','1500','2100']
    '''
    lam = fname[fname.find('lam')+3:fname.find('inc')]
    lam_ref = int(lam)/100
    inc = fname[fname.find('inc')+3:fname.find('.fits')]
    npix = fname[fname.find('npix')+4:fname.find('lam')]
    dir = fname.split('/npix')[0]
    datadir = '/home/ipa/ipaszulagyi/users/xuechen/thesis/RADMC3D/'

    ########## Set sim parameters ##########

    d = 100 # pc

    if lam  == '378':
        cmds = sim.UserCommands('metis_image_LM.config')
        cmds['INST_FILTER_TC'] = 'TC_filter_L.dat'
        cmds['SIM_LAM_MIN'] = 3.38
        cmds['SIM_LAM_MAX'] = 4.31
        filter = 'L'
        metis_scale = 0.00525

    elif lam == '476':
        cmds = sim.UserCommands('metis_image_LM.config')
        cmds['INST_FILTER_TC'] = 'TC_filter_Mp.dat'
        cmds['SIM_LAM_MIN'] = 4.38
        cmds['SIM_LAM_MAX'] = 5.21        
        filter = 'Mp'
        metis_scale = 0.00525

    elif lam == '1020':
        cmds = sim.UserCommands('metis_image_N.config')
        cmds['INST_FILTER_TC'] = 'TC_filter_N2.dat'
        cmds['SIM_LAM_MIN'] = 10.15
        cmds['SIM_LAM_MAX'] = 13.10
        filter = 'N2'
        metis_scale = 0.00679

    
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
    oversample = pixscale / metis_scale

    print("Image pixel scale:", pixscale, "arcsec/pixel, oversample:", oversample)

    # Import spectrum
    spectrum = readSpectrum(datadir+dir+'/sed_star.out'.format(inc))
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
    src = sim.source.source_from_image(sourceimg, lams, flux, units="ph/s/m2", # spectra in [ph/s/m2/bin]
                                        plate_scale=pixscale, # The plate scale of the simulated image
                                        oversample=oversample, conserve_flux=True) 


    ########## Run full simmetis simulation ##########
    # Build the optical train and the detector
    opttrain = sim.OpticalTrain(cmds)
    fpa = sim.Detector(cmds, small_fov=True)

    # Apply optical train to source object
    src.apply_optical_train(opttrain, fpa)

    # Read out detectors, results are FITS HDU objects
    sim_image_2h = fpa.read_out(OBS_DIT=7200)

    ########## Crop final image ##########
    size = sim_image_2h[0].data.shape[0]
    cen = int(size/2)
    r = int(int(npix)/2)

    if cen-r > 0:
        sim_image_2h[0].data = sim_image_2h[0].data[cen-r:cen+r,cen-r:cen+r]

    print("Saving image...")
    if starflag == 'coron':
        sim_image_2h.writeto('./image{}_sim/{}/{}_lam{}inc{}_coron.fits'.format(lam, dir, filter, lam, inc), overwrite=True)
    elif starflag == 'star':
        sim_image_2h.writeto('./image{}_sim/{}/{}_lam{}inc{}.fits'.format(lam, dir, filter, lam, inc), overwrite=True)
    print("ALL DONE.")


if __name__ == "__main__":
    fname = sys.argv[1]
    try:
        starflag = sys.argv[2]  # 'star', 'nostar', 'coron'
    except:
        starflag = 'star'
    sim_wrapper(fname, starflag)

