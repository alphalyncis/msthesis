# FOR NIRISS: convert fits image in ADU unit
import sys
import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy import constants as const


def to_adu(fname,  starflag='star'):
    '''Included:
    - Telescope aperture (D)
    - Filter width (dlam)
    - Sensor gain in e-/ADU (gain)
    - Quantum efficiency (QE)

    NOT included:
    - Filter transmissivity

    fname = '10jup50au/npix84lam378inc0.fits'
    '''
    datadir = '/home/ipa/ipaszulagyi/users/xuechen/thesis/RADMC3D/'
    lam = fname[fname.find('lam')+3:fname.find('inc')]

    photfnu = { '125': 2.773018e-31,   #F115W_n
                '164': 3.123540e-31,   #F150W_n
                '214': 2.879494e-31,   #F200W_n
                '277': 2.827052e-31,   #F277W_n
                '378': 1.295626e-30,   #F380M_n
                '430': 1.613895e-30,   #F430M_n
                '476': 1.453752e-30 }  #F480M_n

    Jy2DN = 1e-23 / photfnu[lam]

    # convert fits files
    with fits.open(datadir+fname) as inhdul:
        img = inhdul[0].data[0]
        hdr = inhdul[0].header

    if starflag == 'star':
        pass
    elif starflag == 'nostar':
        a = len(img[0])
        if (a % 2) == 0:
            img[int(a/2-1):int(a/2+1), int(a/2-1):int(a/2+1)] = 0
        else:
            img[int((a-1)/2),int((a-1)/2)] = 0
    elif starflag == 'coron':
        a = len(img[0])
        if (a % 2) == 0:
            img[int(a/2-1):int(a/2+1), int(a/2-1):int(a/2+1)] *= 0.01
        else:
            img[int((a-1)/2),int((a-1)/2)] *= 0.01
    elif starflag == 'staronly':
        a = len(img[0])
        newimg = np.zeros(img.shape)
        if (a % 2) == 0:
            newimg[int(a/2-1):int(a/2+1), int(a/2-1):int(a/2+1)] = img[int(a/2-1):int(a/2+1), int(a/2-1):int(a/2+1)]
        else:
            newimg[int((a-1)/2),int((a-1)/2)] = img[int((a-1)/2),int((a-1)/2)]

    if starflag == 'staronly':
        imgADU = newimg * Jy2DN
    else:
        imgADU = img * Jy2DN

    hdu = fits.PrimaryHDU(imgADU)
    hdulist = fits.HDUList([hdu])
    hdulist[0].header = hdr

    if starflag =='star':
        print('writing to '+fname[:-5]+'ADU.fits ...')
        hdu.writeto('./image'+lam+'_data/'+fname[:-5]+'ADU.fits', overwrite=True)
    else:
        print('writing to '+fname[:-5]+'ADU_'+starflag+'.fits ...')
        hdu.writeto('./image'+lam+'_data/'+fname[:-5]+'ADU_'+starflag+'.fits', overwrite=True)



if __name__ == "__main__":
    fname = sys.argv[1]
    try:
        starflag = sys.argv[2]  # 'star', 'nostar', 'coron', 'staronly'
    except:
        starflag = 'star'
    to_adu(fname, starflag)