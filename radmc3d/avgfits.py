import numpy as np
from astropy.io import fits
import os
import sys

num = ['01','02','03','04','05']
dirs = ['10jup50au/2/','5jup50au/2/','10jup30au/2/','5jup30au/2/',
        '10jup50au/1/','5jup50au/1/','10jup30au/1/','5jup30au/1/',
        '10jup50au/0/','5jup50au/0/','10jup30au/0/','5jup30au/0/',
        '10jup50au/5/','5jup50au/5/','10jup30au/5/','5jup30au/5/']

def avg_fits(lam, cube):
    if cube == False:
        for dir in dirs:
            lam = str(lam)
            # find npix
            fnames = os.listdir('./'+dir)
            npixs = []
            for name in fnames:
                if ('npix' in name) and ('lam'+lam in name) and ('cube' not in name) and ('_0' in name):
                    npix = name[name.find('npix')+4:name.find('lam')]
                    if npix not in npixs:
                        npixs.append(npix)
            npixs.sort(reverse=True)

            # average images for npix_n
            for n in npixs:
                for inc in ['0','30','60']:
                    img = []
                    fname = './'+dir+'npix'+n+'lam'+lam+'inc'+inc
                    for i in range(5):
                        if os.path.exists(fname+'_'+num[i]+'.fits'):
                            with fits.open(fname+'_'+num[i]+'.fits') as inhdul:
                                print('reading '+fname+'_'+num[i]+'.fits ...')
                                hdr = inhdul[0].header
                                img.append(inhdul[0].data[0])
                    img = np.array(img)
                    imgav = np.median(img, axis=0)
                    print('writing to '+fname+'.fits ...')
                    hdu = fits.PrimaryHDU([imgav])
                    hdulist = fits.HDUList([hdu])
                    hdulist[0].header = hdr
                    hdu.writeto(fname+'.fits', overwrite=True)

    else:
        for dir in dirs:
            lam = str(lam)
            # find npix
            fnames = os.listdir('./'+dir)
            npixs = []
            for name in fnames:
                if ('npix' in name) and ('lam'+lam in name) and ('cube5' in name) and ('_0' in name):
                    npix = name[name.find('npix')+4:name.find('lam')]
                    if npix not in npixs:
                        npixs.append(npix)
            npixs.sort(reverse=True)

            # average images for npix_n
            for n in npixs:
                for inc in ['0','30','60']:
                    img = []
                    fname = './'+dir+'cube5npix'+n+'lam'+lam+'inc'+inc
                    for i in range(5):
                        if os.path.exists(fname+'_'+num[i]+'.fits'):
                            with fits.open(fname+'_'+num[i]+'.fits') as inhdul:
                                print('reading '+fname+'_'+num[i]+'.fits ...')
                                hdr = inhdul[0].header
                                img.append(inhdul[0].data)
                    img = np.array(img)
                    imgav = np.median(img, axis=0)
                    print('writing to '+fname+'.fits ...')
                    hdu = fits.PrimaryHDU(imgav)
                    hdulist = fits.HDUList([hdu])
                    hdulist[0].header = hdr
                    hdu.writeto(fname+'.fits', overwrite=True)

if __name__ == "__main__":
    lam = sys.argv[1]
    try:
        ifcube = sys.argv[2]
        if ifcube == 'cube':
            cube = True
    except:
        cube = False
    avg_fits(lam, cube)