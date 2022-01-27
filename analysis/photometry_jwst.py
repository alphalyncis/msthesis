# Calculating planet magnitudes
import sys
import numpy as np
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm      
import matplotlib.cbook

from astropy.io import fits
from astropy import units as u
from astropy import constants as const
import matplotlib.patches as patches

import csv


def calc_mag(fname, instru='nrc', casedir=''):
    '''
    fname = '10jup50au/F200W_a1_lam214inc0_coron.fits'
    instrus = ['nrc','nrs','miri']
    '''
    #--------------------- common data -----------------------
    print(f'Opening {instru} {fname}...')

    datadir = '/home/ipa/ipaszulagyi/users/xuechen/thesis/simtools/'
    lam = fname[fname.find('lam')+3:fname.find('inc')]
    lam_ref = int(lam)/100
    model = fname[:fname.find('/')]
    if instru == 'nrc':
        if lam_ref < 3:
            instru = 'nrcshort'
        else:
            instru = 'nrclong'

    pixscale = {'nrcshort': 0.031,
                'nrclong': 0.063,
                'nrs': 0.065,
                'nrsami': 0.065,
                'miri': 0.11} # arcsec/pix
    '''
    pixscale = {'nrc': 0.031/0.063,
                'nrs': 0.065,
                'miri': 0.11}          # arcsec/pix
    '''


    psf_fwhm = {'nrclong':{'125': 1.298 ,'164': 1.628, '214': 2.304, '378': 1.901, '476': 2.574}, # pixel
                'nrcshort':{'125': 1.298 ,'164': 1.628, '214': 2.304, '378': 1.901, '476': 2.574}, # pixel
                'nrs':{'125': 0.627, '164': 0.770, '214': 0.991, '378': 1.847, '476': 2.319}, # pixel
                'nrsami':{'277':1.318, '378':1.847, '430':2.067, '476':2.319}, # pixel
                'miri':{'1020': 2.888, '1500': 4.354, '2100': 5.989}} # pixel
    
    snr_crit = {'nrclong': 1.9,
                'nrcshort': 1.9,
                'nrs': 1.5,
                'nrsami': 1.5,
                'miri': 1.78}


    #--------------------- open fits files -----------------------

    if instru == 'nrclong' or instru == 'nrcshort':
        imgdir = f'mirage/NIRCAM/image{lam}_sim/'

    elif instru == 'nrs' or instru == 'nrsami':
        imgdir = f'mirage/NIRISS/image{lam}_sim/'

    elif instru == 'miri':
        imgdir = f'miri/image{lam}_sim/'

    with fits.open(datadir+imgdir+fname) as hdulist:
        img = hdulist['SCI'].data
        hdr = hdulist[0].header
        if instru == 'nrclong':
            img = img[1:,1:]

    #--------------------- convert unit to Jy -----------------------

    imgJy = img * 1e6 / (4.25e10 / (pixscale[instru] ** 2)) # from MJy/sr to Jy/pix
    size = img.shape[0]
    bgrate = max(0, imgJy.min())
    print('bg rate is:', bgrate)

    # for miri an additional cutting is needed
    if instru == 'miri':
        imgJy = imgJy[int(size/2-11):int(size/2+11),int(size/2-11):int(size/2+11)]
        size = imgJy.shape[0]

    #--------------------- define aperture on planet -----------------------

    # estimate aperture size from psf fwhm
    if instru == 'nrclong': # for nrclong
        r = psf_fwhm['nrclong'][lam] * 0.95
    else:
        r = psf_fwhm[instru][lam]
    print('aperture size:', 2 * r)

    # get the planet position
    if instru == 'nrcshort': # for nrcshort
        xcen, ycen = 0.61 * size/2, size/2
        xacen = 1.39 * size/2
    else:
        xcen, ycen = 0.6 * size/2, size/2
        xacen = 1.4 * size/2

    # save fig with patch marked
    fig, ax = plt.subplots(figsize=(5,5))
    rec = patches.Rectangle((xcen-r-0.5, ycen-r-0.5), 2*r, 2*r, linewidth=0.5,edgecolor='r',facecolor='none')
    ax.add_patch(rec)
    im = ax.imshow(imgJy, norm=LogNorm(vmin=1e-8, vmax=1e-2),origin='lower')
    fig.colorbar(im, ax=ax).set_label('Jy')
    plt.savefig(f"photom_plots/{instru}{casedir}_{model}_{fname[fname.find('/')+3:-5]}")

    plt.close()

    # CPD
    cpd = imgJy[int(ycen-r+0.5):int(ycen+r+0.5), int(xcen-r+0.5):int(xcen+r+0.5)]

    # replace saturated pixel
    flag = ' '
    if 'nrc' in instru:
        if casedir == 0:
            pass
        elif casedir == '5' or casedir =='1':
            if ('10jup' in fname) or ('5jup' in fname):
                if 0 in cpd:
                    flag = 'sat'
                    cpd[cpd==0] = imgJy.max()
        elif casedir == '2':
            if '30au' in fname:
                if 0 in cpd:
                    flag = 'sat'
                    cpd[cpd==0] = imgJy.max()
            elif '50au' in fname:
                if ('10jup' in fname) or ('5jup' in fname) or ('1jup' in fname):
                    if 0 in cpd:
                        flag = 'sat'
                        cpd[cpd==0] = imgJy.max()  
    elif instru == 'nrsami':
        if casedir == '2':
            if ('10jup' in fname) or ('5jup' in fname):
                if 0 in cpd:
                    flag = 'sat'
                    cpd[cpd==0] = imgJy.max()                
    elif instru == 'miri':
        if casedir == '0':
            pass
        elif casedir == '1' or casedir =='0':
            if '10jup30au' in fname:
                if 0 in cpd:
                    flag = 'sat'
                    cpd[cpd==0] = imgJy.max()
        elif casedir == '2':
            if ('10jup' in fname) or ('5jup' in fname):
                if 0 in cpd:
                    flag = 'sat'
                    cpd[cpd==0] = imgJy.max()                                    


    # antiCPD
    acpd = imgJy[int(ycen-r+0.5):int(ycen+r+0.5), int(xacen-r+0.5):int(xacen+r+0.5)]

    # bg
    bg = bgrate * np.ones(cpd.shape)

    # integrate over aperture
    f_cpd_tot = cpd.sum()
    f_cpd_minus_bg = cpd.sum() - bg.sum()
    f_acpd_minus_bg = acpd.sum() - bg.sum()
    f_cpd_minus_acpd = cpd.sum() - abs(acpd.sum())
    f_bg = bg.sum()
    snr_cpd = f_cpd_minus_bg / abs(f_acpd_minus_bg)

    # mark each image with flag
    if flag != 'sat':
        if f_cpd_minus_acpd < 0 or snr_cpd < 1:
            flag = 'nondetect'
        elif snr_cpd > snr_crit[instru]:
            flag = 'detect'
        else:
            if 'nrc' in instru:
                flag = 'nondetect'
            elif instru == 'miri' or instru == 'nrsami':
                flag = 'asymm' if snr_cpd > 1.1 else 'nondetect'
                

    # to magnitudes
    zps = {'125':1587, '164':1074, '214':653,'277':427, '378':253, '430':194, '476':150, '1020': 34.9,'1500':18,'2100':8}
    m_cpd_minus_bg = 2.5 * np.log10(zps[lam]/f_cpd_minus_bg)
    m_acpd_minus_bg = 2.5 * np.log10(zps[lam]/f_acpd_minus_bg)
    m_cpd_minus_acpd = 2.5 * np.log10(zps[lam]/f_cpd_minus_acpd)


    print('CPD tot:',f_cpd_tot,'Background:',f_bg, 'CPD-bg:', f_cpd_minus_bg, 'antiCPD-bg:', f_acpd_minus_bg, 'CPD-antiCPD:',f_cpd_minus_acpd,
           'mag_CPD-bg:', m_cpd_minus_bg, 'mag_CPD-antiCPD:', m_cpd_minus_acpd, 'SNR:', snr_cpd, 'flag:', flag)

    return f_cpd_tot, f_bg, f_cpd_minus_bg, f_acpd_minus_bg, f_cpd_minus_acpd, m_cpd_minus_bg, m_cpd_minus_acpd, snr_cpd, flag

def make_table(instru, casedir, inc):  

    models = ['10jup50au','5jup50au','1jup50au','1sat50au','10jup30au','5jup30au','1jup30au']
    filts = {'nrc': ['F115W','F150W','F210M','F360M','F480M'],
             'nrs': ['F115W','F150W','F200W','F380M','F480M'],
             'nrsami':['F277W','F380M','F430M','F480M'],
             'miri': ['F1000W','F1500W','F2100W']}

    lams = {'nrc': ['125','164','214','378','476'],
            'nrs': ['125','164','214','378','476'],
            'nrsami':['277','378','430','476'],
            'miri': ['1020','1500','2100']}

    table = []
    for model in models:
        for filt, lam in zip(filts[instru], lams[instru]):
            if instru == 'nrsami':
                fname = f'{model}/{casedir}/{filt}_lam{lam}inc{inc}_ami_cal.fits'
            else:
                fname = f'{model}/{casedir}/{filt}_lam{lam}inc{inc}_coron_cal.fits'
            mags = list(calc_mag(fname, instru, casedir=casedir))
            mags = ['{:.2e}'.format(i) for i in mags[:-4]] + ['{:.2f}'.format(i) for i in mags[-4:-1]] + [mags[-1]]

            # negative fluxes are ignored
            # for i in range(len(mags)):
            #     if mags[i] == 'nan' or mags[i] == 'inf' or (mags[i] != 'sat' and eval(mags[i]) < 0):
            #         mags[i] = '/'

            # snr < 1 means no detection
            # if mags[-2] == '/' or eval(mags[-2]) <= 1:
            #     mags[-2] = '/'
            line = [model, filt] + mags
            table.append(line)
    tabname = f"flux_tables/table_{instru}{casedir}_inc{inc}_pipl.csv"

    with open(tabname,"w+") as f:
        csvWriter = csv.writer(f,delimiter=',')
        csvWriter.writerow(['model','filter','CPD tot','Background','CPD-bg', 'antiCPD-bg', 'CPD-antiCPD','mag CPD-bg','mag CPD-antiCPD','SNR_CPD','flag'])
        csvWriter.writerows(table)
    

def test():
    fnames = ['5jup50au/F200W_lam214inc60_coron.fits', # nrcshort
              '5jup50au/F356W_lam378inc60_coron.fits', # nrclong
              '5jup50au/F200W_lam214inc60_coron.fits',        # nrs
              '5jup50au/F1000W_lam1020inc60_coron.fits',   # miri
              '5jup50au/Ks_lam214inc60_coron.fits',             # micado
              '5jup50au/L_lam378inc60_coron.fits',              # metisLM
              '5jup50au/N2_lam1020inc60_coron.fits',            # metisN
              '5jup50au/hK_lam214inc60_coron.fits']             # gmt
    instrus = ['nrcshort','nrclong','nrs','miri','micado','metisLM','metisN','gmt']
    for fname, instru in zip(fnames,instrus):
        calc_mag(fname, instru)



if __name__ == "__main__":
    instru = sys.argv[1]
    casedir = sys.argv[2]
    inc = sys.argv[3]

    make_table(instru, casedir, inc)