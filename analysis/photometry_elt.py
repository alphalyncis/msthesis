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


def calc_mag(fname, instru='nrcshort'):
    '''
    fname = '10jup50au/0/F200W_a1_lam214inc0_coron.fits'
    instrus = ['nrcshort','nrclong','nrs','miri','micado','metisLM','metisN','gmt']
    '''
    #--------------------- common data -----------------------
    datadir = '/home/ipa/ipaszulagyi/users/xuechen/thesis/simtools/'
    lam = fname[fname.find('lam')+3:fname.find('inc')]
    lam_ref = int(lam)/100
    model = fname[:fname.find('/')]
    casedir = fname[fname.find('/')+1:fname.find('/')+2]
    '''
    pixscale = {'micado': 0.0015,
                'metis': 0.00525/0.00679
                'gmt': 0.005}          # arcsec/pix
    '''

    psf_fwhm = {'micado':{'125': 0.006/0.0015, '164':0.008/0.0015 ,'214':0.011/0.0015},
                'metis':{'378': 0.02/0.00525,'476': 0.025/0.00525,'1020':0.06/0.00679}}

    QE = {'micado': {'125':0.88, '164':0.9, '214':0.85}, # simcado/data/TC_detector_H2RG.dat
          'metis': {'378':0.86, '476':0.71, '1020':0.8}} # SimMETIS-0.3/data/TC_detector_METIS_LM.dat(TC_detector_GeoSnap_N.dat)                                    # set by self

    gain = {'micado': 1.0,                               # fits header
            'metis': 1.0}                                # fits header

    D = {'micado': 39 * u.m,
         'metis': 39 * u.m}

    snr_crit = {'micado': 1.5,
                'metis': 1.5}
    
    Jy2phot = (1 * u.Jy).to(u.ph * u.s**-1 * u.m**-2 * u.um**-1, equivalencies=u.spectral_density(lam_ref * u.um))

    #--------------------- open fits files -----------------------
    print(f'Opening {instru} {fname}...')

    if instru == 'micado':
        imgdir = f'simcado/image{lam}_sim/'
        with fits.open(datadir+imgdir+fname) as hdulist:
            img = hdulist[0].data
            hdr = hdulist[0].header
        t_int = hdr['HIERARCH OBS_DIT']

        dlam = {'125':(1.345-1.15) * u.um,'164':(1.78-1.49) * u.um,'214':(2.32-1.97) * u.um}
        Jy2DN = (Jy2phot * np.pi*(D[instru]/2)**2 * dlam[lam] * QE[instru][lam] * (1 / gain[instru])).value

    elif instru == 'metis':
        imgdir = f'simmetis/image{lam}_sim/'
        with fits.open(datadir+imgdir+fname) as hdulist:
            img = hdulist[0].data
            hdr = hdulist[0].header
        t_int = hdr['HIERARCH OBS_DIT']        

        dlam = {'378': 0.6 * u.um, '476': 0.5 * u.um, '1020':(13.10-10.15) * u.um}
        Jy2DN = (Jy2phot * np.pi*(D[instru]/2)**2 * dlam[lam] * QE[instru][lam] * (1 / gain[instru])).value

    #--------------------- convert unit to Jy -----------------------

    imgJy = img / t_int / Jy2DN
    size = img.shape[0]

    # calculate bg rate by a patch on the corner
    if instru == 'micado':
        bgrate = np.median(imgJy[:5,:5])
    elif instru == 'metis':
        bgrate = np.min(imgJy[:5,:5])
    print('bg rate is:', bgrate)

    #--------------------- define aperture on planet -----------------------

    # estimate aperture size from psf fwhm
    if instru == 'micado':
        r = psf_fwhm[instru][lam] * 3
    elif instru == 'metis':
        r = psf_fwhm[instru][lam] * 1.5
    print('aperture size:', 2 * r)

    # get the planet position
    if instru == 'micado' and '50au' in fname:
        xcen, ycen = 0.616 * size/2, size/2
        xacen = 1.384 * size/2
    elif instru == 'metis' and int(lam) < 1000: # for metisLM
        if '50au' in fname:
            xcen, ycen = 0.63 * size/2, size/2
            xacen = 1.37 * size/2
        else:
            xcen, ycen = 0.615 * size/2, size/2
            xacen = 1.385 * size/2
    else:
        xcen, ycen = 0.6 * size/2, size/2
        xacen = 1.4 * size/2

    # save fig with patch marked
    fig, ax = plt.subplots(figsize=(5,5))
    rec = patches.Rectangle((xcen-r-0.5, ycen-r-0.5), 2*r, 2*r, linewidth=0.5,edgecolor='r',facecolor='none')
    ax.add_patch(rec)
    imgJysub = imgJy - bgrate*np.ones((size,size))
    ax.imshow(imgJysub, norm=LogNorm(),origin='lower')
    plt.savefig(f"photom_plots/{instru}{casedir}_{model}_{fname[fname.find('/')+3:-5]}")
    plt.close()

    # CPD
    cpd = imgJy[int(ycen-r+0.5):int(ycen+r+0.5), int(xcen-r+0.5):int(xcen+r+0.5)]

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
    flag = ' '
    if instru == 'micado':
        if f_cpd_minus_acpd < 0 or snr_cpd < 1.1:
            flag = 'nondetect'
        elif 1.1 < snr_cpd < snr_crit[instru]:
            flag = 'asymm'
        elif snr_cpd > snr_crit[instru]:
            flag = 'detect'
    elif instru == 'metis':
        if f_cpd_minus_acpd < 0 or snr_cpd <= 1.05:
            flag = 'nondetect'
        elif 1.05 < snr_cpd < snr_crit[instru]:
            flag = 'asymm'
        elif snr_cpd >= snr_crit[instru]:
            flag = 'detect'

    # to magnitudes
    zps = {'125':1587, '164':1074, '214':653, '378':253, '476':150, '1020': 34.9,'1500':18,'2100':8}
    m_cpd_minus_bg = 2.5 * np.log10(zps[lam]/f_cpd_minus_bg)
    m_acpd_minus_bg = 2.5 * np.log10(zps[lam]/f_acpd_minus_bg)
    m_cpd_minus_acpd = 2.5 * np.log10(zps[lam]/f_cpd_minus_acpd)


    print('CPD tot:',f_cpd_tot,'Background:',f_bg, 'CPD-bg:', f_cpd_minus_bg, 'antiCPD-bg:', f_acpd_minus_bg, 'CPD-antiCPD:',f_cpd_minus_acpd,
           'mag_CPD-bg:', m_cpd_minus_bg, 'mag_CPD-antiCPD', m_cpd_minus_acpd, 'SNR', snr_cpd, 'flag:', flag)

    return f_cpd_tot, f_bg, f_cpd_minus_bg, f_acpd_minus_bg, f_cpd_minus_acpd, m_cpd_minus_bg, m_cpd_minus_acpd, snr_cpd, flag

def make_table(instru, casedir, inc):  

    models = ['10jup50au','5jup50au','1jup50au','1sat50au','10jup30au','5jup30au','1jup30au']
    filts = {'micado':['J','H','Ks'],
             'metis': ['L','Mp','N2'],
             'gmt': ['hJ','hH','hK']}
    lams = {'micado':['125','164','214'],
            'metis': ['378','476','1020'],
            'gmt': ['125','164','214']}

    table = []
    for model in models:
        for filt, lam in zip(filts[instru], lams[instru]):
            fname = f'{model}/{casedir}/{filt}_lam{lam}inc{inc}_coron.fits'
            mags = list(calc_mag(fname, instru))
            mags = ['{:.2e}'.format(i) for i in mags[:-4]] + ['{:.2f}'.format(i) for i in mags[-4:-1]] + [mags[-1]]
            # # negative fluxes are ignored
            # for i in range(len(mags)):
            #     if mags[i] == 'nan' or eval(mags[i]) < 0:
            #         mags[i] = '/'
            # # snr < 1 means no detection
            # if mags[-1] == '/' or eval(mags[-1]) <= 1:
            #     mags[-1] = '/'
            line = [model, filt] + mags
            table.append(line)
    with open(f"flux_tables/table_{instru}{casedir}_inc{inc}.csv","w+") as f:
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