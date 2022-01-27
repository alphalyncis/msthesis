# Plot image grids
import sys
import numpy as np
import matplotlib
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm  
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

import matplotlib.cbook
from matplotlib import colors,cm
cmap = plt.cm.inferno
cmap.set_bad(color=cmap(0))
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib import patches
matplotlib.rcParams.update({'font.size': 12})
from pylab import MaxNLocator
import matplotlib.ticker as ticker

from astropy.io import fits
from astropy import units as u
from astropy import constants as const

def norm1(img):
    min = np.median(img) + 0.5 * img.std()
    max = img.max() 
    norm = colors.LogNorm(min, max, clip=True)
    return norm

def norm2(img):
    min = np.median(img)
    max = img.max()
    norm = colors.LogNorm(min, max)
    return norm


def plot_lams(instru, casedir, inc='0', sep='50', vmin=2, vmax=7, sub=True, starflag='coron'):
    '''
    fname = '10jup50au/F200W_a1_lam214inc0_coron.fits'
    instru in ['nrcshort','nrclong','nrs','miri','micado','metisLM','metisN','gmt']
    lam in ['125','164','214','378','476','1020','1500','2100']
    '''
    #--------------------- common data -----------------------
    datadir = '/home/ipa/ipaszulagyi/users/xuechen/thesis/simtools/'
    models = {'50':['10jup50au','5jup50au','1jup50au','1sat50au'],'30':['10jup30au','5jup30au','1jup30au']}
    incs = ['0','30','60']
    text = ['10 $M_{jup}$','5 $M_{jup}$', '1 $M_{jup}$', '1 $M_{sat}$']
    extent = {'50':[-1.25, 1.25, -1.25 , 1.25], '30':[-0.75, 0.75, -0.75, 0.75]}

    QE = {'micado': {'125':0.88, '164':0.9, '214':0.85}, # simcado/data/TC_detector_H2RG.dat
          'metis': {'378':0.86, '476':0.71, '1020':0.8}} # SimMETIS-0.3/data/TC_detector_METIS_LM.dat(TC_detector_GeoSnap_N.dat)                                    # set by self

    gain = {'micado': 1.0,                               # fits header
            'metis': 1.0}                                # fits header

    D = {'micado': 39 * u.m,
         'metis': 39 * u.m}

    pixscale = {'nrcshort': 0.031,
                'nrclong': 0.063,
                'nrs': 0.065,
                'nrsami': 0.065,
                'miri': 0.11,
                'micado': 0.0015,
                'metisL': 0.00525,
                'metisN': 0.00679} # arcsec/pix

    if instru == 'nrc':
        insdir = 'mirage/NIRCAM'
        filts = ['F115W','F150W','F210M','F360M','F480M']
        lams = ['125','164','214','378','476']

    elif instru == 'nrs':
        insdir = 'mirage/NIRISS'
        filts = ['F115W','F150W','F200W','F356W','F444W']
        lams = ['125','164','214','378','476']

    elif instru == 'nrsami' or instru == 'nrsamiref':
        insdir = 'mirage/NIRISS'
        filts = ['F277W','F380M','F430M','F480M']
        lams = ['277','378','430', '476']

    elif instru == 'miri':
        insdir = 'miri'
        filts = ['F1000W','F1500W','F2100W']
        lams = ['1020','1500','2100']

    elif instru == 'micado':
        insdir = 'simcado'
        filts = ['J','H','Ks']
        lams = ['125','164','214']

    elif instru == 'metis':
        insdir = 'simmetis'
        filts = ['L','Mp','N2']
        lams = ['378','476','1020']

    elif instru == 'gmt':
        insdir = 'gmt'
        filts = ['hJ','hH','hK']
        lams = ['125','164','214']

    with plt.rc_context({'axes.edgecolor':'white', 'figure.facecolor':'white'}):
        if instru == 'miri' or instru == 'metis':
            axes_pad = 0.05
            cbar_mode = 'each'
            cbar_pad = 0
            cbar_size = '5%'
        else:
            axes_pad = 0
            cbar_mode = 'single'
            cbar_pad = 0.1
            cbar_size = '2%'

        fig = plt.figure(figsize=(17.5, 2.5*len(filts)))
        grid = ImageGrid(fig, 111, nrows_ncols=(len(filts), len(models[sep])), axes_pad=axes_pad, label_mode="L", share_all=True, 
                        cbar_location="right", cbar_mode=cbar_mode, cbar_size=cbar_size, cbar_pad=cbar_pad)

        if instru == 'micado' or instru == 'metis':
            suffix = 'coron'
        elif instru == 'nrsami':
            suffix = 'ami_cal'
        elif instru == 'nrsamiref':
            suffix = 'amiref_cal'
        else:
            suffix = 'coron_cal'
        for i in range(len(filts)):
            for j in range(len(models[sep])):
                ax = grid[i*len(models[sep])+j]

                # open file
                fname = f'{insdir}/image{lams[i]}_sim/{models[sep][j]}/{casedir}/{filts[i]}_lam{lams[i]}inc{inc}_{suffix}.fits'
                print('opening '+fname)
                with fits.open(fname) as inhdul:
                    if instru == 'gmt':
                        img = inhdul[0].data[0]
                    elif instru == 'nrc':
                        img = inhdul['SCI'].data
                        if lams[i] == '378' or lams[i] == '476':
                            if sep == '50':
                                img[13,32] = img[13,33]
                            img = img[1:,1:]
                    elif instru == 'nrs' or instru == 'nrsami' or instru == 'miri':
                        img = inhdul['SCI'].data
                    elif instru == 'nrsamiref':
                        img = inhdul['FIT'].data
                    else:
                        img = inhdul[0].data
                        hdr = inhdul[0].header
                        t_int = hdr['HIERARCH OBS_DIT']

                # subtract by background if needed
                size = img.shape[0]
                if sub:
                    img = img - (np.ones((size,size))*img.min() - 1)

                # crop image if needed
                if instru == 'nrs' and sep == '30':
                    img = img[:-1,:-1]
                elif instru == 'nrsami':
                    img = img[1:,1:]
                elif instru == 'nrsamiref':
                    if sep =='50':
                        img = img[int(size/2-20):int(size/2+20),int(size/2-20):int(size/2+20)]
                    elif sep == '30':
                        img = img[int(size/2-11):int(size/2+11),int(size/2-11):int(size/2+11)]               
                elif instru == 'miri':
                    if sep == '50':
                        img = img[int(size/2-11):int(size/2+11),int(size/2-11):int(size/2+11)]
                    elif sep == '30':
                        img = img[int(size/2-7):int(size/2+7),int(size/2-7):int(size/2+7)]               
                size = img.shape[0]

                # fix saturated pixel for entire sys plots
                xcen, ycen = 0.6 * size/2, size/2
                r = 5
                if 'nrc' in instru and casedir == '5':
                    if ('10jup' in fname) or ('5jup' in fname):
                        img[int(ycen-r+0.5):int(ycen+r+0.5), int(xcen-r+0.5):int(xcen+r+0.5)][img[int(ycen-r+0.5):int(ycen+r+0.5), int(xcen-r+0.5):int(xcen+r+0.5)] == 0] = img.max()

                # convert from MJy to Jy if jwst pipl
                if instru in ['nrsami']:
                    img = img/1e6
                # convert from DN to Jy/sr if needed
                if instru in ['metis']:
                    dlam = {'378': 0.6 * u.um, '476': 0.5 * u.um, '1020':(13.10-10.15) * u.um}
                    lam = fname[fname.find('lam')+3:fname.find('inc')]
                    lam_ref = int(lam)/100
                    Jy2phot = (1 * u.Jy).to(u.ph * u.s**-1 * u.m**-2 * u.um**-1, equivalencies=u.spectral_density(lam_ref * u.um))
                    Jy2DN = (Jy2phot * np.pi*(D[instru]/2)**2 * dlam[lam] * QE[instru][lam] * (1 / gain[instru])).value
                    imgJypix = img / t_int / Jy2DN 
                    # from Jy/pix to Jy/sr
                    pixscale = {'378': 0.00525, '476': 0.00525, '1020':0.00679}
                    imgJysr = imgJypix * (4.25e10 / (pixscale[lam] ** 2))
                    img = imgJysr

                # plot file
                if instru == 'miri':
                    im = ax.imshow(np.log10(img),vmin=np.log10(np.median(img)), cmap=cmap, extent=extent[sep], origin='lower')
                
                elif instru == 'metis':
                    img[int(size/2-5):int(size/2+5),int(size/2-5):int(size/2+5)]*=0
                    im = ax.imshow(np.log10(img),vmin=np.log10(np.median(img))-0.2, cmap=cmap, extent=extent[sep], origin='lower')
                
                elif instru == 'nrsami' or instru == 'nrsamiref':
                    im = ax.imshow(img, vmin=vmin, vmax=vmax, cmap=cmap, extent=extent[sep], origin='lower')
                
                # if instru == 'miri' and i < 1.1:
                #     im = ax.imshow(np.log10(img), vmin=vmin+2.7, vmax=vmax,cmap=cmap, extent=extent[sep], origin='lower')
                #     if sep == '50':
                #         ax.text(0.43, -1.15, 'vmin='+str(vmin+2.7), color='white',fontsize=12)
                #     elif sep == '30':
                #         ax.text(0.23, -0.65, 'vmin='+str(vmin+2.7), color='white',fontsize=12)
                else:
                    im = ax.imshow(np.log10(img), vmin=vmin, vmax=vmax,cmap=cmap, extent=extent[sep], origin='lower')


                # add a mask at center
                if instru=='nrsamiref' or instru=='nrsami':
                    pass
                else:
                    if sep == '50':
                        r = 0.17
                    elif sep == '30':
                        r = 0.1
                    if instru == 'nrs' and sep == '30':
                        mask = patches.Circle((0.04,0.03), radius=r, linewidth=0.5, edgecolor='white',facecolor=cmap(0.3) ,zorder=10)
                    else:
                        mask = patches.Circle((0.01,0), radius=r, linewidth=0.5, edgecolor='white',facecolor=cmap(0.3) ,zorder=10)
                    ax.add_patch(mask)

                ax.tick_params(axis='x', direction='in', color='white')
                ax.tick_params(axis='y', direction='in', color='white')
                ax.set_xlabel('x(arcsec)')
                ax.set_ylabel('y(arcsec)')

                if instru in ['nrc','nrsami','miri']:
                    filtext = filts[i]
                else:
                    filtext = f'{filts[i]} band'
                if cbar_mode == 'single':
                    fsize = 11 if instru == 'micado' else 12
                    if sep == '50':
                        ax.text(0.5, 1, text[j], color='white',fontsize=fsize)
                        ax.text(0.5, 0.76, filtext, color='white',fontsize=fsize)
                    elif sep == '30':
                        ax.text(0.3, 0.62, text[j], color='white',fontsize=fsize)
                        ax.text(0.3, 0.47, filtext, color='white',fontsize=fsize)
                elif cbar_mode == 'each':
                    fsize = 12
                    if sep == '50':
                        ax.text(0.15, 1, text[j], color='white',fontsize=fsize)
                        ax.text(0.15, 0.75, filtext, color='white',fontsize=fsize)
                    elif sep == '30':
                        ax.text(0.1, 0.6, text[j], color='white',fontsize=fsize)
                        ax.text(0.1, 0.45, filtext, color='white',fontsize=fsize)                    
                # if j==0 and sep == '50':
                #     if 'cal' in fname:
                #         ax.text(-2.3,-0.3, filts[i],rotation=90)
                #     else:
                #         ax.text(-2.3,-0.1, filts[i],rotation=90)
                
                if cbar_mode == 'each':
                    cb = grid.cbar_axes[i*len(models[sep])+j].colorbar(im) #, ticks=ticker.MultipleLocator(0.1)
                    cb.ax.tick_params(labelsize=7,left=True,right=False,length=2,color='white',pad=0.1,labelleft=True,labelright=False,labelcolor='white')
            if cbar_mode == 'each':
                if instru == 'micado' or instru == 'nrsami':
                    cb.set_label_text('log$_{10}$DN')
                if instru == 'metis':
                    cb.set_label_text('log$_{10}F_{\lambda}$ (Jy/sr)') 
                else:
                    cb.set_label_text('log$_{10}F_{\lambda}$ (MJy/sr)')
        if cbar_mode == 'single':
            #cb = grid.cbar_axes[0].colorbar(im)
            if instru == 'micado':
                cb = grid.cbar_axes[0].colorbar(im)
                cb.set_label_text('log$_{10}$DN')
            if instru == 'metis':
                cb = grid.cbar_axes[0].colorbar(im)
                cb.set_label_text('log$_{10}F_{\lambda}$ (Jy/sr)')            
            elif instru == 'nrsami' or instru == 'nrsamiref':
                cb = grid.cbar_axes[0].colorbar(im, format='%.2f')
                cb.set_label_text('$F_{\lambda}$ (Jy/sr)')
            else:
                cb = grid.cbar_axes[0].colorbar(im)
                cb.set_label_text('log$_{10}F_{\lambda}$ (MJy/sr)')

        if sep == '50':
            grid.axes_llc.set_xticks([-1, 0, 1])
            grid.axes_llc.set_yticks([-1, 0, 1])
        elif sep == '30':
            grid.axes_llc.set_xticks([-0.5, 0, 0.5])
            grid.axes_llc.set_yticks([-0.5, 0, 0.5]) 

    fig.suptitle(sep+'au')
    fig.tight_layout()
    
    if (not sub) and ('cal' in fname): # nrc nrsami miri
        plt.savefig('paper/'+instru+casedir+'_'+sep+'au_inc'+inc+'_nosub_cal', bbox_inches='tight', dpi=200)
    elif sub and ('cal' in fname):
        print('sub')
        plt.savefig('paper/'+instru+casedir+'_'+sep+'au_inc'+inc+'_sub_cal', bbox_inches='tight', dpi=200)
    elif (not sub) and ('cal' not in fname):
        plt.savefig('paper/'+instru+casedir+'_'+sep+'au_inc'+inc+'_nosub', bbox_inches='tight', dpi=200)
    elif sub and ('cal' not in fname): # metis micado
        print('sub')
        plt.savefig('paper/'+instru+casedir+'_'+sep+'au_inc'+inc+'_sub', bbox_inches='tight', dpi=150)
    
#plot_lams(instru='nrc', inc='0', sep='30',vmin=2,vmax=7)


if __name__ == "__main__":
    instru = sys.argv[1]
    casedir = sys.argv[2]
    inc = sys.argv[3]
    sep = sys.argv[4]
    vmin = eval(sys.argv[5])
    vmax = eval(sys.argv[6])
    try:
        ifsub = sys.argv[7]
        if ifsub == 'nosub':
            sub = False
    except:
        sub = True
    plot_lams(instru=instru, casedir=casedir, inc=inc, sep=sep,vmin=vmin,vmax=vmax,sub=sub)
