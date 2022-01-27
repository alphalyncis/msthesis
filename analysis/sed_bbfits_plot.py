# Remove photon noise from spectrum
import sys
import numpy as np
from astropy.io import fits
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm      
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
import astropy.units as u
from astropy.constants import h, c, k_B, b_wien
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter1d

fnames = ['/5/sed_nostar.out','/0/sed_nostar.out','/1/sed_nostar.out','/2/sed_nostar.out','/3/sed_nostar.out','/4/sed_nostar.out']
labels = ['entire system (CSD+CPD+planet)','CSD only','CPD+planet only','planet only','CPD only - with perfect absorber CSD','planet only - with perfect absorber CPD']
colors = ['#1f77b4','#8c564b', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

def read_spectrum(fname=''):
    with open(fname, 'r') as rfile:
        dum = rfile.readline()        # Read the format number
        nwav = int(rfile.readline())  # Read the number of wavelengths
        dum = rfile.readline()        # Read a blank line

        res = np.zeros([nwav, 2], dtype=np.float64)
        for iwav in range(nwav):
            dum = rfile.readline().split()
            res[iwav, 0] = float(dum[0])
            res[iwav, 1] = float(dum[1])
    return res

def plot_spectrum(arr, option='smooth', dpc=100, **kwargs):
    ''' option = 'full', 'clean', 'smooth' or both '''
    lumfact = 1e+23  # from cgs to Jy
    distfact = 1.0 / (dpc ** 2)
    x_coord = arr[:,0]
    y_coord = arr[:,1] * lumfact * distfact

    if 'full' in option:
        plt.plot(x_coord, y_coord, **kwargs)

    for i in range(2, len(y_coord)-2):
        if (y_coord[i] > 1.5*y_coord[i-1] and y_coord[i] > 1.5*y_coord[i-1]):
            y_coord[i] = (y_coord[i-2] + y_coord[i+2]) /2

    y_clean = np.array([y_coord[0]]+[np.min(y_coord[i-1:i+2]) for i in range(1, len(y_coord)-1)]+[y_coord[-1]])
    y_smooth = gaussian_filter1d(y_clean, sigma=5)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\lambda [\mu\mathrm{m}]$')
    plt.ylabel(r'$F_{\nu} [Jy]$')
    plt.ylim(1e-10,10)

    if 'clean' in option:
        plt.plot(x_coord, y_clean, **kwargs)
    if 'smooth' in option:
        plt.plot(x_coord, y_smooth, **kwargs)
    return np.array([x_coord, y_clean]).T

def blackbody(nu, T, amp):
    """ Blackbody as a function of wavelength (um) and temperature (K).
    returns in units of 
    """
    nu = nu * u.Hz
    T = T * u.K
    B_nu = amp * (2*h*nu**3 / c**2 / (np.exp(h*nu / (k_B*T)) - 1) /u.sr).to(u.Jy/u.sr) # W m-2 Hz-1 sr-1
    return B_nu.value

def fit_bb(specs_set, dir):
    plt.figure(figsize=(10,5))
    efftemps=[]
    for i in range(len(fnames)):
        plt.subplot(2,3,i+1)
        if i == 5 and dir in ['1jup30au', '1jup50au','1sat50au','5jup30au']:
            start = 200
        else:
            start = 20
        x = c.value/(specs_set[i][start:,0] * 1e-6)  # Hz
        y = specs_set[i][start:,1]
        coeff, cov = curve_fit(blackbody, x, y, p0=(200,1e-10), maxfev=2000)
        efftemps.append(coeff)
        #print(int(coeff[0]),coeff[1])
        xlam = 1e6*c.value/x
        plt.plot(xlam, blackbody(x, coeff[0],coeff[1]), zorder=10)
        plt.plot(xlam,y)
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(1e-10,10)
        plt.xlim(1,2e5)
        if i == 3 or i == 4 or i == 5:
            plt.xlabel(r'$\lambda [\mu\mathrm{m}]$')
        plt.title(labels[i], fontsize=10)
        if i == 0 or i == 3:
            plt.ylabel(r'$F_{\nu} [Jy]$')
        plt.text(1000,0.1,'T_eff = {} K'.format(int(coeff[0]+0.5)), color='#1f77b4')
    plt.tight_layout()
    plt.savefig('bbfit_plots/'+dir+'_fit.png')
    return efftemps

def plot_seds(dir):
    specs_full = []
    specs_clean = []
    for i in range(len(fnames)):
        specs_full.append(read_spectrum(dir+fnames[i]))
        specs_clean.append(plot_spectrum(specs_full[i]))

    # get eff temp of seds by bb fitting
    efftemps = fit_bb(specs_clean, dir)

    # FINAL PLOT
    plt.figure(figsize=(6,4.5))
    plt.ylim(1e-10,10)
    plt.xlim(0.5,2e4)
    specs_clean = []
    labels = ['entire system (CSD+CPD+planet)','CSD only','CPD+planet only','planet only','CPD only - with perfect absorber CSD','planet only - with perfect absorber CPD']
    colors = ['#1f77b4','#8c564b', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    for i in range(len(fnames)):
        if i == 0:
            zorder = 10
        else:
            zorder = 1
        specs_clean.append(plot_spectrum(specs_full[i], label='{} - {}K'.format(labels[i], int(efftemps[i][0]+0.5)), color=colors[i], zorder=zorder))

    plt.legend(loc=3, prop={'size': 9})
    #plt.title('SED '+dir)

    plt.gca().axvspan(5, 28, ymax=0.987, alpha=0.04, color='darkred')
    plt.text(5.2, 4, 'MIRI', alpha=0.3, color='darkred', fontsize=7.5)

    plt.gca().axvspan(3, 13, ymax=0.95, alpha=0.07, color='darkolivegreen')
    plt.text(3, 1.45, 'METIS', alpha=0.3, color='darkolivegreen', fontsize=7.5)

    plt.gca().axvspan(0.8, 2.45, ymax=0.95, alpha=0.06, color='chocolate')
    plt.text(0.8, 1.45, 'MICADO', alpha=0.4, color='chocolate', fontsize=7.5)

    #plt.gca().axvspan(0.8, 4.9, ymax=0.91, alpha=0.05, color='midnightblue')
    #plt.text(0.8, 0.55, 'NIRISS', alpha=0.3, color='midnightblue', fontsize=7.5)

    plt.gca().axvspan(0.6, 4.8, ymax=0.91, alpha=0.04, color='midnightblue')
    plt.text(0.6, 0.55, 'NIRCam', alpha=0.3, color='midnightblue', fontsize=7.5)

    plt.savefig('bbfit_plots/'+dir+'_seds.png')

if __name__ == "__main__":
    dir = sys.argv[1]
    plot_seds(dir)