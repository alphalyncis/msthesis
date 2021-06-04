import numpy as np
from astropy.io import fits
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm      
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

from tofits import readImage

fname = 'image200.out'

# ------- Plot image.out with matplolib -------

im = readImage(fname)
cmap = plt.cm.viridis
cmap.set_bad(color=cmap(0))
img = im.imageJyppix.reshape((im.nx, im.ny)).swapaxes(1, 0)
plt.imshow(img, cmap=cmap, norm=LogNorm(vmin=1e-3,vmax=1e3), origin='lower')  
cbar = plt.colorbar()
cbar.set_label('Jy/pixel')

'''
# ------- Plot fits image with matplolib -------

with fits.open(fname) as inhdul:
    sourceimg = inhdul[0].data[0]
plt.imshow(sourceimg, cmap=cmap, norm=LogNorm(vmin=1e-10, vmax=1e2), origin='lower')
cbar = plt.colorbar()
cbar.set_label('Jy/pixel')
'''