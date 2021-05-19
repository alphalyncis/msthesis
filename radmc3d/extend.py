import numpy as np
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit
import shutil

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


class Radmc3dModel():

    def __init__(self):
        self.grid = []             # 1d grid info
        self.nrlayers = 0
        self.nx = []               # cell numbers for each layer
        self.ny = []
        self.nynew = []
        self.nz = []
        self.ix = []               # cell vertices for base layer
        self.iy = []
        self.iynew = []
        self.iz = []
        self.amrinfo = []          # array containing layer sizes
        self.nrcellsnew = 0
        self.gas_dens_hdr = []     # file header
        self.gas_dens_data = []    # 1d original data 
        self.gas_dens = []         # reshaped array
        self.gas_dens_extended = []# 1d extended data
        self.gas_temp_hdr = []     # file header
        self.gas_temp_data = []    # 1d original data
        self.gas_temp_extended = []# 1d extended data

        ''' Example:
            model_1 = Radmc3dModel()
            model_1.read_grid()

        '''

    def read_grid(self, fname='amr_grid.inp'):

        with open(fname,'r') as gfile:
            self.grid = np.fromfile(gfile, count=-1, sep=' ')
            # grid[0] = iformat
            # grid[1] = grid style
            # grid[2] = coordsystem
            # grid[3] = gridinfo
            # grid[4:7] = include x,y,z
            # grid[7:10] = nx,ny,nz
            # grid[10] = nrlevels
            # grid[11] = nrlayers
            # grid[-7*nrlayers:] = amr layers info

        self.nrlayers = int(self.grid[11]) + 1 # including 0th layer

        self.nx = np.zeros(self.nrlayers, dtype=int)
        self.ny = np.zeros(self.nrlayers, dtype=int)
        self.nz = np.zeros(self.nrlayers, dtype=int)
        self.nx[0], self.ny[0], self.nz[0] = self.grid[7:10]

        self.amrinfo = np.reshape([int(i) for i in self.grid[-7*(self.nrlayers-1):]], (self.nrlayers-1,7))
        for n in range(self.nrlayers-1):
            self.nx[n+1], self.ny[n+1], self.nz[n+1] = 2 * self.amrinfo[n,4], 2 * self.amrinfo[n,5], 2 * self.amrinfo[n,6]

        self.ix, self.iy, self.iz = self.grid[12:13+self.nx[0]], self.grid[13+self.nx[0]: 14+self.nx[0]+self.ny[0]], self.grid[14+self.nx[0]+self.ny[0]:15+self.nx[0]+self.ny[0]+self.nz[0]]  # cell vertices

    def extend_grid(self):
        '''Write new amr_grid_suffix.inp file with suffix'''

        ystep = self.iy[1] - self.iy[0]
        ystart = self.iy[0] - 30 * ystep
        self.iynew = np.linspace(ystart, ystart+100*ystep, 101)
        self.nynew = np.concatenate((np.array([100]), self.ny[1:]))
        self.amrinfo[0][2] += 30

        with open('amr_grid_'+suffix+'.inp','w') as wfile:
            wfile.write('%d\n' % self.grid[0])  # iformat
            wfile.write('%d\n' % self.grid[1])  # grid style
            wfile.write('%d\n' % self.grid[2])  # coordsystem
            wfile.write('%d\n' % self.grid[3])  # grid info
            wfile.write('%d %d %d \n' % (self.grid[4], self.grid[5], self.grid[6]))
            wfile.write('%d %d %d \n' % (self.nx[0], self.nynew[0], self.nz[0]))
            wfile.write('%d %d \n' % (self.grid[10], self.grid[11]))
            for i in range(self.nx[0]+1):
                wfile.write('%.5e  ' % self.ix[i])
            wfile.write('\n')
            for i in range(self.nynew[0]+1):
                if i == 50:
                    wfile.write('1.57079632679489661923132169164  ')
                else:
                    wfile.write('%.5f  ' % self.iynew[i])
            wfile.write('\n')
            for i in range(self.nz[0]+1):
                wfile.write('%.5f  ' % self.iz[i])
            wfile.write('\n')
            for row in self.amrinfo:
                line = ' '.join([str(elem) for elem in row])
                wfile.write(line+'\n')
        wfile.close()

    def read_gas_dens(self, fname='gas_density.binp'):

        with open(fname,'r') as rfile:
            self.gas_dens_hdr = np.fromfile(rfile, count=4, dtype=np.int64)
            # hdr[0] = iformat
            # hdr[1] = precision
            # hdr[2] = nrcells
            # hdr[3] = nrspecs
            self.gas_dens_data = np.fromfile(rfile, count=self.gas_dens_hdr[2], dtype=np.float64)

        dens = []
        start = 0
        for n in range(self.nrlayers):
            if n == 0:
                layer = self.gas_dens_data[:self.nx[n] * self.ny[n] * self.nz[n]]
            else:
                start += len(layer)
                layer = self.gas_dens_data[start:start + self.nx[n] * self.ny[n] * self.nz[n]]
            #print(len(layer))
            layer3d = np.reshape(layer, (self.nx[n], self.ny[n], self.nz[n]), order="F")
            dens.append(layer3d)

        self.gas_dens = dens

    def extend_gas_dens(self):
        '''Write new extended 'gas_density_suffix.binp' file with suffix'''

        print('Extending gas density by Gaussian extrapolation...')

        # ------- Gaussian fit -------

        def g(x, amp, mean, stddev):
            z = (x - mean) / stddev
            y = amp * np.exp(-z**2 / 2)
            return y

        yold = self.iy[:40]
        ynew = self.iynew

        dens0new = np.zeros((214, 100, 680))

        # fixing the boundaries to be continous
        sigma =np.ones(len(yold)) 
        sigma[0] = 0.01
        sigma[-1] = 0.01  

        # do fitting and extrapolation 
        for x in range(214):
            for z in range(680):
                d = self.gas_dens[0][x,:,z]
                p_init = [d.max(), 1.57, 0.01]
                coeff, cov = curve_fit(g, yold, d, p0=p_init, sigma=sigma, maxfev=3000)
                dnew = g(ynew, coeff[0], coeff[1], coeff[2])

                for y in range(100):
                    if ynew[y] > yold.max() or ynew[y] < yold.min():
                        dens0new[x,y,z] = dnew[y]
                    else:
                        dens0new[x,y,z] = d[y-30]

        ''' Plot density slices:
        cmap = plt.cm.viridis
        cmap.set_bad(color=cmap(0))
        plt.imshow(dens0new[:,:,341].T,cmap=cmap, norm=LogNorm(vmin=5e-20, vmax=5e-10))
        plt.colorbar(orientation="horizontal")
        '''

        # ------- Change gas_density.binp -------

        dens0new1d = dens0new.flatten("F")
        self.gas_dens_extended = np.concatenate((dens0new1d, self.gas_dens_data[self.nx[0]*self.ny[0]*self.nz[0]:]))

        nrcellsnew = 0
        for n in range(self.nrlayers):
            if n == 0:
                nrcellsnew += self.nx[n] * self.nynew[n] * self.nz[n]
            else:
                nrcellsnew += self.nx[n] * self.ny[n] * self.nz[n]

        self.nrcellsnew = nrcellsnew

        with open('gas_density_'+suffix+'.binp','w') as wfile:
            hdrnew = np.array([self.gas_dens_hdr[0], self.gas_dens_hdr[1], self.nrcellsnew, self.gas_dens_hdr[3]], dtype=np.int64)
            hdrnew.tofile(wfile)
            self.gas_dens_extended.tofile(wfile)

    def read_gas_temp(self, fname='gas_temperature.binp'):
        with open(fname,'r') as rfile:
            self.gas_temp_hdr = np.fromfile(rfile, count=4, dtype=np.int64)
            # hdr[0] = iformat
            # hdr[1] = precision
            # hdr[2] = nrcells
            # hdr[3] = nrspecs
            self.gas_temp_data = np.fromfile(rfile, count=self.gas_temp_hdr[2], dtype=np.float64)

    def extend_gas_temp(self):
        '''Write new extended 'gas_temperature_suffix.binp' file with suffix'''

        print('Extending gas temperature...')

        temp0 = self.gas_temp_data[:self.nx[0] * self.ny[0] * self.nz[0]]
        temp03d = np.reshape(temp0, (self.nx[0], self.ny[0], self.nz[0]), order="F")

        tempnew = np.zeros((self.nx[0], self.nynew[0], self.nz[0]))
        for x in range(self.nx[0]):
            for z in range(self.nz[0]):
                temp_top = temp03d[x,0,z]
                temp_bot = temp03d[x,39,z]
                for y in range(0,30):
                    tempnew[x,y,z] = temp_top
                for y in range(30,70):
                    tempnew[x,y,z] = temp03d[x,y-30,z]
                for y in range(70,100):
                    tempnew[x,y,z] = temp_bot


        # write temperature files
        tempnew1d = tempnew.flatten("F")
        self.gas_temp_extended = np.concatenate((tempnew1d, self.gas_dens_data[self.nx[0]*self.ny[0]*self.nz[0]:]))

        with open('gas_temperature_'+suffix+'.binp','w') as wfile:
            hdr = np.array([self.gas_temp_hdr[0], self.gas_temp_hdr[1], self.nrcellsnew, self.gas_temp_hdr[3]], dtype=np.int64)
            hdr.tofile(wfile)
            self. gas_temp_extended.tofile(wfile)
        wfile.close()

    def write_dust_temp(self):
        '''Write new extended 'dust_temperature_suffix.binp' file with suffix'''
        shutil.copy('gas_temperature_'+suffix+'.binp', 'dust_temperature_'+suffix+'.bdat')

    def write_dust_dens(self):
        with open('gas_density_'+suffix+'.binp','r') as rfile:
            hdr = np.fromfile(rfile, count=4, dtype=np.int64)
            # hdr[0] = iformat
            # hdr[1] = precision
            # hdr[2] = nrcells
            # hdr[3] = nrspecs
            gas_dens = np.fromfile(rfile, count=hdr[2], dtype=np.float64)

        dust_dens = 0.01 * gas_dens

        with open('dust_density_'+suffix+'.binp','w') as wfile:
            hdr.tofile(wfile)
            dust_dens.tofile(wfile)
        wfile.close()

    def include_evap(self):
        '''Write new extended 'gas_temperature_suffix.binp' file with suffix'''
        print('Including dust evaporation...')

        with open('dust_density_'+suffix+'.binp','r') as rfile:
            hdr = np.fromfile(rfile, count=4, dtype=np.int64)
            # hdr[0] = iformat
            # hdr[1] = precision
            # hdr[2] = nrcells
            # hdr[3] = nrspecs
            dust_dens = np.fromfile(rfile, count=hdr[2], dtype=np.float64)

        # dust evaporation > 1500 K
        change = 0
        for i in range(len(dust_dens)):
            if self.gas_temp_extended[i] > 1500:
                dust_dens[i] = 0
                change += 1

        with open('dust_density_'+suffix+'_evap.binp','w') as wfile:
            hdr.tofile(wfile)
            dust_dens.tofile(wfile)
        wfile.close()
        print(change, 'cells dust density set to 0.')

    def extend_model(self):
        self.extend_grid()
        self.extend_gas_dens()
        self.extend_gas_temp()
        self.write_dust_dens()
        self.write_dust_temp()
        print('Extension done.')

def read_in_model():
    print('Reading in model...')
    model = Radmc3dModel()
    model.read_grid()
    model.read_gas_dens()
    model.read_gas_temp()
    return model

# Example:
suffix = 'test'
testmodel = read_in_model()
testmodel.extend_model()
testmodel.include_evap()