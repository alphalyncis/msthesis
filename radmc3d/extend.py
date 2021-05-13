import numpy as np
from astropy.modeling import models, fitting
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm



# ------- Read grids ---------

with open('amr_grid.inp','r') as gfile:
    grid = np.fromfile(gfile, count=-1, sep=' ')
    # grid[0] = iformat
    # grid[1] = grid style
    # grid[2] = coordsystem
    # grid[3] = gridinfo
    # grid[4:7] = include x,y,z
    # grid[7:10] = nx,ny,nz
    # grid[10] = nrlevels
    # grid[11] = nrlayers
    # grid[-7*nrlayers:] = amr layers info

nrlayers = int(grid[11]) + 1 # including 0th layer

nx = np.zeros(nrlayers, dtype=int)
ny = np.zeros(nrlayers, dtype=int)
nz = np.zeros(nrlayers, dtype=int)
nx[0], ny[0], nz[0] = grid[7:10]

amrinfo = np.reshape([int(i) for i in grid[-7*(nrlayers-1):]], (nrlayers-1,7))
for n in range(nrlayers-1):
    nx[n+1], ny[n+1], nz[n+1] = 2 * amrinfo[n,4], 2 * amrinfo[n,5], 2 * amrinfo[n,6]

ix, iy, iz = grid[12:13+nx[0]], grid[13+nx[0]: 14+nx[0]+ny[0]], grid[14+nx[0]+ny[0]:15+nx[0]+ny[0]+nz[0]]  # cell vertices


# ------- Change amr_grid.inp -------

ystep = iy[1] - iy[0]
ystart = iy[0] - 30 * ystep
iynew = np.linspace(ystart, ystart+100*ystep, 101)
nynew = 100
amrinfo[0][2] += 30

with open('amr_grid_ext.inp','w') as wfile:
    wfile.write('%d\n' % grid[0])  # iformat
    wfile.write('%d\n' % grid[1])  # grid style
    wfile.write('%d\n' % grid[2])  # coordsystem
    wfile.write('%d\n' % grid[3])  # grid info
    wfile.write('%d %d %d \n' % (grid[4], grid[5], grid[6]))
    wfile.write('%d %d %d \n' % (nx[0], nynew, nz[0]))
    wfile.write('%d %d \n' % (grid[10], grid[11]))
    for i in range(nx[0]+1):
        wfile.write('%.5e  ' % ix[i])
    wfile.write('\n')
    for i in range(ny[0]+1):
        if i == 50:
            wfile.write('1.57079632679489661923132169164  ')
        else:
            wfile.write('%.5f  ' % iynew[i])
    wfile.write('\n')
    for i in range(nz[0]+1):
        wfile.write('%.5f  ' % iz[i])
    wfile.write('\n')
    for row in amrinfo:
        line = ' '.join([str(elem) for elem in row])
        wfile.write(line+'\n')
wfile.close()


# -------- Read densities -------

with open('gas_density.binp','r') as rfile:
    hdr = np.fromfile(rfile, count=4, dtype=np.int64)
    # hdr[0] = iformat
    # hdr[1] = precision
    # hdr[2] = nrcells
    # hdr[3] = nrspecs
    data = np.fromfile(rfile, count=hdr[2], dtype=np.float64)

dens = []
start = 0
for n in range(nrlayers):
    if n == 0:
        layer = data[:nx[n] * ny[n] * nz[n]]
    else:
        start += len(layer)
        layer = data[start:start + nx[n] * ny[n] * nz[n]]
    #print(len(layer))
    layer3d = np.reshape(layer, (nx[n], ny[n], nz[n]), order="F")
    dens.append(layer3d)


# ------- Gaussian fit -------

def g(x, amp, mean, stddev):
    z = (x - mean) / stddev
    y = amp * np.exp(-z**2 / 2)
    return y

yold = iy[:40]
ynew = iynew

dens0new = np.zeros((214, 100, 680))

# fixing the boundaries to be continous
sigma =np.ones(len(yold)) 
sigma[0] = 0.01
sigma[-1] = 0.01  

# do fitting and extrapolation 
for x in range(214):
    for z in range(680):
        d = dens[0][x,:,z]
        p_init = [d.max(), 1.57, 0.01]
        coeff, cov = curve_fit(g, yold, d, p0=p_init, sigma=sigma, maxfev=3000)
        dnew = g(ynew, coeff[0], coeff[1], coeff[2])

        for y in range(100):
            if ynew[y] > yold.max() or ynew[y] < yold.min():
                dens0new[x,y,z] = dnew[y]
            else:
                dens0new[x,y,z] = d[y-30]

# plot density slices
cmap = plt.cm.viridis
cmap.set_bad(color=cmap(0))
plt.imshow(dens0new[:,:,341].T,cmap=cmap, norm=LogNorm(vmin=5e-20, vmax=5e-10))
plt.colorbar(orientation="horizontal")



# ------- Change gas_density.binp -------

dens0new1d = dens0new.flatten("F")
datanew = np.concatenate((dens0new1d, data[nx[0]*ny[0]*nz[0]:]))

nynew = 100
nrcellsnew = 0
for n in range(nrlayers):
    if n == 0:
        nrcellsnew += nx[n] * nynew * nz[n]
    else:
        nrcellsnew += nx[n] * ny[n] * nz[n]

with open('gas_density_ext.binp','w') as wfile:
    hdrnew = np.array([hdr[0], hdr[1], nrcellsnew, hdr[3]], dtype=np.int64)
    hdrnew.tofile(wfile)
    datanew.tofile(wfile)



# ------- Change gas_temperature.binp -------

with open('gas_temperature.binp','r') as rfile:
    hdr = np.fromfile(rfile, count=4, dtype=np.int64)
    # hdr[0] = iformat
    # hdr[1] = precision
    # hdr[2] = nrcells
    # hdr[3] = nrspecs
    gas_temp = np.fromfile(rfile, count=hdr[2], dtype=np.float64)

temp0 = gas_temp[:nx[0] * ny[0] * nz[0]]
temp03d = np.reshape(temp0, (nx[0], ny[0], nz[0]), order="F")

tempnew = np.zeros((nx[0], nynew, nz[0]))
for x in range(nx[0]):
    for z in range(nz[0]):
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
gas_temp_new = np.concatenate((tempnew1d, gas_temp[nx[0]*ny[0]*nz[0]:]))
nrcellsnew = 0
for n in range(nrlayers):
    if n == 0:
        nrcellsnew += nx[n] * nynew * nz[n]
    else:
        nrcellsnew += nx[n] * ny[n] * nz[n]

with open('gas_temperature_ext.binp','w') as wfile:
    hdr = np.array([hdr[0], hdr[1], nrcellsnew, hdr[3]], dtype=np.int64)
    hdr.tofile(wfile)
    gas_temp_new.tofile(wfile)
wfile.close()



# ------- Change dust_temperature.bdat -------

with open('dust_temperature_ext.bdat','w') as wfile:
    hdr = np.array([hdr[0], hdr[1], nrcellsnew, hdr[3]], dtype=np.int64)
    hdr.tofile(wfile)
    gas_temp_new.tofile(wfile)
wfile.close()



# -------- Change dust_density.binp -------

with open('gas_density_ext.binp','r') as rfile:
    hdr = np.fromfile(rfile, count=4, dtype=np.int64)
    # hdr[0] = iformat
    # hdr[1] = precision
    # hdr[2] = nrcells
    # hdr[3] = nrspecs
    gas_dens = np.fromfile(rfile, count=hdr[2], dtype=np.float64)

dust_dens = 0.01 * gas_dens

with open('dust_density_ext.binp','w') as wfile:
    hdr.tofile(wfile)
    dust_dens.tofile(wfile)
wfile.close()
