# Plotting hydro fields from radmc3d models
import numpy as np
from astropy.io import fits
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm      
import matplotlib.cbook
from matplotlib import patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.cbook
import warnings
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
cmap = plt.cm.viridis
cmap.set_bad(color=cmap(0))


from radmcutils import *

def incl_amr(model, field, vmin=1e-22, vmax=1e-12, axis=True):
    # inset axes....
    for i in range(model.nrlayers-1):
        if i == 0:
            xmin = (model.amrinfo[0][1]-1)/nx 
            ymin = (model.amrinfo[0][2]-1)/ny
            xlen = model.amrinfo[0][4]/nx 
            ylen = model.amrinfo[0][5]/ny 
        else: 
            t1 = xlen * (model.amrinfo[i][1]-1)/(2*model.amrinfo[i-1][4])
            t2 = ylen * (model.amrinfo[i][2]-1)/(2*model.amrinfo[i-1][5])
            xmin += t1
            ymin += t2
            xlen *= model.amrinfo[i][4]/(model.amrinfo[i-1][4]*2)
            ylen *= model.amrinfo[i][5]/(model.amrinfo[i-1][5]*2)
        bound = [xmin, ymin, xlen, ylen]

        axins = ax.inset_axes(bound)
        extent = (0, model.amrinfo[i][4]*2, 0, model.amrinfo[i][5]*2)
        z = int(model.nz[i+1]/2)
        axins.imshow(field[i+1][:,:,z].T, extent=extent, origin="lower", 
                    norm=LogNorm(vmin=vmin, vmax=vmax))
        if axis == False:
            axins.axis('off')
        else:
            axins.xaxis.set_ticklabels([])
            axins.xaxis.set_ticks_position('none')
            axins.yaxis.set_ticklabels([])
            axins.yaxis.set_ticks_position('none')

def plot_slice(model, field, viewaxis='z', vmin=1e-20, vmax=1e-10, axis=True):
    if viewaxis == 'z':
        fig, ax = plt.subplots(figsize=[nx/10, ny/10], tight_layout={'pad':0})
        y = int(model.nz[0]/2)
        ax.imshow(field[0][:,:,z].T, origin="lower", norm=LogNorm(vmin=vmin, vmax=vmax))
    elif viewaxis == 'y':
        fig, ax = plt.subplots(figsize=[nx/10, nz/10], tight_layout={'pad':0})
        y = int(model.ny[0]/2)
        ax.imshow(field[0][:,y,:].T, origin="lower", norm=LogNorm(vmin=vmin, vmax=vmax))

    for i in range(model.nrlayers-1):
        if i == 0:
            xmin = (model.amrinfo[0][1]-1)/nx 
            ymin = (model.amrinfo[0][2]-1)/ny
            zmin = (model.amrinfo[0][3]-1)/nz
            xlen = model.amrinfo[0][4]/nx 
            ylen = model.amrinfo[0][5]/ny
            zlen = model.amrinfo[0][6]/nz
        else: 
            t1 = xlen * (model.amrinfo[i][1]-1)/(2*model.amrinfo[i-1][4])
            t2 = ylen * (model.amrinfo[i][2]-1)/(2*model.amrinfo[i-1][5])
            t3 = zlen * (model.amrinfo[i][3]-1)/(2*model.amrinfo[i-1][6])
            xmin += t1
            ymin += t2
            zmin += t3
            xlen *= model.amrinfo[i][4]/(model.amrinfo[i-1][4]*2)
            ylen *= model.amrinfo[i][5]/(model.amrinfo[i-1][5]*2)
            zlen *= model.amrinfo[i][6]/(model.amrinfo[i-1][6]*2)

        if viewaxis == 'z':
            index = int(model.nz[i+1]/2)
            bound = [xmin, ymin, xlen, ylen]
            extent = (0, model.amrinfo[i][4]*2, 0, model.amrinfo[i][5]*2)
            axins = ax.inset_axes(bound)
            axins.imshow(field[i+1][:,:,index].T, extent=extent, origin="lower", 
                    norm=LogNorm(vmin=vmin, vmax=vmax))

        elif viewaxis == 'y':
            index = int(model.ny[i+1]/2)
            bound = [xmin, zmin, xlen, zlen]
            extent = (0, model.amrinfo[i][4]*2, 0, model.amrinfo[i][6]*2)
            axins = ax.inset_axes(bound)
            axins.imshow(field[i+1][:,index,:].T, extent=extent, origin="lower", 
                        norm=LogNorm(vmin=vmin, vmax=vmax))

        if axis == False:
            axins.axis('off')
        else:
            axins.xaxis.set_ticklabels([])
            axins.xaxis.set_ticks_position('none')
            axins.yaxis.set_ticklabels([])
            axins.yaxis.set_ticks_position('none')

for dir in ['10jup5au/', '5jup5au/', '1jup5au/', '1sat5au/', '1jup30au/']:
    model = Radmc3dModel()
    model.read_grid(dir=dir)
    model.read_file('dust_density.binp')
    model.read_file('dust_temperature.bdat')
    print(model.amrinfo)
    nx, ny, nz = model.nx[0], model.ny[0], model.nz[0]

    # plot 1: separate dens layers, to see if evap
    plt.figure(figsize=(12,8))
    for i in range(model.nrlayers):
        plt.subplot(2,3,i+1)
        zindex = int(model.nz[i]/2)
        plt.imshow(model.dust_dens[i][:,:,zindex].T,cmap=cmap, norm=LogNorm(vmin=1e-25, vmax=1e-10))
        plt.title("z="+str(zindex))
    plt.colorbar()
    plt.savefig('./hydrofield_plots/'+dir[:-1]+'_dens_sep.png')
    plt.close()

    # plot 2: dens z
    fig, ax = plt.subplots(figsize=[nx/15, ny/15], tight_layout={'pad':0})
    z = int(model.nz[0]/2)
    im = ax.imshow(model.dust_dens[0][:,:,z].T, origin="lower", norm=LogNorm(vmin=1e-22, vmax=1e-12))
    #fig.colorbar(im, ax=ax)
    incl_amr(model, field=model.dust_dens, axis=True)
    plt.savefig('./hydrofield_plots/'+dir[:-1]+'_dens_z.png')
    plt.close()

    # plot 3: temp z
    fig, ax = plt.subplots(figsize=[nx/15, ny/15], tight_layout={'pad':0})
    im = ax.imshow(model.dust_temp[0][:,:,341].T, origin="lower", norm=LogNorm(vmin=10, vmax=1e3))
    #fig.colorbar(im, ax=ax)
    incl_amr(model, model.dust_temp, vmin=10, vmax=1e3, axis=True)
    plt.savefig('./hydrofield_plots/'+dir[:-1]+'_temp_z.png')
    plt.close()

    # plot 4: dens y
    plot_slice(model, field=model.dust_dens, viewaxis='y')
    plt.savefig('./hydrofield_plots/'+dir[:-1]+'_dens_y.png')
    plt.close()

# plot 5: dens colorbar
fig, ax = plt.subplots(figsize=[nx/15, ny/15], tight_layout={'pad':0})
z = int(model.nz[0]/2)
im = ax.imshow(model.dust_dens[0][:,:,z].T, origin="lower", norm=LogNorm(vmin=1e-22, vmax=1e-12))
fig.colorbar(im, ax=ax)
incl_amr(model, field=model.dust_dens, axis=True)
plt.savefig('./hydrofield_plots/dens_colorbar.png')
plt.close()

# plot 3: temp colorbar
fig, ax = plt.subplots(figsize=[nx/15, ny/15], tight_layout={'pad':0})
im = ax.imshow(model.dust_temp[0][:,:,341].T, origin="lower", norm=LogNorm(vmin=10, vmax=1e3))
fig.colorbar(im, ax=ax)
incl_amr(model, model.dust_temp, vmin=10, vmax=1e3, axis=True)
plt.savefig('./hydrofield_plots/temp_colorbar.png')
plt.close()
