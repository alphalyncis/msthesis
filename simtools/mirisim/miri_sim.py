# import the configuration file parsers so they can be written to file
from mirisim.config_parser import SimConfig, SimulatorConfig, SceneConfig

# import scene component generators
from mirisim.skysim import Background
from mirisim.skysim import externalsources

from mirisim import MiriSimulation

#other things to be used
import numpy as np
import glob                 # glob is used to find the output directory
import os                   # for listing directory contents
import sys
import shutil
from astropy.io import fits # for reading FITS file contents

#jwst pipeline imports
from jwst.pipeline import Detector1Pipeline
from jwst.pipeline import Image2Pipeline


def sim_wrapper(fname, pipl):
    '''
    e.g.
    fname = '10jup50au/cube5npix23lam1020inc0.fits'
    lams = ['1020','1500','2100']
    '''
    lam = fname[fname.find('lam')+3:fname.find('inc')]
    if lam == '1020':
        filt = '1000'
    else:
        filt = lam
    inc = fname[fname.find('inc')+3:fname.find('.fits')]
    dir = fname.split('/cube')[0]
    datadir = '/home/ipa/ipaszulagyi/users/xuechen/thesis/RADMC3D/'
    
    ########## Create a scene ##########
    bg = Background(level = 'low', gradient = 5., pa = 45.)

    skycube = externalsources.Skycube(datadir+fname)
    scene_config = SceneConfig.makeScene(name="sky0", targets=[skycube], background=bg, loglevel=1)


    ########## Config simulation parameters ##########
    filters = {'1020':'F1000W', '1500':'F1500W', '2100':'F2100W'}
    filter = filters[lam]

    sim_config = SimConfig.makeSim(
    name = 'CPD_sim',    # name given to simulation
    scene = 'scene.ini', # name of scene file to input
    rel_obsdate = 0.0,          # relative observation date (0 = launch, 1 = end of 5 yrs)
    POP = 'IMA',                # Component on which to center (Imager or MRS)
    ConfigPath = 'IMA_SUB256',  # Configure the Optical path (MRS sub-band)
    Dither = False,             # Don't Dither
    StartInd = 1,               # start index for dither pattern [NOT USED HERE]
    NDither = 2,                # number of dither positions [NOT USED HERE]
    DitherPat = 'ima_recommended_dither.dat', # dither pattern to use [NOT USED HERE]
    disperser = 'SHORT',        # [NOT USED HERE]
    detector = 'SW',            # [NOT USED HERE]
    mrs_mode = 'SLOW',          # [NOT USED HERE]
    mrs_exposures = 2,          # [NOT USED HERE]
    mrs_integrations = 3,       # [NOT USED HERE]
    mrs_frames = 5,             # [NOT USED HERE]
    ima_exposures = 1,          # number of exposures
    ima_integrations = 10,      # number of integrations
    ima_frames = 50,            # number of groups (for MIRI, # Groups = # Frames)
    ima_mode = 'FAST',          # Imager read mode (default is FAST ~ 2.3 s)
    filter = filter,          # Imager Filter to use
    readDetect = 'SUB256'         # Portion of detector to read out
    )


    ########## Run mirage simulation ##########

    simulator_config = SimulatorConfig.from_default()
    
    mysim = MiriSimulation(sim_config, scene_config, simulator_config)
    mysim.run()

    outputdir = sorted(glob.glob('*_*_mirisim'),key=os.path.getmtime)[-1]
    infits = glob.glob('{}/{}/*.fits'.format(outputdir, 'det_images'))[0]

    ########## Data reduction by jwst pipeline ##########
    if pipl:
        result1 = Detector1Pipeline.call(infits, steps={"dark_current":{"skip": True}}, save_results=True, output_dir=outputdir)
        result2 = Image2Pipeline.call(result1, save_results=True, output_dir=outputdir)
        infits = outputdir+'/det_image_seq1_MIRIMAGE_F'+filt+'Wexp1_cal.fits'


    ########## Crop final image ##########

    print("Cropping final image...")
    r = 30
    cen = int(256/2)
    xmin, xmax, ymin, ymax = cen-r, cen+r, cen-r, cen+r

    if not pipl:
        with fits.open(infits) as hdul:
            lin_data = hdul['SCI'].data
            integ,frames,nx,ny = hdul['SCI'].data.shape

            print('rewriting'+infits+'...')
            cropped = np.flip(lin_data[:,frames-1, xmin:xmax, ymin:ymax], axis=2)
            hdul['SCI'].data = cropped
            hdul.writeto('./image{}_sim/{}/{}_lam{}inc{}.fits'.format(lam, dir, filter, lam, inc), overwrite=True)
        
        newdir = './image{}_sim/{}20int_{}_lam{}inc{}'.format(lam, filter, dir, lam, inc)
        if os.path.exists(newdir):
            shutil.rmtree(newdir, ignore_errors=True)
        os.rename(outputdir, './image{}_sim/{}20int_{}_lam{}inc{}'.format(lam, filter, dir, lam, inc))

    else:
        with fits.open(infits) as hdul:
            lin_data = hdul['SCI'].data

            print('rewriting'+infits+'...')
            cropped = np.flip(lin_data[xmin:xmax, ymin:ymax], axis=1)
            hdul['SCI'].data = cropped
            hdul.writeto('./image{}_sim/{}/{}_lam{}inc{}_cal.fits'.format(lam, dir, filter, lam, inc), overwrite=True)

        os.rename(outputdir, 'outputdir')

    print("ALL DONE.")



if __name__ == "__main__":
    fname = sys.argv[1]
    try:
        pipl = sys.argv[2]
        if pipl == 'pipl':
            pipl = True
    except:
        pipl = False
    sim_wrapper(fname, pipl)
