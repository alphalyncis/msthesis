# FOR NIRISS: mirage imaging simlation wrapper
import os
import sys

# For processing outputs
from astropy.io import fits

# mirage imports
from mirage import imaging_simulator
from mirage.seed_image import catalog_seed_image
from mirage.yaml import yaml_generator

#jwst pipeline imports
from jwst.pipeline import Detector1Pipeline, Image2Pipeline, Ami3Pipeline


def sim_wrapper(fname, pipl, ref):
    lam = fname[fname.find('lam')+3:fname.find('inc')]

    # Set environment variables
    os.environ["MIRAGE_DATA"] = "/home/ipa/ipaszulagyi/users/xuechen/thesis/simtools/mirage/mirage/mirage_data"
    os.environ["CRDS_DATA"] = "/home/ipa/ipaszulagyi/users/xuechen/thesis/simtools/mirage/mirage/mirage_data/crds_cache"
    os.environ["CRDS_SERVER_URL"] = "https://jwst-crds.stsci.edu"
    
    ########## Rewrite source_catalog_file.cat ##########
    print('Rewriting source cat file...')

    with open(f'./image{lam}_data/source_catalog_file.cat', 'r+') as catfile:
        lines = catfile.readlines()
        lines[-1] = f'1  47.52000  -10.09000  None  None  ./image{lam}_data/{fname}'
        catfile.seek(0)
        catfile.truncate()
        catfile.writelines(lines)



    ########## Config simulation ##########

    # Specify the xml and pointing files exported from APT
    xml_file = f'./image{lam}_data/nrsami_grp50int10.xml'
    pointing_file = f'./image{lam}_data/nrsami.pointing'

    # Source catalogs to be used.
    cat_dict = {'extended':f'image{lam}_data/source_catalog_file.cat'}

    # Set reference file values. 
    # Setting to 'crds_full_name' will search for and download needed calibration reference files (commonly referred to as CRDS reference files) when the yaml_generator is run. 
    # Setting to 'crds' will put placeholders in the yaml files and save the downloading for when the simulated images are created.
    reffile_defaults = 'crds'

    # Optionally set the cosmic ray library and rate
    cosmic_rays = {'library': 'SUNMIN', 'scale': 0}

    # Optionally set the background signal rates to be used
    background = 'medium'

    # Optionally set the telescope roll angle (PAV3) for the observations
    pav3 = 0

    # Optionally set the observation date to use for the data. Note that this information
    # is placed in the headers of the output files, but not used by Mirage in any way.
    dates = '2022-10-31'

    ghosts = False
    convolve_ghosts = False

    # Set the directory into which the yaml files and simulation output will be written
    yaml_dir = f'./image{lam}_data/'
    simulation_dir = f'./image{lam}_sim/'

    # specify the data reduction state of the Mirage outputs
    datatype = 'linear'


    ########## Generate yaml file ##########

    # For simu for the 1st time, run the yaml generator 
    yam = yaml_generator.SimInput(input_xml=xml_file, pointing_file=pointing_file,
                                catalogs=cat_dict, cosmic_rays=cosmic_rays,
                                background=background, roll_angle=pav3,
                                dates=dates, reffile_defaults=reffile_defaults,
                                add_ghosts=ghosts, convolve_ghosts_with_psf=convolve_ghosts,
                                verbose=True, output_dir=yaml_dir,
                                simdata_output_dir=simulation_dir,
                                datatype=datatype)
    yam.create_inputs()


    ########## Run mirage simulation ##########

    # set file dir
    yfile = f'./image{lam}_data/jw00005001001_01101_00005_nis.yaml'
    linfile = f'./image{lam}_sim/jw00005001001_01101_00005_nis_linear.fits'

    # for finding central position later
    cat = catalog_seed_image.Catalog_seed()
    cat.paramfile = yfile
    cat.make_seed()

    # run simulation
    img_sim = imaging_simulator.ImgSim()
    img_sim.paramfile = yfile
    img_sim.create()


    ########## Data reduction by jwst pipeline ##########
    if pipl:
        result1 = Detector1Pipeline.call(linfile, steps={"dark_current":{"skip": True}, "persistence":{"skip":True}}, save_results=True, output_dir=f'./image{lam}_sim/')
        result2 = Image2Pipeline.call(result1, steps={"photom":{"skip": False}}, save_results=True, output_dir=f'./image{lam}_sim/')
        if not ref:
            os.chdir(f'image{lam}_sim')
            result3 = Ami3Pipeline.call('niriss_asn.json')
            os.chdir('..')


    ########## Crop final image ##########
    print("Cropping final image...")
    if '50au' in fname:
        r = 20
    else:
        r = 12
    ycen, xcen = int(get_center_pos(cat)[0]+0.5), int(get_center_pos(cat)[1]+0.5)
    ymin, ymax, xmin, xmax = ycen-r, ycen+r, xcen-r, xcen+r
    with fits.open(f'./image{lam}_sim/jw00005001001_01101_00005_nis_linear_cal.fits', 'update') as hdul:
        lin_data = hdul['SCI'].data
        print('rewriting'+linfile+'...')
        cropped = lin_data[ymin:ymax, xmin:xmax]
        hdul['SCI'].data = cropped
        hdul.flush()

    # if pipl:
    #     if ref:
    #         pass
    #     else:    
    #         with fits.open(f'./image{lam}_sim/combined_aminorm.fits', 'update') as hdul:
    #             lin_data = hdul['FIT'].data
    #             print('rewriting'+linfile+'...')
    #             cropped = lin_data[ymin:ymax, xmin:xmax]
    #             hdul['FIT'].data = cropped
    #             hdul.flush()
    # print("ALL DONE.")


def get_center_pos(cat):
    filename = cat.params['simSignals']['extended']
    lines, pixelflag, magsys = cat.read_point_source_file(filename)
    indexes = lines['index']
    indexes, lines = cat.remove_outside_fov_sources(indexes, lines, pixelflag, 4096)
    for indexnum, values in zip(indexes, lines):
        pixelx, pixely, ra, dec, ra_str, dec_str = cat.get_positions(values['x_or_RA'], values['y_or_Dec'], pixelflag, 4096)
    return pixely, pixelx


if __name__ == "__main__":
    fname = sys.argv[1]
    try:
        option = sys.argv[2]
        if option == 'pipl':
            pipl = True
            ref = False
        elif option == 'piplref':
            pipl = True
            ref = True
    except:
        pipl = False
        ref = False
    sim_wrapper(fname, pipl, ref)
