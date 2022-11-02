import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

import numpy as np
import scipy.interpolate
import os, sys, glob
import gc
import argparse
import pandas
import time
import yt
yt.enable_parallelism()
import logging
import pickle

import athena_input as ai
import athena_output as ao
import helpers
helpers.logging_setup()

import cooling
import yt_custom
from clump_analysis import compute_surface_area

def setup_argparse():
    parser = argparse.ArgumentParser(description="Compute surface areas.",
                                     formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("out_prefix", help = "The prefix to use for the output.")
    parser.add_argument("vtk_filelst", help = "List of vtk files (expands `*`)", nargs = "+")
    parser.add_argument("-nstep", help = "Only analyze every `nstep` files.", default = 1,
                        type = int)

    return parser.parse_args()

    
def main():
    args = setup_argparse()
    if '*' in args.vtk_filelst:
        fnlst = glob.glob(args.vtk_filelst)
        if yt.is_root(): logging.info("Using path.")
    elif not os.path.isfile(args.vtk_filelst[0]):
        fnlst = glob.glob(os.path.abspath(args.vtk_filelst[0]) + "/*.vtk")
        if yt.is_root(): logging.info("Using all files in folder.")
    else:
        fnlst = args.vtk_filelst
        if yt.is_root(): logging.info("Using files given.")

    assert not '.vtk' in args.out_prefix, "vtk file given as out prefix...?"


    fnlst = sorted([ i for i in fnlst if os.path.isfile(i) ])[::args.nstep]

    assert len(fnlst) > 0, "No vtk files found...?"

    # General prep
    cfg = ai.read(ai.get_filename())

    tcool_func =  cooling.Cooling().tcool
    def field_tcool(field, data):
        return yt_custom.field_tcool(field, data, tcool_func) / ai.get_tcrush(cfg)

    storage = {}

    for sto, cfn in yt.parallel_objects(fnlst, storage = storage):
        logging.info(">>>> Computing quantities for %s", cfn)
        helpers.logging_flush()

        ds = yt.load(cfn)
        ds.add_field(('gas','temp'),function=yt_custom.field_temp, units='cm*dyne/g')
        yt_custom.add_units(ds, cfg)
        ad = ds.all_data()

        sto.result_id = cfn
        sto.result = {}

        # Add cooling time
        ds.add_field(('gas','t_cool'),function=field_tcool, units='')

        # Add additional fields to store here
        ncpus_area = 6
        sto.result['time'] = float(ds.current_time.value)
        th = [1.1,2.,3.,5.,7.,9.]
        cA = compute_surface_area(ds, th, ncpus = ncpus_area)
        for i in range(len(th)):
            sto.result['area_temp_%g' %(th[i])] = cA[i]
        logging.info("Done temp.")

        th = [0.2, 0.5, 0.7, 0.9]
        cA = compute_surface_area(ds, th, mode = 'density', ncpus = ncpus_area)
        for i in range(len(th)):
            sto.result['area_density_%g' %(th[i])] = cA[i]
        logging.info("Done density.")

        th = [0.25, 0.5, 1., 2., 4.]
        cA = compute_surface_area(ds, th, mode = 't_cool', ncpus = ncpus_area)
        for i in range(len(th)):
            sto.result['area_t_cool_%g' %(th[i])] = cA[i]
        logging.info("Done t_cool.")
        

        del ad
        del ds
        gc.collect()
        logging.info(">>>> Done with %s", cfn)

    logging.info(">>> Done with loop.")
    helpers.logging_flush()
    
    if yt.is_root():
        outpre = args.out_prefix
        logging.info("Gathering all data and save to %s*", outpre)
        
        df = pandas.DataFrame(storage).T
        # sort columns
        colnames = ['time'] + [ i for i in df.keys() if i != 'time']
        df = df[colnames]
        df.to_hdf(outpre + ".hdf5", 'dat')
        hdr = ['# ' + df.columns[0]] + [i for i in df.columns[1:]]
        df.to_csv(outpre + ".dat", sep = "\t", index=False, header=hdr, na_rep = 'nan')

        logging.info(">>> Done with everything.")
    return storage

if __name__ == "__main__":
    dat = main()

    
