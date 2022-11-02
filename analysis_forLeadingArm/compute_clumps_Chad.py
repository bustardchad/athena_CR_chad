import matplotlib
#matplotlib.use('Agg')

import matplotlib.pyplot as plt

import numpy as np
import scipy.interpolate
import pandas
import time
import yt
yt.enable_parallelism()

from clump_analysis import clump_analysis


areas = []
times = []
ts = yt.load('Mach6/cells16/*.vtk',parallel=False)
for ds in ts.piter():
        ad = ds.all_data()

        # Add additional fields to store here
        ncpus_area = 1
        times.append(float(ds.current_time.value))
       # th = [1.1,2.,3.,5.,7.,9.]
        th = [1,1,1]
       # cA = clump_analysis(ds,threshold = 0.001333334,thresh_hot=0.01,outdir='Mach15/cells16')
        cA = clump_analysis(ds,threshold = 9.3258e-12,thresh_hot=6.994e-11,outdir='Mach6/cells16')
       # cA = clump_analysis(ds,threshold = 2.0,thresh_hot=15.0,thresh_unit="T_cl", outdir='Mach15/cells16')

"""
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
"""
