import matplotlib
#matplotlib.use('Agg')

import matplotlib.pyplot as plt

import numpy as np
import scipy.interpolate
#import pandas
import time
import yt
yt.enable_parallelism()

from clump_analysis_Chad import compute_surface_area


areas = []
times = []
ts = yt.load('*.vtk',parallel=False)
for ds in ts.piter():
        ad = ds.all_data()

        # Add additional fields to store here
        ncpus_area = 1
        times.append(float(ds.current_time.value))
       # th = [1.1,2.,3.,5.,7.,9.]
        th = [1,1,1]
        cA = compute_surface_area(ds, th, ncpus = ncpus_area)
        print("Area: %g" %(cA))
        areas.append(cA)

print("times: ")
print(times)
print("areas: ")
print(areas)
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
