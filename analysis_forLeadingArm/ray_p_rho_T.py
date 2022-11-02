# plots slices of density for various times

import yt
from yt.units import dimensions
from yt.units import pc
import matplotlib.pyplot as plt
import numpy as np

yt.enable_parallelism()

chi = 100
T0 = 1.2e-12
P0 = 0.0167
i = 0
ts = yt.DatasetSeries('../merged/*.vtk')
for ds in ts:
 ad = ds.all_data()
 cold_ad = ad.cut_region(["obj['density'] > 33.3"])

 print(ad.quantities.extrema('pressure'))
 print(ad.quantities.extrema('density'))
 print(ad.quantities.extrema('temperature'))
 ray = ds.ortho_ray(0, (0, 0))                      # define the line along which the profiles are taken
 xray = ray['x']                                 # x coordinate
 pray = ray['pressure']/P0                        # E_cr
 dray = ray["density"]/chi
 tray = ray["temperature"]/T0
 time = ds.current_time/7.5
 plt.semilogy(xray, pray, label = r"$P/P_{cl}$")
 plt.semilogy(xray, dray, label = r"$\rho/\rho_{cl}$")
 plt.semilogy(xray, tray, label = r"$T/T_{cl}$")
 plt.ylim(0.005,150)
 plt.xlim(-1,1)
 plt.xlabel('x')
 plt.legend()
 plt.tight_layout()
 plt.title(r"t = %3.2f t$_{cc}$" % time)
 plt.savefig('ray_full_'+str(i).zfill(3)+'.png')
 plt.close()
 i+=1

