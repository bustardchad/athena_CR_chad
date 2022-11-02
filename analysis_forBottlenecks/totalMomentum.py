# plots 1D profile of CR energy density E_cr for various times

import matplotlib.pyplot as plt
import yt
from yt.units import erg, pc
from yt.units.yt_array import YTQuantity
import numpy as np

# conversion factors
edenstocgs = 6.54e-11
denstocgs = 6.85e-27
Myr = 1.
kpc = 1.

yt.enable_parallelism()

base="cr.out2."

labels = []
times  = range(0,100,1)   # use this to determine which times to plot

plt.clf()

totalMomentumArr = []


def momentumX(field,data): #was defined as negative in Athena++, so needs a minus sign
      return data['vel1']*1.e8*data['density']*denstocgs*data['dx']*data['dy']*(3.0856e21)*(3.0856e21)

# conversion factors
yt.add_field(("gas","momentumX"), function=momentumX,units='code_mass/code_time',display_name="X-Momentum")


ts = yt.DatasetSeries('../cr.out1*',parallel=False)
i = 0
for ds in ts.piter():
#for i in times:
# ds = yt.load(dir+base+str(i).zfill(5)+'.athdf')    # load the data
 dd = ds.all_data()
# print(dd.quantities.extrema('user_out_var7'))
 time = ds.current_time.v/Myr                       # store the time (with units removed)
 
 totalMomentum = dd.quantities.total_quantity("momentumX")
# totalColl = sum(collisions)*3.155e13
# totalDamp = sum(damping)*3.155e13
# totalColl = -dd.sum(["user_out_var2"])
 totalMomentumArr.append(totalMomentum)
 i = i+1

 
print("momentumXTotal")
print(totalMomentumArr)

plt.plot(times,totalMomentumArr[0:len(times)])
plt.xlabel("Time (Myrs)")
plt.ylabel(r"Total X-Momentum (g cm$^{-2}$ s$^{-1}$)")
#plt.ylim(1.e-13,3.e-11)
#plt.title("$v_{st} = v_{A}^{ion}$ + Ion-Neutral Damping")
plt.title("$v_{st} = v_{A}$")
#plt.legend()
plt.tight_layout()
plt.savefig('TotalXMomentum.png')
plt.close()
