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

#dir="./noCooling/GMC_lowerDensity/30MyrPulse_higherCRPres/Gammas_IonAlfven_Damping/"
#dir="./noCooling/Everett2011Cloud/30MyrPulse/widerInterface/higherB_lowCRPres/"
#dir="./noCooling/Everett2011Cloud/30MyrPulse/widerInterface_higherMaxVst/ionfunction2/min_fion1e-6/maxVst10/lowB_lowCRPres/"
#dir="./noCooling/lowerDensCloud/higherIonFraction/ionAlfven_Damping/fraction1e-2/"
#dir="./noCooling/lowerDensCloud/higherIonFraction_higherVMax/ionAlfven_Damping/fraction1e-1/"
#dir="./noCooling/lowerDensCloud/higherIonFraction/ionAlfven_Damping/fraction1e-1/moreuservars/"
#dir="./noCooling/lowerDensCloud_rad10pc/noIonAlfven_noDamping/res4/interface_5pc/"
dir="./noCooling/lowerDensCloud_rad10pc/ionAlfven_Damping/res1/interface_5pc/"
#dir="./noCooling/lowerDensCloud/higherIonFraction/ionAlfven_Damping/fraction1e-2/moreuservars/"
base="cr.out2."


def collisionalLum(field,data): #was defined as negative in Athena++, so needs a minus sign
          return -data['user_out_var2']*data['dx']*data['dy']*(3.0856e21)**2.0*YTQuantity(1,"erg/s/cm**3")

def collisionlessLum(field,data): #was defined as negative in Athena++, so needs a minus sign
          return -data['user_out_var7']*data['dx']*data['dy']*(3.0856e21)**2.0*YTQuantity(1,"erg/s/cm**3")*edenstocgs/3.155e13

yt.add_field(("gas","collisionalLum"), function=collisionalLum,display_name="Collisional",units="erg/cm/s")
yt.add_field(("gas","collisionlessLum"), function=collisionlessLum,display_name="Collisionless",units="erg/cm/s")


labels = []
times  = range(0,100,1)   # use this to determine which times to plot

plt.clf()

totalCollArr = []
totalDampArr = []
totalCRLossArr = []
ts = yt.DatasetSeries('../cr.out2*',parallel=False)
i = 0
for ds in ts.piter():
#for i in times:
# ds = yt.load(dir+base+str(i).zfill(5)+'.athdf')    # load the data
 dd = ds.all_data()
# print(dd.quantities.extrema('user_out_var7'))
 time = ds.current_time.v/Myr                       # store the time (with units removed)
 ray = ds.ortho_ray(0, (0, 0))                      # define the line along which the profiles are taken
 xpc = ray['x']/kpc                                 # x coordinate
 collisions = -ray['user_out_var2']                        # E_cr
 damping= -ray["user_out_var7"]*edenstocgs/3.155e13
 
# totalColl = -dd.quantities.total_quantity("user_out_var2")
 totalColl = dd.quantities.total_quantity("collisionalLum")
# totalColl = sum(collisions)*3.155e13
# totalDamp = sum(damping)*3.155e13
# totalColl = -dd.sum(["user_out_var2"])
# totalDamp = -dd.quantities.total_quantity("user_out_var7")*edenstocgs/3.155e13
 totalDamp = dd.quantities.total_quantity("collisionlessLum")
 totalCRLoss = totalColl+totalDamp
# print(totalColl)
 totalCollArr.append(totalColl)
 totalDampArr.append(totalDamp)
 totalCRLossArr.append(totalCRLoss)
 i = i+1

 
print("Collisional")
print(totalCollArr)
print("Collisionless")
print(totalDampArr)
print("Total")
print(totalCRLossArr)

plt.plot(times,totalCollArr[0:len(times)],label="Collisional")
plt.plot(times,totalDampArr[0:len(times)],label="Collisionless")
plt.plot(times,totalCRLossArr[0:len(times)],label="Total")
plt.xlabel("Time (Myrs)")
plt.ylabel(r"CR Energy Loss (erg cm$^{-3}$ s$^{-1}$)")
#plt.ylim(1.e-13,3.e-11)
plt.title("$v_{st} = v_{A}^{ion}$ + Ion-Neutral Damping")
#plt.title("$v_{st} = v_{A}$)
plt.legend()
plt.tight_layout()
plt.savefig('TotalCRLoss_lum.png')
plt.close()
