# plots slices of gamma ray emission for various times
import yt
from yt.units import dimensions
from yt.units import erg
from yt.units.yt_array import YTQuantity
import numpy as np
import matplotlib.pyplot as plt

yt.enable_parallelism()

edenstocgs = 6.54e-11
denstocgs = 6.85e-27
gammatocgs = edenstocgs*denstocgs/1.67e-24 # ecr * rho/1.67e-24

Myr = 1.
kpc = 1.
def gammaRay_emission(field,data): #was defined as negative in Athena++, so needs a minus sign
      return -0.3333*0.7*0.728*data['user_out_var2']*YTQuantity(1,'erg/cm**3/s')

def gammaRay_lum(field,data): #was defined as negative in Athena++, so needs a minus sign
      return -0.3333*0.7*0.728*((3.0856e21)**(3.0))*data['user_out_var2']*YTQuantity(1,'erg/cm**3/s')*data['cell_volume']

yt.add_field(("gas","gammaRay_emission"), function=gammaRay_emission,units="erg/cm**3/s")

# conversion factors
yt.add_field(("gas","gammaRay_lum"), function=gammaRay_lum,units="erg/s")

# conversion factors

#dir="./noCooling/Everett2011Cloud/temp1e6/30MyrPulse/higherCRPres/Gammas_IonAlfven_Damping/"
#dir="./noCooling/Everett2011Cloud/temp1e6/60MyrPulse/fiducialCRPres/Gammas_noIonAlfven_noDamping/"
#dir="./noCooling/Everett2011Cloud/temp1e6/60MyrPulse/higherCRPres/Gammas_IonAlfven_Damping_Ambipolar/"
#dir="./noCooling/GMC_lowerDensity/30MyrPulse/Gammas_IonAlfven_Damping/"
#dir="./noCooling/Everett2011Cloud/temp1e6/multipleClouds/30MyrPulse/higherCRPres/setup2/"
dir="./noClumps/"  # L = 10, alpha = 0.5
base="cr.out2."

timeArr = []
gammaArr_L10_alpha01 = []
times = range(0,100,1)

ts = yt.DatasetSeries(dir+'cr.out2.00*',parallel=False)
for ds in ts.piter():
 dd = ds.all_data()
 time = ds.current_time.v/Myr
 totalLum = dd.quantities.total_quantity("gammaRay_lum")
 timeArr.append(time)
 gammaArr_L10_alpha01.append(totalLum)

print("noClumps:")
print(gammaArr_L10_alpha01)


dir="./noIonAlfven_noDamping/"  # L = 10, alpha = 0.5
base="cr.out2."

gammaArr_L10_alpha01_no = []
times = range(0,100,1)

ts = yt.DatasetSeries(dir+'cr.out2.00*',parallel=False)
for ds in ts.piter():
 dd = ds.all_data()
 time = ds.current_time.v/Myr
 totalLum = dd.quantities.total_quantity("gammaRay_lum")
 gammaArr_L10_alpha01_no.append(totalLum)

print("noIonAlfven_noDamping")
print(gammaArr_L10_alpha01_no)

dir="./ionAlfven_Damping/"  # L = 10, alpha = 0.5
base="cr.out2."

gammaArr_L10_alpha01_no = []
times = range(0,100,1)

ts = yt.DatasetSeries(dir+'cr.out2.00*',parallel=False)
for ds in ts.piter():
 dd = ds.all_data()
 time = ds.current_time.v/Myr
 totalLum = dd.quantities.total_quantity("gammaRay_lum")
 gammaArr_L10_alpha01_no.append(totalLum)

print("ionAlfven_Damping")
print(gammaArr_L10_alpha01_no)

"""
dir="./L2_alpha01/"  # L = 10, alpha = 0.5
base="cr.out2."

gammaArr_L2_alpha01 = []
times = range(0,100,1)

ts = yt.DatasetSeries(dir+'cr.out2.0*',parallel=False)
for ds in ts.piter():
 dd = ds.all_data()
 time = ds.current_time.v/Myr
 totalLum = dd.quantities.total_quantity("gammaRay_lum")
 gammaArr_L2_alpha01.append(totalLum)

print("L = 2, alpha = 0.1")
print(gammaArr_L2_alpha01)

"""
"""
plt.semilogy(timeArr,gammaArr_L10_alpha01,label="ion Alfven + Damping")
plt.semilogy(timeArr,gammaArr_L10_alpha01_no,label="no ion Alfven and no damping")
#plt.semilogy(timeArr,gammaArr_L2_alpha01,label="L = 2, alpha = 0.1")
plt.xlabel('Time (Myrs)',fontsize=20)
plt.ylabel('Gamma-Ray Luminosity (ergs/s)',fontsize=20)
plt.ylim(1.e34,1e36)
#plt.title("2D Cloud, Everett+ 2011 Parameters")
#plt.title("2D Cloud")
plt.legend()
plt.tight_layout()
plt.savefig("gammaRayLuminosity_L10_alpha01.pdf")
plt.close()
"""
