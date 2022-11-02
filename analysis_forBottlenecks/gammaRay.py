# plots slices of gamma ray emission for various times
import yt
from yt.units import dimensions
from yt.units import erg
from yt.units.yt_array import YTQuantity
import numpy as np
import matplotlib.pyplot as plt

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
dir="./lognormal/medCRPres/ionAlfven_Damping/L10_alpha05/"
base="cr.out2."

plotvar = 'gammaRay_emission'
varmin = 1.e-30
varmax = 1.e-26
#varmin = 1.e-34
#varmax = 1.e-30


timeArr = []
gammaArr = []
times = range(0,40,1)

for i in times:
 ds = yt.load(dir+base+str(i).zfill(5)+'.athdf')

 
 dd = ds.all_data()
 time = ds.current_time.v/Myr

 slc = yt.SlicePlot(ds, 'z', plotvar)
 slc.set_zlim(plotvar, varmin, varmax)
 slc.set_cmap(field=plotvar, cmap='plasma')
 slc.set_colorbar_label(plotvar,r"Gamma Ray Emission (erg cm$^{-3}$ s$^{-1}$)")
# slc.set_width((1.4*kpc, 0.5*kpc))
# slc.set_width((0.25*kpc, 0.25*kpc))
 #slc.set_width((0.5*kpc, 0.5*kpc))
 slc.set_log(plotvar, True)
 slc.annotate_title("t = %3.0f Myr" % time)
 slc.save(dir+'plots/'+plotvar+'_slice'+str(i).zfill(3)+'.png')

 totalLum = dd.quantities.total_quantity("gammaRay_lum")
 timeArr.append(time)
 gammaArr.append(totalLum)


plt.semilogy(timeArr,gammaArr)
plt.xlabel('Time (Myrs)',fontsize=20)
plt.ylabel('Gamma-Ray Luminosity (ergs/s)',fontsize=20)
plt.ylim(1.e34,1e37)
#plt.title("2D Cloud, Everett+ 2011 Parameters")
#plt.title("2D Cloud")
plt.tight_layout()
plt.savefig(dir+"plots/gammaRayLuminosity.pdf")
plt.close()

