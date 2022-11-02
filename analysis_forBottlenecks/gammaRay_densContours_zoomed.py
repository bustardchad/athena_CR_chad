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

def _ecr(field, data):
  return data['Ec']*edenstocgs*YTQuantity(1,"erg/cm**3")

def _d(field, data):
  return data['density']*denstocgs

def gammaRay_emission(field,data): #was defined as negative in Athena++, so needs a minus sign
          return 0.3333*0.7*0.728*(1.022e-15)*(data['d']/1.67e-24) * (data['Ec']*edenstocgs) * YTQuantity(1,"erg/s/g")


yt.add_field(('gas', u'd'), function = _d, units="g/cm**3",display_name=r"Density", dimensions=dimensions.density)
yt.add_field(('gas', u'ecr'), function = _ecr, units="g/(cm*s**2)",display_name=r"E$_{cr}$", dimensions=dimensions.pressure)
yt.add_field(("gas","gammaRay_emission"), function=gammaRay_emission,display_name="Gamma-Ray Emission",units="erg/cm**3/s")



# conversion factors

#dir="./noCooling/Everett2011Cloud/temp1e6/30MyrPulse/higherCRPres/Gammas_IonAlfven_Damping/"
#dir="./noCooling/Everett2011Cloud/temp1e6/60MyrPulse/fiducialCRPres/Gammas_noIonAlfven_noDamping/"
#dir="./noCooling/Everett2011Cloud/temp1e6/60MyrPulse/higherCRPres/Gammas_IonAlfven_Damping_Ambipolar/"
#dir="./noCooling/GMC_lowerDensity/30MyrPulse/Gammas_IonAlfven_Damping/"
#dir="./noCooling/Everett2011Cloud/temp1e6/multipleClouds/30MyrPulse/higherCRPres/setup2/"
dir="./../"
base="cr.out1."

plotvar = 'gammaRay_emission'
varmin = 1.e-31
varmax = 1.e-27
#varmin = 1.e-34
#varmax = 1.e-30


timeArr = []
gammaArr = []
times = range(0,40,1)

ts = yt.DatasetSeries('../cr.out1*',parallel=10)
for ds in ts.piter():
 dd = ds.all_data()
 time = ds.current_time.v/Myr

 slc = yt.SlicePlot(ds, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc],fontsize=23)
 slc.set_zlim(plotvar, varmin, varmax)
 slc.set_xlabel('x (kpc)')
 slc.set_ylabel('y (kpc)')
 slc.set_cmap(field=plotvar, cmap='plasma')
 slc.set_colorbar_label(plotvar,r"Gamma Ray Emission (erg cm$^{-3}$ s$^{-1}$)")
 slc.set_background_color(plotvar)
# slc.annotate_contour('d',ncont = 1,clim=[2e-25,2e-25],take_log = True, label = False,
#                plot_args={"colors": "white", "linewidths": 0.5})
 slc.annotate_contour('d',ncont = 1,clim=[2e-24,2e-24],take_log = True, label = False,
                plot_args={"colors": "white", "linewidths": 1.0})
 slc.annotate_contour('d',ncont = 1,clim=[2e-23,2e-23],take_log = True, label = False,
                plot_args={"colors": "cyan", "linewidths": 1.0})
 slc.annotate_contour('d',ncont = 1,clim=[2e-22,2e-22],take_log = True, label = False,
                plot_args={"colors": "green", "linewidths": 1.0})
# slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=2)
 slc.set_width((1.0*kpc, 2.0*kpc))
# slc.set_width((0.25*kpc, 0.25*kpc))
 #slc.set_width((0.5*kpc, 0.5*kpc))
 slc.set_log(plotvar, True)
 slc.annotate_title("t = %3.0f Myr" % time)
 slc.save()

