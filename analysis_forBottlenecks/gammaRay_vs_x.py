import yt
from yt.units.yt_array import YTQuantity
from yt.units import dimensions
import numpy as np
from numpy import *
import matplotlib as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import h5py
from yt.fields import interpolated_fields
from yt.fields.field_detector import \
    FieldDetector
from yt.utilities.linear_interpolators import \
    BilinearFieldInterpolator

yt.enable_parallelism()

pUnit = YTQuantity(1, 'cm**2/s**2')

# conversion factors
denstocgs = 6.85e-27
edenstocgs = 6.54e-11
Myr = 1.
kpc = 1.



def _ecr(field, data):
  return data['Ec']*edenstocgs*YTQuantity(1,"erg/cm**3")

def _fcr(field, data):
  return data['Fc1']*edenstocgs*(3.0856e21/3.155e13) * YTQuantity(1,"erg/cm**2/s")

def _d(field, data):
  return data['density']*denstocgs

def gammaRay_lum(field,data): #was defined as negative in Athena++, so needs a minus sign
          return 0.3333*0.7*0.728*((3.0856e21)**(3.0))*(1.022e-15)*(data['d']/1.67e-24) * (data['Ec']*edenstocgs) * YTQuantity(1,"erg/s/g")


yt.add_field(('gas', u'd'), function = _d, units="g/cm**3",display_name=r"Density", dimensions=dimensions.density)
yt.add_field(('gas', u'ecr'), function = _ecr, units="g/(cm*s**2)",display_name=r"E$_{cr}$", dimensions=dimensions.pressure)
yt.add_field(('gas', u'fcr'), function = _fcr, units="erg/cm**2/s",display_name=r"CR Energy Flux")
yt.add_field(("gas","gammaRay_lum"), function=gammaRay_lum,display_name="Gamma-Ray Emission",units="erg/cm**3/s")

ts = yt.DatasetSeries("../cr.out1.0*",parallel=10)
#ts = yt.DatasetSeries("lognormal/lowCRPres/noIonAlfven_noDamping/lowCRPres/L10_alpha05/cr.out1.00040*")
rho_ex = []
times = []
for ds in ts.piter():
 #       time = ds.current_time.in_units("Kyrs")
        dd = ds.all_data()
       
        time = str(ds.current_time.in_units('Myr'))
        time = (time[:3]) if len(time) > 3 else time
        t = "t = {} Myrs".format(str(time))
        plot = yt.ProfilePlot(dd, "x", ["gammaRay_lum"],fractional=False, accumulation=False)
       # plot.set_xlim(1.E-27,1.E-21)
        plot.set_log("x", False)
       # plot.set_ylim("gammaRay_lum",1.E34,5.E36)
        plot.set_ylim("gammaRay_lum",1.E35,5.E37)
       # plot.annotate_title(r"$v_{st} = v_{A}^{ion}$ + Ion-Neutral Damping")
        plot.annotate_title(r"$v_{st} = v_{A}$")
        plot.save(suffix='png')
        
