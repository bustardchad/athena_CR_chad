import yt
from yt.units.yt_array import YTQuantity
from yt.units import dimensions
from yt.units import pc
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

def gammaRay_lum(field,data): #was defined as negative in Athena++, so needs a minus sign
      return -0.3333*0.7*0.728*((3.0856e21)**(3.0))*data['user_out_var2']*YTQuantity(1,'erg/s')

def diffusionLength(field,data): #was defined as negative in Athena++, so needs a minus sign
      return 3.e29*data['user_out_var6']*YTQuantity(1,'cm')/3.e10

# conversion factors
yt.add_field(("gas","gammaRay_lum"), function=gammaRay_lum,units="erg/s")
yt.add_field(("gas","diffusionLength"), function=diffusionLength,units="pc")

ts = yt.DatasetSeries("../cr.out2.0*",parallel=10)
#ts = yt.DatasetSeries("lognormal/lowCRPres/noIonAlfven_noDamping/lowCRPres/L10_alpha05/cr.out1.00040*")
rho_ex = []
times = []
for ds in ts.piter():
 #       time = ds.current_time.in_units("Kyrs")
        dd = ds.all_data()
       
        time = str(ds.current_time.in_units('Myr'))
        time = (time[:3]) if len(time) > 3 else time
        t = "t = {} Myrs".format(str(time))
        plot = yt.ProfilePlot(dd, "diffusionLength", ["gammaRay_lum"],fractional=True, accumulation=False,weight_field=None)
        plot.set_xlim(1.E-3,1.E2)
        plot.set_ylim("gammaRay_lum",1.E-3,1.E0)
        plot.annotate_title(r"$v_{st} = v_{A}^{ion}$ + Ion-Neutral Damping")
       # plot.annotate_title(r"$v_{st} = v_{A}$")
        plot.save(suffix='png')
        
