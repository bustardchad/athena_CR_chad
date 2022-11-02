# plots slices of density for various times
import yt
from yt.units.yt_array import YTQuantity
from yt.units import dimensions
from yt.units import pc, Myr
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


def diffusionLength(field,data): #was defined as negative in Athena++, so needs a minus sign
      return 3.e29*data['user_out_var6']*YTQuantity(1,'cm')/3.e10

# conversion factors
yt.add_field(("gas","diffusionLength"), function=diffusionLength,units="pc",display_name="Diffusion Length")

#dir="./ionneutral/pres1/test_March20_vmax100/"
plotvar = 'diffusionLength'
varmin = 1.e-3
varmax = 1.e2

ts = yt.load('../cr.out2*',parallel=10)

for ds in ts.piter():
 dd = ds.all_data()
 time = ds.current_time.v/Myr

 slc = yt.SlicePlot(ds, 'z', plotvar)
 slc.set_zlim(plotvar, varmin, varmax)
 slc.set_cmap(field=plotvar, cmap='plasma')
 slc.set_log(plotvar, True)
 slc.annotate_title("t = %3.0f Myr" % time)
 slc.save()
