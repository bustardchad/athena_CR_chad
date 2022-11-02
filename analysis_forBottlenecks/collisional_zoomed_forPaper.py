# plots slices of density for various times

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

def collisions(field,data): #was defined as negative in Athena++, so needs a minus sign
        #  return data['user_out_var7']*gammatocgs*YTQuantity(1,"erg/cm**3/s")
          return -data['user_out_var2']*YTQuantity(1,"erg/cm**3/s")


yt.add_field(("gas","collisions"), function=collisions,display_name="Cosmic-Ray Collisional Losses",units="erg/cm**3/s")


plotvar = 'collisions'
varmin = 1.e-28
varmax = 1.e-25

times = range(0,46,1)



ts = yt.DatasetSeries('../cr.out2.00030*',parallel=10)
for ds in ts.piter():
 dd = ds.all_data()
 print(ds.field_list)
 time = ds.current_time.v/Myr
# slc = yt.SlicePlot(ds, 'z', plotvar, origin='native', center=[1.5*kpc, 0, 0])
# slc.set_width((0.6*kpc, .25*kpc))

 slc = yt.SlicePlot(ds, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc],fontsize=20)
 slc.set_zlim(plotvar, varmin, varmax)
 slc.set_cmap(field=plotvar, cmap='plasma')
 slc.set_colorbar_label(plotvar,r"Cosmic-Ray Collisional Losses (erg cm$^{-3}$ s$^{-1}$)")
 slc.set_xlabel('x (kpc)')
 slc.set_ylabel('y (kpc)')
 slc.set_width((1.0*kpc, 2.0*kpc))
 slc.set_background_color(plotvar)
# plot = slc.plots[plotvar]
# colorbar=plot.cb
# slc._setup_plots()
# colorbar.set_ticks([1e-28])
# colorbar.set_ticklabels(['$10^{-28}$']
 slc.set_log(plotvar, True)
# slc.annotate_streamlines('Bcc1','Bcc2')
# slc.annotate_quiver('vel1','vel2')
 slc.annotate_title("t = %3.0f Myr" % time)
 slc.save(suffix='pdf')
