# plots slices of density for various times

import yt
from yt.units import dimensions
from yt.units import pc

yt.enable_parallelism()

# conversion factors
denstocgs = 6.85e-27
Myr = 1.
kpc = 1.

def _d(field, data):
 return data['density']*denstocgs

yt.add_field(('gas', u'd'), function = _d, units="g/cm**3",display_name=r"Density", dimensions=dimensions.density)


plotvar = 'd'
#varmax = 1.e-24
varmax = 1.e-23
#varmax = 2.3e-27
varmin = 1.e-25

ts = yt.DatasetSeries('../cr.out1*',parallel=10)
for ds in ts.piter():
 time = ds.current_time.v/Myr
 slc = yt.SlicePlot(ds, 'z', plotvar,origin='native', center=[1.0*kpc, 0*pc, 0*pc],fontsize=20)
 slc.set_zlim(plotvar, varmin, varmax)
 slc.set_cmap(field=plotvar, cmap='Blues')
 slc.set_xlabel('x (kpc)')
 slc.set_ylabel('y (kpc)')
 slc.set_width((2.0*kpc, 2.0*kpc))
# slc.set_log(plotvar, False)
# slc.annotate_streamlines('Bcc1','Bcc2')
# slc.annotate_quiver('vel1','vel2')
 slc.annotate_title("t = %3.0f Myr" % time)
 slc.save()
