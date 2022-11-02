# plots slices of CR energy density for various times

import yt
from yt.units import dimensions


yt.enable_parallelism()

# conversion factors
denstocgs = 6.85e-27
edenstocgs = 6.54e-11
Myr = 1.
kpc = 1.

# define a new field for yt to use
# this might be more than is necessary, but this ensures the quantity plotted is
# in cgs units and also HAS cgs units
def _ecr(field, data):
 return ds.arr(data['Ec'].v*edenstocgs, 'g/(cm*s**2)')

yt.add_field(('gas', u'ecr'), function = _ecr, units="g/(cm*s**2)", dimensions=dimensions.pressure)

dir="./../"
base="cr.out1."

plotvar = 'ecr'
varmax = 1.e-11									# color bar max value
varmin = 1.e-14									# color bar min value

times  = range(0,90,1)								# pick which times to plot

ts = yt.DatasetSeries('../cr.out1*',parallel=10)				# load data
for ds in ts.piter():
 time = ds.current_time.v/Myr							# store time
# slc = yt.SlicePlot(ds, 'z', plotvar, origin='native', center=[1.1*kpc, 0, 0])	# make slice plot (with defined center)
 slc = yt.SlicePlot(ds, 'z', plotvar)	# make slice plot (with defined center)
# slc.set_width((0.6*kpc, .25*kpc))						# set width of domain to be plotted
 slc.set_zlim(plotvar, varmin, varmax)						# set color bar max and min
 slc.set_cmap(field=plotvar, cmap='Reds')					# set color bar color map
 slc.set_log(plotvar, True)							# plot on linear scale
# slc.annotate_streamlines('Bcc1','Bcc2')					# add magnetic field streamlines
# slc.annotate_quiver('vel1','vel2')						# add velocity arrows
 slc.annotate_title("t = %3.0f Myr" % time)					# add title with the time
 slc.save()			# save figure
