# plots slices of total pressure and pressure components at various times

import yt
from yt.units import pc
from yt.units import dimensions
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

# conversion factors
denstocgs = 6.85e-27
edenstocgs = 6.54e-11
kpc    = 1.
Myr    = 1.

scale = 1.e-13			# scale chosen to keep color bar labels clean

def _pcr(field, data):
 return ds.arr(data['Ec'].v*edenstocgs/(3.*scale), 'g/(cm*s**2)')
def _pth(field, data):
 return ds.arr(data['press'].v*edenstocgs/scale, 'g/(cm*s**2)')
def _pmag(field, data):
 return ds.arr(((data['Bcc1']**2+data['Bcc2']**2)/2.).v*edenstocgs/scale, 'g/(cm*s**2)')
def _ptot(field, data):
 return (data['pmag']+data['pth']+data['pcr'])*scale	# keep original nonscaled values for the total pressure

yt.add_field(('gas', u'pcr'), function = _pcr, units="g/(cm*s**2)", dimensions=dimensions.pressure)
yt.add_field(('gas', u'pth'), function = _pth, units="g/(cm*s**2)", dimensions=dimensions.pressure)
yt.add_field(('gas', u'pmag'), function = _pmag, units="g/(cm*s**2)", dimensions=dimensions.pressure)
yt.add_field(('gas', u'ptot'), function = _ptot, units="g/(cm*s**2)", dimensions=dimensions.pressure)

#dir="./Everett2011Cloud_noDamping/"
dir="./exactCooling/noCRSource/setup1/"
base="cr.out1."

fig = plt.figure()
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (4,1),
                direction = "column",
                axes_pad = 0.0,
                label_mode = "1",
                #share_all = True,
                cbar_location = "right",
                cbar_mode = "each",
                cbar_size = "2%",
                cbar_pad = "0%")

times = range(0,46,1)

for i in times:
 ds = yt.load(dir+base+str(i).zfill(5)+'.athdf')
 time = ds.current_time.v/Myr

 plotvars = ['ptot','pth','pcr','pmag']
 varmax = [5.2e-13, 4.8, 1.1, 4.8]											# max for each pressure (note ptot is not scaled)
 varmin = [2.6e-13, 0.01, 0.01, 0.01]											# min for each pressure
 colors = ['Greens', 'Blues', 'Reds', 'Purples']									# color schemes 
 clabels = ['$P_{tot}$ $\left(\\frac{\\mathrm{dyne}}{\\mathrm{cm}^2}\\right)$', '$P_{th}$', '$P_{cr}$', '$P_B$']	# color bar labels

 for j in range(4):
  var = plotvars[j]
  vmax = varmax[j]
  vmin = varmin[j]
  vcol = colors[j]

  slc = yt.SlicePlot(ds, 'z', var, origin='native', center=[0.9*kpc, 0*pc, 0*pc])
  slc.set_cmap(field=var, cmap=vcol)
  slc.set_zlim(var, vmin, vmax)
  slc.set_log(var, False)
  slc.set_width((1.4*kpc, 0.25*kpc))
  slc.set_colorbar_label(var, clabels[j])
  #slc.set_font({'size': 24})
  if j==0:
   slc.annotate_title("t = %3.0f Myr" % time)

  plot = slc.plots[var]
  plot.figure = fig
  plot.axes = grid[j].axes
  plot.cax = grid.cbar_axes[j]
 
  slc._setup_plots()
  fig.set_size_inches(12,8)

 plt.savefig(dir+'plots/pres_'+str(i).zfill(2)+'.png')
