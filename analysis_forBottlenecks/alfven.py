# plots slices of total pressure and pressure components at various times

import yt
from yt.units import pc
from yt.units import dimensions
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import numpy as np
# conversion factors
denstocgs = 6.85e-27
edenstocgs = 6.54e-11
kpc    = 1.
Myr    = 1.

scale = 1.e-13			# scale chosen to keep color bar labels clean

def _pmag(field, data):
 return ds.arr((data['Bcc1']**2+data['Bcc2']**2)/2., 'g/(cm*s**2)')

def _va(field, data):
 return ds.arr((data['pmag'])/np.sqrt(data['density']), 'cm/s')	# keep original nonscaled values for the total pressure

yt.add_field(('gas', u'pmag'), function = _pmag,units='g/(cm*s**2)')
yt.add_field(('gas', u'va'), function = _va,units='cm/s')

dir="./ionneutral/diffuseFlux/lowerDensCloud/maxNeut1/"
base="cr.out1."

plotvar = 'alfven_speed'

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

times = range(0,146,1)

for i in times:
 ds = yt.load(dir+base+str(i).zfill(5)+'.athdf')
 time = ds.current_time.v/Myr

 slc = yt.SlicePlot(ds, 'z', plotvar)
 slc.set_cmap(field=plotvar,cmap='plasma')
 slc.set_zlim(plotvar, 1.e-3, 1.e0)
 slc.set_log(plotvar, True)
 slc.set_width((1.4*kpc, 0.25*kpc))
 slc.save('plots/va_'+str(i).zfill(2)+'.png')
