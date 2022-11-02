# plots slices of density for various times

import yt
from yt.units import dimensions

# conversion factors
denstocgs = 6.85e-27
Myr = 1.
kpc = 1.

#dir="./ionneutral/pres1/test_March20_vmax100/"
dir="./ionneutral/pres1/noAdaptive/"
base="cr.out2."

plotvar = 'user_out_var5'
varmin = 1.e-2
varmax = 1.e1

times = range(0,46,1)

for i in times:
 ds = yt.load(dir+base+str(i).zfill(5)+'.athdf')
 dd = ds.all_data()
 print(ds.field_list)
 time = ds.current_time.v/Myr
# slc = yt.SlicePlot(ds, 'z', plotvar, origin='native', center=[1.5*kpc, 0, 0])
# slc.set_width((0.6*kpc, .25*kpc))

 print(dd.quantities.extrema("user_out_var5"))

 slc = yt.SlicePlot(ds, 'z', plotvar)
 slc.set_zlim(plotvar, varmin, varmax)
 slc.set_cmap(field=plotvar, cmap='plasma')
 slc.set_colorbar_label(plotvar,r"Streaming Speed (scaled units)")
 slc.set_width((1.4*kpc, 0.25*kpc))
# plot = slc.plots[plotvar]
# colorbar=plot.cb
# slc._setup_plots()
# colorbar.set_ticks([1e-28])
# colorbar.set_ticklabels(['$10^{-28}$']
 slc.set_log(plotvar, True)
# slc.annotate_streamlines('Bcc1','Bcc2')
# slc.annotate_quiver('vel1','vel2')
 slc.annotate_title("t = %3.0f Myr" % time)
 slc.save('plots/'+plotvar+'_slice'+str(i).zfill(2)+'.png')
