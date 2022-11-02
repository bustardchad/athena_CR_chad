# plots slices of density for various times

import yt
from yt.units import dimensions
from yt.units import pc

yt.enable_parallelism()

mach = 4.5

units_override = {"length_unit":(120.0,"pc"),
                  "mass_unit":(5.076e61,"g"),
                  "time_unit":(0.3,"Myr")}

def column(field,data): #was defined as negative in Athena++, so needs a minus sign
         # return data['density']*6.171e16*mach
          return data['density']*1e-4


yt.add_field(('gas', 'column'), function = column, units='g/cm**3')



plotvar = 'density'

ts = yt.DatasetSeries('../merged/cloud.009*',parallel=10,units_override=units_override)
for ds in ts.piter():
 ad = ds.all_data()
# print(ds.field_list)
# ds.define_unit("myUnit", (900, "cm"))
# print(ad.quantities.extrema("density"))
 cold_ad = ad.cut_region(["obj['density'] > 33.3"])
# time = ds.current_time.v/Myr
 slc = yt.ProjectionPlot(ds, 'z', plotvar,fontsize=24,weight_field="ones") #,data_source=cold_ad)
# slc.set_zlim(plotvar, 1.e17, 1.e21)
 slc.set_zlim(plotvar, 1, 1.e2)
 slc.set_cmap(field=plotvar, cmap='kamae')
# slc.set_axes_unit('pc')
 slc.set_xlabel("Distance in kpc ( x n$_{-4}$)")
 slc.set_ylabel("")
# slc.set_colorbar_label(plotvar,r"N$_{tot}$ (cm$^{-2}$)")
 slc.set_colorbar_label(plotvar,r"Overdensity $\chi$")
 slc.set_background_color(plotvar)
 tintcc = ds.current_time.v*(0.3/15) * 1.0
 time = ds.current_time.v*0.3
 slc.hide_axes()
# slc.annotate_title(r"Time = %2.1ft$_{cc}$ = %2.1f/n$_{-4}$ Myrs" % (tintcc, time))
# slc.annotate_contour('column',ncont = 1,clim=[1e18,1e18],take_log = True, label = False,
#               plot_args={"colors": "white", "linewidths": 2.0})
# slc.set_width((150, 30))
# slc.set_log(plotvar, False)
# slc.annotate_streamlines('Bcc1','Bcc2')
# slc.annotate_quiver('vel1','vel2')
# slc.annotate_title("t = %3.0f Myr" % time)
 slc.save()
