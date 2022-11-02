import yt
from yt.units import dimensions
from yt.units import erg
from yt.units.yt_array import YTQuantity
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

yt.enable_parallelism()

edenstocgs = 6.54e-11
denstocgs = 6.85e-27
gammatocgs = edenstocgs*denstocgs/1.67e-24 # ecr * rho/1.67e-24

Myr = 1.
kpc = 1.

def _ecr(field, data):
  return data['Ec']*edenstocgs*YTQuantity(1,"erg/cm**3")
def _d(field, data):
  return data['density']*denstocgs


def collisional(field,data): #was defined as negative in Athena++, so needs a minus sign
          return (1.022e-15)*(data['density']*denstocgs/1.67e-24) * (data['Ec']*edenstocgs) * YTQuantity(1,"erg/s/g")

def heating(field,data): #was defined as negative in Athena++, so needs a minus sign
        #  return data['user_out_var7']*gammatocgs*YTQuantity(1,"erg/cm**3/s")
          return abs(data['user_out_var7'])*(edenstocgs/3.155e13)*YTQuantity(1,"erg/cm**3/s")

def ionAlfven(field,data): #was defined as negative in Athena++, so needs a minus sign
          return data['user_out_var5'] * 3.0856e21/3.155e13 * YTQuantity(1,"cm/s")





fn = "../cr.out1.00030.athdf"
ds = yt.load(fn) # load data
fn2 = "../cr.out2.00030.athdf"
ds2 = yt.load(fn2) # load data

ds.add_field(('gas', u'd'), function = _d, units="g/cm**3",display_name=r"Density", dimensions=dimensions.density)

ds2.add_field(("gas","ionAlfven"), function=ionAlfven,display_name=r"$v_{A}^{ion}$",units="cm/s")

ds2.add_field(("gas","heating"), function=heating,display_name="Collisionless Loss",units="erg/cm**3/s")

ds.add_field(('gas', u'ecr'), function = _ecr, units="g/(cm*s**2)",display_name=r"E$_{cr}$", dimensions=dimensions.pressure)
ds.add_field(("gas","collisional"), function=collisional,display_name="Collisional Loss",units="erg/cm**3/s")

fig = plt.figure()

# See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
# These choices of keyword arguments produce a four panel plot that includes
# four narrow colorbars, one for each plot.  Axes labels are only drawn on the
# bottom left hand plot to avoid repeating information and make the plot less
# cluttered.
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (1, 4),
                axes_pad = 0.5,
                label_mode = "1",
                share_all = True,
                cbar_location="right",
                cbar_mode="each",
                cbar_size="3%",
                cbar_pad="0%",
		direction='column')

fields = ['ecr', 'ionAlfven', 'heating', 'collisional']

for j in range(4):
  plotvar = fields[j]
  if (j==0):
    varmin = 1.e-14
    varmax = 1.e-11
    slc = yt.SlicePlot(ds, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc])
    slc.set_zlim(plotvar, varmin, varmax)
    slc.set_xlabel('x (kpc)')
    slc.set_ylabel('y (kpc)')
    slc.set_cmap(field=plotvar, cmap='viridis')
   # slc.set_colorbar_label(plotvar,r"Cosmic Ray Energy Density (erg cm$^{-3}$)")
    slc.set_background_color(plotvar)
    # slc.annotate_contour('d',ncont = 1,clim=[2e-25,2e-25],take_log = True, label = False,
    #                plot_args={"colors": "white", "linewidths": 0.5})
    slc.annotate_contour('d',ncont = 1,clim=[2e-24,2e-24],take_log = True, label = False,
                   plot_args={"colors": "white", "linewidths": 1.0})
    slc.annotate_contour('d',ncont = 1,clim=[2e-23,2e-23],take_log = True, label = False,
                   plot_args={"colors": "yellow", "linewidths": 1.0})
    slc.annotate_contour('d',ncont = 1,clim=[2e-22,2e-22],take_log = True, label = False,
                   plot_args={"colors": "red", "linewidths": 1.0})
    slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=6,plot_args={'linewidth':0.3})
    slc.set_width((1.0*kpc, 2.0*kpc))
    # slc.set_width((0.25*kpc, 0.25*kpc))
    #slc.set_width((0.5*kpc, 0.5*kpc))
    slc.set_log(plotvar, True)
    plot = slc.plots[plotvar]
    plot.figure = fig
    plot.axes = grid[j].axes
    plot.cax = grid.cbar_axes[j]
    slc._setup_plots()
    fig.set_size_inches(14,8)
  if (j==1):
    varmin = 1.e5
    varmax = 1.e8
    slc = yt.SlicePlot(ds2, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc])
    slc.set_zlim(plotvar, varmin, varmax)
    slc.set_cmap(field=plotvar, cmap='afmhot')
    slc.set_xlabel('x (kpc)')
    slc.set_ylabel('y (kpc)')
    slc.set_width((1.0*kpc,2.0*kpc))
    slc.set_background_color(plotvar)
    # plot = slc.plots[plotvar]
    # colorbar=plot.cb
    # slc._setup_plots()
    # colorbar.set_ticks([1e-28])
    # colorbar.set_ticklabels(['$10^{-28}$']
    slc.set_log(plotvar, True)
    plot = slc.plots[plotvar]
    plot.figure = fig
    plot.axes = grid[j].axes
    plot.cax = grid.cbar_axes[j]
    slc._setup_plots()
    fig.set_size_inches(14,8)
  if (j==2):
    varmin = 1.e-28
    varmax = 1.e-25
    slc = yt.SlicePlot(ds2, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc])
    slc.set_zlim(plotvar, varmin, varmax)
    slc.set_cmap(field=plotvar, cmap='plasma')
   # slc.set_colorbar_label(plotvar,r"Collisionless Loss (erg cm$^{-3}$ s$^{-1}$)")
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
    plot = slc.plots[plotvar]
    plot.figure = fig
    plot.axes = grid[j].axes
    plot.cax = grid.cbar_axes[j]
    slc._setup_plots()
    fig.set_size_inches(14,8)
  if (j==3):
    varmin = 1.e-28
    varmax = 1.e-25
    slc = yt.SlicePlot(ds, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc])
    slc.set_zlim(plotvar, varmin, varmax)
    slc.set_xlabel('x (kpc)')
    slc.set_ylabel('y (kpc)')
    slc.set_cmap(field=plotvar, cmap='plasma')
   # slc.set_colorbar_label(plotvar,r"Collisional Loss (erg cm$^{-3}$ s$^{-1}$)")
    slc.set_background_color(plotvar)
    # slc.annotate_contour('d',ncont = 1,clim=[2e-25,2e-25],take_log = True, label = False,
    #                plot_args={"colors": "white", "linewidths": 0.5})
    slc.annotate_contour('d',ncont = 1,clim=[2e-24,2e-24],take_log = True, label = False,
                   plot_args={"colors": "white", "linewidths": 1.0})
    slc.annotate_contour('d',ncont = 1,clim=[2e-23,2e-23],take_log = True, label = False,
                   plot_args={"colors": "cyan", "linewidths": 1.0})
    slc.annotate_contour('d',ncont = 1,clim=[2e-22,2e-22],take_log = True, label = False,
                  plot_args={"colors": "green", "linewidths": 1.0})
    # slc.annotate_streamlines('magnetic_field_x','magnetic_field_y',density=2)
    slc.set_width((1.0*kpc, 2.0*kpc))
    # slc.set_width((0.25*kpc, 0.25*kpc))
    #slc.set_width((0.5*kpc, 0.5*kpc))
    slc.set_log(plotvar, True)

    plot = slc.plots[plotvar]
    plot.figure = fig
    plot.axes = grid[j].axes
    plot.cax = grid.cbar_axes[j]
   # cbar=fig.colorbar(plot.cax,orientation='horizontal')
   # plot.cax.colorbar(orientation='horizontal')
    slc._setup_plots()
    fig.set_size_inches(14,8)


# For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
# axes.
#for i, field in enumerate(fields):
#    plot = slc.plots[field]
#    plot.figure = fig
#    plot.axes = grid[i].axes
#    plot.cax = grid.cbar_axes[i]

# Finally, redraw the plot on the AxesGrid axes.

plt.savefig('multiplot_damping.png')
