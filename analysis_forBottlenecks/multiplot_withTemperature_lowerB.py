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
prestocgs = edenstocgs
temptocgs = prestocgs/denstocgs


Myr = 1.
kpc = 1.

def _ecr(field, data):
  return data['Ec']*edenstocgs*YTQuantity(1,"erg/cm**3")
def _d(field, data):
  return data['density']*denstocgs

def _tmp(field, data):
 return data['temperature']*temptocgs




def collisional(field,data): #was defined as negative in Athena++, so needs a minus sign
          return (1.022e-15)*(data['density']*denstocgs/1.67e-24) * (data['Ec']*edenstocgs) * YTQuantity(1,"erg/s/g")

def heating(field,data): #was defined as negative in Athena++, so needs a minus sign
        #  return data['user_out_var7']*gammatocgs*YTQuantity(1,"erg/cm**3/s")
          return abs(data['user_out_var7'])*(edenstocgs/3.155e13)*YTQuantity(1,"erg/cm**3/s")

def ionAlfven(field,data): #was defined as negative in Athena++, so needs a minus sign
          return data['user_out_var5'] * 3.0856e21/3.155e13 * YTQuantity(1,"cm/s")





fn = "noIonAlfven_noDamping/noClumpsNearBoundary/L5/alpha1_5_lowerB/cr.out1.00040.athdf"
ds = yt.load(fn) # load data

fn2 = "noIonAlfven_noDamping/noClumpsNearBoundary/L5/alpha1_5_lowerB/cr.out2.00040.athdf"
ds2 = yt.load(fn2) # load data

fn3 = "ionAlfven_Damping/noClumpsNearBoundary/L5_alpha1_5/lowerB/cr.out1.00040.athdf"
ds3 = yt.load(fn3) # load data

fn4 = "ionAlfven_Damping/noClumpsNearBoundary/L5_alpha1_5/lowerB/cr.out2.00040.athdf"
ds4 = yt.load(fn4) # load data

ds.coordinates.x_axis[2] = 1
ds.coordinates.x_axis['z'] = 1
ds.coordinates.y_axis[2] = 0
ds.coordinates.y_axis['z'] = 0
ds2.coordinates.x_axis[2] = 1
ds2.coordinates.x_axis['z'] = 1
ds2.coordinates.y_axis[2] = 0
ds2.coordinates.y_axis['z'] = 0
ds3.coordinates.x_axis[2] = 1
ds3.coordinates.x_axis['z'] = 1
ds3.coordinates.y_axis[2] = 0
ds3.coordinates.y_axis['z'] = 0
ds4.coordinates.x_axis[2] = 1
ds4.coordinates.x_axis['z'] = 1
ds4.coordinates.y_axis[2] = 0
ds4.coordinates.y_axis['z'] = 0

ds.add_field(('gas', u'tmp'), function = _tmp, units="K",display_name=r"T", dimensions=dimensions.temperature)
ds.add_field(('gas', u'd'), function = _d, units="g/cm**3",display_name=r"Density", dimensions=dimensions.density)
ds2.add_field(("gas","ionAlfven"), function=ionAlfven,display_name=r"$v_{A}^{ion}$",units="cm/s")
ds2.add_field(("gas","heating"), function=heating,display_name="Collisionless Loss",units="erg/cm**3/s")
ds.add_field(('gas', u'ecr'), function = _ecr, units="g/(cm*s**2)",display_name=r"E$_{cr}$", dimensions=dimensions.pressure)
ds.add_field(("gas","collisional"), function=collisional,display_name="Collisional Loss",units="erg/cm**3/s")

ds3.add_field(('gas', u'tmp'), function = _tmp, units="K",display_name=r"T", dimensions=dimensions.temperature)
ds3.add_field(('gas', u'd'), function = _d, units="g/cm**3",display_name=r"Density", dimensions=dimensions.density)
ds4.add_field(("gas","ionAlfven"), function=ionAlfven,display_name=r"$v_{A}^{ion}$",units="cm/s")
ds4.add_field(("gas","heating"), function=heating,display_name="Collisionless",units="erg/cm**3/s")
ds3.add_field(('gas', u'ecr'), function = _ecr, units="g/(cm*s**2)",display_name=r"E$_{cr}$", dimensions=dimensions.pressure)
ds3.add_field(("gas","collisional"), function=collisional,display_name="Collisional",units="erg/cm**3/s")

fig = plt.figure()

# See http://matplotlib.org/mpl_toolkits/axes_grid/api/axes_grid_api.html
# These choices of keyword arguments produce a four panel plot that includes
# four narrow colorbars, one for each plot.  Axes labels are only drawn on the
# bottom left hand plot to avoid repeating information and make the plot less
# cluttered.
grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (5, 2),
                axes_pad = 0.5,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="edge",
                cbar_size="3%",
                cbar_pad="1%",
		direction='row')


fsize = 18

fields = ['ecr', 'ionAlfven', 'tmp', 'heating', 'collisional']

for j in range(5):
  plotvar = fields[j]
  if (j==0):
    varmin = 1.e-14
    varmax = 1.e-11
    slc = yt.SlicePlot(ds, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc],fontsize=fsize)
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
    #slc.annotate_streamlines('magnetic_field_y','magnetic_field_x',density=6,plot_args={'linewidth':0.3})
    slc.set_width((2.0*kpc, 1.0*kpc))
    slc.annotate_title("$f_{ion}^{min} = 1.0$")
    # slc.set_width((0.25*kpc, 0.25*kpc))
    #slc.set_width((0.5*kpc, 0.5*kpc))
    slc.set_log(plotvar, True)
    plot = slc.plots[plotvar]
    plot.figure = fig
    plot.axes = grid[0].axes
    plot.cax = grid.cbar_axes[0]
    slc._setup_plots()

    slc = yt.SlicePlot(ds3, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc],fontsize=fsize)
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
   # slc.annotate_streamlines('magnetic_field_y','magnetic_field_x',density=6,plot_args={'linewidth':0.3})
    slc.set_width((2.0*kpc, 1.0*kpc))
    slc.annotate_title("$f_{ion}^{min} = 10^{-4}$")
    # slc.set_width((0.25*kpc, 0.25*kpc))
    #slc.set_width((0.5*kpc, 0.5*kpc))
    slc.set_log(plotvar, True)
    plot = slc.plots[plotvar]
    plot.figure = fig
    plot.axes = grid[1].axes
    plot.cax = grid.cbar_axes[j]
    slc._setup_plots()
  if (j==1):
    varmin = 1.e5
    varmax = 1.e8
   # varmin = 1.e3
   # varmax = 1.e6
    slc = yt.SlicePlot(ds2, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc],fontsize=fsize)
    slc.set_zlim(plotvar, varmin, varmax)
    slc.set_cmap(field=plotvar, cmap='afmhot')
   # slc.set_cmap(field=plotvar, cmap='coolwarm')
    slc.set_xlabel('x (kpc)')
    slc.set_ylabel('y (kpc)')
    slc.set_width((2.0*kpc,1.0*kpc))
    slc.set_background_color(plotvar)
    # plot = slc.plots[plotvar]
    # colorbar=plot.cb
    # slc._setup_plots()
    # colorbar.set_ticks([1e-28])
    # colorbar.set_ticklabels(['$10^{-28}$']
    slc.set_log(plotvar, True)
    plot = slc.plots[plotvar]
    plot.figure = fig
    plot.axes = grid[2].axes
    plot.cax = grid.cbar_axes[j]
    slc._setup_plots()

    slc = yt.SlicePlot(ds4, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc],fontsize=fsize)
    slc.set_zlim(plotvar, varmin, varmax)
    slc.set_cmap(field=plotvar, cmap='afmhot')
    #slc.set_cmap(field=plotvar, cmap='coolwarm')
    slc.set_xlabel('x (kpc)')
    slc.set_ylabel('y (kpc)')
    slc.set_width((2.0*kpc,1.0*kpc))
    slc.set_background_color(plotvar)
    # plot = slc.plots[plotvar]
    # colorbar=plot.cb
    # slc._setup_plots()
    # colorbar.set_ticks([1e-28])
    # colorbar.set_ticklabels(['$10^{-28}$']
    slc.set_log(plotvar, True)
    plot = slc.plots[plotvar]
    plot.figure = fig
    plot.axes = grid[3].axes
    plot.cax = grid.cbar_axes[j]
    slc._setup_plots()
  if (j==2):
   # varmin = 1.e5
   # varmax = 1.e8
    varmin = 1.e3
    varmax = 1.e6
    slc = yt.SlicePlot(ds, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc],fontsize=fsize)
    slc.set_zlim(plotvar, varmin, varmax)
   # slc.set_cmap(field=plotvar, cmap='afmhot')
    slc.set_cmap(field=plotvar, cmap='coolwarm')
    slc.set_xlabel('x (kpc)')
    slc.set_ylabel('y (kpc)')
    slc.set_width((2.0*kpc,1.0*kpc))
    slc.set_background_color(plotvar)
    # plot = slc.plots[plotvar]
    # colorbar=plot.cb
    # slc._setup_plots()
    # colorbar.set_ticks([1e-28])
    # colorbar.set_ticklabels(['$10^{-28}$']
    slc.set_log(plotvar, True)
    plot = slc.plots[plotvar]
    plot.figure = fig
    plot.axes = grid[4].axes
    plot.cax = grid.cbar_axes[j]
    slc._setup_plots()

    slc = yt.SlicePlot(ds3, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc],fontsize=fsize)
    slc.set_zlim(plotvar, varmin, varmax)
   # slc.set_cmap(field=plotvar, cmap='afmhot')
    slc.set_cmap(field=plotvar, cmap='coolwarm')
    slc.set_xlabel('x (kpc)')
    slc.set_ylabel('y (kpc)')
    slc.set_width((2.0*kpc,1.0*kpc))
    slc.set_background_color(plotvar)
    # plot = slc.plots[plotvar]
    # colorbar=plot.cb
    # slc._setup_plots()
    # colorbar.set_ticks([1e-28])
    # colorbar.set_ticklabels(['$10^{-28}$']
    slc.set_log(plotvar, True)
    plot = slc.plots[plotvar]
    plot.figure = fig
    plot.axes = grid[5].axes
    plot.cax = grid.cbar_axes[j]
    slc._setup_plots()
  if (j==3):
    varmin = 1.e-28
    varmax = 1.e-25
    slc = yt.SlicePlot(ds2, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc],fontsize=fsize)
    slc.set_zlim(plotvar, varmin, varmax)
    slc.set_cmap(field=plotvar, cmap='plasma')
    slc.set_xlabel('x (kpc)')
    slc.set_ylabel('y (kpc)')
    slc.set_width((2.0*kpc, 1.0*kpc))
    slc.set_background_color(plotvar)
    slc.set_log(plotvar, True)
    plot = slc.plots[plotvar]
    plot.figure = fig
    plot.axes = grid[6].axes
    plot.cax = grid.cbar_axes[j]
    slc._setup_plots()
    slc = yt.SlicePlot(ds4, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc],fontsize=fsize)
    slc.set_zlim(plotvar, varmin, varmax)
    slc.set_cmap(field=plotvar, cmap='plasma')
    slc.set_xlabel('x (kpc)')
    slc.set_ylabel('y (kpc)')
    slc.set_width((2.0*kpc, 1.0*kpc))
    slc.set_background_color(plotvar)
    slc.set_log(plotvar, True)
    plot = slc.plots[plotvar]
    plot.figure = fig
    plot.axes = grid[7].axes
    plot.cax = grid.cbar_axes[j]
    slc._setup_plots()
  if (j==4):
    varmin = 1.e-28
    varmax = 1.e-25
    slc = yt.SlicePlot(ds, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc],fontsize=fsize)
    slc.set_zlim(plotvar, varmin, varmax)
    slc.set_xlabel('x (kpc)')
    slc.set_ylabel('y (kpc)')
    slc.set_cmap(field=plotvar, cmap='plasma')
    slc.set_background_color(plotvar)
    slc.annotate_contour('d',ncont = 1,clim=[2e-24,2e-24],take_log = True, label = False,
                   plot_args={"colors": "white", "linewidths": 1.0})
    slc.annotate_contour('d',ncont = 1,clim=[2e-23,2e-23],take_log = True, label = False,
                   plot_args={"colors": "cyan", "linewidths": 1.0})
    slc.annotate_contour('d',ncont = 1,clim=[2e-22,2e-22],take_log = True, label = False,
                  plot_args={"colors": "green", "linewidths": 1.0})
    slc.set_width((2.0*kpc, 1.0*kpc))
    slc.set_log(plotvar, True)

    plot = slc.plots[plotvar]
    plot.figure = fig
    plot.axes = grid[8].axes
    plot.cax = grid.cbar_axes[j]
    slc._setup_plots()
    slc = yt.SlicePlot(ds3, 'z', plotvar,origin='native', center=[0.5*kpc, 0*kpc, 0*kpc],fontsize=fsize)
    slc.set_zlim(plotvar, varmin, varmax)
    slc.set_xlabel('x (kpc)')
    slc.set_ylabel('y (kpc)')
    slc.set_cmap(field=plotvar, cmap='plasma')
    slc.set_background_color(plotvar)
    slc.annotate_contour('d',ncont = 1,clim=[2e-24,2e-24],take_log = True, label = False,
                   plot_args={"colors": "white", "linewidths": 1.0})
    slc.annotate_contour('d',ncont = 1,clim=[2e-23,2e-23],take_log = True, label = False,
                   plot_args={"colors": "cyan", "linewidths": 1.0})
    slc.annotate_contour('d',ncont = 1,clim=[2e-22,2e-22],take_log = True, label = False,
                  plot_args={"colors": "green", "linewidths": 1.0})
    slc.set_width((2.0*kpc, 1.0*kpc))
    slc.set_log(plotvar, True)

    plot = slc.plots[plotvar]

    plot.figure = fig
    plot.axes = grid[9].axes
    plot.cax = grid.cbar_axes[j]
    slc._setup_plots()
    #fig.set_size_inches(14,14)
    fig.set_size_inches(14,16)


# For each plotted field, force the SlicePlot to redraw itself onto the AxesGrid
# axes.
#for i, field in enumerate(fields):
#    plot = slc.plots[field]
#    plot.figure = fig
#    plot.axes = grid[i].axes
#    plot.cax = grid.cbar_axes[i]

# Finally, redraw the plot on the AxesGrid axes.

plt.savefig('multiplot_damping_nodamping_withTemperature_lowerB.pdf')
