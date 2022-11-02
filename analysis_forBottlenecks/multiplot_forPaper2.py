import yt
import numpy as np
from yt.visualization.base_plot_types import get_multi_plot
import matplotlib.colorbar as cb
from matplotlib.colors import LogNorm
from yt.units import dimensions
from yt.units import erg
from yt.units.yt_array import YTQuantity


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



fn = "../cr.out1.00030.athdf" # dataset to load
fn2 = "../cr.out2.00030.athdf" # dataset to load
orient = 'horizontal'

ds = yt.load(fn) # load data
ds2 = yt.load(fn2) # load data

ds.add_field(('gas', u'd'), function = _d, units="g/cm**3",display_name=r"Density", dimensions=dimensions.density)

ds2.add_field(("gas","ionAlfven"), function=ionAlfven,display_name=r"$v_{A}^{ion}$",units="cm/s")

ds2.add_field(("gas","heating"), function=heating,display_name="Collisionless Loss",units="erg/cm**3/s")

ds.add_field(('gas', u'ecr'), function = _ecr, units="g/(cm*s**2)",display_name=r"E$_{cr}$", dimensions=dimensions.pressure)
ds.add_field(("gas","collisional"), function=collisional,display_name="Collisional Loss",units="erg/cm**3/s")


width = (1.0, 'unitary')
res = [256, 256]

# There's a lot in here:
#   From this we get a containing figure, a list-of-lists of axes into which we
#   can place plots, and some axes that we'll put colorbars.
# We feed it:
#   Number of plots on the x-axis, number of plots on the y-axis, and how we
#   want our colorbars oriented.  (This governs where they will go, too.
#   bw is the base-width in inches, but 4 is about right for most cases.
fig, axes, colorbars = get_multi_plot(4, 1, colorbar=orient, bw = 4)
#
slc = ds.proj(0,0.5)
proj = ds2.slice(0,0.5)

#slc = yt.SlicePlot(ds, 'z', fields=["ecr","collisional"])
#proj = yt.SlicePlot(ds2, 'z', fields=["ionAlfven","heating"])

#slc_frb = slc.data_source.to_frb((1.0,2.0), 512)
slc_frb = slc.to_frb(width,res)
#proj_frb = proj.data_source.to_frb(1.0, 512)
proj_frb = proj.to_frb(width, res)

ecr_axes = [axes[0][0], axes[0][0]]
va_axes = [axes[0][1], axes[0][1]]
collision_axes = [axes[0][2], axes[0][2]]
collisionless_axes = [axes[0][3], axes[0][3]]

for eax, vax, collax, collessax in zip(ecr_axes, va_axes, collision_axes,collisionless_axes) :

    eax.xaxis.set_visible(False)
    eax.yaxis.set_visible(False)
    vax.xaxis.set_visible(False)
    vax.yaxis.set_visible(False)
    collax.xaxis.set_visible(False)
    collax.yaxis.set_visible(False)
    collessax.xaxis.set_visible(False)
    collessax.yaxis.set_visible(False)

# Converting our Fixed Resolution Buffers to numpy arrays so that matplotlib
# can render them
slc_frb1 = FixedResolutionBuffer(ds.proj(0,"ecr"),(0.0,0.0,1.0,2.0),(512,512))
slc_frb2 = FixedResolutionBuffer(ds.proj(0,"collisional"),(0.0,0.0,1.0,2.0),(512,512))
proj_frb1 = FixedResolutionBuffer(ds2.proj(0,"ionAlfven"),(0.0,0.0,1.0,2.0),(512,512))
proj_frb2 = FixedResolutionBuffer(ds2.proj(0,"collisionless"),(0.0,0.0,1.0,2.0),(512,512))
slc_ecr = np.array(slc_frb1)
#slc_ecr = np.array(slc_frb['ecr'])
#proj_va = np.array(proj_frb['ionAlfven'])
proj_va = np.array(proj_frb1)
#slc_coll = np.array(slc_frb['collisional'])
slc_coll = np.array(slc_frb2)
#proj_colless = np.array(proj_frb['heating'])
proj_colless = np.array(proj_frb2)

plots = [ecr_axes[0].imshow(slc_ecr,origin='lower', norm=LogNorm()),
         va_axes[0].imshow(proj_va,origin='lower',norm=LogNorm()),
         collision_axes[0].imshow(slc_coll,origin='lower', norm=LogNorm()),
         collisionless_axes[0].imshow(proj_colless,origin='lower', norm=LogNorm())]

plots[0].set_clim((1.0e-14,1.0e-11))
plots[0].set_cmap("viridis")
plots[1].set_clim((1.0e5,1.0e8))
plots[1].set_cmap("afmhot")
plots[2].set_clim((1.0e-28,1.0e-25))
plots[2].set_cmap("plasma")
plots[3].set_clim((1.0e-28,1.0e-25))
plots[3].set_cmap("plasma")

titles=[r'$\mathrm{Cosmic Ray Energy Density}\ (\mathrm{erg\ cm^{-3}})$',
        r'$\mathrm{v_{A}^{ion}}\ (\mathrm{cm/s})$',
        r'$\mathrm{Collisional Energy Loss}\ (\mathrm{erg\ cm^{-3} s^{-1}})$',
        r'$\mathrm{Collisionless Energy Loss}\ (\mathrm{erg\ cm^{-3} s^{-1}})$']

for p, cax, t in zip(plots[0:4:1], colorbars, titles):
    cbar = fig.colorbar(p, cax=cax, orientation=orient)
    cbar.set_label(t)

# And now we're done!
#fig.savefig("%s_3x2" % ds)
fig.savefig("multiplot_damping2.png")
