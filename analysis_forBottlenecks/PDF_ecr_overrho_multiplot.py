import yt
from yt.units.yt_array import YTQuantity
from yt.units import dimensions
import numpy as np
from numpy import *
import matplotlib as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import h5py
from yt.fields import interpolated_fields
from yt.fields.field_detector import \
    FieldDetector
from yt.utilities.linear_interpolators import \
    BilinearFieldInterpolator

yt.enable_parallelism()

pUnit = YTQuantity(1, 'cm**2/s**2')

# conversion factors
denstocgs = 6.85e-27
edenstocgs = 6.54e-11
Myr = 1.
kpc = 1.


def _ecr(field, data):
  return data['Ec']*edenstocgs*YTQuantity(1,"erg/cm**3")

def _d(field, data):
  return data['density']*denstocgs

def _ecrrho(field,data):
  return data['Ec']*edenstocgs*YTQuantity(1,"erg/cm**3")/(data['density']*denstocgs)

def gammaRay_lum(field,data): #was defined as negative in Athena++, so needs a minus sign
          return 0.3333*0.7*0.728*((3.0856e21)**(3.0))*(1.022e-15)*(data['d']/1.67e-24) * (data['Ec']*edenstocgs) * YTQuantity(1,"erg/s/g")

def collision(field,data): #was defined as negative in Athena++, so needs a minus sign
          return 0.3333*0.7*0.728*((3.0856e21)**(3.0))*(1.022e-15)*(data['d']/1.67e-24) * (data['Ec']*edenstocgs) * YTQuantity(1,"erg/s/g")



yt.add_field(('gas', u'd'), function = _d, units="g/cm**3",display_name=r"Density", dimensions=dimensions.density)
yt.add_field(('gas', u'ecr'), function = _ecr, units="g/(cm*s**2)",display_name=r"E$_{cr}$", dimensions=dimensions.pressure)
yt.add_field(('gas', u'ecrrho'), function = _ecrrho, units="erg/g",display_name=r"E$_{cr}/\rho$", dimensions=dimensions.pressure)
yt.add_field(("gas","gammaRay_lum"), function=gammaRay_lum,display_name="Gamma-Ray",units="erg/cm**3/s")
yt.add_field(("gas","collision"), function=collision,display_name="Collisional Loss",units="erg/cm**3/s")

profiles = []
labels = []
plot_specs = []

dslist = []

ds = yt.load('noIonAlfven_noDamping/noClumpsNearBoundary/L5/alpha1_5/cr.out1.00040.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.25"])
       

time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profiles.append(yt.create_profile(ad, "d", ["ecrrho"],weight_field=None, fractional=True, accumulation=False,extrema={'ecrrho': (1e-3, 1e0),
                                      'd': (1e-25, 1e-22)}))

# Add labels
labels.append(r"")
plot_specs.append(dict(linestyle="-",color='k', linewidth=3))
dslist.append(ds)


ds = yt.load('noIonAlfven_noDamping/noClumpsNearBoundary/L5/alpha1_5/cr.out1.00040.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.25"])
       
time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profiles.append(yt.create_profile(ad, "d", ["collision"],weight_field=None, fractional=True, accumulation=False,extrema={'collision': (1e-3, 1e0),
                                      'd': (1e-25, 1e-22)}))

# Add labels
labels.append(r"")
plot_specs.append(dict(linestyle="-",color='k', linewidth=3))

dslist.append(ds)


ds = yt.load('ionAlfven_Damping/noClumpsNearBoundary/L5_alpha1_5/cr.out1.00040.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.25"])
       
time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profiles.append(yt.create_profile(ad, "d", ["ecrrho"],weight_field=None, fractional=True, accumulation=False,extrema={'ecrrho': (1e-3, 1e0),
                                      'd': (1e-25, 1e-22)}))

# Add labels
labels.append(r"")
plot_specs.append(dict(linestyle="-",color='r', linewidth=3))

dslist.append(ds)


ds = yt.load('ionAlfven_Damping/noClumpsNearBoundary/L5_alpha1_5/cr.out1.00040.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.25"])
       
time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profiles.append(yt.create_profile(ad, "d", ["collision"],weight_field=None, fractional=True, accumulation=False,extrema={'collision': (1e-3, 1e0),
                                      'd': (1e-25, 1e-22)}))

# Add labels
labels.append(r"")
plot_specs.append(dict(linestyle="-",color='r', linewidth=3))
dslist.append(ds)

################33
# Create the profile plot from the list of profiles.
plot = yt.ProfilePlot.from_profiles(profiles, labels=labels, plot_specs=plot_specs)

ax = plot.axes['gas','collision']
labels = ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels()
labels += [ax.title, ax.xaxis.label, ax.yaxis.label,
           ax.xaxis.get_offset_text(), ax.yaxis.get_offset_text()]
for label in labels:
    label.set_fontsize(22)
ax = plot.axes['gas','ecrrho']
labels = ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels()
labels += [ax.title, ax.xaxis.label, ax.yaxis.label,
           ax.xaxis.get_offset_text(), ax.yaxis.get_offset_text()]
for label in labels:
    label.set_fontsize(22)
#plot.set_xlim(1.E-27,1.E-21)
#plot.set_xlim(1.E-27,1.E-21)
#plot.set_ylim("gammaRay_lum",1.E-3,1.E0)
#plot.set_ylabel("ecr", r"E$_{CR}$ PDF")
#plot.annotate_title(r"Lognormal  Distributions")
#plot.annotate_title(r"L = 5, $\alpha$ = 1.5")
plot.save(suffix='pdf')
        
