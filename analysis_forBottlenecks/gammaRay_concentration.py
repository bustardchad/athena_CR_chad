import yt
from yt.units.yt_array import YTQuantity
from yt.units import dimensions
import numpy as np
from numpy import *
import matplotlib.pyplot as plt
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
prestocgs = 6.54e-11
Myr = 1.
kpc = 1.
temptocgs = prestocgs/denstocgs


def _ecr(field, data):
  return data['Ec']*edenstocgs*YTQuantity(1,"erg/cm**3")

def _d(field, data):
  return data['density']*denstocgs

def gammaRay_lum(field,data): #was defined as negative in Athena++, so needs a minus sign
          return 0.3333*0.7*0.728*((3.0856e21)**(3.0))*(1.022e-15)*(data['d']/1.67e-24) * (data['Ec']*edenstocgs) * YTQuantity(1,"erg/s/g")
def _tmp(field, data):
 return data['temperature']*temptocgs


yt.add_field(('gas', u'tmp'), function = _tmp, units="K",display_name=r"Temperature", dimensions=dimensions.temperature)
yt.add_field(('gas', u'd'), function = _d, units="g/cm**3",display_name=r"Density", dimensions=dimensions.density)
yt.add_field(('gas', u'ecr'), function = _ecr, units="g/(cm*s**2)",display_name=r"E$_{cr}$", dimensions=dimensions.pressure)
yt.add_field(("gas","gammaRay_lum"), function=gammaRay_lum,display_name="Gamma-Ray",units="erg/cm**3/s")

profiles = []
labels = []
plot_specs = []

ds = yt.load('noIonAlfven_noDamping/noClumpsNearBoundary/L5/alpha1_5/cr.out1.00040.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
ad = dd.cut_region(['(obj["x"] < 0.5) & (obj["x"] > 0.1)'])
       
time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profile_ecr = yt.create_profile(ad, "tmp", ["gammaRay_lum"],weight_field=None, fractional=True, accumulation=False,extrema={'gammaRay_lum': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))


print(profile_ecr.x)
va_L5_alpha15 = profile_ecr["gas", "gammaRay_lum"]
print(profile_ecr["gas", "gammaRay_lum"]) 
cray = np.array(profile_ecr["gas", "gammaRay_lum"])
profile_vol = yt.create_profile(ad, "tmp", ["cell_volume"],weight_field=None, fractional=True, accumulation=False,extrema={'cell_volume': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))

print(profile_vol.x)
print(profile_vol["gas", "cell_volume"]) 
vol = np.array(profile_vol["gas", "cell_volume"])


concentration_va_L5_alpha15 = cray/vol
print(concentration_va_L5_alpha15)



ds = yt.load('ionAlfven_Damping/noClumpsNearBoundary/L5_alpha1_5/cr.out1.00040.athdf')

dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
ad = dd.cut_region(['(obj["x"] < 0.5) & (obj["x"] > 0.1)'])
       
time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profile_ecr = yt.create_profile(ad, "tmp", ["gammaRay_lum"],weight_field=None, fractional=True, accumulation=False,extrema={'gammaRay_lum': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))


print(profile_ecr.x)
vaion_L5_alpha15 = profile_ecr["gas", "gammaRay_lum"]
print(profile_ecr["gas", "gammaRay_lum"]) 
cray = np.array(profile_ecr["gas", "gammaRay_lum"])
profile_vol = yt.create_profile(ad, "tmp", ["cell_volume"],weight_field=None, fractional=True, accumulation=False,extrema={'cell_volume': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))

print(profile_vol.x)
print(profile_vol["gas", "cell_volume"]) 
vol = np.array(profile_vol["gas", "cell_volume"])


concentration_vaion_L5_alpha15 = cray/vol
print(concentration_vaion_L5_alpha15)


ds = yt.load('noIonAlfven_noDamping/noClumpsNearBoundary/L2/alpha1_5/res2/cr.out1.00040.athdf')

dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
ad = dd.cut_region(['(obj["x"] < 0.5) & (obj["x"] > 0.1)'])
       
time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profile_ecr = yt.create_profile(ad, "tmp", ["gammaRay_lum"],weight_field=None, fractional=True, accumulation=False,extrema={'gammaRay_lum': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))


print(profile_ecr.x)
va_L2_alpha15 = profile_ecr["gas", "gammaRay_lum"]
print(profile_ecr["gas", "gammaRay_lum"]) 
cray = np.array(profile_ecr["gas", "gammaRay_lum"])
profile_vol = yt.create_profile(ad, "tmp", ["cell_volume"],weight_field=None, fractional=True, accumulation=False,extrema={'cell_volume': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))

print(profile_vol.x)
print(profile_vol["gas", "cell_volume"]) 
vol = np.array(profile_vol["gas", "cell_volume"])


concentration_va_L2_alpha15 = cray/vol
print(concentration_va_L2_alpha15)

ds = yt.load('ionAlfven_Damping/noClumpsNearBoundary/L2_alpha1_5/res2/cr.out1.00040.athdf')

dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
ad = dd.cut_region(['(obj["x"] < 0.5) & (obj["x"] > 0.1)'])
       
time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profile_ecr = yt.create_profile(ad, "tmp", ["gammaRay_lum"],weight_field=None, fractional=True, accumulation=False,extrema={'gammaRay_lum': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))


print(profile_ecr.x)
vaion_L2_alpha15 = profile_ecr["gas", "gammaRay_lum"]
print(profile_ecr["gas", "gammaRay_lum"]) 
cray = np.array(profile_ecr["gas", "gammaRay_lum"])
profile_vol = yt.create_profile(ad, "tmp", ["cell_volume"],weight_field=None, fractional=True, accumulation=False,extrema={'cell_volume': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))

print(profile_vol.x)
print(profile_vol["gas", "cell_volume"]) 
vol = np.array(profile_vol["gas", "cell_volume"])


concentration_vaion_L2_alpha15 = cray/vol
print(concentration_vaion_L2_alpha15)


ds = yt.load('noIonAlfven_noDamping/noClumpsNearBoundary/L5/alpha1_5_lowerB/cr.out1.00040.athdf')

dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
ad = dd.cut_region(['(obj["x"] < 0.5) & (obj["x"] > 0.1)'])
       
time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profile_ecr = yt.create_profile(ad, "tmp", ["gammaRay_lum"],weight_field=None, fractional=True, accumulation=False,extrema={'gammaRay_lum': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))


print(profile_ecr.x)
va_L5_alpha15_lowerB = profile_ecr["gas", "gammaRay_lum"]
print(profile_ecr["gas", "gammaRay_lum"]) 
cray = np.array(profile_ecr["gas", "gammaRay_lum"])
profile_vol = yt.create_profile(ad, "tmp", ["cell_volume"],weight_field=None, fractional=True, accumulation=False,extrema={'cell_volume': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))

print(profile_vol.x)
print(profile_vol["gas", "cell_volume"]) 
vol = np.array(profile_vol["gas", "cell_volume"])


concentration_va_L5_alpha15_lowerB = cray/vol
print(concentration_va_L5_alpha15_lowerB)


ds = yt.load('ionAlfven_Damping/noClumpsNearBoundary/L5_alpha1_5/lowerB/cr.out1.00040.athdf')

dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
ad = dd.cut_region(['(obj["x"] < 0.5) & (obj["x"] > 0.1)'])
       
time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profile_ecr = yt.create_profile(ad, "tmp", ["gammaRay_lum"],weight_field=None, fractional=True, accumulation=False,extrema={'gammaRay_lum': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))


print(profile_ecr.x)
vaion_L5_alpha15_lowerB = profile_ecr["gas", "gammaRay_lum"]
print(profile_ecr["gas", "gammaRay_lum"]) 
cray = np.array(profile_ecr["gas", "gammaRay_lum"])
profile_vol = yt.create_profile(ad, "tmp", ["cell_volume"],weight_field=None, fractional=True, accumulation=False,extrema={'cell_volume': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))

print(profile_vol.x)
print(profile_vol["gas", "cell_volume"]) 
vol = np.array(profile_vol["gas", "cell_volume"])


concentration_vaion_L5_alpha15_lowerB = cray/vol
print(concentration_vaion_L5_alpha15_lowerB)


ds = yt.load('noIonAlfven_noDamping/noClumpsNearBoundary/L2/alpha1_5_lowerB/res2/cr.out1.00040.athdf')

dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
ad = dd.cut_region(['(obj["x"] < 0.5) & (obj["x"] > 0.1)'])
       
time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profile_ecr = yt.create_profile(ad, "tmp", ["gammaRay_lum"],weight_field=None, fractional=True, accumulation=False,extrema={'gammaRay_lum': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))


print(profile_ecr.x)
va_L2_alpha15_lowerB = profile_ecr["gas", "gammaRay_lum"]
print(profile_ecr["gas", "gammaRay_lum"]) 
cray = np.array(profile_ecr["gas", "gammaRay_lum"])
profile_vol = yt.create_profile(ad, "tmp", ["cell_volume"],weight_field=None, fractional=True, accumulation=False,extrema={'cell_volume': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))

print(profile_vol.x)
print(profile_vol["gas", "cell_volume"]) 
vol = np.array(profile_vol["gas", "cell_volume"])


concentration_va_L2_alpha15_lowerB = cray/vol
print(concentration_va_L2_alpha15_lowerB)

ds = yt.load('ionAlfven_Damping/noClumpsNearBoundary/L2_alpha1_5_lowerB/res2/cr.out1.00040.athdf')

dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
ad = dd.cut_region(['(obj["x"] < 0.5) & (obj["x"] > 0.1)'])
       
time = str(ds.current_time.in_units('Myr'))
time = (time[:3]) if len(time) > 3 else time
t = "t = {} Myrs".format(str(time))
profile_ecr = yt.create_profile(ad, "tmp", ["gammaRay_lum"],weight_field=None, fractional=True, accumulation=False,extrema={'gammaRay_lum': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))


print(profile_ecr.x)
vaion_L2_alpha15_lowerB = profile_ecr["gas", "gammaRay_lum"]
print(profile_ecr["gas", "gammaRay_lum"]) 
cray = np.array(profile_ecr["gas", "gammaRay_lum"])
profile_vol = yt.create_profile(ad, "tmp", ["cell_volume"],weight_field=None, fractional=True, accumulation=False,extrema={'cell_volume': (1e-3, 1e0),
                                      'tmp': (1e3, 5e5)},n_bins=(32))

print(profile_vol.x)
print(profile_vol["gas", "cell_volume"]) 
vol = np.array(profile_vol["gas", "cell_volume"])


concentration_vaion_L2_alpha15_lowerB = cray/vol
print(concentration_vaion_L2_alpha15_lowerB)






plt.semilogx(np.array(profile_ecr.x),concentration_va_L5_alpha15,'k-',label=r'L = 5, $\alpha = 1.5$, B = 5 $\mu$G')
plt.semilogx(np.array(profile_ecr.x),concentration_vaion_L5_alpha15,'k--')
plt.semilogx(np.array(profile_ecr.x),concentration_va_L2_alpha15,'g-',label=r'L = 2, $\alpha = 1.5$, B = 5 $\mu$G')
plt.semilogx(np.array(profile_ecr.x),concentration_vaion_L2_alpha15,'g--',)
plt.semilogx(np.array(profile_ecr.x),concentration_va_L5_alpha15_lowerB,'b-',label=r'L = 5, $\alpha = 1.5$, B = 1 $\mu$G')
plt.semilogx(np.array(profile_ecr.x),concentration_vaion_L5_alpha15_lowerB,'b--')
plt.semilogx(np.array(profile_ecr.x),concentration_va_L2_alpha15_lowerB,'m-',label=r'L = 2, $\alpha = 1.5$, B = 1 $\mu$G')
plt.semilogx(np.array(profile_ecr.x),concentration_vaion_L2_alpha15_lowerB,'m--')
plt.ylim(0.0,5.0)
#plt.legend()
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel("Gamma-Ray Concentration",fontsize=20)
plt.xlabel("Temperature (K)",fontsize=20)
plt.tight_layout()
plt.savefig("gammaRay_concentration.pdf")
plt.close()

plt.loglog(np.array(profile_ecr.x),va_L5_alpha15,'k-',label=r'L = 5, $\alpha = 1.5$, B = 5 $\mu$G')
plt.loglog(np.array(profile_ecr.x),vaion_L5_alpha15,'k--')
plt.loglog(np.array(profile_ecr.x),va_L2_alpha15,'g-',label=r'L = 2, $\alpha = 1.5$, B = 5 $\mu$G')
plt.loglog(np.array(profile_ecr.x),vaion_L2_alpha15,'g--',)
plt.loglog(np.array(profile_ecr.x),va_L5_alpha15_lowerB,'b-',label=r'L = 5, $\alpha = 1.5$, B = 1 $\mu$G')
plt.loglog(np.array(profile_ecr.x),vaion_L5_alpha15_lowerB,'b--')
plt.loglog(np.array(profile_ecr.x),va_L2_alpha15_lowerB,'m-',label=r'L = 2, $\alpha = 1.5$, B = 1 $\mu$G')
plt.loglog(np.array(profile_ecr.x),vaion_L2_alpha15_lowerB,'m--')
plt.ylim(5.E-3,5.E-1)
#plt.legend()
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.ylabel("Gamma-Ray Emission PDF",fontsize=20)
plt.xlabel("Temperature (K)",fontsize=20)
plt.tight_layout()
plt.savefig("gammaRay_PDF.pdf")
plt.close()
