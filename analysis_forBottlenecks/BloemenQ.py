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
prestocgs = 6.54e-11
#temptocgs = 0.6*1.67e-24*prestocgs/(denstocgs*1.38e-16)
temptocgs = prestocgs/(denstocgs)

def _ecr(field, data):
  return data['Ec']*edenstocgs*YTQuantity(1,"erg/cm**3")

def _d(field, data):
  return data['density']*denstocgs

def gammaRay_lum(field,data): #was defined as negative in Athena++, so needs a minus sign
          return 0.3333*0.7*0.728*((3.0856e21)**(3.0))*(1.022e-15)*(data['d']/1.67e-24) * (data['Ec']*edenstocgs) * YTQuantity(1,"erg/s/g")

def collision(field,data): #was defined as negative in Athena++, so needs a minus sign
          return 0.3333*0.7*0.728*((3.0856e21)**(3.0))*(1.022e-15)*(data['d']/1.67e-24) * (data['Ec']*edenstocgs) * YTQuantity(1,"erg/s/g")

def _tmp(field, data):
 return data['temperature']*temptocgs

def qparam(field,data): #was defined as negative in Athena++, so needs a minus sign
          return 0.3333*0.7*0.728*((3.0856e21)**(3.0))*(1.022e-15) * (data['Ec']*edenstocgs) * YTQuantity(1,"erg/s/g")

yt.add_field(('gas', u'tmp'), function = _tmp, units="K",display_name=r"Temperature", dimensions=dimensions.temperature)
yt.add_field(('gas', u'd'), function = _d, units="g/cm**3",display_name=r"Density", dimensions=dimensions.density)
yt.add_field(('gas', u'ecr'), function = _ecr, units="g/(cm*s**2)",display_name=r"E$_{cr}$", dimensions=dimensions.pressure)
yt.add_field(("gas","gammaRay_lum"), function=gammaRay_lum,display_name="Gamma-Ray",units="erg/cm**3/s")
yt.add_field(("gas","qparam"), function=qparam,display_name="Gamma-Ray Emissivity",units="erg/g/s")
yt.add_field(("gas","collision"), function=collision,display_name="Collisional Loss",units="erg/cm**3/s")

profiles = []
labels = []
plot_specs = []

dslist = []

tcut = 5e3

qarr = []
qH2arr = []
volFracarr = []

ds = yt.load('noIonAlfven_noDamping/noClumpsNearBoundary/L5/alpha1_5/cr.out1.00010.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.25"])
adcold = dd.cut_region(['(obj["tmp"] < 1e4) & (obj["x"] < 0.25)'])
       
q = ad.quantities.total_quantity("qparam")
qH2 = adcold.quantities.total_quantity("qparam")
vol1 = ad.quantities.total_quantity("cell_volume")
vol2 = adcold.quantities.total_quantity("cell_volume")


qarr.append(q)
qH2arr.append(qH2/q)
volFracarr.append(vol2/vol1)

ds = yt.load('ionAlfven_Damping/noClumpsNearBoundary/L5_alpha1_5/cr.out1.00010.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.25"])
       
adcold = dd.cut_region(['(obj["tmp"] < 1e4) & (obj["x"] < 0.25)'])
q = ad.quantities.total_quantity("qparam")
qH2 = adcold.quantities.total_quantity("qparam")
vol1 = ad.quantities.total_quantity("cell_volume")
vol2 = adcold.quantities.total_quantity("cell_volume")
volFracarr.append(vol2/vol1)
qarr.append(q)
qH2arr.append(qH2/q)


ds = yt.load('noIonAlfven_noDamping/noClumpsNearBoundary/L2/alpha1_5/res2/cr.out1.00010.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.25"])

adcold = dd.cut_region(['(obj["tmp"] < 1e4) & (obj["x"] < 0.25)'])
q = ad.quantities.total_quantity("qparam")
qH2 = adcold.quantities.total_quantity("qparam")
qarr.append(q)
qH2arr.append(qH2/q)
vol1 = ad.quantities.total_quantity("cell_volume")
vol2 = adcold.quantities.total_quantity("cell_volume")
volFracarr.append(vol2/vol1)

        
ds = yt.load('ionAlfven_Damping/noClumpsNearBoundary/L2_alpha1_5/res2/cr.out1.00010.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.25"])

adcold = dd.cut_region(['(obj["tmp"] < 1e4) & (obj["x"] < 0.25)'])
q = ad.quantities.total_quantity("qparam")
qH2 = adcold.quantities.total_quantity("qparam")
qarr.append(q)
qH2arr.append(qH2/q)
vol1 = ad.quantities.total_quantity("cell_volume")
vol2 = adcold.quantities.total_quantity("cell_volume")
volFracarr.append(vol2/vol1)

"""
ds = yt.load('noIonAlfven_noDamping/noClumpsNearBoundary/L5/alpha1_5_lowerB/cr.out1.00040.athdf')
dd = ds.all_data()
adcold = dd.cut_region(["obj['tmp'] < 1e4"])
       
q = dd.quantities.total_quantity("qparam")
qH2 = adcold.quantities.total_quantity("qparam")

qarr.append(q)
qH2arr.append(qH2/q)



ds = yt.load('ionAlfven_Damping/noClumpsNearBoundary/L5_alpha1_5/lowerB/cr.out1.00040.athdf')
dd = ds.all_data()
adcold = dd.cut_region(["obj['tmp'] < 1e4"])
       
q = dd.quantities.total_quantity("qparam")
qH2 = adcold.quantities.total_quantity("qparam")

qarr.append(q)
qH2arr.append(qH2/q)

"""
ds = yt.load('noIonAlfven_noDamping/noClumpsNearBoundary/L2/alpha1_5_lowerB/res2/cr.out1.00010.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.25"])

adcold = dd.cut_region(['(obj["tmp"] < 1e4) & (obj["x"] < 0.25)'])
q = dd.quantities.total_quantity("qparam")
qH2 = adcold.quantities.total_quantity("qparam")
qarr.append(q)
qH2arr.append(qH2/q)
vol1 = ad.quantities.total_quantity("cell_volume")
vol2 = adcold.quantities.total_quantity("cell_volume")
volFracarr.append(vol2/vol1)


ds = yt.load('ionAlfven_Damping/noClumpsNearBoundary/L2_alpha1_5_lowerB/res2/cr.out1.00010.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.25"])

adcold = dd.cut_region(['(obj["tmp"] < 1e4) & (obj["x"] < 0.25)'])
q = ad.quantities.total_quantity("qparam")
qH2 = adcold.quantities.total_quantity("qparam")
qarr.append(q)
qH2arr.append(qH2/q)
vol1 = ad.quantities.total_quantity("cell_volume")
vol2 = adcold.quantities.total_quantity("cell_volume")
volFracarr.append(vol2/vol1)


print(qH2arr)
print(volFracarr)
