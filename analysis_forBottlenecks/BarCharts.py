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

Cold1arr = []
Cold2arr = []
Cold3arr = []
Hotarr = []
volCold1arr = []

ds = yt.load('noIonAlfven_noDamping/noClumpsNearBoundary/L5/alpha1_5/cr.out1.00040.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
adcold1 = dd.cut_region(['(obj["tmp"] < 5e3) & (obj["x"] < 0.5)'])
adcold2 = dd.cut_region(['(5e3 < obj["tmp"]) & (obj["tmp"] < 8e3) & (obj["x"] < 0.5)'])
adcold3 = dd.cut_region(['(8e3 < obj["tmp"]) & (obj["tmp"] < 2e4) & (obj["x"] < 0.5)'])
adhot = dd.cut_region(['(obj["tmp"] > 2e4) & (obj["x"] < 0.5)'])
       
vol = ad.quantities.total_quantity("cell_volume")
volCold1 = adcold1.quantities.total_quantity("cell_volume")/vol
volCold2 = adcold2.quantities.total_quantity("cell_volume")/vol
volCold3 = adcold3.quantities.total_quantity("cell_volume")/vol
volHot = adhot.quantities.total_quantity("cell_volume")/vol
q = ad.quantities.total_quantity("qparam")
qCold1 = adcold1.quantities.total_quantity("qparam")/q/volCold1/1.
qCold2 = adcold2.quantities.total_quantity("qparam")/q/volCold2/1.
qCold3 = adcold3.quantities.total_quantity("qparam")/q/volCold3/1.
qHot = adhot.quantities.total_quantity("qparam")/q/volHot/1.

Cold1arr.append(qCold1)
Cold2arr.append(qCold2)
Cold3arr.append(qCold3)
Hotarr.append(qHot)


ds = yt.load('ionAlfven_Damping/noClumpsNearBoundary/L5_alpha1_5/cr.out1.00040.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
adcold1 = dd.cut_region(['(obj["tmp"] < 5e3) & (obj["x"] < 0.5)'])
adcold2 = dd.cut_region(['(5e3 < obj["tmp"]) & (obj["tmp"] < 8e3) & (obj["x"] < 0.5)'])
adcold3 = dd.cut_region(['(8e3 < obj["tmp"]) & (obj["tmp"] < 2e4) & (obj["x"] < 0.5)'])
adhot = dd.cut_region(['(obj["tmp"] > 2e4) & (obj["x"] < 0.5)'])
       
vol = ad.quantities.total_quantity("cell_volume")
volCold1 = adcold1.quantities.total_quantity("cell_volume")/vol
volCold2 = adcold2.quantities.total_quantity("cell_volume")/vol
volCold3 = adcold3.quantities.total_quantity("cell_volume")/vol
volHot = adhot.quantities.total_quantity("cell_volume")/vol
q = ad.quantities.total_quantity("qparam")
qCold1 = adcold1.quantities.total_quantity("qparam")/q/volCold1/1.
qCold2 = adcold2.quantities.total_quantity("qparam")/q/volCold2/1.
qCold3 = adcold3.quantities.total_quantity("qparam")/q/volCold3/1.
qHot = adhot.quantities.total_quantity("qparam")/q/volHot/1.

Cold1arr.append(qCold1)
Cold2arr.append(qCold2)
Cold3arr.append(qCold3)
Hotarr.append(qHot)



ds = yt.load('noIonAlfven_noDamping/noClumpsNearBoundary/L2/alpha1_5/res2/cr.out1.00040.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
adcold1 = dd.cut_region(['(obj["tmp"] < 5e3) & (obj["x"] < 0.5)'])
adcold2 = dd.cut_region(['(5e3 < obj["tmp"]) & (obj["tmp"] < 8e3) & (obj["x"] < 0.5)'])
adcold3 = dd.cut_region(['(8e3 < obj["tmp"]) & (obj["tmp"] < 2e4) & (obj["x"] < 0.5)'])
adhot = dd.cut_region(['(obj["tmp"] > 2e4) & (obj["x"] < 0.5)'])
       
vol = ad.quantities.total_quantity("cell_volume")
volCold1 = adcold1.quantities.total_quantity("cell_volume")/vol
volCold2 = adcold2.quantities.total_quantity("cell_volume")/vol
volCold3 = adcold3.quantities.total_quantity("cell_volume")/vol
volHot = adhot.quantities.total_quantity("cell_volume")/vol
q = ad.quantities.total_quantity("qparam")
qCold1 = adcold1.quantities.total_quantity("qparam")/q/volCold1/1.
qCold2 = adcold2.quantities.total_quantity("qparam")/q/volCold2/1.
qCold3 = adcold3.quantities.total_quantity("qparam")/q/volCold3/1.
qHot = adhot.quantities.total_quantity("qparam")/q/volHot/1.

Cold1arr.append(qCold1)
Cold2arr.append(qCold2)
Cold3arr.append(qCold3)
Hotarr.append(qHot)


        
ds = yt.load('ionAlfven_Damping/noClumpsNearBoundary/L2_alpha1_5/res2/cr.out1.00040.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
adcold1 = dd.cut_region(['(obj["tmp"] < 5e3) & (obj["x"] < 0.5)'])
adcold2 = dd.cut_region(['(5e3 < obj["tmp"]) & (obj["tmp"] < 8e3) & (obj["x"] < 0.5)'])
adcold3 = dd.cut_region(['(8e3 < obj["tmp"]) & (obj["tmp"] < 2e4) & (obj["x"] < 0.5)'])
adhot = dd.cut_region(['(obj["tmp"] > 2e4) & (obj["x"] < 0.5)'])
       
vol = ad.quantities.total_quantity("cell_volume")
volCold1 = adcold1.quantities.total_quantity("cell_volume")/vol
volCold2 = adcold2.quantities.total_quantity("cell_volume")/vol
volCold3 = adcold3.quantities.total_quantity("cell_volume")/vol
volHot = adhot.quantities.total_quantity("cell_volume")/vol
q = ad.quantities.total_quantity("qparam")
qCold1 = adcold1.quantities.total_quantity("qparam")/q/volCold1/1.
qCold2 = adcold2.quantities.total_quantity("qparam")/q/volCold2/1.
qCold3 = adcold3.quantities.total_quantity("qparam")/q/volCold3/1.
qHot = adhot.quantities.total_quantity("qparam")/q/volHot/1.

Cold1arr.append(qCold1)
Cold2arr.append(qCold2)
Cold3arr.append(qCold3)
Hotarr.append(qHot)


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
ds = yt.load('noIonAlfven_noDamping/noClumpsNearBoundary/L2/alpha1_5_lowerB/res2/cr.out1.00040.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
adcold1 = dd.cut_region(['(obj["tmp"] < 5e3) & (obj["x"] < 0.5)'])
adcold2 = dd.cut_region(['(5e3 < obj["tmp"]) & (obj["tmp"] < 8e3) & (obj["x"] < 0.5)'])
adcold3 = dd.cut_region(['(8e3 < obj["tmp"]) & (obj["tmp"] < 2e4) & (obj["x"] < 0.5)'])
adhot = dd.cut_region(['(obj["tmp"] > 2e4) & (obj["x"] < 0.5)'])
       
vol = ad.quantities.total_quantity("cell_volume")
volCold1 = adcold1.quantities.total_quantity("cell_volume")/vol
volCold2 = adcold2.quantities.total_quantity("cell_volume")/vol
volCold3 = adcold3.quantities.total_quantity("cell_volume")/vol
volHot = adhot.quantities.total_quantity("cell_volume")/vol
q = ad.quantities.total_quantity("qparam")
qCold1 = adcold1.quantities.total_quantity("qparam")/q/volCold1/1.
qCold2 = adcold2.quantities.total_quantity("qparam")/q/volCold2/1.
qCold3 = adcold3.quantities.total_quantity("qparam")/q/volCold3/1.
qHot = adhot.quantities.total_quantity("qparam")/q/volHot/1.

Cold1arr.append(qCold1)
Cold2arr.append(qCold2)
Cold3arr.append(qCold3)
Hotarr.append(qHot)




ds = yt.load('ionAlfven_Damping/noClumpsNearBoundary/L2_alpha1_5_lowerB/res2/cr.out1.00040.athdf')
dd = ds.all_data()
ad = dd.cut_region(["obj['x'] < 0.5"])
adcold1 = dd.cut_region(['(obj["tmp"] < 5e3) & (obj["x"] < 0.5)'])
adcold2 = dd.cut_region(['(5e3 < obj["tmp"]) & (obj["tmp"] < 8e3) & (obj["x"] < 0.5)'])
adcold3 = dd.cut_region(['(8e3 < obj["tmp"]) & (obj["tmp"] < 2e4) & (obj["x"] < 0.5)'])
adhot = dd.cut_region(['(obj["tmp"] > 2e4) & (obj["x"] < 0.5)'])
       
vol = ad.quantities.total_quantity("cell_volume")
volCold1 = adcold1.quantities.total_quantity("cell_volume")/vol
volCold2 = adcold2.quantities.total_quantity("cell_volume")/vol
volCold3 = adcold3.quantities.total_quantity("cell_volume")/vol
volHot = adhot.quantities.total_quantity("cell_volume")/vol
q = ad.quantities.total_quantity("qparam")
qCold1 = adcold1.quantities.total_quantity("qparam")/q/volCold1/1.
qCold2 = adcold2.quantities.total_quantity("qparam")/q/volCold2/1.
qCold3 = adcold3.quantities.total_quantity("qparam")/q/volCold3/1.
qHot = adhot.quantities.total_quantity("qparam")/q/volHot/1.

Cold1arr.append(qCold1)
Cold2arr.append(qCold2)
Cold3arr.append(qCold3)
Hotarr.append(qHot)





# Plotting the bar chart
width = 0.15       # the width of the bars: can also be len(x) sequence

ind = np.arange(6)
ind2 = [0.75,2.75,4.75]
p1 = plt.bar(ind, Cold1arr, width)
p2 = plt.bar(ind + width, Cold2arr,width) #, width, bottom = np.array(Cold1arr))
p3 = plt.bar(ind+width+width, Cold3arr,width) #, width, bottom = np.array(Cold1arr) + np.array(Cold2arr))
p4 = plt.bar(ind+width+width+width, Hotarr,width) #, width, bottom = np.array(Cold1arr) + np.array(Cold2arr) + np.array(Cold3arr))
#p1 = plt.bar(ind, Cold1arr, width)
#p2 = plt.bar(ind, Cold2arr, width, bottom = np.array(Cold1arr))
#p3 = plt.bar(ind, Cold3arr, width, bottom = np.array(Cold1arr) + np.array(Cold2arr))
#p4 = plt.bar(ind, Hotarr, width, bottom = np.array(Cold1arr) + np.array(Cold2arr) + np.array(Cold3arr))
plt.title(r'E$_{CR}$ Concentration')
#plt.xticks(ind, (r'L = 5, $\alpha = 1.5$' '\n' 'B = 5 $\mu$G',r'L = 5, $\alpha = 1.5$' '\n' 'B = 5 $\mu$G', r'L = 2, $\alpha = 1.5$' '\n' 'B = 5 $\mu$G', r'L = 2, $\alpha = 1.5$' '\n' 'B = 5 $\mu$G', r'L = 2, $\alpha = 1.5$' '\n' 'B = 5 $\mu$G', r'L = 2, $\alpha = 1.5$' '\n' 'B = 1 $\mu$G'),fontsize=10)
plt.xticks(ind2, (r'L = 5, $\alpha = 1.5$' '\n' 'B = 5 $\mu$G', r'L = 2, $\alpha = 1.5$' '\n' 'B = 5 $\mu$G', r'L = 2, $\alpha = 1.5$' '\n' 'B = 1 $\mu$G'),fontsize=12)
plt.yticks(np.arange(0, 2.0, 0.2))
plt.legend((p1[0], p2[0], p3[0], p4[0]), (r'T < 5 x 10$^{3}$ K', r'5 x 10$^{3}$ K < T < 8 x 10$^{3}$ K', r'8 x 10$^{3}$ K < T < 2 x 10$^{4}$ K', r'T  > 2 x 10$^{4}$ K'),ncol=2)

plt.tight_layout()

plt.savefig('BarPlot2.pdf')
