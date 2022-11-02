import matplotlib.pyplot as plt
import yt
from yt.units import erg, pc
from yt.units.yt_array import YTQuantity
import numpy as np
# conversion factors
edenstocgs = 6.54e-11
denstocgs = 6.85e-27
Myr = 1.
kpc = 1.
 
yt.enable_parallelism()

#Streaming energy losses
def crscale(field,data):
          return np.abs(data['user_out_var8'])*YTQuantity(1,"cm")
 
 
yt.add_field(("gas","crscale"), function=crscale,display_name="CR Scale Length along B",units="cm") 
 
times = []
totalCollArr = []  # array holding the collisional energy loss
totalDampArr = []  # array holding the streaming loss
totalDampDensArr = []  # array holding the streaming loss
totalCRLossArr = [] # array holding the total (collisional + streaming) loss
ts = yt.DatasetSeries('../cr.out2*',parallel=False)
i = 0

for ds in ts.piter():
  ad = ds.all_data()
  time = ds.current_time.v/Myr                       # store the time (with units removed)
  times.append(time)
  ray = ds.ortho_ray(0,(0,0))
  print(ray['crscale'])
