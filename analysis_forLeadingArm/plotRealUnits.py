import matplotlib.pyplot as plt
import numpy as np

# Athena history dump for level=0 domain=0 volume=1.350000e+05
#   [1]=time      [2]=dt         [3]=mass       [4]=total E    [5]=x1 Mom.    [6]=x2 Mom.    [7]=x3 Mom.    [8]=x1-KE      [9]=x2-KE      [10]=x3-KE     [11]=scalar 0  [12]=m13  [13]=m110  [14]=mT2  [15]=Mx13  [16]=Erad  [17]=<c>  [18]=<c^2>  [19]=<c * E>  [20]=<c * x1>  [21]=<c * Vx>  [22]=<c * Vy>  [23]=<c * Vz>  [24]=<(c * Vx)^2>  [25]=<(c * Vy)^2>  [26]=<(c * Vz)^2>  [27]=dye entropy  [28]=x_shift  [29]=v_flow
#

time, m13, m110, mT2, x_shift, vflow = np.loadtxt("../out/id0/cloud.hst",skiprows=3,usecols=(0,11,12,13,27,28),unpack=True)


m13_new = np.array(m13)/m13[0]
m110_new = np.array(m110)/m110[0]
mT2_new = np.array(mT2)/mT2[0]


# take inputs of tcc, rcl, M, nhalo
nhalo = input("Enter halo density in units of 1e-5 cm^-3: ")
stri = str(nhalo).zfill(3)
tccinput = input("Enter t_cc in code units: ")
Minput = input("Enter Mach number in code units: ")
rclinput = input("Enter cloud radius in code units: ")

timeRatio = 4.5/tccinput

rclpc = 2.0*(tccinput/4.5)*Minput/(1e3*nhalo) # r_cl

tccMyrs = 0.075*rclpc/Minput # crushing times in Myrs
tMyrs = 0.075*rclpc/(tccinput*Minput) # time unit to multiply to time array
dkpc = 0.001*rclpc/rclinput # distance unit in kpc to multiply to distance array

timeMyrs = time*tMyrs
distance = x_shift*dkpc
vel = vflow*(distance/timeMyrs)*3.0856e16/3.155e13 # in km/s

num = "ColdMassIncrease_LA_{}.png".format(stri)
plt.plot(timeMyrs,m13_new,'k-',label='m13')
plt.plot(timeMyrs,m110_new,'k--',label='m110')
plt.plot(timeMyrs, mT2_new,'k-.',label='mT2')
plt.title("Leading Arm Test")
plt.legend()
plt.xlabel("Time (Myrs)")
plt.ylabel("Cold Gas Mass Increase")
plt.ylim(0,20)
plt.savefig(num,format='png')
plt.close()


num = "Distance_LA_{}.png".format(stri)
plt.plot(timeMyrs,distance,'k-')
plt.title("Leading Arm Test")
plt.xlabel("Time (Myrs)")
plt.ylabel("Distance (kpc)")
plt.savefig(num,format='png')
plt.close()

num = "VFlow_LA_{}.png".format(stri)
plt.plot(timeMyrs,vel,'k-')
plt.title("Leading Arm Test")
plt.xlabel("Time (Myrs)")
plt.ylabel("Flow Velocity (km/s)")
plt.savefig(num,format='png')
plt.close()
