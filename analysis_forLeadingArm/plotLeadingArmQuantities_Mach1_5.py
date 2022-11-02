import matplotlib.pyplot as plt
import numpy as np

# Athena history dump for level=0 domain=0 volume=1.350000e+05
#   [1]=time      [2]=dt         [3]=mass       [4]=total E    [5]=x1 Mom.    [6]=x2 Mom.    [7]=x3 Mom.    [8]=x1-KE      [9]=x2-KE      [10]=x3-KE     [11]=scalar 0  [12]=m13  [13]=m110  [14]=mT2  [15]=Mx13  [16]=Erad  [17]=<c>  [18]=<c^2>  [19]=<c * E>  [20]=<c * x1>  [21]=<c * Vx>  [22]=<c * Vy>  [23]=<c * Vz>  [24]=<(c * Vx)^2>  [25]=<(c * Vy)^2>  [26]=<(c * Vz)^2>  [27]=dye entropy  [28]=x_shift  [29]=v_flow
#

time, m13, m110, mT2, x_shift, vflow = np.loadtxt("../out/id0/cloud.hst",skiprows=3,usecols=(0,11,12,13,27,28),unpack=True)

m13_new = np.array(m13)/m13[0]
m110_new = np.array(m110)/m110[0]
mT2_new = np.array(mT2)/mT2[0]

plt.plot(time,m13,'k-',label='m13')
plt.plot(time,m110,'k--',label='m110')
plt.plot(time, mT2,'k-.',label='mT2')
plt.title("Leading Arm Test")
plt.legend()
plt.xlabel("Time (code units)")
plt.ylabel("Cold Gas Mass (code units)")
plt.savefig("ColdMass_LA.png")
plt.close()



plt.plot(time,m13_new,'k-',label="m13")
plt.plot(time,m110_new,'k--',label="m110")
plt.plot(time,mT2_new,'k-.',label="mT2")
plt.title("Leading Arm Test")
plt.ylim(0.0,8.0)
plt.legend()
plt.xlabel("Time (code units)")
plt.ylabel("Cold Gas Mass Increase")
plt.savefig("ColdMassIncrease_LA.png")
plt.close()

plt.plot(time,x_shift,'k-')
plt.title("Leading Arm Test")
plt.xlabel("Time (code units)")
plt.ylabel("Distance")
plt.savefig("Distance_LA.png")
plt.close()

plt.plot(time,vflow,'k-')
plt.title("Leading Arm Test")
plt.xlabel("Time (code units)")
plt.ylabel("Flow Velocity")
plt.savefig("vflow_LA.png")
plt.close()
