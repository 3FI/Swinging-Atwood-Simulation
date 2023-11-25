import OneSwinging_2D_Explicit
import TwoSwinging_2D_Explicit

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np

tf = 100
dt = 0.01
maxU = 25
U = np.linspace(1.1,maxU+1,200)
x0 = np.array([ 3, 0, np.pi/2, 0 ])

DET = np.zeros(len(U))
t = np.zeros(len(U))

for i,u in enumerate(U):
    r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0)
    DET[i] = max(abs(dET))
    t[i] = len(r)*dt
plt.plot(U,DET,label="Maximal Fractional Energy Loss")
plt.plot(U,t/max(t),label="Fractional Time Length Of The Simulation")

plt.xlabel("Mu")
plt.ylabel("Maximal dE (%)")

plt.savefig('../results/dE_By_MU-1S-Explicit.png')
plt.show()


plt.plot(U[t/max(t)>0.95],DET[t/max(t)>0.95],label="Maximal Fractional Energy Loss")
plt.xlabel("Mu")
plt.ylabel("Maximal dE (%)")
plt.savefig('../results/dE_By_MU-Masked-1S-Explicit.png')
plt.show()


x0 = np.array([ 3, 0, np.pi/2, 0 , np.pi/8 , 0 ])
r0M = 3

DET = np.zeros(len(U))
t = np.zeros(len(U))

for i,u in enumerate(U):
    r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0,r0M=r0M)
    DET[i] = min(100,max(abs(dET)))
    t[i] = len(r)*dt
plt.plot(U,DET,label="Clamped Maximal Fractional Energy Loss")
plt.plot(U,t/max(t),label="Fractional Time Length Of The Simulation")

plt.xlabel("Mu")
plt.ylabel("Maximal dE (%)")

plt.savefig('../results/dE_By_MU-2S-Explicit.png')
plt.show()


plt.plot(U[t/max(t)>0.95],DET[t/max(t)>0.95],label="Maximal Fractional Energy Loss")
plt.xlabel("Mu")
plt.ylabel("Maximal dE (%)")
plt.savefig('../results/dE_By_MU-Masked-2S-Explicit.png')
plt.show()