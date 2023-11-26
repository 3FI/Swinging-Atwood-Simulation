import OneSwinging_2D_Explicit
import TwoSwinging_2D_Explicit

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np
import scipy

tf = 2*np.pi
max_nStep = 10**5
STEPS = np.linspace(10**1,max_nStep,20) 
DT = tf/STEPS
u=3
x0 = np.array([ 3, 0, np.pi/2, 0 ])

DET = np.zeros(len(DT))

for i,dt in enumerate(DT):
    r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0)
    E = pr**2 / (2*(u+1)) + ptheta**2 / (2*1*r**2) + u*scipy.constants.g*r - 1*scipy.constants.g*r*np.cos(theta)
    dET = E - E[0] 
    DET[i] = abs(dET[-1])
plt.plot(STEPS,DET,label="Maximal Fractional Energy Loss")

n_vec = 10**np.linspace(1.5,3.7,10)
plt.plot(n_vec, 1/n_vec,":", label=r'$1/N$')
plt.plot(n_vec, 1/n_vec**2,":", label=r'$1/N^2$')
plt.plot(n_vec, 1/n_vec**3,":", label=r'$1/N^3$')
plt.plot(n_vec, 1/n_vec**4,":", label=r'$1/N^4$')
n_vec = 10**np.linspace(4.5,5,10)
plt.plot(n_vec, 1e-17*n_vec**0.5,":", label=r'$N^{1/2}$')

plt.xscale("log")
plt.yscale("log")

plt.xlabel("Number of Step")
plt.ylabel("Maximal dE")

plt.savefig('../results/dE_By_dt-1S-Explicit.png')
plt.show()





x0 = np.array([ 3, 0, np.pi/2, 0 , np.pi/8 , 0 ])
r0M = 3
ropeLength = r0M + x0[0]

DET = np.zeros(len(DT))

for i,dt in enumerate(DT):
    r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0,r0M=r0M)
    E = pr**2 / (2*(u+1)) + ptheta_M**2 / (2*u*(ropeLength-r)**2) + ptheta_m**2 / (2*1*r**2) + u*scipy.constants.g*(ropeLength - (ropeLength-r)*np.cos(theta_M)) + 1*scipy.constants.g*(ropeLength-r*np.cos(theta_m))
    dET = E - E[0] 
    DET[i] = abs(dET[-1])
plt.plot(STEPS,DET,label="Clamped Maximal Fractional Energy Loss")

n_vec = 10**np.linspace(1.5,3.7,10)
plt.plot(n_vec, 1/n_vec,":", label=r'$1/N$')
plt.plot(n_vec, 1/n_vec**2,":", label=r'$1/N^2$')
plt.plot(n_vec, 1/n_vec**3,":", label=r'$1/N^3$')
plt.plot(n_vec, 1/n_vec**4,":", label=r'$1/N^4$')
n_vec = 10**np.linspace(4.5,5,10)
plt.plot(n_vec, 1e-17*n_vec**0.5,":", label=r'$N^{1/2}$')

plt.xscale("log")
plt.yscale("log")

plt.xlabel("Number of Step")
plt.ylabel("Maximal dE")

plt.savefig('../results/dE_By_dt-2S-Explicit.png')
plt.show()