import OneSwinging_2D_Explicit
import TwoSwinging_2D_Explicit

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np
import scipy

#This is the same as dE_by_dt but in the periodic case instead of quasiperiodic

tf = np.pi
max_nStep = 10**5
STEPS = np.linspace(10**1,max_nStep,20) 
DT = tf/STEPS

#The only difference with the other file : the mu is now in the periodic case
u=1.665
x0 = np.array([ 3, 0, np.pi/2, 0 ])

DET = np.zeros(len(DT))

for i,dt in enumerate(DT):
    r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0)
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
plt.ylabel("Final dE")
plt.title("Energy Error at the End Dependant on the Number of Steps Taken")
plt.legend()
plt.savefig('../results/dE_By_dt-Periodic-1S-Explicit.png')
plt.show()