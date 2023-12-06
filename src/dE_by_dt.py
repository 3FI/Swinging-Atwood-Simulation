import OneSwinging_2D_Explicit
import TwoSwinging_2D_Explicit

import OneSwinging_2D_Symplectic
import TwoSwinging_2D_Symplectic

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np
import scipy



#This is the plot of the energy error at the final time by the number of steps used in the integration

#The time of the quasi-period
tf = np.pi

#The array of number of steps we want to test for
max_nStep = 10**5
STEPS = np.linspace(10**1,max_nStep,20) 
#Converting it in an array of time interval so our integrator can use it
DT = tf/STEPS

#The mu and initial parameters
u=3
x0 = np.array([ 3, 0, np.pi/2, 0 ])

#The array to stock the energy error of each system
DET = np.zeros(len(DT))

#Do the integration for each dt
for i,dt in enumerate(DT):
    r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0)
    #We stock the absolute value of the error at the final time to look at the conservation over the entire quasiperiod
    DET[i] = abs(dET[-1])
plt.plot(STEPS,DET,label="Maximal Fractional Energy Loss")

#Plot comparison lines
n_vec = 10**np.linspace(1.5,3.7,10)
plt.plot(n_vec, 1/n_vec,":", label=r'$1/N$')
plt.plot(n_vec, 1/n_vec**2,":", label=r'$1/N^2$')
plt.plot(n_vec, 1/n_vec**3,":", label=r'$1/N^3$')
plt.plot(n_vec, 1/n_vec**4,":", label=r'$1/N^4$')
n_vec = 10**np.linspace(4.5,5,10)
plt.plot(n_vec, 1e-17*n_vec**0.5,":", label=r'$N^{1/2}$')

plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.xlabel("Number of Step")
plt.ylabel("Final dE")
plt.title("Energy Error at the End Dependant on the Number of Steps Taken")

plt.savefig('../results/dE_By_dt-1S-Explicit.png')
plt.show()


# THE 2 SWINGING CASE


#Initial parameters for 2 swinging mass
x0 = np.array([ 3, 0, np.pi/2, 0 , np.pi/8 , 0 ])
#The initial R of M
r0M = 3

#Refresh the array
DET = np.zeros(len(DT))

for i,dt in enumerate(DT):
    r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0,r0M=r0M)
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
plt.legend()
plt.xlabel("Number of Step")
plt.ylabel("Final dE")
plt.title("Energy Error at the End Dependant on the Number of Steps Taken")

plt.savefig('../results/dE_By_dt-2S-Explicit.png')
plt.show()

########################################################################
# repeat the same plots for the leapfrog integrator
########################################################################

#The time of the quasi-period
tf = np.pi

#The array of number of steps we want to test for
max_nStep = 10**5
STEPS = np.linspace(10**1,max_nStep,20) 
#Converting it in an array of time interval so our integrator can use it
DT = tf/STEPS

#The mu and initial parameters
u=3
x0 = np.array([ 3, 0, np.pi/2, 0 ])

#The array to stock the energy error of each system
DET = np.zeros(len(DT))

#Do the integration for each dt
for i,dt in enumerate(DT):
    r,pr,theta,ptheta,dET = OneSwinging_2D_Symplectic.integrate(u=u, tf=tf, dt=dt, x0=x0)
    #We stock the absolute value of the error at the final time to look at the conservation over the entire quasiperiod
    DET[i] = abs(dET[-1])
plt.plot(STEPS,DET,label="Maximal Fractional Energy Loss")

#Plot comparison lines
n_vec = 10**np.linspace(1.5,3.7,10)
plt.plot(n_vec, 1/n_vec,":", label=r'$1/N$')
plt.plot(n_vec, 1/n_vec**2,":", label=r'$1/N^2$')
plt.plot(n_vec, 1/n_vec**3,":", label=r'$1/N^3$')
plt.plot(n_vec, 1/n_vec**4,":", label=r'$1/N^4$')
n_vec = 10**np.linspace(4.5,5,10)
plt.plot(n_vec, 1e-17*n_vec**0.5,":", label=r'$N^{1/2}$')

plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.xlabel("Number of Step")
plt.ylabel("Final dE")
plt.title("Energy Error at the End Dependant on the Number of Steps Taken")

plt.savefig('../results/dE_By_dt-1S-Symplectic.png')
plt.show()


# THE 2 SWINGING CASE


#Initial parameters for 2 swinging mass
x0 = np.array([ 3, 0, np.pi/2, 0 , np.pi/8 , 0 ])
#The initial R of M
r0M = 3

#Refresh the array
DET = np.zeros(len(DT))

for i,dt in enumerate(DT):
    r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Symplectic.integrate(u=u, tf=tf, dt=dt, x0=x0,r0M=r0M)
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
plt.legend()
plt.xlabel("Number of Step")
plt.ylabel("Final dE")
plt.title("Energy Error at the End Dependant on the Number of Steps Taken")

plt.savefig('../results/dE_By_dt-2S-Symplectic.png')
plt.show()