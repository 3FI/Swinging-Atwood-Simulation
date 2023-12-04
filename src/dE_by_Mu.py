import OneSwinging_2D_Explicit
import TwoSwinging_2D_Explicit

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np

#This is the plot of the energy error by the value of mu and it's correlation with the lifetime by mu

tf = 500 #The final time of the simulation
dt = 0.01 #The time interval

#The range of mu we are covering.
maxU = 25
U = np.linspace(1.1,maxU+1,200)

#The initial parameters of the simulation
x0 = np.array([ 3, 0, np.pi/2, 0 ])

#The array which will hold the error in energy
DET = np.zeros(len(U))
#The array which will hold the lifetime
t = np.zeros(len(U))

for i,u in enumerate(U):
    #Does the integration
    r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0)

    #Stock the values for the current mu
    DET[i] = max(abs(dET))
    t[i] = len(r)*dt

#Plot the data
plt.plot(U,DET,label="Maximal Fractional Energy Loss")
plt.plot(U,5*t/max(t),label="5 * Fractional Time Length Of The Simulation")
plt.xlabel("Mu")
plt.ylabel("Maximal dE (%)")
plt.legend()
plt.savefig('../results/dE_By_MU-1S-Explicit.png')
plt.show()

#How correlated are the two line for the singular (lifetime < final time)
print("The correlation between the data is : " + str(np.corrcoef(t[t/max(t)<1]/max(t),DET[t/max(t)<1]/max(DET))))

#Plot the energy correlation for the non-singular only by using the lifetime as a mask
#This didn't lead to any concrete results.
plt.plot(U[t/max(t)==1],DET[t/max(t)==1],label="Maximal Fractional Energy Loss")
plt.xlabel("Mu")
plt.ylabel("Maximal dE (%)")
plt.legend()
plt.savefig('../results/dE_By_MU-Masked-1S-Explicit.png')
plt.show()



#Now let's do it for the 2 swinging mass case

#The final time is low because we want to compare with high time to show it's meaningless
tf = 100
#The new initial parameters
x0 = np.array([ 3, 0, np.pi/2, 0 , np.pi/8 , 0 ])
#The initial R of M
r0M = 3

#Refresh the arrays
DET = np.zeros(len(U))
t = np.zeros(len(U))

for i,u in enumerate(U):
    r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0,r0M=r0M)
    DET[i] = min(3,max(abs(dET)))
    t[i] = len(r)*dt

plt.plot(U,DET,label="Maximal Fractional Energy Loss")
plt.plot(U,1*t/max(t),label="1 * Fractional Time Length Of The Simulation")
plt.xlabel("Mu")
plt.ylabel("Maximal dE (%)")
plt.legend()
plt.savefig('../results/dE_By_MU-2S-Explicit-100s.png')
plt.show()


#We now change the value of the final time to compare
tf = 1000
DET = np.zeros(len(U))
t = np.zeros(len(U))

for i,u in enumerate(U):
    r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0,r0M=r0M)
    DET[i] = min(3,max(abs(dET)))
    t[i] = len(r)*dt

plt.plot(U,DET,label="Maximal Fractional Energy Loss")
plt.plot(U,1*t/max(t),label="1 * Fractional Time Length Of The Simulation")
plt.xlabel("Mu")
plt.ylabel("Maximal dE (%)")
plt.legend()
plt.savefig('../results/dE_By_MU-2S-Explicit-1000s.png')
plt.show()
