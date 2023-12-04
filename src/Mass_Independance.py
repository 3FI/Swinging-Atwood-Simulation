import OneSwinging_2D_Explicit
import TwoSwinging_2D_Explicit

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np


#This is the plot of the independance of the motion under multiplication of both mass

fig,ax = plt.subplots(1,2,figsize=(8,4))

#Set the mass ratio which we will always keep constant 
u = 3
#The final time of the integration
tf = 10
#The time interval
dt = 0.01

#The time array of the integration (for drawing purpose)
t = np.arange(0,tf,dt)

#The initial parameters
x0 = np.array([ 3, 0, np.pi/4, 0  ])

#We integrate by using the variation of integrate which takes both mass instead of a singlu mu
m = 1
r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrateMass(m=m, M=u*m, tf=tf, dt=dt, x0=x0)
ax[0].plot(t,r,label="m=1, M=3")
ax[1].plot(t,theta,label="m=1, M=3")

#Do the same but with the m=2
m = 2
r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrateMass(m=m, M=u*m, tf=tf, dt=dt, x0=x0)
ax[0].plot(t,r,label="m=2, M=6")
ax[1].plot(t,theta,label="m=2, M=6")

#Do the same but with the m=3
m = 3
r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrateMass(m=m, M=u*m, tf=tf, dt=dt, x0=x0)
ax[0].plot(t,r,label="m=3, M=9")
ax[1].plot(t,theta,label="m=3, M=9")

ax[0].set(xlabel='t', ylabel='r')
ax[1].set(xlabel='t', ylabel='$\\theta$')
plt.legend()
plt.suptitle("Independance of the motion under mass multiplication")
plt.savefig('../results/Mass-Independance-1S-Explicit.png')
plt.show()

#Verify that they are indeed the same by differentiating them
r1,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrateMass(m=m, M=u*m, tf=tf, dt=dt, x0=x0)
r2,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrateMass(m=m, M=u*m, tf=tf, dt=dt, x0=x0)
print(np.mean(r1-r2))