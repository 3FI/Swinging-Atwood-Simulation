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
t = np.zeros(len(U))
for i,u in enumerate(U):
    r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0)
    t[i] = len(r)*dt
plt.plot(U,t)
plt.savefig('../results/Time_Until_Singular-1S-Explicit.png')
plt.show()


tf = 100
dt = 0.01
maxU = 25
U = np.linspace(1.1,maxU+1,200)
x0 = np.array([ 3, 0, np.pi/2, 0 , np.pi/8 , 0 ])
r0M = 3
t = np.zeros(len(U))
for i,u in enumerate(U):
    r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0,r0M=r0M)
    t[i] = len(r)*dt
plt.plot(U,t)
plt.savefig('../results/Time_Until_Singular-2S-Explicit.png')
plt.show()