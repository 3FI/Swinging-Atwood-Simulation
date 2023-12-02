import OneSwinging_2D_Explicit
import TwoSwinging_2D_Explicit

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np

fig,ax = plt.subplots(1,2,figsize=(8,4))

u = 3
m = 1
tf = 10
dt = 0.01
x0_M = -10
r0_M = 3
t = np.arange(0,tf,dt)

x0 = np.array([ 3, 0, np.pi/4, 0  ])
ropeLength = r0_M + x0[0]

r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrateMass(m=m, M=u*m, tf=tf, dt=dt, x0=x0)

ax[0].plot(t,r,label="m=1, M=3")
ax[1].plot(t,theta,label="m=1, M=3")

m = 2
r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrateMass(m=m, M=u*m, tf=tf, dt=dt, x0=x0)

ax[0].plot(t,r,label="m=2, M=6")
ax[1].plot(t,theta,label="m=2, M=6")


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


r1,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrateMass(m=m, M=u*m, tf=tf, dt=dt, x0=x0)
r2,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrateMass(m=m, M=u*m, tf=tf, dt=dt, x0=x0)

print(np.mean(r1-r2))