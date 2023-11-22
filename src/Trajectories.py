import OneSwinging_2D_Explicit
import TwoSwinging_2D_Explicit

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np


tf = 100
dt = 0.01
x0 = np.array([ 3, 0, np.pi/2, 0  ])
r0M = 3

fig,ax = plt.subplots(2,4,figsize=(40,10))
plt.subplots_adjust(wspace = 0.3, hspace = 0.3)
i=0
for u in [2,3,5,6,16,19,21,24]:    
    r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0)
    X = r*np.sin(theta); Y = -r*np.cos(theta)
    sc = ax[int(i/4)][i%4].scatter(X,Y,s=0.5,c=dET,alpha=0.8,)
    ax[int(i/4)][i%4].set_aspect('equal','datalim')
    ax[int(i/4)][i%4].set_title('mu = ' + str(u))

    (plt.colorbar(sc,ax=ax[int(i/4)][i%4])).formatter.set_powerlimits((0, 0))

    i+=1

plt.suptitle('Trajectories for diverse mass ratio (Theta_0 = pi/2  &  No initial velocity)')
plt.savefig('../results/Trajectories-1S-Explicit.png')
plt.show()




x0 = np.array([ 3, 0, np.pi/2, 0 , np.pi/8 , 0 ])

fig,ax = plt.subplots(2,4,figsize=(40,10))
plt.subplots_adjust(wspace = 0.2)
i=0
for u in [2,3,5,6,16,19,21,24]:    
    r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0,r0M=r0M)
    X = r*np.sin(theta_m); Y = -r*np.cos(theta_m)
    sc = ax[int(i/4)][i%4].scatter(X,Y,s=0.5,c=dET,alpha=0.8,)
    ax[int(i/4)][i%4].set_aspect('equal','datalim')
    ax[int(i/4)][i%4].set_title('mu = ' + str(u))

    (plt.colorbar(sc,ax=ax[int(i/4)][i%4])).formatter.set_powerlimits((0, 0))

    i+=1

plt.suptitle('Trajectories for diverse mass ratio (Theta_0 = pi/2  &  No initial velocity)')
plt.savefig('../results/Trajectories-2S-Explicit.png')
plt.show()


fig,ax = plt.subplots(2,4,figsize=(40,10))
plt.subplots_adjust(wspace = 0.2)
i=0
for u in [2,3,5,6,16,19,21,24]:    
    r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0,r0M=r0M)
    X_M = (r0M + x0[0]-r)*np.sin(theta_M); Y_M = -(r0M + x0[0]-r)*np.cos(theta_M)
    sc = ax[int(i/4)][i%4].scatter(X,Y,s=0.5,c=dET,alpha=0.8,)
    ax[int(i/4)][i%4].set_aspect('equal','datalim')
    ax[int(i/4)][i%4].set_title('mu = ' + str(u))

    (plt.colorbar(sc,ax=ax[int(i/4)][i%4])).formatter.set_powerlimits((0, 0))

    i+=1

plt.suptitle('Trajectories for diverse mass ratio (Theta_0 = pi/2  &  No initial velocity)')
plt.savefig('../results/Trajectories-2S-Explicit_M.png')
plt.show()