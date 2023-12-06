import OneSwinging_2D_Explicit
import TwoSwinging_2D_Explicit

import OneSwinging_2D_Symplectic
import TwoSwinging_2D_Symplectic

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np

#This is the plot of the trajectories for multiple selected mu 

#The total time over which we carry the integration
tf = 100
#The time interval used in the integration
dt = 0.01
#The initial parameters
x0 = np.array([ 3, 0, np.pi/2, 0  ])

fig,ax = plt.subplots(2,4,figsize=(20,10))
plt.subplots_adjust(wspace = 0.3, hspace = 0.3)
i=0
#For each selected mu, carry the integration and do it's plot
for u in [2,3,5,6,16,19,21,24]:    
    #Integrate for the current mu
    r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0)
    #Get back to cartesian coordinates
    X = r*np.sin(theta); Y = -r*np.cos(theta)
    #Draw the scatter of the trajectory
    sc = ax[int(i/4)][i%4].scatter(X,Y,s=0.5,c=dET,alpha=0.8,)

    #Set the title
    ax[int(i/4)][i%4].set_aspect('equal','datalim')
    ax[int(i/4)][i%4].set_title('$\\mu$ = ' + str(u))

    #Make the colorbar for the energy
    (plt.colorbar(sc,ax=ax[int(i/4)][i%4])).formatter.set_powerlimits((0, 0))

    i+=1

plt.suptitle('Trajectories for diverse mass ratio ($\\theta_0$ = $\\pi$/2  &  No initial velocity)')
plt.savefig('../results/Trajectories-1S-Explicit.png')
plt.show()

# do the same simulation with the leapfrog integrator
fig,ax = plt.subplots(2,4,figsize=(20,10))
plt.subplots_adjust(wspace = 0.3, hspace = 0.3)
i=0
#For each selected mu, carry the integration and do it's plot
for u in [2,3,5,6,16,19,21,24]:    
    #Integrate for the current mu
    r,pr,theta,ptheta,dET = OneSwinging_2D_Symplectic.integrate(u=u, tf=tf, dt=dt, x0=x0)
    #Get back to cartesian coordinates
    X = r*np.sin(theta); Y = -r*np.cos(theta)
    #Draw the scatter of the trajectory
    sc = ax[int(i/4)][i%4].scatter(X,Y,s=0.5,c=dET,alpha=0.8,)

    #Set the title
    ax[int(i/4)][i%4].set_aspect('equal','datalim')
    ax[int(i/4)][i%4].set_title('$\\mu$ = ' + str(u))

    #Make the colorbar for the energy
    (plt.colorbar(sc,ax=ax[int(i/4)][i%4])).formatter.set_powerlimits((0, 0))

    i+=1

plt.suptitle('Trajectories for diverse mass ratio ($\\theta_0$ = $\\pi$/2  &  No initial velocity)')
plt.savefig('../results/Trajectories-1S-Symplectic.png')
plt.show()


# 2 SWINGING CASE

#The initial parameters
x0 = np.array([ 3, 0, np.pi/2, 0 , np.pi/8 , 0 ])
#The inital R of M
r0M = 3

#We do the exact same as before but with the 2 swinging
fig,ax = plt.subplots(2,4,figsize=(20,10))
plt.subplots_adjust(wspace = 0.2)
i=0
for u in [2,3,5,6,16,19,21,24]:    
    r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0,r0M=r0M)
    X = r*np.sin(theta_m); Y = -r*np.cos(theta_m)
    sc = ax[int(i/4)][i%4].scatter(X,Y,s=0.5,c=dET,alpha=0.8,)
    ax[int(i/4)][i%4].set_aspect('equal','datalim')
    ax[int(i/4)][i%4].set_title('$\\mu$ = ' + str(u))

    (plt.colorbar(sc,ax=ax[int(i/4)][i%4])).formatter.set_powerlimits((0, 0))

    i+=1

plt.suptitle('Trajectories of m for diverse mass ratio ($\\theta_0$ = $\\pi$  &  No initial velocity)')
plt.savefig('../results/Trajectories-2S-Explicit.png')
plt.show()


#But this was only the trajectories of m. Now M also has some interesting trajectories so we do it again.
fig,ax = plt.subplots(2,4,figsize=(20,10))
plt.subplots_adjust(wspace = 0.2)
i=0
for u in [2,3,5,6,16,19,21,24]:    
    r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0,r0M=r0M)
    X_M = (r0M + x0[0]-r)*np.sin(theta_M); Y_M = -(r0M + x0[0]-r)*np.cos(theta_M)
    sc = ax[int(i/4)][i%4].scatter(X_M,Y_M,s=0.5,c=dET,alpha=0.8,)
    ax[int(i/4)][i%4].set_aspect('equal','datalim')
    ax[int(i/4)][i%4].set_title('$\\mu$ = ' + str(u))

    (plt.colorbar(sc,ax=ax[int(i/4)][i%4])).formatter.set_powerlimits((0, 0))

    i+=1

plt.suptitle('Trajectories of M for diverse mass ratio ($\\theta_0$ = $\\pi$  &  No initial velocity)')
plt.savefig('../results/Trajectories-2S-Explicit_M.png')
plt.show()



# Integrate for the smaller mass m with the leapfrog symplectic
fig,ax = plt.subplots(2,4,figsize=(20,10))
plt.subplots_adjust(wspace = 0.2)
i=0
for u in [2,3,5,6,16,19,21,24]:    
    r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Symplectic.integrate(u=u, tf=tf, dt=dt, x0=x0,r0M=r0M)
    X = r*np.sin(theta_m); Y = -r*np.cos(theta_m)
    sc = ax[int(i/4)][i%4].scatter(X,Y,s=0.5,c=dET,alpha=0.8,)
    ax[int(i/4)][i%4].set_aspect('equal','datalim')
    ax[int(i/4)][i%4].set_title('$\\mu$ = ' + str(u))

    (plt.colorbar(sc,ax=ax[int(i/4)][i%4])).formatter.set_powerlimits((0, 0))

    i+=1

plt.suptitle('Trajectories of m for diverse mass ratio ($\\theta_0$ = $\\pi$  &  No initial velocity)')
plt.savefig('../results/Trajectories-2S-Symplectic.png')
plt.show()

# integrate for the motion of the larger mass M
fig,ax = plt.subplots(2,4,figsize=(20,10))
plt.subplots_adjust(wspace = 0.2)
i=0
for u in [2,3,5,6,16,19,21,24]:    
    r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Symplectic.integrate(u=u, tf=tf, dt=dt, x0=x0,r0M=r0M)
    X_M = (r0M + x0[0]-r)*np.sin(theta_M); Y_M = -(r0M + x0[0]-r)*np.cos(theta_M)
    sc = ax[int(i/4)][i%4].scatter(X_M,Y_M,s=0.5,c=dET,alpha=0.8,)
    ax[int(i/4)][i%4].set_aspect('equal','datalim')
    ax[int(i/4)][i%4].set_title('$\\mu$ = ' + str(u))

    (plt.colorbar(sc,ax=ax[int(i/4)][i%4])).formatter.set_powerlimits((0, 0))

    i+=1

plt.suptitle('Trajectories of M for diverse mass ratio ($\\theta_0$ = $\\pi$  &  No initial velocity)')
plt.savefig('../results/Trajectories-2S-Symplectic_M.png')
plt.show()
