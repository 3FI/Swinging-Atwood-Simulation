import libraries.custom_ODE_Solver as custom_ODE_Solver

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np
import scipy
import scipy.constants

#See my notes in ref
def derivs(t,x):
    r,pr,theta_m,p_theta_m,theta_M,p_theta_M = x
    dr = pr/(M+m)
    dpr = -p_theta_M**2 / (M*(ropeLength-r)**3) + p_theta_m**2 / (m*r**3) - M*scipy.constants.g*np.cos(theta_M) + m*scipy.constants.g*np.cos(theta_m)
    dtheta_m = p_theta_m / (m*r**2)
    dp_theta_m = -m*scipy.constants.g*r*np.sin(theta_m)
    dtheta_M = p_theta_M / (M*(ropeLength-r)**2)
    dp_theta_M = -M*scipy.constants.g*(ropeLength-r)*np.sin(theta_M)
    
    return np.array([dr,dpr,dtheta_m,dp_theta_m,dtheta_M,dp_theta_M])

def derivs2(t, x):
    r,pr,theta_m,p_theta_m,theta_M,p_theta_M = x
    
    # first derivatives
    dr = pr/(M+m)
    dpr = -p_theta_M**2 / (M*(ropeLength-r)**3) \
        + p_theta_m**2 / (m*r**3) - M*scipy.constants.g*np.cos(theta_M) \
            + m*scipy.constants.g*np.cos(theta_m)
    dtheta_m = p_theta_m / (m*r**2)
    dp_theta_m = -m*scipy.constants.g*r*np.sin(theta_m)
    dtheta_M = p_theta_M / (M*(ropeLength-r)**2)
    dp_theta_M = -M*scipy.constants.g*(ropeLength-r)*np.sin(theta_M)

    dx = np.array([dr, dpr, dtheta_m, dp_theta_m, dtheta_M, dp_theta_M])

    # second derivatives
    d2r = dpr/(M+m)
    d2pr = -2*p_theta_M*dp_theta_M/(M*(ropeLength-r)**3) - 3*p_theta_M**2*dr/(M*(ropeLength-r)**4) \
        + 2*p_theta_m*dp_theta_m/(m*r**3) - 3*p_theta_m**2*dr/(m*r**4) \
            + scipy.constants.g*(M*np.sin(theta_M)*dtheta_M - m*np.sin(theta_m)*dtheta_m)
    d2theta_m = dp_theta_m/(m*r**2) - 2*p_theta_m*dr/(m*r**3)
    d2p_theta_m = -m*scipy.constants.g*(np.sin(theta_m)*dr - r*np.cos(theta_m)*dtheta_m)
    d2theta_M = dp_theta_M/(M*(ropeLength-r)**2) + 2*p_theta_M*dr/(M*(ropeLength-r)**3)
    d2p_theta_M = M*scipy.constants.g*(np.sin(theta_M)*dr - (ropeLength-r)*np.cos(theta_M)*dtheta_M)
    
    d2x = np.array([d2r, d2pr, d2theta_m, d2p_theta_m, d2theta_M, d2p_theta_M])

    return dx, d2x


fig,ax = plt.subplots(1,2,figsize=(9,6))

#Initial parameters
tf = 10         #Length of the simulation
dt = 0.02       #Time interval of the simulation
M =  2          #Mass of the larger mass
m = 1          #Mass of the smaller mass
r0_M = 3        #Initial radius of mass M
dx_m_M = -10    #Spacing between m and M on the drawing

#Initial : (
#           ropelength of m, 
#           momentum of ropelength of m,
#           angle of m
#           rotational momentum of m
#           angle of M
#           rotational momentum of M
#          )
x0 = np.array([ 3, 0, np.pi/4, 0 , np.pi/8 , 0 ])
ropeLength = r0_M + x0[0] #Total length on BOTH sides (m & M)

#SOLVE THE ODE
path = custom_ODE_Solver.euler_cromer(int(tf/dt),dt,x0,derivs2)

r = path[:,0]
pr = path[:,1]
theta_m = path[:,2]
ptheta_m  = path[:,3]
theta_M = path[:,4]
ptheta_M  = path[:,5]

#Transfer to planar coordinates
X = r*np.sin(theta_m); Y = -r*np.cos(theta_m)
X_M = (ropeLength-r)*np.sin(theta_M) + dx_m_M; Y_M = -(ropeLength-r)*np.cos(theta_M)

#Calculate the energy for each frame (equation of the hamiltonian)
ET = pr**2 / (2*(M+m)) + ptheta_M**2 / (2*M*(ropeLength-r)**2) + ptheta_m**2 / (2*m*r**2) - M*scipy.constants.g*(ropeLength-r)*np.cos(theta_M) - m*scipy.constants.g*r*np.cos(theta_m)
ET = 100* (ET - ET[0])/ET[0]

#Dotted trace of the mass m
scat = ax[0].scatter(X[0],Y[0],s=2, label='Euler-Cromer')
#TODO : ADD TRACE for mass M
mPoint = ax[0].scatter(X[0],Y[0] ,s=40*m, c='k')
MPoint = ax[0].scatter(X_M[0],Y_M[0] ,s=40*M, c='k')
#The drawing of the rope starting from M then going to the anchor point of both M and m and finally reaching m
Rope = ax[0].add_line(mlines.Line2D([X_M[0],dx_m_M,0,X[0]],[Y_M[0],0,0,Y[0]],c='k'))

#TODO : UPDATE SO IT SCALES
ax[0].set(xlim=[-12, 6], ylim=[-6, 4], xlabel='x', ylabel='y')

#Subplot 2 is the energy
energy = (ax[1].plot(np.arange(0,tf,dt)[0],ET[0]))[0]
ax[1].set(xlim=[0,tf], ylim=[min(ET), max(ET)], xlabel='t', ylabel='dE (%)')

#Update the drawing every 10 frames
def update(num, x, y, x_M, y_M, scat, energy, mPoint, MPoint, Rope):
    #Update the trace
    scat.set_offsets(np.stack([x[:num], y[:num]]).T)
    scat.set_facecolor(cm.rainbow(np.arange(x[:num].size)%300))

    #Update the mass positions
    mPoint.set_offsets(np.stack([x[num], y[num]]).T)
    MPoint.set_offsets(np.stack([x_M[num], y_M[num]]).T)
    #Update the rope drawing
    Rope.set_data([x_M[num],dx_m_M,0,x[num]],[y_M[num],0,0,y[num]])

    #Update the energy plot
    energy.set_data(np.arange(0,tf,dt)[:num],ET[:num])
    return scat,energy, mPoint, MPoint, Rope

test = animation.FuncAnimation(fig, update, frames=len(X), interval=10, fargs=[X, Y, X_M, Y_M, scat, energy, mPoint, MPoint, Rope])
plt.show()
#test.save('../results/test.gif', fps=60,writer="pillow")


fig,ax = plt.subplots(2,4)
i=0
for M in [2,3,5,6,16,19,21,24]:
    m=1
    tf = 100
    x0 = np.array([ 3, 0, np.pi/2, 0 , 0, 0 ])
    path = custom_ODE_Solver.euler_cromer(int(tf/dt),dt,x0,derivs2)
    r = path[:,0]
    theta_m = path[:,2]
    X = r*np.sin(theta_m); Y = -r*np.cos(theta_m)
    ax[int(i/4)][i%4].plot(X,Y,'k-',lw=1,alpha=0.8)
    ax[int(i/4)][i%4].set_aspect('equal','datalim')
    ax[int(i/4)][i%4].set_title('mu = ' + str(M))
    i+=1

plt.suptitle('Trajectories for diverse mass ratio (Theta_0 = pi/2  &  No initial velocity)')
#\plt.show()