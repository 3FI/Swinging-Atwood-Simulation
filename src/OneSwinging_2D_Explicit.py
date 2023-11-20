import libraries.custom_ODE_Solver as custom_ODE_Solver

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np
import scipy
import scipy.constants

#Taken from https://en.wikipedia.org/wiki/Swinging_Atwood%27s_machine
def derivs(t,x):
    r,pr,theta,p_theta = x
    dr = pr/(M+m)
    dpr = p_theta**2 / (m*r**3) - M*scipy.constants.g + m*scipy.constants.g*np.cos(theta)
    dtheta = p_theta / (m*r**2)
    dp_theta = -m*scipy.constants.g*r*np.sin(theta)
    return np.array([dr,dpr,dtheta,dp_theta])

fig,ax = plt.subplots(1,2,figsize=(9,6))


tf = 10
dt = 0.02
M = 3
m = 2
y0_M = 3
x0_M = -10
x0 = np.array([ 3, 0, np.pi/4, 0  ])
ropeLength = y0_M + x0[0]

path = custom_ODE_Solver.RK4(int(tf/dt),dt,x0,derivs)
r = path[:,0]
pr = path[:,1]
theta = path[:,2]
ptheta  = path[:,3]
X = r*np.sin(theta); Y = -r*np.cos(theta)
ET = pr**2 / (2*(M+m)) + ptheta**2 / (2*m*r**2) + M*scipy.constants.g*r - m*scipy.constants.g*r*np.cos(theta)
ET = 100* (ET - ET[0])/ET[0]

scat = ax[0].scatter(X[0],Y[0],s=2, label='RK4')
mPoint = ax[0].scatter(X[0],Y[0] ,s=40*m, c='k')
Rope = ax[0].add_line(mlines.Line2D([x0_M,x0_M,0,X[0]],[-y0_M,0,0,Y[0]],c='k'))
MPoint = ax[0].scatter(x0_M,-y0_M ,s=40*M, c='k')
ax[0].set(xlim=[-12, 6], ylim=[-6, 4], xlabel='x', ylabel='y')

energy = (ax[1].plot(np.arange(0,tf,dt)[0],ET[0]))[0]
ax[1].set(xlim=[0,tf], ylim=[min(ET), max(ET)], xlabel='t', ylabel='dE (%)')

def update(num, x, y, scat, energy, mPoint, MPoint, Rope):
    #Update the trace
    scat.set_offsets(np.stack([x[:num], y[:num]]).T)
    scat.set_facecolor(cm.rainbow(np.arange(x[:num].size)%300))

    #Update the mass positions
    mPoint.set_offsets(np.stack([x[num], y[num]]).T)
    MPoint.set_offsets(np.stack([x0_M, r[num]-ropeLength]).T)
    #Update the rope drawing
    Rope.set_data([x0_M,x0_M,0,X[num]],[r[num]-ropeLength,0,0,Y[num]])

    #Update the energy plot
    energy.set_data(np.arange(0,tf,dt)[:num],ET[:num])
    return scat,energy, mPoint, MPoint, Rope

test = animation.FuncAnimation(fig, update, frames=len(X), interval=10, fargs=[X, Y, scat, energy, mPoint, MPoint, Rope])
plt.show()
#test.save('../results/test.gif', fps=60,writer="pillow")


fig,ax = plt.subplots(2,4)
i=0
for M in [2,3,5,6,16,19,21,24]:
    m=1
    tf = 100
    x0 = np.array([ 3, 0, np.pi/2, 0  ])
    path = custom_ODE_Solver.RK4(int(tf/dt),dt,x0,derivs)
    r = path[:,0]
    theta = path[:,2]
    X = r*np.sin(theta); Y = -r*np.cos(theta)
    ax[int(i/4)][i%4].plot(X,Y,'k-',lw=1,alpha=0.8)
    ax[int(i/4)][i%4].set_aspect('equal','datalim')
    ax[int(i/4)][i%4].set_title('mu = ' + str(M))
    i+=1

plt.suptitle('Trajectories for diverse mass ratio (Theta_0 = pi/2  &  No initial velocity)')
plt.show()