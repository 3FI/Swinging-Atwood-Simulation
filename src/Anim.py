
import OneSwinging_2D_Explicit
import TwoSwinging_2D_Explicit

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np



def plot_anim(X,Y,X_M,Y_M,tf,dt,dET,u,x0_M,saveName="noName"):
    fig,ax = plt.subplots(1,2,figsize=(9,6))

    scat = ax[0].scatter(X[0],Y[0],s=2, label='RK4')
    mPoint = ax[0].scatter(X[0],Y[0] ,s=40, c='k')
    Rope = ax[0].add_line(mlines.Line2D([X_M[0],x0_M,0,X[0]],[Y_M[0],0,0,Y[0]],c='k'))
    MPoint = ax[0].scatter(X_M[0],Y_M[0] ,s=40*u, c='k')
    ax[0].set(xlim=[-12, 6], ylim=[-6, 4], xlabel='x', ylabel='y')
    
    energy = (ax[1].plot(np.arange(0,tf,dt)[0],dET[0]))[0]
    ax[1].set(xlim=[0,tf], ylim=[min(dET), max(dET)], xlabel='t', ylabel='dE (%)')
    
    def update(num, x, y, scat, energy, mPoint, MPoint, Rope):
        #Update the trace
        scat.set_offsets(np.stack([x[:num], y[:num]]).T)
        scat.set_facecolor(cm.rainbow(np.arange(x[:num].size)%300))
    
        #Update the mass positions
        mPoint.set_offsets(np.stack([x[num], y[num]]).T)
        MPoint.set_offsets(np.stack([X_M[num], Y_M[num]]).T)
        #Update the rope drawing
        Rope.set_data([X_M[num],x0_M,0,X[num]],[Y_M[num],0,0,Y[num]])
    
        #Update the energy plot
        energy.set_data(np.arange(0,tf,dt)[:num],dET[:num])
        return scat,energy, mPoint, MPoint, Rope
    global anim
    anim = animation.FuncAnimation(fig, update, frames=len(X), interval=10, fargs=[X, Y, scat, energy, mPoint, MPoint, Rope])
    anim.save('../results/'+saveName, fps=60,writer="pillow")
    plt.show()

u = 3
tf = 10
dt = 0.01
x0_M = -10
r0_M = 3

x0 = np.array([ 3, 0, np.pi/4, 0  ])
ropeLength = r0_M + x0[0]

r,pr,theta,ptheta,dET = OneSwinging_2D_Explicit.integrate(u=u, tf=tf, dt=dt, x0=x0)
X = r*np.sin(theta); Y = -r*np.cos(theta)
X_M = np.full(len(r),x0_M); Y_M = r - ropeLength
plot_anim(X,Y,X_M,Y_M,tf,dt,dET,u,x0_M,"Anim_u=3-1S-Explicit.gif")


#Initial : (
#           ropelength of m, 
#           momentum of ropelength of m,
#           angle of m
#           rotational momentum of m
#           angle of M
#           rotational momentum of M
#          )
x0 = np.array([ 3, 0, np.pi/2, 0 , np.pi/8 , 0 ])
ropeLength = r0_M + x0[0]

r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET = TwoSwinging_2D_Explicit.integrate(u=3, tf=10, dt=0.01, x0=x0,r0M=r0_M)
#Transfer to planar coordinates
X = r*np.sin(theta_m); Y = -r*np.cos(theta_m)
X_M = (ropeLength-r)*np.sin(theta_M) + x0_M; Y_M = -(ropeLength-r)*np.cos(theta_M)
plot_anim(X,Y,X_M,Y_M,tf,dt,dET,u,x0_M,"Anim_u=3-2S-Explicit.gif")
