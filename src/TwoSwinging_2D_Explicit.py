import libraries.custom_ODE_Solver as custom_ODE_Solver

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np
import scipy
import scipy.constants


#x0 : (
#           ropelength of m, 
#           momentum of ropelength of m,
#           angle of m
#           rotational momentum of m
#           angle of M
#           rotational momentum of M
#          )
def integrate(u, tf, dt, x0, r0M) : 

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
    
    m = 1
    M = u
    ropeLength = r0M + x0[0] #Total length on BOTH sides (m & M)
    
    #SOLVE THE ODE
    path, i = custom_ODE_Solver.RK4(int(tf/dt),dt,x0,derivs)

    r = path[:i,0]
    pr = path[:i,1]
    theta_m = path[:i,2]
    ptheta_m  = path[:i,3]
    theta_M = path[:i,4]
    ptheta_M  = path[:i,5]

    #Calculate the energy for each frame (equation of the hamiltonian)
    ET = pr**2 / (2*(M+m)) + ptheta_M**2 / (2*M*(ropeLength-r)**2) + ptheta_m**2 / (2*m*r**2) - M*scipy.constants.g*(ropeLength-r)*np.cos(theta_M) - m*scipy.constants.g*r*np.cos(theta_m)
    dET = 100* (ET - ET[0])/ET[0]

    return r,pr,theta_m,ptheta_m,theta_M,ptheta_M,dET