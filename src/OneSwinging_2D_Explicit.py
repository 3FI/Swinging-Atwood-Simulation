import libraries.custom_ODE_Solver as custom_ODE_Solver

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np
import scipy
import scipy.constants

def integrate(u, tf, dt, x0) : 

    #Taken from https://en.wikipedia.org/wiki/Swinging_Atwood%27s_machine
    def derivs(t,x):
        r,pr,theta,p_theta = x
        dr = pr/(M+m)
        dpr = p_theta**2 / (m*r**3) - M*scipy.constants.g + m*scipy.constants.g*np.cos(theta)
        dtheta = p_theta / (m*r**2)
        dp_theta = -m*scipy.constants.g*r*np.sin(theta)
        return np.array([dr,dpr,dtheta,dp_theta])
    
    m = 1
    M = u

    path, i = custom_ODE_Solver.RK4(int(tf/dt),dt,x0,derivs)
    r = path[:i,0]
    pr = path[:i,1]
    theta = path[:i,2]
    ptheta  = path[:i,3]

    ET = pr**2 / (2*(M+m)) + ptheta**2 / (2*m*r**2) + M*scipy.constants.g*r - m*scipy.constants.g*r*np.cos(theta)
    dET = 100* (ET - ET[0])/ET[0]

    return r,pr,theta,ptheta,dET