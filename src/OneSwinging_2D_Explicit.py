import libraries.custom_ODE_Solver as custom_ODE_Solver

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np
import scipy
import scipy.constants

def integrate(u, tf, dt, x0) : 
    """
    The integration function for the 1 swinging case

    Parameters
    ----------
    u : float
        The mass ratio
    tf : float
        The time up to which we want to carry the integration if we can
    dt : float
        The time interval between each step
    x0 : floats[4]
        The inital parameters of the form :(
            radius of m, 

            momentum of radius of m,

            angle of m,

            angular momentum of m
        )

    Returns
    ----------
    r : float[]
        The radius for each time
    pr : float[]
        The radius momentum for each time
    theta : float[]
        The angle of m for each time
    ptheta : float[]
        The angular momentum of m for each time
    dET : float[]
        The energy error in % for each time
    """
    #See the derivations in the report
    def derivs(t,x):
        r,pr,theta,p_theta = x
        dr = pr/(M+m)
        dpr = p_theta**2 / (m*r**3) - M*scipy.constants.g + m*scipy.constants.g*np.cos(theta)
        dtheta = p_theta / (m*r**2)
        dp_theta = -m*scipy.constants.g*r*np.sin(theta)
        return np.array([dr,dpr,dtheta,dp_theta])
    
    #Because of independance under mass multiplication we set m=1
    m = 1
    M = u

    #We solve the ODE
    path, i = custom_ODE_Solver.RK4(int(tf/dt),dt,x0,derivs)
    r = path[:i,0]
    pr = path[:i,1]
    theta = path[:i,2]
    ptheta  = path[:i,3]

    #Calculate the energy using the hamiltonian found in the report
    ET = pr**2 / (2*(M+m)) + ptheta**2 / (2*m*r**2) + M*scipy.constants.g*r - m*scipy.constants.g*r*np.cos(theta)
    dET = 100* (ET - ET[0])/ET[0]

    return r,pr,theta,ptheta,dET


# This is an alteration of the function above which is only used in the mass independance analysis. Instead of taking
# a mu it takes both masses explicitely.
def integrateMass(m, M, tf, dt, x0) : 
    def derivs(t,x):
        r,pr,theta,p_theta = x
        dr = pr/(M+m)
        dpr = p_theta**2 / (m*r**3) - M*scipy.constants.g + m*scipy.constants.g*np.cos(theta)
        dtheta = p_theta / (m*r**2)
        dp_theta = -m*scipy.constants.g*r*np.sin(theta)
        return np.array([dr,dpr,dtheta,dp_theta])
    path, i = custom_ODE_Solver.RK4(int(tf/dt),dt,x0,derivs)
    r = path[:i,0]
    pr = path[:i,1]
    theta = path[:i,2]
    ptheta  = path[:i,3]
    ET = pr**2 / (2*(M+m)) + ptheta**2 / (2*m*r**2) + M*scipy.constants.g*r - m*scipy.constants.g*r*np.cos(theta)
    dET = 100* (ET - ET[0])/ET[0]
    return r,pr,theta,ptheta,dET