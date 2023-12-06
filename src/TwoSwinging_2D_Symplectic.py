import libraries.custom_ODE_Solver as custom_ODE_Solver

import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.lines as mlines
from matplotlib import cm
import numpy as np
import scipy
import scipy.constants

def integrate(u, tf, dt, x0, r0M) : 
    """
     The integration function for the 2 swinging case

    Parameters
    ----------
    u : float
        The mass ratio
    tf : float
        The time up to which we want to carry the integration if we can
    dt : float
        The time interval between each step
    x0 : floats[6]
        The inital parameters of the form :(
            radius of m, 

            momentum of radius of m,

            angle of m,

            angular momentum of m,

            angle of M,

            angular momentum of M,
        )
    r0M: float
        the initial R of M

    Returns
    ----------
    r : float[]
        The radius for each time
    pr : float[]
        The radius momentum for each time
    theta_m : float[]
        The angle of m for each time
    ptheta_m : float[]
        The angular momentum of m for each time
    theta_M : float[]
        The angle of M for each time
    ptheta_M : float[]
        The angular momentum for each time
    dET : float[]
        The energy error in % for each time
    """
    #See the derivations in the report
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
    
    #Because of independance under mass multiplication we set m=1
    m = 1
    M = u
    ropeLength = r0M + x0[0] #Total length on BOTH sides (m & M)

    #We solve the ODE
    path, i = custom_ODE_Solver.leapfrog(int(tf/dt),dt,x0,derivs2)
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