import numpy as np

def euler(nsteps, dt, x0, derivs):
    x = np.zeros((nsteps, len(x0)))
    x[0] = x0
    for i in range(1,nsteps):    
        f = derivs(i*dt, x[i-1])
        x[i] = x[i-1] + f*dt
    return x

def midpoint(nsteps, dt, x0, derivs):
    x = np.zeros((nsteps, len(x0)))
    x[0] = x0
    for i in range(1,nsteps):    
        f = derivs(i*dt, x[i-1])
        x1 = x[i-1] + f*dt/2
        f1 = derivs(i*dt + dt/2, x1)
        x[i] = x[i-1] + f1*dt
    return x

def RK4(nsteps, dt, x0, derivs):
    x = np.zeros((nsteps, len(x0)))
    x[0] = x0
    for i in range(1,nsteps):    
        f = derivs(i*dt, x[i-1])
        x1 = x[i-1] + f*dt/2
        f1 = derivs(i*dt + dt/2, x1)
        x2 = x[i-1] + f1*dt/2
        f2 = derivs(i*dt + dt/2, x2)
        x3 = x[i-1] + f2*dt
        f3 = derivs(i*dt + dt, x3)
        x[i] = x[i-1] + dt/6 * (f+2*f1+2*f2+f3)
    return x
