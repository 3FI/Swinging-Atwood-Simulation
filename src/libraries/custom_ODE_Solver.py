import numpy as np

def euler(nsteps, dt, x0, derivs):
    x = np.zeros((nsteps, len(x0)))
    x[0] = x0
    for i in range(1,nsteps):
        f = derivs((i-1)*dt, x[i-1])
        x[i] = x[i-1] + f*dt
    return x

def midpoint(nsteps, dt, x0, derivs):
    x = np.zeros((nsteps, len(x0)))
    x[0] = x0
    for i in range(1,nsteps):
        f = derivs((i-1)*dt, x[i-1])
        f1 = derivs((i-1)*dt + dt/2, x[i-1] + f*dt/2)
        x[i] = x[i-1] + f1*dt   
    return x

def RK4(nsteps, dt, x0, derivs):
    x = np.zeros((nsteps, len(x0)))
    x[0] = x0
    for i in range(1,nsteps):
        f = derivs((i-1)*dt, x[i-1])
        f1 = derivs((i-1)*dt + dt/2, x[i-1] + f*dt/2)
        f2 = derivs((i-1)*dt + dt/2, x[i-1] + f1*dt/2)
        f3 = derivs(i*dt, x[i-1] + f2*dt)
        x[i] = x[i-1] + dt*(f + 2*f1 + 2*f2 + f3)/6
    return x

def leapfrog(nsteps, dt, x0, v0):
    x = np.zeros((nsteps, 2))
    v = np.zeros((nsteps, 2))
    x[0] = x0
    v[0] = v0
    for i in range(1, nsteps):
        r = (x[i-1,0]**2 + x[i-1,1]**2)**0.5
        a = -x[i-1] / r**3
        v12 = v[i-1] + a * dt/2
        x[i] = x[i-1] + v12 * dt
        r = (x[i,0]**2 + x[i,1]**2)**0.5
        a = -x[i] / r**3
        v[i] = v12 + a * dt/2
    return x[:,0], x[:,1], v[:,0], v[:,1]