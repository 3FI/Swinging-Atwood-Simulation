import numpy as np

#This is a file containing various ODE solver


def RK4(nsteps, dt, x0, derivs):
    x = np.zeros((nsteps, len(x0)))
    x[0] = x0
    break_i = nsteps
    for i in range(1,nsteps):
        
        f = derivs((i-1)*dt, x[i-1])
        f1 = derivs((i-1)*dt + dt/2, x[i-1] + f*dt/2)
        f2 = derivs((i-1)*dt + dt/2, x[i-1] + f1*dt/2)
        f3 = derivs(i*dt, x[i-1] + f2*dt)
        x[i] = x[i-1] + dt*(f + 2*f1 + 2*f2 + f3)/6

        if x[i][0] <= 0.1 :
            break_i = i
            break
        
    return x, break_i



#These are the unsuccessful symplectic integrators

def euler_cromer(nsteps, dt, x0, derivs2):
    x = np.zeros((nsteps, len(x0)))
    x[0] = x0
    break_i = nsteps
    for i in range(1, nsteps):
        dx, d2x = derivs2(dt, x[i-1])
        dx_new = dx + d2x*dt
        x[i] =x[i-1] + dx_new*dt        
        if x[i][0] <= 0.1 : 
            break_i = i

    return x, break_i


def leapfrog(nsteps, dt, x0, derivs2):

    x = np.zeros((nsteps, len(x0)))
    x[0] = x0
    break_i = nsteps

    for i in range(1, nsteps):
        dx, d2x = derivs2((i-1)*dt, x[i-1])
        x12 = x[i-1] + 0.5*dt*dx
        dx_12, d2x_12 = derivs2((i+0.5)*dt, x12)
        dx_new = dx + dt*d2x_12
        x[i] = x12 + 0.5*dt*dx_new

        if x[i][0] <= 0.1 :
            break_i = i
            break

    return x, break_i


#These are test integrators

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