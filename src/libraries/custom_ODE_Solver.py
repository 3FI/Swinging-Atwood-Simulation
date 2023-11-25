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

def leapfrog(nsteps, dt, x0, derivs2):
    
    r = x0[0]
    pr = x0[1]
    theta_m = x0[2]
    ptheta_m  = x0[3]
    theta_M = x0[4]
    ptheta_M  = x0[5]

    x = np.zeros((nsteps, len(x0)))
    break_i = nsteps

    for i in range(1, nsteps):
        dx, d2x = derivs2((i-1)*dt, x[i-1])
        v_12 = x[i-1][1::2] + d2x * dt/2
        x[i][0::2] = x[i-1][0::2] + v_12 * dt
        dx, d2x = derivs2(i*dt, x[i])
        x[i][1::2] = v_12 + d2x * dt/2

        if x[i][0] <= 0.1 :
            break_i = i
            break
    return x, break_i


def euler_cromer(nsteps, dt, x0, derivs2):
    x = np.zeros((nsteps, len(x0)))
    x[0] = x0
    v = np.zeros((nsteps, int(len(x0)/2)))
    dx0, d2x0 = derivs2(dt, x0)
    v[0] = dx0[0::2]
    break_i = nsteps
    for i in range(1, nsteps):
        dx, d2x = derivs2(dt, x[i-1])
        v[i] = v[i-1] + dt*d2x[0::2]
        x[i][1::2] = x[i-1][1::2] + dt*d2x[1::2]
        x[i][0::2] = x[i-1][0::2] + dt*v[i]

        if x[i][0] <= 0.1 : 
            break_i = i
        '''r_new, _, theta_new, _ = x[i]
        dr_new, _, dtheta_new, _ = dx_new
        pr_new = (M+m)*dr_new
        ptheta_new = m*r_new**2*dtheta_new
        x[i][1] = pr_new
        x[i][3] = ptheta_new'''

    return x, break_i