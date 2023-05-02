import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# making a function for the set of differential equations
def system(state, t):
    """
    Parameters
    ----------
    state : array with 3 values
    t : array for time steps

    Returns
    -------
    du/dt, dTe/dt, dTw/dt : array of of values for given state and t
    """

    u, T_e, T_w = state

    # constants, to the right of each constant i noted its original value
    A = 1.0 # 1.0
    B = 663.0 # 663.0
    C = 3.0 # 3.0
    T_star = 12.0 # 12.0
    u_star = -14.2 # -14.2
    delta_x = 7.5 # 7.5

    # setting up the 3 differential equations
    du = B / delta_x * (T_e - T_w) - C * (u - u_star) #  = du/dt
    dTe = u * T_w / (2 * delta_x) - A * (T_e - T_star) #  = dTe/dt
    dTw = -u * T_e / (2 * delta_x) - A * (T_w - T_star) #  = dTw/dt

    return du, dTe, dTw

# calculates the norms/distances between two systems
def norm(funcs1, funcs2, t):
    """
    Parameters
    ----------
    funcs1 : odeint output for first system
    funcs2 : odeint output for second system
    t : array for time steps
    
    Returns
    -------
    returns array of distances bewteen the fist and second system
    """
    x1 = funcs1[:,0]
    y1 = funcs1[:,1]
    z1 = funcs1[:,2]

    x2 = funcs2[:,0]
    y2 = funcs2[:,1]
    z2 = funcs2[:,2]

    norms = []
    for i in range(len(t)):
        norms.append( np.sqrt( (x1[i] - x2[i])**2 + (y1[i] - y2[i])**2 + (z1[i] - z2[i])**2) )

    return norms


dt = .001 # resolution
t_start = 0 # this will allways be 0
t_end = 50 # the maximum value of time

t = np.arange(t_start, t_end, dt)

state_0 = [10, 10, 14] # initial conditions
y = odeint(system, state_0, t) # y is a list with elemets that are list with 3 entries containing u, Te, Tw for each time step
u = y[:,0]
Te = y[:,1]
Tw = y[:,2]

# this bit is to check how the distance between the system and the same system with slight chaged initial conditions
# changes over time.
state_0_plus = [10.001, 10, 14] # original inital conditions [10, 10, 14]
y_plus = odeint(system, state_0_plus, t) # justify why small change in initial conditions is small enough
u_plus = y[:,0]
Te_plus = y[:,1]
Tw_plus = y[:,2]

plt.plot(t, norm(y, y_plus, t))
plt.xlabel("t")
plt.ylabel("distance")
plt.show()
