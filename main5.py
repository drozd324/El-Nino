import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# in the first half, we show all of the functions we will use in the program because it think it will be easier to read.
# Then we will have all of the plots and results we want to obtain in the second half.

""" FIRST HALF: DECLARING FUNCTIONS """

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


dt = .01 # resolution
t_start = 0 # this will allways be 0
t_end = 100 # the maximum value of time

t = np.arange(t_start, t_end, dt)
state_0 = [10, 10, 14] # initial conditions

y = odeint(system, state_0, t) # y is a list with elemets that are list with 3 entries containing u, Te, Tw for each time step
# this corresponds to the functions u, Te, Tw which are function of time
u = y[:,0]
Te = y[:,1]
Tw = y[:,2]

# plot of the fractal 8 figure thing
plt.plot(u, Te - Tw)
#plt.title("Difference in Temperature against current velocity (fractal)".format(t_end))
plt.xlabel("u")
plt.ylabel(r'$T_e - T_w$')
plt.ylim((-30,30))
plt.savefig("Difference in Temperature against current velocity (fractal thing) for {} years".format(t_end))
plt.show()