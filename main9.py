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