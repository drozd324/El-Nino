import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


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


#this is the final modified system with u* = u * (1 + *( 3*sin(2*pi*t)) )
def mod_system(state, t):
    """
    Parameters
    ----------
    state : array with 3 values
    t : array for time steps
    Returns
    -------
    (du/dt, dTe/dt, dTw/dt) : array of of values for given state and t
    """

    u, T_e, T_w = state

    # constants, to the right of each constant i noted its original value
    A = 1.0 # 1.0
    #if you double A, current velocity flattens after about 7 years
    B = 663.0 # 663.0
    C = 3.0 # 3.0
    T_star = 12.0 # 12.0
    u_star = -14.2 # -14.2
    delta_x = 7.5 # 7.5

    # setting up the 3 differential equations
    du = B / delta_x * (T_e - T_w) - C * (u - u_star * ( 1 + 3 * np.sin( 2 * np.pi * t ))) #  = du/dt
    dTe = u * T_w / (2 * delta_x) - A * (T_e - T_star) #  = dTe/dt
    dTw = -u * T_e / (2 * delta_x) - A * (T_w - T_star) #  = dTw/dt

    return du, dTe, dTw

# function for taking the derivative of a function (symmetric derivative)
def diff(func, x, dx):
    """
    Parameters
    ----------
    func : function in array form (from odeint)
    x : numpy array of changing valiable
    dx : step size
    Returns
    -------
    derivative of function in array form
    """

    dfs = np.zeros(len(x) - 2)
    for i in range(1, len(x) - 1):
        dfs[i-1] = (func[i + 1] - func[i - 1]) / (2 * dx)
    return dfs

# this will turn the numpy array into a something that works like a mathematical function
def make_func(list_func, x):
    """
    Parameters
    ----------
    list_func : function in array form (from odeint)
    x : value at which the function will be evaluated at
    Returns
    -------
    Value of function at x
    """
    # first bit makes sure that we wont round into an index outside the reach of the list
    if len(list_func) - 1 < int(round(x/dt)):
        return list_func[-1]
    else:
        return list_func[int(round(x/dt))]

# function for root finding (bisection method)
def find_root(a, b, func):
    """
    Parameters
    ----------
    a : float number for left bracket
    b : float number for right bracket
    func : function in array form (from odeint)
    Returns
    -------
    root of function between on the interval [a,b], if no root is found then
    reutuns a string
    """

    x1 = 0
    x2 = 0
    x3 = 0
    if make_func(func, a) * make_func(func, b) > 0:
        x3 = "None" # im makeing the output a string if there is no root
        run = False
    if make_func(func, a) < 0:
        x1 = a
        x2 = b
    else:
        x1 = b
        x2 = a

    run = True
    times_run = 0
    while run:
        times_run += 1
        x3 = (x1 + x2)/2
        if make_func(func, x3) < 0:
            x1 = x3
        else:
            x2 = x3
        if np.abs(x1 - x2) < 10**(-3): # precison
            run = False
        elif times_run > 100:
            run = False
            x3 = "None"  # im makeing the output a string if there is no root

    return x3

# sexy root finder
def find_all_maxima(f, df, start_step, end_step, step_size, y_shift):
    """
    Parameters
    ----------
    f : function in array form (from odeint)
    df : derivative of function in array form
    start_step : start of interval on which the function is defined from
    end_step : end of interval on which the function is defined on
    step_size : step size
    y_shift : shifts the function along the y axis by given value
    Returns
    -------
    return the values of highest peaks of the function
    """
    roots = []
    brackets = []
    a = start_step
    b = a + step_size
    brackets_start = 0
    for i in range(int(end_step/step_size)):
        if (make_func(f, a) - y_shift) * (make_func(f, b) - y_shift) < 0:
            if len(brackets) == 0:
                if (make_func(f, a) - y_shift) > 0:
                    brackets_start = 1
            brackets.append(a)
        a += step_size
        b += step_size

    for i in range(int(len(brackets) / 2)):
        j = 2*i + brackets_start
        a = brackets[j]
        b = brackets[j + 1]
        root = find_root(a, b, df)
        roots.append(root)

    return roots

dt = .001 # resolution
t_start = 0 # this will allways be 0
t_end = 500 # the maximum value of time

t = np.arange(t_start, t_end, dt)
state_0 = [10, 10, 14] # initial conditions

y = odeint(system, state_0, t) # y is a list with elemets that are list with 3 entries containing u, Te, Tw for each time step
# this corresponds to the functions u, Te, Tw which are function of time
u = y[:,0]
Te = y[:,1]
Tw = y[:,2]

# calculating u, Te, Tw for modified system
mod_y = odeint(mod_system, state_0, t)
mod_u = mod_y[:,0]
mod_Te = mod_y[:,1]
mod_Tw = mod_y[:,2]

mod_du = diff(mod_u, t, dt)

mod_roots = find_all_maxima(mod_u, mod_du, t_start, t_end, .1, 100)
ENSO_months = []
for i in range(len(mod_roots)):
    ENSO_months.append(mod_roots[i] % 1)

# histogram plot for in which months ENSO events occur
plt.hist(ENSO_months, bins = 12, edgecolor="white")
#plt.title("Times of the year an ENSO event occurs over {} years".format(t_end))
plt.ylabel("ENSO event")
plt.xlabel("Months")
plt.savefig("Times of the year an ENSO event occurs over {} years".format(t_end))
plt.show()