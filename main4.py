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
t_end = 20000 # the maximum value of time

t = np.arange(t_start, t_end, dt)
state_0 = [10, 10, 14] # initial conditions

y = odeint(system, state_0, t) # y is a list with elemets that are list with 3 entries containing u, Te, Tw for each time step
# this corresponds to the functions u, Te, Tw which are function of time
u = y[:,0]
Te = y[:,1]
Tw = y[:,2]

du = diff(u, t, dt) # calculating derivative of u

roots = find_all_maxima(u, du, t_start, t_end, .1, 100) # finding all the

# number of El-Nino events in the time run
num_of_elNinos = len(roots)

# calculating the periods of time that elapses between the El-Nino events
times_between_ENSO = [] # periods of time
for i in range(num_of_elNinos - 1):
    times_between_ENSO.append(roots[i+1] - roots[i]) # roots is sorted so dont need abs() function


# plot of ti against ti+10, ignoring the first 10 events
plt.scatter(times_between_ENSO[10:-10], times_between_ENSO[20:] , s = 5)
#plt.title(r"Scatter plot of $T_i$ vs $T_{(i+10)}$ for {} years")
plt.xlabel(r'$T_{(i+10)}$')
plt.ylabel(r'$T_i$')
plt.savefig("Scatter plot of T_i vs T_(i+10) for {} years".format(t_end))
plt.show()
# if uncorrelated, graph would have no distinguishable pattern