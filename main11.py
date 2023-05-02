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

# intagrates a function on its entire interval
def integrate(func, dx):
    """
    Parameters
    ----------
    func : function in array form (from odeint)
    dx : step size

    Returns
    -------
    integral of function for where its defined
    """
    result = 0
    for i in range(len(func) - 1):
        result += func[i] + func[i + 1]
    return result * (dx/2)

# calculates the fourier series of a function in array form
def fourier_series(func, x, dx, end, n, total_series):

    a = []
    for j in range(n):
        cos = []
        for i in range(len(x)):
            cos.append(np.cos((j * np.pi * x[i]) / end))
        cos = np.array(cos)
        a_n = integrate(func * cos, dx)
        a.append(a_n)
    a[0] = a[0]/2

    b = []
    for j in range(n):
        sin = []
        for i in range(len(x)):
            sin.append(np.sin(j * np.pi * x[i]) / end)
        sin = np.array(sin)
        b_n = integrate(func * sin, dx)
        b.append(b_n)

    fourier_list = []
    if total_series == True:
        for i in range(len(x)):
            slot = 0
            for j in range(n):
                slot += (2/end) * (a[j] * np.cos((j * np.pi * x[i]) / end) + b[j] * np.sin((j * np.pi * x[i]) / end))
            fourier_list.append(slot)
    else:
        for i in range(len(x)):
            slot = []
            for j in range(n):
                slot.append((2/end) * (a[j] * np.cos((j * np.pi * x[i]) / end) + b[j] * np.sin((j * np.pi * x[i]) / end)))
            fourier_list.append(slot)

    return fourier_list

dt = .001 # resolution
t_start = 0 # this will allways be 0
t_end = 20 # the maximum value of time

t = np.arange(t_start, t_end, dt)
state_0 = [10, 10, 14] # initial conditions

y = odeint(system, state_0, t) # y is a list with elemets that are list with 3 entries containing u, Te, Tw for each time step
# this corresponds to the functions u, Te, Tw which are function of time
u = y[:,0]
Te = y[:,1]
Tw = y[:,2]


# fourier series stuff ( fun but i think ultimately useless unless you see some patterns which i didnt notice)
series_degree = 40
sift = 1
# this bit is the fourier series up to series degree we chose
u_fourier = fourier_series(u, t, dt, t_end, series_degree, True)
plt.plot(t, u_fourier, "--")
plt.plot(t, u) #, label="u")
# this bit docomposes the fourier series into sin and cosine waves
"""un_fourier = fourier_series(u, t, dt, t_end, series_degree, False)
for i in range(int(series_degree/sift)):
    n_term = np.array(un_fourier)[:,sift*i]
    plt.plot(t, n_term, label = "{}".format(i))
plt.legend()"""
plt.show()