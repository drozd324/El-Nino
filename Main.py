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

# calculates the average of a list of numbers
def mean(mylist):
    """
    Parameters
    ----------
    mylist : list of values

    Returns
    -------
    average/mean of mylist
    """
    list_lenght = len(mylist)
    if list_lenght == 0:
        return 0
    else:
        my_sum = 0
        for i in range(list_lenght):
            my_sum += mylist[i]
        return my_sum / list_lenght

# calculates the standard deviation of a list of numbers
def standard_dev(mylist):
    """
    Parameters
    ----------
    mylist : list of values

    Returns
    -------
    standard deviation of mylist
    """
    list_length = len(mylist)
    if list_length == 0:
        return 0
    else:
        my_sum = 0
        average = mean(mylist)
        for i in range(list_length):
            my_sum += (mylist[i] - average)**2
        return np.sqrt(my_sum/list_length)

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

def norm_untill(funcs1, funcs2, t, max_dist):
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

""" SECOND HALF: OBTAINING RESULTS """

dt = .001 # resolution
t_start = 0 # this will allways be 0
t_end = 50 # the maximum value of time

t = np.arange(t_start, t_end, dt)
state_0 = [10, 10, 14] # initial conditions

y = odeint(system, state_0, t) # y is a list with elemets that are list with 3 entries containing u, Te, Tw for each time step
# this corresponds to the functions u, Te, Tw which are function of time
u = y[:,0]
Te = y[:,1]
Tw = y[:,2]


# this bit is to check how the distance between the system and the same system with slight chaged initial conditions
# changes over time.
"""state_0_plus = [10, 10.01, 14] # original inital conditions [10, 10, 14]
y_plus = odeint(system, state_0_plus, t) # justify why small change in initial conditions is small enough
u_plus = y[:,0]
Te_plus = y[:,1]
Tw_plus = y[:,2]
plt.plot(t, norm(y, y_plus, t))
plt.show()"""


# first plots
"""plt.plot(t, u) #, label="u")
plt.title("Current velocity against time ({} years)".format(t_end))
plt.xlabel("Time t (years)")
plt.ylabel("Current Velocity u ") # 10^3km / year  units?
plt.ylim((-400,400))
plt.xlim((0, t_end))
#plt.legend()
#plt.grid(axis = 'y')
plt.show()"""

"""plt.plot(t, Te - Tw)
plt.title("Difference in Temp against time ({} years)".format(t_end))
plt.xlabel("Time t (years)")
plt.ylabel(r'$T_e - T_w$')
plt.ylim((-30,30))
plt.show()"""


du = diff(u, t, dt) # calculating derivative of u
# this will plot the derivative of u
#plt.plot(t[1:-1], du, "--")

roots = find_all_maxima(u, du, t_start, t_end, .1, 100) # finding all the
# this will mark all the maxima in the list of roots on the graph
"""for i in range(len(roots)):
    plt.scatter(roots[i], make_func(u, roots[i]))
plt.show() """

# number of El-Nino events in the time run
num_of_elNinos = len(roots)
print("The number of ENSO events in {} years is: {}".format(t_end, num_of_elNinos))

# calculating the periods of time that elapses between the El-Nino events
times_between_ENSO = [] # periods of time
for i in range(num_of_elNinos - 1):
    times_between_ENSO.append(roots[i+1] - roots[i]) # roots is sorted so dont need abs() function
#times_between_ENSO.sort()

# finding the mean, and standard deviations (igonring the first 10 ENSO events)
average_T = mean(times_between_ENSO[10:])
dev_T = standard_dev(times_between_ENSO[10:])
print("The mean time between ENSO events (ignoring the first 10) is: {} years".format(average_T))
print("The standard deviation of ENSO events (ignoring the first 10) is: {}".format(dev_T))


# histogram plot for times between ENSO events
"""plt.hist(times_between_ENSO, bins = 260, edgecolor="white", range=[1.5, 7])
plt.title("Histogram plot of times between ENSO events for {} years".format(t_end))
plt.ylabel("Number of ENSO events")
plt.xlabel("Time beteen ENSO events (years)")
plt.show()"""

# plot of the fractal 8 figure thing
"""plt.plot(u, Te - Tw)
plt.title("Difference in Temperature against current velocity".format(t_end))
plt.xlabel("u")
plt.ylabel(r'$T_e - T_w$')
plt.ylim((-30,30))
plt.show()"""

# plot of ti against ti+10, ignoring the first 10 events
"""plt.scatter(times_between_ENSO[10:-10], times_between_ENSO[20:] , s = 5)
plt.title(r"Scatter plot of $T_i$ vs $T_{(i+10)}$ for {} years")
plt.xlabel(r'$T_{(i+10)}$')
plt.ylabel(r'$T_i$')
plt.show()"""
# if uncorrelated, graph would have no distinguishable pattern

# calculating u, Te, Tw for modified system
"""mod_y = odeint(mod_system, state_0, t)
mod_u = mod_y[:,0]
mod_Te = mod_y[:,1]
mod_Tw = mod_y[:,2]
#plt.plot(t, mod_u, label="mod_u")
mod_roots = find_all_maxima(mod_u, du, t_start, t_end, .1, 100)

ENSO_months = []
for i in range(len(mod_roots)):
    ENSO_months.append(mod_roots[i] % 1)

# histogram plot for in which months ENSO events occur
plt.hist(ENSO_months, bins = 12, edgecolor="white")
plt.title("Times of the year an ENSO event occurs over {} years".format(t_end))
plt.ylabel("ENSO event")
plt.xlabel("Months")
plt.show()"""

# 3d plot
"""fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(u, Te, Tw)
ax.set_xlabel(r'$u$')
ax.set_ylabel(r'$T_e$')
ax.set_zlabel(r'$T_w$')
ax.set_title("Strange Attractor")
plt.show()"""

# fourier series stuff ( fun but i think ultimately useless unless you see some patterns which i didnt notice)
series_degree = 10
sift = 1
# this bit is the fourier series up to series degree we chose
"""u_fourier = fourier_series(u, t, dt, t_end, series_degree, True)
plt.plot(t, u_fourier, "--")"""
# this bit docomposes the fourier series into sin and cosine waves
"""un_fourier = fourier_series(u, t, dt, t_end, series_degree, False)
for i in range(int(series_degree/sift)):
    n_term = np.array(un_fourier)[:,sift*i]
    plt.plot(t, n_term, label = "{}".format(i))
plt.legend()"""
#plt.show()
