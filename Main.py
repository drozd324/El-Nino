import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# declaring constants

A = 1.0
B = 663.0
C = 3.0
T_star = 12.0
u_star = -14.2
delta_x = 7.5

# making the vectorised differential equations

def f(state, t):
    u, T_e, T_w = state

    # setting up the 3 differential equations

    du = B / delta_x * (T_e - T_w) - C * (u - u_star) #  = du/dt
    dTe = u * T_w / (2 * delta_x) - A * (T_e - T_star) #  = dTe/dt
    dTw = -u * T_e / (2 * delta_x) - A * (T_w - T_star) #  = dTw/dt

    return du, dTe, dTw 

dt = .00001 # resolution

t_start = 0 # this will allways be 0
t_end = 100 # this will be for how many years we want to see this run # check 40 and 50

t = np.arange(t_start, t_end, dt)
state_0 = [10, 10, 14] # initial conditions
y = odeint(f, state_0, t) # y is a list with elemets that are list with 3 entries containing u, Te, Tw for each time step

# this corresponds to the functions u, Te, Tw which are function of time
u = y[:,0]
Te = y[:,1]
Tw = y[:,2]

# first plots

plt.plot(t, u)
plt.title("Current velocity against time ({} years)".format(t_end))
plt.xlabel("Time t (years)")
plt.ylabel("Current Velocity u ") # 10^3km / year  units?
plt.ylim((-400,400))
plt.grid(axis = 'y')
#plt.show()

"""
plt.plot(t, Te - Tw)
plt.title("Difference in Temp against time ({} years)".format(t_end))
plt.xlabel("Time t (years)")
plt.ylabel("T_e - T_w ")
plt.ylim((-30,30))
plt.show()
"""

# function for taking the derivative of a function (symmetric derivative)
def diff(func, x, dx): 
    dfs = np.zeros(len(x) - 2)
    for i in range(1, len(x) - 1):
        dfs[i-1] = (func[i + 1] - func[i - 1]) / (2 * dx)
    return dfs

def make_func(list_func, x):
    return list_func[int(round(x/dt))]

# fix single root finder (maybe doesnt sometimes work for decreasing derivatives)
# function for root finding (bisection method) 
def find_root(a, b, func): 
    x1 = 0
    x2 = 0
    x3 = 0
    if a * b > 0:
        x3 = "None"
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

du = diff(u, t, dt)

def func_u(x):
    return make_func(u, x) # remove the / (sigfig * dt) possibly

# idea: run the rootfinder of increments of .5 years in t to find roots, there will 
# be spots in which there is no root but rootfinder should stop after a while
# also notice that it looks like all maxima we care about are above y = 100

# this function will run the bisection method on intervals of time to find all roots and the pick out the ones we want
# currently this function below consistently decides to miss out on 2 to 6 roots for some reason and they are allways in the end
# funny fix idea: run programm for more than you need and then just look at the values you care about
def iterate_find_root(dx ,start_step, end_step, increment, min_y, num_iterations):
    roots = []
    for j in range(num_iterations):
        a = start_step + (j*increment)/num_iterations
        b = a + increment
        for i in range(int(end_step / increment)):
            root = find_root(a, b, dx)
            a += increment
            if type(root) == float: 
                if func_u(root) > min_y:
                    roots.append(root)

    # removes entries in roots which are very close to eachother in and makes a new set
    if num_iterations > 0:
        roots.sort()
        good_roots = roots
        for i, root_1 in enumerate(roots):
            for root_2 in roots[i+1:]:
                if root_2 - root_1 < .1:
                    good_roots.remove(root_2)

    return good_roots

roots = iterate_find_root(du, t_start, t_end, .01 , 100, 2)
roots.sort()
#print(roots)

for i in range(len(roots)):
    plt.scatter(roots[i], func_u(roots[i]))
   
# code for single root finder

"""rot = find_root(83, 84, du)
print(rot)
print(func_u(rot))
plt.scatter(rot, func_u(rot))
"""

# number of El-Nino events in the time run
num_of_elNinos = len(roots) 
print("The number of ENSO events in {} years is: {}".format(t_end, num_of_elNinos))

# calculating the periods of time that elapses between the El-Nino events
times_between_ENSO = [] # periods of time
for i in range(num_of_elNinos - 1):
    times_between_ENSO.append(roots[i+1] - roots[i]) # roots is sorted so dont need abs() function

# calculates the average of a list of numbers
def mean(mylist):
    list_lenght = len(mylist)
    my_sum = 0
    for i in range(list_lenght):
        my_sum += mylist[i]
    return my_sum / list_lenght

# calculates the mean of a list of numbers
def standard_dev(mylist):
    list_length = len(mylist)
    my_sum = 0
    average = mean(mylist)
    for i in range(list_length):
        my_sum += (mylist[i] - average)**2
    return np.sqrt(my_sum/list_length)

# to ignore the first 10 events repace times_between_ENSO with times_between_ENSO[11:]
average_T = mean(times_between_ENSO)
dev_T = standard_dev(times_between_ENSO)
print("The mean time between ENSO events is: {} years".format(average_T))
print("The standard deviation of ENSO events is: {}".format(dev_T))

# this will plot the derivative of u
#plt.plot(t[1:-1], du, "--", label = "dfs")
plt.show()

# histgram plot # make nicer
"""
plt.hist(times_between_ENSO, bins= 30)
plt.ylabel("Number of ENSO events")
plt.xlabel("Time beteen ENSO events")
plt.show()
"""

# butterfly
"""plt.plot(Te - Tw, u)
plt.show()
"""

"""
The number of ENSO events in 1000 years is: 248
The mean time between ENSO events is: 3.9895605733842947 years
The standard deviation of ENSO events is: 2.0507974941721994
"""
"""
The number of ENSO events in 10000 years is: 2424
The mean time between ENSO events is: 4.115222184869135 years
The standard deviation of ENSO events is: 2.564259102078737
"""
