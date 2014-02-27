#!/usr/bin/env python

# Wed Feb 26 21:12:40 EET 2014, nickkouk

"""
Problem 1, Question 3
for 2 pairs of phi-psi values for a certain isentropic efficiency rate 
draw the non-thick blades (rotor & stator)

Reminder: This is the plotting function, 
for the calculation of velocities, angles etc see problem_data.dat
"""



# Initializations first
import numpy as np
import scipy as sp
import matplotlib.pylab as plt

# first pair
stator_1 = lambda x: 1.2022189 * x ** 2 + 1.08905 * x
rotor_1  = lambda x: -1.202146 * (x - 1.5) ** 2 - 1.50555 * (x - 1.5) + 2.3 

# second pair
stator_2 = lambda x: 0.244975 * x ** 2 + 1.42078 * x
rotor_2 = lambda x: -2.1410 * (x - 1.5) ** 2 - 1.0791 * (x - 1.5) + 1.6

pairs = [(rotor_1, stator_1, (0.2, 0.481, 0.95)), (rotor_2, stator_2, (0.4, 0.22, 0.9))]

for pair in pairs:
    # range of the phi-axis
    x1 = np.linspace(0, 1, 100)
    x2 = np.linspace(1.5, 2.5, 100)

    # psi-values for stator and rotor
    y_stator = pair[1](x1)
    y_rotor = pair[0](x2)

    # Initiate the graph
    plt.figure()
    plt.hold('on')

    # setting the graph parameters
    fig = plt.gcf()
    fig.suptitle("Stator-Rotor Stage\nPhi = {0}, Psi = {1}, htt = {2}".format(pair[2][0], pair[2][1], pair[2][2]), fontsize=14)
    frame = plt.gca()
    plt.grid()

    # plotting the stator, rotor lines
    plt.plot(x1, y_stator, 'b', label = 'stator', linewidth = 2) # Plotting the stator line 
    plt.plot(x2, y_rotor, 'r', label = 'rotor', linewidth = 2) # Plotting the rotor line

    plt.legend()

# show both graphs on the screen
plt.show()
