#!/usr/bin/env python

# Wed Feb 26 20:08:19 EET 2014, nickkouk
# Problem 1, Question 3

"""
problem1_3
Calculation of velocities and angles for 2 fixed pairs of phi-psi values 
for htt = 0.95, 0.90 respectively
for the calculations the equations of p. 6 of the course's handbook
"""
# Imports first
from __future__ import division
import numpy as np
from numpy import pi

# function initialization
V1_U = lambda x, y, r: np.sqrt(x**2 + (1 - r - y / 2) ** 2)
V2_U = lambda x, y, r: np.sqrt(x**2 + (1 - r + y / 2) ** 2)
W3_U = lambda x, y, r: np.sqrt(x ** 2 + (r + y / 2) ** 2)
W2_U = lambda x, y, r: np.sqrt(x ** 2 + (r - y / 2) ** 2)

a1 = lambda x, y, r: np.arctan(1 / x  * (1 - r - y / 2))
a2 = lambda x, y, r: np.arctan(1 / x  * (1 - r + y / 2))
b2 = lambda x, y, r: np.arctan(1 / x  * (- r + y / 2))
b3 = lambda x, y, r: np.arctan(1 / x  * (- r - y / 2))

# Data input
r = 13 / 24.
pair_1 = [0.2, 0.481, r,  0.95]
pair_2 = [0.4, 0.22, r,  0.90]
data = [pair_1, pair_2]

# printing out the data in a nice form

print "r = {}".format(np.round(r, 4))
for pair in data:
    print "\n" + "*" * 30
    print "htt = {h}\nPair: Phi = {phi},\
                            Psi = {psi}".format(phi = pair[0], 
                                                psi = pair[1], 
                                                h = pair[3])
    print "\n" + "*" * 30
    print "V1/U = {0}\tV2/U = {1}\nW3/U = {2}\
           \tW2/U = {3}".format(V1_U(pair[0], pair[1], pair[2]), 
                                V2_U(pair[0], pair[1], pair[2]),
                                W3_U(pair[0], pair[1], pair[2]),
                                W2_U(pair[0], pair[1], pair[2]))
    print "\n" + "*" * 30
    print "In Radians"
    print "a1 ={0}\ta2 = {1}\nb2 = {2}\
           \tb3 = {3}".format(a1(pair[0], pair[1], pair[2]), 
                              a2(pair[0], pair[1], pair[2]), 
                              b2(pair[0], pair[1], pair[2]), 
                              b3(pair[0], pair[1], pair[2]))
    print "\n" + "*" * 30
    print "In Degrees"
    print "a1 = {0}\ta2 = {1}\nb2 = {2}\
           \tb3 = {3}".format(180 / pi * a1(pair[0], pair[1], pair[2]), 
                              180 / pi * a2(pair[0], pair[1], pair[2]), 
                              180 / pi * b2(pair[0], pair[1], pair[2]), 
                              180 / pi * b3(pair[0], pair[1], pair[2]))
