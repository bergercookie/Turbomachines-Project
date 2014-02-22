#!/usr/bin/env python

# Tue Feb 18 14:21:07 EET 2014 nickkouk

# Initialization
from __future__ import division
import numpy as np
import scipy as sp
from scipy.optimize import newton
from numpy import pi
import sys

# Initial values
R = 287.
Rh_in = 22.5e-2
Rf_in = 11.3e-2
Rh_out = 19.0e-2
Rf_out = 14.95e-2
Cp = 1004.5
N = 15025
w = 2 * pi * N / 60.
m = 18.98
gamma = 1.4

DRf = (Rf_out - Rf_in) / 15. # 14 blade rows => 15 positions
DRh = (Rh_out - Rh_in) / 15.

# Printing arguement
np.set_printoptions(precision=3)

# Radius changes by a constant portion after every blade row
Rm_1 = []
Rf_1 = []
Rh_1 = []
Rm = []
Rh = []
Rf = []
for i in range(14):
    Rh_1.append(Rh_in + i * DRh)
    Rf_1.append(Rf_in + i * DRf)
    Rm_1.append(np.mean((Rh_1[i], Rf_1[i])))
for i in range(7):
    Rm.append((Rm_1[2 * i], Rm_1[2 * i + 1], Rm_1[2 * i + 1])) # append the second element twice, 3(vathm) == 1(vathm + 1)
    Rh.append((Rh_1[2 * i], Rh_1[2 * i + 1], Rh_1[2 * i + 1])) # append the second element twice, 3(vathm) == 1(vathm + 1)
    Rf.append((Rf_1[2 * i], Rf_1[2 * i + 1], Rf_1[2 * i + 1])) # append the second element twice, 3(vathm) == 1(vathm + 1)
    
Rm = np.array(Rm)
Rh = np.array(Rh)
Rf = np.array(Rf)

# a-angles
a = np.zeros((7,3))
a[0,0] = 2
for vathm in range(1, 7):
    a[vathm, 0] = 1 + (1 + vathm)
    a[vathm - 1, 2 ] = a[vathm, 0] 

a = pi * a / 180
#print a

U = Rm * w

#print "w =\n", w
#print "U =\n", U
#print "Rm =\n", Rm
#print len(Rm)

# Computing the Area in a given position
def area(vathm, thesi):
    return pi * (Rh[vathm, thesi] ** 2 - Rf[vathm, thesi] ** 2)

#print area(1, 2) 

# Computing the velocity v  
def v_func(v, Tt, Pt, a, area):
    G = m*R*Tt/(Pt*area*np.cos(a))
    return v - np.sqrt(2.*Cp*Tt*(1 - (v/G)**(1-gamma)))

def v_func_2(V, Tt, Pt, area, Vu):
    equation = V * Pt * area * np.sqrt(1 - (Vu / V) ** 2 ) - m* R* Tt * (1 - V**2 / (2 * Cp * Tt)) ** (-1 / (gamma - 1))
    return equation
     

