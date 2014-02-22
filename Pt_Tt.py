#!/usr/bin/env python

# Mon Feb 17 21:42:45 EET 2014 nickkouk

"""
function for computing Pressures and Temperatures in every blade row 
of the 7-stage compressor
"""

# Imports first
from __future__ import division
import numpy as np
from history_1 import DTt
from synarthseis import *

# Initializations
Pt = np.zeros((7, 3))
Tt = np.zeros((7, 3))
Pt37 = 4.24
tol = 0.01 # Given tolerance
n = 2.1 # Initial guess of the polytropic exponent
neo = 0
step = 0.5
diafora = 100

# Python has 0-based indexing
# first stage => 0th row

# Tt values
Tt[0, 0] = 300 
Tt[0, 1] = Tt[0, 2] = Tt[0, 0] + DTt

for vathm in range(1, 7):
    Tt[vathm, 0] = Tt[vathm - 1, 2]
    Tt[vathm, 2] = Tt[vathm, 1] = Tt[vathm, 0] + DTt

# Pentry value
Pt[0, 0]  = 1.


while tol <= diafora:
    if Pt[vathm, 2] > Pt37: # greater value than the wanted one, increase n
        n+= step
    else:
        n -= step
    for vathm in range(0, 7):
        Pt[vathm, 1] = Pt[vathm, 0] * (Tt[vathm, 1] / Tt[vathm, 0])**(n/(n-1))
        Pt[vathm, 2] = (Pt[vathm, 1] + 0.12 * Pt[vathm, 0]) / 1.12
        if vathm != 6:
            Pt[vathm + 1, 0] = Pt[vathm, 2]
        else:
            diafora = np.abs(Pt[vathm, 2] - Pt37)
    
    step *= 0.8

Pt = Pt * 1.013 * 10 ** 5 # in Pa 

print "\n" + "*" * 30 + "\n"
print "Finding the polytropic exponent of every Rotor:\nDesired tolerance = {0}\nCurrent tolerance = {1}\nPt37 = {2}\nn = {3}".format(tol, diafora, Pt[vathm, 2], n)

ht = Tt * Cp

V = np.zeros((7,3))
for vathm in range(7):

    # Arithmetic Computation of the velocity v in the entrance point of every stage
    for v_guess in range(100, 700, 5):
        try :
            v_calc = newton(lambda x : v_func(x, Tt[vathm, 0], Pt[vathm, 0], 
                                              a[vathm, 0], area(vathm, 0)), 
                            v_guess)
            #if vathm != 0 and v_calc <= V[vathm - 1, 0]:
                #print "Can't have the velocity V decrease in a compressor"
                #continue
            break
        except RuntimeError as e:
            pass
            #print "RuntimeError: {0}".format(e.message)
        except:
            print "Unexpected error:", sys.exc_info()[0]
            raise 
    V[vathm, 0] = v_calc
    if vathm != 0:
        V[vathm - 1, 2] = V[vathm, 0]

#print "V = \n{}\n".format(V)

Va = np.zeros((7,3))
Vu = np.zeros((7,3))

for vathm in range(7):
    Va[vathm, 0] = V[vathm, 0] * np.cos(a[vathm, 0])
    Vu[vathm, 0] = V[vathm, 0] * np.sin(a[vathm, 0])
    if vathm != 0:
       Va[vathm - 1, 2] = Va[vathm, 0]
       Vu[vathm -1, 2] = Vu[vathm, 0]


# Euler 's Theorem for every rotor ==> Vu2(vathm)

for vathm in range(7):
    Vu[vathm, 1] = (U[vathm, 0] * Vu[vathm, 0] + (ht[vathm, 1] - ht[vathm, 0])) / U[vathm, 1]

# Compute the Velocity V for every stage after the rotor
for vathm in range(7):
    for v_guess in range(int(Vu[vathm, 1]) + 1, 1300, 5): # V larger than Vu
        try :
            v_calc1 = newton(lambda v:v_func_2(v, Tt[vathm, 1], Pt[vathm, 1], area(vathm, 1), Vu[vathm, 1]), v_guess)
            break           
        except RuntimeError as e:
            pass
            print "RuntimeError: {0}".format(e.message)
        except:
            print "Unexpected error:", sys.exc_info()[0]
    
    V[vathm, 1] = v_calc1
    a[vathm, 1] = np.arcsin(Vu[vathm, 1] / v_guess)

# Finishing the Vu, Va matrices
for vathm in range(7):
    Va[vathm, 1] = V[vathm, 1] * np.cos(a[vathm, 1])

T_stat = Tt - 1 / (2 * Cp) * V**2
P_stat = Pt * (T_stat / Tt) ** (gamma / (gamma - 1))

# W velocities
Wa = Va
Wu = Vu - U
W = np.sqrt(Wa ** 2 + Wu ** 2)
beta = np.arctan(Wu / Wa)

# Energy Loss due to Entropy
DS = np.zeros((7, 2))
for vathm in range(7):
    DS[vathm, 0] = Cp * np.log(T_stat[vathm, 1] / T_stat[vathm, 0]) - R * np.log(P_stat[vathm, 1] / P_stat[vathm, 0])
    DS[vathm, 1] = Cp * np.log(T_stat[vathm, 2] / T_stat[vathm, 1]) - R * np.log(P_stat[vathm, 2] / P_stat[vathm, 1])

# Polytropic & Isentropic Process Computations for every stage
polytropic_vathm = np.zeros((7, 2))
for vathm in range(7):
    kati = np.log(Pt[vathm, 2] / Pt[vathm, 0]) / np.log(Tt[vathm, 2] / Tt[vathm, 0])
    polytropic_vathm[vathm, 0] = -kati / (1 - kati) # Exponent of the polytropic process
    n = polytropic_vathm[vathm, 0]
    polytropic_vathm[vathm, 1] = ((gamma - 1) * n) / (gamma * (n - 1))

isentropic_vathm = np.zeros(7)
for vathm in range(7):
    pc = Pt[vathm, 2] / Pt[vathm, 0]
    n = polytropic_vathm[vathm, 0]
    isentropic_vathm[vathm] = (pc ** ((gamma - 1) / gamma) - 1) / (pc **  ((n - 1) / n) - 1)


print "\n" + "*" * 30 + "\n"
print "Angles:\na =\n{0}\nb =\n{1}".format(a, beta) 
print "\n" + "*" * 30 + "\n"
print "Velocities:\nV =\n{0}\nVa =\n{1}\nVu =\n{2}".format(V, Va, Vu)
print "W =\n{0}\nWa =\n{1}\nWu =\n{2}".format(W, Wa, Wu)
print "U =\n{0}".format(U)
print "\n" + "*" * 30 + "\n"
print "Thermodynamic Properties:\nTt =\n{0}\nPt =\n{1}".format(Tt, Pt)
print "T_stat =\n{0}\nP_stat =\n{1}".format(T_stat, P_stat)
print "\n" + "*" * 30 + "\n"
print "Entropy Losses= \n{}".format(DS)
print "\n" + "*" * 30 + "\n"
print "Polytropic Process (per stage):\nExponential of process (n) = \n{0}\nPolytropic Process Efficiency =\n{1}".format(polytropic_vathm[:, 0], polytropic_vathm[:, 1]) 
print "\n" + "*" * 30 + "\n"
print "Isentropic Process Efficiency (per stage):\n{0}\n".format(isentropic_vathm) 
