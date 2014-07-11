#!/usr/bin/env python

# Sun Feb 23 15:59:44 EET 2014, nickkouk

# Initializations first
from exceptions import *
from sympy import atan
import numpy as np
import sympy as sm
from scipy.optimize import newton, fsolve, newton_krylov, curve_fit
from termcolor import cprint
import sys
from math import isnan
import matplotlib.pyplot as plt
import sys



# Needed values
#ON = 13
ON = 1
r = ON / 24.

# Symbolic Manipulations
phi, psi = sm.symbols('phi psi')

# function definitions
"""
Breaking it into pieces, 
defining en, er, zn, zr as different functions 
and then summing them up to create teh efficiency rate
"""
en = atan(-psi * phi / (phi ** 2 + (1 - r) ** 2 - (psi ** 2) / 4))
er = atan(psi * phi / (phi ** 2 + r ** 2 - psi ** 2 / 4))

zn = 0.04 + 0.06 * (en / 100.) ** 2
zr = 0.04 + 0.06 * (er / 100.) ** 2

sm_htt = 1 / (1  + (phi ** 2 * (zr + zn) + zr * (r + psi / 2) ** 2 + zn * (1 - r + psi ** 2 / 2)) / (2 * psi))

# Derivative handling
sm_htt_phi = sm.diff(sm_htt, phi)
sm_htt_phi2 = sm.diff(sm_htt, phi, 2)
sm_htt_psi = sm.diff(sm_htt, psi)
sm_htt_psi2 = sm.diff(sm_htt, psi, 2)
sm_htt_phipsi = sm.diff(sm_htt, phi, psi)

# Turn into numerical expressions
# numpy label, otherwise atan will not be translated to arctan
htt = sm.lambdify((phi, psi), sm_htt, "numpy")
htt_phi = sm.lambdify((phi, psi), sm_htt_phi, "numpy")
htt_phi2 = sm.lambdify((phi, psi), sm_htt_phi2, "numpy")
htt_psi = sm.lambdify((phi, psi), sm_htt_psi, "numpy")
htt_psi2 = sm.lambdify((phi, psi), sm_htt_psi2, "numpy")
htt_phipsi = sm.lambdify((phi, psi), sm_htt_phipsi, "numpy")

# Problem 1

#Constant_initialization
h_range = (0.900, 0.925, 0.950)
color = ['red', 'blue', 'purple', 'black']
phi_range =np.arange(0.1, 0.5, 0.05)
guess_range = np.delete(np.linspace(0, 2, 100), 0)   # needs a damn good range that newton thing!


# Compute the required points
points_total = []
for h in h_range: # Looping for each one of the wanted efficiency rates
    points_1 = []
    for phi_new in phi_range: # for every phi_value
        for guess in guess_range:
            try:
                # Computing the new psi-value
                psi_new =  newton(lambda y: htt(round(phi_new, 5), y)- h, guess)
                if psi_new <= 5 and psi_new >= 0.1:
                    points_1.append((round(phi_new, 5), psi_new, h))
                    break
            except RuntimeError:
                pass
    # finally, this is a valid point after all, so add it up to the total_point array
    points_total.append(points_1)


points_total = np.array(points_total)

# Problem 2

def equations(vec):
    """
    function for solving the 2-equation system of the first derivatives of htt
    takes as input a vector of phi-psi points and then computes and returns 
    as another vector the first derivatives with regards to phi and psi
    """
    (phi, psi) = vec
    return htt_phi(phi, psi), htt_psi(phi, psi)

# Compute the required points for the maximum htt computed above
"""
Compute the htt_max for one solution of the first-derivatives system of equations
and then proceed to the previous method selecting phi and computing psi in order to get the needed
(phi, psi) pairs
"""
result = fsolve(equations, (0.4, 0.4))
h = np.round(htt(result[0], result[1]))

phi_range =np.arange(0.1, 0.5, 0.05)
guess_range = np.delete(np.linspace(0.1, 2, 100), 0)
max_points = []
for phi_new in phi_range:
    for guess in guess_range:
        try:
            psi_new = fsolve(lambda y: htt(round(phi_new, 5), y)- h, guess)
            if psi_new <= 5 and psi_new >= 0.1:
                max_points.append((round(phi_new, 5), psi_new[0], np.round(h, 4)))
                break
        except RuntimeError:
            pass

# turning it into a numpy_array
max_points = np.array(max_points)

# Print out the points 
plt.hold('True')
plt.grid('on')
print "Problem 1"
for one_h in points_total:
    cprint("htt={}".format(one_h[0, 2]), "blue")
    for one_row in one_h:
        print "\t{0} ==> {1}".format(one_row[0], np.round(one_row[1], 4))

print "Problem 2"
print "htt=h_max"
for row in max_points:
    print "\t{0} ==> {1}".format(row[0], np.round(row[1], 4))

# Plotting everything
line_num = 0
for one_h in points_total:
    plt.plot(one_h[:, 0], one_h[:, 1], color[line_num], linewidth = 2, label='htt = {}'.format(one_h[0, 2]))
    line_num += 1
plt.plot(max_points[:, 0], max_points[:, 1], color[line_num], linewidth = 2, label = 'htt = h_max')

plt.legend(loc='upper left')
plt.title("Phi-Psi")
plt.show()
