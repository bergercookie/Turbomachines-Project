#!usr/bin/env python

# Tue Feb 18 23:51:10 EET 2014, nickkouk

"""
The goal of this module is to compute DTt
Was designed seperately from synarthseis.py 
for simplicity reasons
"""

# Imports first
from __future__ import division
import sys
from scipy.optimize import newton

Tt11 = 300.
Cp = 1004.5

def paronomastis(n_pol):
    return logos_p**((n_pol -1)/n_pol)-1
def arithmitis(n_pol):
    return (logos_p**((gamma-1)/gamma) - 1)
n_func = lambda n_pol: hisentropic - arithmitis(n_pol) / paronomastis(n_pol)
hisentropic = 0.85
logos_p = 4.24
gamma = 1.4
try:
    n_pol = newton(n_func, 1.5)
except:
    print "Unexpected error:", sys.exc_info()[0]
    raise 
        
hpol = (gamma - 1)/(gamma*(1 - 1/n_pol))

# DTt

Tt37 = Tt11 * logos_p ** ((n_pol-1)/n_pol)
DTt = (Tt37 - Tt11) / 7.
