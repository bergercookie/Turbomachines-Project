#!/usr/bin/env python

from __future__ import division
from scipy.optimize import newton, fsolve
import numpy as np
from numpy import pi
from history_1 import DTt


m = 18.98
R = 287.
Tt = 300.
Pt = 1.013 * 10 ** 5
Rf = 11.3 * 10 ** -2
Rh = 22.5 * 10 ** -2
Drf = +0.26
Drh = -0.25

Area = pi * (Rh ** 2 - Rf ** 2)
l = 1.
a = (1 + l)*pi / 180.
Cp =  1004.5
gamma = 1.4
N = 1529.
w = N * pi/30.
Rm = (Rf + Rh) / 2.
Rm2 = (Rf + Drf + Rh + Drh) / 2.
G = m*R*Tt/(Pt*Area*np.cos(a))
v_func = lambda v: v - np.sqrt(2.*Cp*Tt*(1 - (v/G)**(1-gamma)))
V = newton(v_func, 150)
Vu1 = V * np.sin(a)
Va1 = V * np.cos(a)
U1 = w * Rm
U2 = w * Rm2
Wu1 = Vu1 - U1
Wa1 = Va1
Vu2 = (+Cp * DTt + U1 * Vu1) / U2
Wu2 = Vu2 - U2
