from __future__ import division
from sympy import atan
import numpy as np
import sympy as sm
from scipy.optimize import fsolve
from termcolor import cprint
import sys
from math import isnan
import matplotlib.pyplot as plt

phi, psi, r = sm.var('phi, psi, r')


def equations(vec):
    (phi, psi) = vec
    return num_htt_x(phi, psi), num_htt_y(phi, psi)

def determinant(point):
    return num_htt_xx(point[0], point[1])*num_htt_yy(point[0], point[1]) - num_htt_xy(point[0], point[1])**2
    
zn = 0.04 + 0.06*(atan((-psi*phi)/(phi**2 + (1-r)**2-(psi**2)/4))/100)**2
zr = 0.04 + 0.06*(atan((psi*phi)/(phi**2 + r**2-(psi**2)/4))/100)**2
htt = 1/(1+ (zr*(r + psi/2)**2+zn*(1-r+psi/2)**2+(zr+zn)*phi**2)/(2*psi))



# Compute the 1st, 2nd derivatives regarding phi, psi

htt_x = sm.diff(htt, phi)
htt_y = sm.diff(htt, psi)
htt_xx = sm.diff(htt_x, phi)
htt_yy = sm.diff(htt_y, psi)
htt_xy = sm.diff(htt_y, phi)

"""
*** Method 2 ***
Solve the non-linear system htt_x = 0 & htt_y = 0 for phi, psi
"""
phi_range =np.arange(0.1, 1, 0.1)
guess_range = np.arange(0.1, 1, 0.05)

phi_range_guesses = np.linspace(1e-3, 1, 5)
psi_range_guesses = np.linspace(1e-3, 2, 5)
r_range = [13/24., 0.67]

for r_num in r_range:
    #cprint("r = {0}".format(round(r_num, 5)), "red")
    zeros_1 = []
    
    htt_wo_r = htt.subs(r, r_num)
    htt_x_wo_r = htt_x.subs(r, r_num)
    htt_y_wo_r = htt_y.subs(r, r_num)
    htt_xx_wo_r = htt_xx.subs(r, r_num)
    htt_yy_wo_r = htt_yy.subs(r, r_num)
    htt_xy_wo_r = htt_xy.subs(r, r_num)
    
    # Scipy expressions 
    num_htt = sm.lambdify((phi, psi), htt_wo_r, "numpy")
    num_htt_x = sm.lambdify((phi, psi), htt_x_wo_r, "numpy")
    num_htt_y = sm.lambdify((phi, psi), htt_y_wo_r, "numpy")
    num_htt_xx = sm.lambdify((phi, psi), htt_xx_wo_r, "numpy")
    num_htt_yy = sm.lambdify((phi, psi), htt_yy_wo_r, "numpy")
    num_htt_xy = sm.lambdify((phi, psi), htt_xy_wo_r, "numpy")
    
    # Solve the non-linear system [num_htt_x, num_htt_y]
    for i in phi_range_guesses:
        for j in psi_range_guesses:
            try:
                zeros_1.append(fsolve(equations, (i, j)))
            except:
                print "Unexpected error:", sys.exc_info()[0]  
                          
    # find the ones which satisfy the second derivatives criteria
    
    topika_megista = []
    asymvata = []
    overflowing = []
    for krisimo in zeros_1:
        x = num_htt_xx(krisimo[0], krisimo[1])
        if isnan(x) == False:
            if x < 0:
                if determinant(krisimo)>0:
                    if krisimo[0] > 0 and krisimo[1] > 0:
                        topika_megista.append(krisimo)
                    else:
                        asymvata.append(krisimo)
        
        else:
            overflowing.append(krisimo)

    """
    *** Drawable 'Phi-Psi' points
    find some proper (turbomachinery) points for Phi, Psi
    Stroggylopoiw ton vathmo apodoshs sta 3 dekadika. 
    kai twra ypologizw shmeia auths ths kampylhs ta opoia mporw na valw panw
    sto diagramma
    """
    
    max_vathmos = round(num_htt(topika_megista[0][0], topika_megista[0][1]), 3)
    points = []
    cprint("htt={}".format(max_vathmos), "blue")
    for phi_new in phi_range:
        for guess in guess_range:
            try:
                psi_new = fsolve(lambda y: num_htt(round(phi_new, 5), y)- max_vathmos, guess)
                cprint("num_htt(phi_new, psi_new):: {0}, phi_new:: {1}, psi_new: {2}".format(num_htt(phi_new, psi_new[0]), phi_new, psi_new[0]), "green")
                points.append((round(phi_new, 5), psi_new[0]))
                break
            except:
                pass
    points_1 = np.array(points)
    plt.figure(r_range.index(r_num) + 1)
    plt.plot(points_1[:, 0], points_1[:, 1], 'g', label='htt_max~{}'.format(max_vathmos))
    plt.legend(loc='upper left')
    plt.show()


