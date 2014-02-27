from __future__ import division
from sympy import atan
import sympy as sm
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pylab as plt
from termcolor import cprint
#import sys

phi, psi, r = sm.symbols('phi psi r')

zn = 0.04 + 0.06*(atan((-psi*phi)/(phi**2 + (1-r)**2-(psi**2)/4))/100)**2
zr = 0.04 + 0.06*(atan((psi*phi)/(phi**2 + r**2-(psi**2)/4))/100)**2
htt=1/(1+ (zr*(r + psi/2)**2+zn*(1-r+psi/2)**2+(zr+zn)*phi**2)/(2*psi))

r_range = [13/24., 0.67]
h_range = (0.900, 0.925, 0.950)
color = ['red', 'blue', 'purple', 'grey']
phi_range =np.arange(0.1, 1, 0.1)
guess_range = np.arange(0.1, 1, 0.05)

for r_num in r_range:
    cprint("r = {0}".format(round(r_num, 5)), "red")
    htt_wo_r = htt.subs(r, r_num)
    num_htt = sm.lambdify((phi, psi), htt_wo_r, "numpy")
    plt.figure()
    plt.hold('True')
    for h in h_range:
        points = []
        cprint("htt={}".format(h), "blue")
        for phi_new in phi_range:
            for guess in guess_range:
                try:
                    psi_new =  fsolve(lambda y: num_htt(round(phi_new, 5), y)- h, guess)
                    points.append((round(phi_new, 5), psi_new[0]))
                    break
                except:
                    pass
        points_1 = np.array(points)
        plt.plot(points_1[:, 0], points_1[:, 1], color[int('{}'.format(h_range.index(h)))], label='htt={}'.format(h))
        
        if __name__ == '__main__':
            for point in points:
                print "{0} ==> {1}".format(point[0], round(point[1], 5))
        plt.legend(loc='upper left')
    
    # Plotting everything
    plt.title("Phi-Psi gia r={}".format(round(r_num, 2)))
    plt.show()
