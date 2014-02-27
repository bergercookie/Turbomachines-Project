# REMINDER: if this is to work you have to integrate it into thema_1.py module

"""
*** First attempt ***
Solve the non-linear system htt_x = 0 & htt_y = 0 for phi, psi
Check Fermat's Theorem to see if it's a local max
"""
guesses = 10
phi_range =np.arange(0.1, 0.5, 0.05)
psi_range = np.arange(0.1, 1, 0.05)
# Solve the non-linear system [htt_x, htt_y]

zeros_1 = []
for psi_1 in psi_guesses:
    for phi_1 in phi_guesses:
        try:
            result = newton_krylov(equations, (phi_1, psi_1))
            if result[0] <= 5 or result[0] >= 0.01:
                zeros_1.append((result[0], result[1]))
        except:
            pass

# Second Derivative Criterium
print "zeros =\n{}".format(zeros_1)
topika_megista = []
asymvata = []
overflowing = []
for krisimo in zeros_1:
    x = htt_phi2(krisimo[0], krisimo[1])
    if not np.isnan(np.sum(x)):
        if x < 0:
            if determinant(krisimo)>0:
                # Both Phi, psi must be over zero
                if krisimo[0] > 0 and krisimo[1] > 0:
                    topika_megista.append(krisimo)
                else:
                    asymvata.append(krisimo)
    
    else:
        overflowing.append(krisimo)

topika_megista = np.array(topika_megista)

"""
*** Second attempt ***
Interpolate through the available from the numerical solutions 
to get a decent curve
"""
# interpolate the data from topika_megista to get the curve

def func(x, a, b, c):
    return a + b*x + c*x*x

x0    = np.array([0.0, 0.0, 0.0])
x_data = topika_megista[:, 0]
y_data = topika_megista[:, 1]
a, b, c =  curve_fit(func, x_data, y_data, x0)[0]

x_points = np.linspace(0, 0.5, 100)
y_points = a * x_points ** 2 + b * x_points + c
max_points = zip(x_points, y_points)

print max_points


