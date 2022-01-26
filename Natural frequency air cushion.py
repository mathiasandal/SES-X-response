import numpy as np
from air_cushion import air_cushion_area, interpolate_fan_characteristics, read_fan_characteristics

# Main parameters of conceptual SES with rectangular side hulls
L = 20  # [m] length of the vessel
B = 10  # [m] total beam of the vessel
b_G = 1.5  # [m] beam of side hulls (gondolas)
b_c = 7  # [m] beam of the air cushion/beam between the side hulls
h_hull = 2  # [m] distance from bottom of side hull to roof of air cushion
d = 0.5  # [m] draught without excess pressure in the air cushion
p_0 = 3500  # [Pa] excess pressure in the air cushion
gamma = 1.4  # [-] Ratio of specific heat for air

# physical parameters
rho = 1025  # [kg/m^3] density of sea water
p_a = 101325  # [Pa] atmospheric pressure
g = 9.81  # [m/s^2] acceleration of gravity


volume_displacement = 2 * L * d * b_G  # [m^3] volume displacement without excess air cushion pressure
M = volume_displacement * rho  # [kg] total mass of the vessel

h = p_0/rho/g  # [m] Difference in water level inside and outside the air cushion

h_b = h_hull + h - d  # [m] height of the air cushion after the water level inside the cushion i suppressed.

A_b = L * b_c  # [m] Area of the air cushion


# Read in fan characteristics
Q, P, rpm_dummy = read_fan_characteristics('Input files/fan characteristics/fan characteristics.csv', '1800rpm')

# Interpolates values
Q_0, dQdp_0 = interpolate_fan_characteristics(p_0, P, Q)

# Q_0 [m^3/s] Mean volumetric flow of air fan
# dQdp_0 [(m^3s^-1)/(Pa)] Slope of fan characteristic curve estimated from Figure 5.6 in Faltinsen High-Speed

'''
Equivalent mass, damping and spring coefficients
m*ddx + c*dx + k*x = F(t)
'''
m = M*h_b/gamma/(p_0 + p_a)  # [m^2 s^2] Mass equivalent coefficient
c = M/A_b/p_0 * (0.5*Q_0 - dQdp_0*p_0)  # [m^2 s] Damping equivalent coefficient
k = A_b  # [m^2] Stiffness equivalent coefficient

'''Calculating natural frequency and damping ratio'''

omega_0 = np.sqrt(k/m)  # [rad/s] Natural frequency of uniform pressure oscillations
c_cr = 2 * np.sqrt(k*m)  # [m^2 s] Critical damping of the system
zeta = c/c_cr  # [-] Damping ratio of the system

'''Print results'''

print('Natural frequency:\t', round(omega_0, 2), '\t[rad/s]')
print('Natural frequency:\t', round(omega_0/2/np.pi, 2), '\t[Hz]')
print('Damping ratio:\t\t', round(zeta, 3), '\t[-]')

