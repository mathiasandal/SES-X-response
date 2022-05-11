"""
This script calculates hydrostatics of a conceptual SES with rectangular side hulls.
"""


# Main parameters of conceptual SES with rectangular side hulls
L = 35  # [m] length of the vessel
L_air_cushion = 28  # [m]
b_G = 0.5  # [m] beam of side hulls (gondolas)
b_c = 8  # [m] beam of the air cushion/beam between the side hulls
h_hull = 2  # [m] distance from bottom of side hull to roof of air cushion
m = 140e3  # [kg] total mass of the vessel
p_0 = 500.  # [mmWc] Mean cushion pressure

B = 2 * b_G + b_c  # [m] total beam of the vessel

# physical parameters
rho = 1025  # [kg/m^3] density of sea water
p_a = 101325  # [Pa] atmospheric pressure
g = 9.81  # [m/s^2] acceleration of gravity

p_0 = rho*g*p_0*1e-3  # [Pa] Mean cushion pressure

h = p_0 / rho / g  # [m] Difference in water level inside and outside the air cushion


A_b = L_air_cushion * b_c  # [m] Area of the air cushion

d = (m*g - A_b*p_0)/(rho*2*L*b_G*g)  # [m] draught with excess pressure in the air cushion

h_b = h_hull + h - d  # [m] height of the air cushion after the water level inside the cushion is suppressed

print('The draught under the mean water level is ' + str(d) + '[m] with an excess pressure in the air cushion of '
      + str(p_0) + '[Pa].')