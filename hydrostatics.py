"""
This script calculates hydrostatics of a conceptual SES with rectangular side hulls.
"""

# Main parameters of conceptual SES with rectangular side hulls
L = 20  # [m] length of the vessel
B = 10  # [m] total beam of the vessel
b_G = 1.5  # [m] beam of side hulls (gondolas)
b_c = 7  # [m] beam of the air cushion/beam between the side hulls
h_hull = 2  # [m] distance from bottom of side hull to roof of air cushion
d = 0.5  # [m] draught without excess pressure in the air cushion
p_0 = 3500  # [Pa] excess pressure in the air cushion

# physical parameters
rho = 1025  # [kg/m^3] density of sea water
p_a = 101325  # [Pa] atmospheric pressure
g = 9.81  # [m/s^2] acceleration of gravity


volume_displacement = 2 * L * d * b_G  # [m^3] volume displacement without excess air cushion pressure
M = volume_displacement * rho  # [kg] total mass of the vessel

h = p_0/rho/g  # [m] Difference in water level inside and outside the air cushion

A_b = L * b_c  # [m] Area of the air cushion

V_b = M / rho - p_0 * A_b / rho / g  # [m^3] hull volume under the mean free surface with a cushion excess pressure.
M_hulls = V_b * rho  # [kg] mass displacement of the side hulls

decimal_precision = 2

air_cushion_volume = A_b * h  # [m^3] volume of the air cushion # TODO: need to check if this is the correct definition

print('h =', h, '[m]')
print('Displacement\t\t\t\t', round(M / 1000, decimal_precision), '[tonnes]', sep='\t')
print('Displacement of the side hulls', round(M_hulls / 1000, decimal_precision), '[tonnes]', sep='\t')
print('Lift from the air cushion\t', round(p_0 * A_b / g / 1000, decimal_precision), '[tonnes]', sep='\t')
print('Volume of the air cushion\t', round(air_cushion_volume, decimal_precision), '[m^3]', sep='\t')
