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

h = p_0 / rho / g  # [m] Difference in water level inside and outside the air cushion

h_b = h_hull + h - d  # [m] height of the air cushion after the water level inside the cushion is suppressed

A_b = L * b_c  # [m] Area of the air cushion

volume_displacement = 2 * L * d * b_G  # [m^3] volume displacement without excess air cushion pressure

lift_hull = volume_displacement * rho * g  # [N] hydrostatic force from the side hulls

lift_air_cushion = p_0 * A_b  # [N]  lift force provided from the air cushion

M = lift_hull / g + lift_air_cushion / g  # [kg] total mass of the vessel

decimal_precision = 2

air_cushion_volume = A_b * (h_hull - d)  # [m^3] volume of the air cushion

print('h =', round(h, 2), '[m]')
print('Displacement\t\t\t\t', round(M / 1000, decimal_precision), '[tonnes]', sep='\t')
print('Displacement of the side hulls', round(lift_hull / g / 1000, decimal_precision), '[tonnes]', sep='\t')
print('Lift from the air cushion\t', round(lift_air_cushion / g / 1000, decimal_precision), '[tonnes]', sep='\t')
print('Volume of the air cushion\t', round(air_cushion_volume, decimal_precision), '[m^3]', sep='\t')

print('Lift fraction provided by the air cushion: ', round(lift_air_cushion / (M * g) * 100, decimal_precision), '%', sep="")
