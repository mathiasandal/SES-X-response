"""Implementation of analysis with spatially varying pressure with the same input as in p. 35 in Steen and Faltinsen (1995). 'Cobblestone
    Oscillations of an SES with Flexible Bag Aft Seal'
"""

# ***** Define input parameters *****

# Physical constants
c = 343  # [m/s] Speed of sound in air
g = 9.81  # [m/s^2] Acceleration of gravity
p_a = 101325  # [Pa] Atmospheric pressure
rho_0 = 0  # [kg/m^3] Density of air at mean cushion pressure  # TODO: What is this value?
rho_a = 1.225  # [kg/m^3] Density of air at atmospheric pressure # TODO: What is this value

# SES main dimensions
L_oa = 35  # [m] Length overall
L = 28  # [m] Air cushion length
b = 8  # [m] Air cushion beam
m = 140  # [tons] Vessel total mass
h_0 = 2.0  # [m] Cushion height
U = 50  # [knots] Velocity

# Fan characteristics
p_0 = 500  # [mmWc] Mean cushion pressure
Q_0 = 150  # [m^3/s] Mean fan flow rate
dQdp_0 = -140  # [m^2/2] Linear fan slope
lcg_fan = 5.6  # [m] Longitudinal fan position (from CG)


# Loading condition
