import numpy as np

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

# Derived parameters
A_c = L*b  # [m^2] Air cushion area  # TODO: Might want to use expression for rectangular cushion shape with triangle at the front
x_cp = 1  # [m] longitudinal centroid of air cushion

# Read hydrodynamic coefficients
# TODO: Temporarily set to zero. Need to read them in correctly
A_33 = 0
B_33 = 0
C_33 = 0

A_35 = 0
B_35 = 0
C_35 = 0

A_53 = 0
B_53 = 0
C_53 = 0

I_55 = 0
A_55 = 0
B_55 = 0
C_55 = 0

omega_e = np.linspace(0, 10, 10)

n_freq = len(omega_e)
# Heave equation

j_max = 4  # number of acoustic modes to include in the calculations

a = np.array([[omega_e, omega_e - 10],
              [np.multiply(omega_e, 2), np.power(omega_e, 2)]])

A_mat = np.zeros([3, 3, n_freq], dtype=complex)

A_mat[0, 0, :] = -(m + A_33)*np.power(omega_e, 2) + B_33 * 1j * omega_e + C_33
A_mat[0, 1, :] = -A_35*np.power(omega_e, 2) + B_35 * 1j * omega_e + C_35
A_mat[0, 2, :] = -A_c * p_0

A_mat[1, 0, :] = -A_53*np.power(omega_e, 2) + B_53 * 1j * omega_e + C_53
A_mat[1, 1, :] = -(I_55+A_55)*np.power(omega_e, 2) + B_55 * 1j * omega_e + C_55
A_mat[1, 2, :] = A_c * p_0 * x_cp

# TODO: Fill in for equation of dynamic pressure, i.e. eq. (3.74) in 'Cobblestone effect on SES'
A_mat[2, 0, :] = 0
A_mat[2, 1, :] = 0
A_mat[2, 2, :] = 0

print(A_mat[0, 0, :])

print('Hei')
