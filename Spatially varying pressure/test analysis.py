import numpy as np
from utilities import K_1, K_2, K_3, K_4, Xi_j, Omega_j, solve_mean_value_relation, A_0_AP, A_0_FP, A_0j, A_3j, A_5j, \
    A_7j, B_0j, B_3j, B_5j, B_7j,  r_j
"""Implementation of analysis with spatially varying pressure with the same input as in p. 35 in Steen and Faltinsen (1995). 'Cobblestone
    Oscillations of an SES with Flexible Bag Aft Seal'
"""

# ***** Define input parameters *****

# Physical constants
c = 343  # [m/s] Speed of sound in air
g = 9.81  # [m/s^2] Acceleration of gravity
p_a = 101325  # [Pa] Atmospheric pressure
rho_0 = 1.0  # [kg/m^3] Density of air at mean cushion pressure  # TODO: What is this value?
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
x_F = 0  # [m]

# Other parameters
zeta_a = 1  # [m] wave amplitude
h_s_AP = 0  # [m] aft seal submergence
h_s_FP = 0  # [m] bow seal submergence
x_g_AP = -L/2  # [m] position of leakage at AP relative to center of gravity
x_g_FP = L/2  # [m] position of leakage at FP relative to center of gravity

# Derived parameters
A_c = L*b  # [m^2] Air cushion area  # TODO: Might want to use expression for rectangular cushion shape with triangle at the front
x_cp = -10  # [m] longitudinal centroid of air cushion

# ***** Read hydrodynamic coefficients *****
# TODO: Temporarily set to 10. Need to read them in correctly
A_33 = 10
B_33 = 10
C_33 = 10

A_35 = 10
B_35 = 10
C_35 = 10

A_53 = 10
B_53 = 10
C_53 = 10

I_55 = 10
A_55 = 10
B_55 = 10
C_55 = 10

omega_0 = np.linspace(1, 3, 10)
k = np.power(omega_0, 2)/g  # wave number of water waves
omega_e = omega_0 + np.power(omega_0, 2)/g*U*0.514444  # encounter frequencies
n_freq = len(omega_e)

# ***** Compute wave pumping *****
F_wp = rho_0 * A_c * np.multiply(omega_e, np.divide(np.sin(k*L/2), k*L/2)) * zeta_a

# ***** Compute constants *****
k_1 = K_1(rho_0, p_0, h_0, A_c)  # [kg] eq. (83) in Steen and Faltinsen (1995)

k_2_AP = K_2(p_0, 1.0)  # [m/s] linearized equivalent outflow velocity constant at stern
k_2_FP = K_2(p_0, 0.61)  # [m/s] linearized equivalent outflow velocity constant at bow

k_3 = K_3(rho_0, p_0, Q_0, dQdp_0)  # air cushion flow airflow constant

# ***** Start of iteration loop to fit N_R and N_B *****
# Guess an initial value for bias and gain value of non-linear leakage
n_R_AP = 0.5
n_R_FP = 0.5
n_B_AP = 0.5
n_B_FP = 0.5

# Solve mean value relation
eta_3m, eta_5m = solve_mean_value_relation(n_B_AP, n_B_FP, L, b, x_cp, A_c, p_0, k_2_AP, k_2_FP, k_3, h_s_AP, h_s_FP, C_33, C_55, C_35, C_53)
# Compute mean leakage area et AP and FP
a_0_AP = A_0_AP(L, b, n_B_AP, eta_3m, eta_5m, h_s_AP)
a_0_FP = A_0_FP(L, b, n_B_FP, eta_3m, eta_5m, h_s_FP)

j_max = 4  # number of acoustic modes to include in the calculations

A_mat = np.zeros([3, 3, n_freq], dtype=complex)  # initialize coefficient matrix for linear system of eq.
f_vec = np.zeros([3, n_freq], dtype=complex)  # initialize column vector on the right hand side of the equation

# Assign values for coefficient matrix in solving A_mat*x = f_vec, with exception of terms due to spat. varying pressure
# Heave equation, i.e. eq. (87) Steen and Faltinsen (1995)
A_mat[0, 0, :] = -(m + A_33)*np.power(omega_e, 2) + B_33 * 1j * omega_e + C_33
A_mat[0, 1, :] = -A_35*np.power(omega_e, 2) + B_35 * 1j * omega_e + C_35
A_mat[0, 2, :] = -A_c * p_0
f_vec[0, :] = 0  # TODO: to be filled in

# Pitch equation, i.e. eq. (88) Steen and Faltinsen (1995)
A_mat[1, 0, :] = -A_53*np.power(omega_e, 2) + B_53 * 1j * omega_e + C_53
A_mat[1, 1, :] = -(I_55+A_55)*np.power(omega_e, 2) + B_55 * 1j * omega_e + C_55
A_mat[1, 2, :] = A_c * p_0 * x_cp
f_vec[1, :] = 0  # TODO: to be filled in

# TODO: Fill in for equation of dynamic pressure, i.e. eq. (3.74) in 'Cobblestone effect on SES'
# Equation of dynamic uniform pressure, i.e. eq. (82) Steen and Faltinsen (1995)
A_mat[2, 0, :] = rho_a * b * (k_2_AP*n_R_AP + k_2_FP*n_R_FP) + rho_0 * A_c * 1j * omega_e
A_mat[2, 1, :] = rho_a * b * L/2 * (k_2_AP*n_R_AP - k_2_FP*n_R_FP) - rho_0 * A_c * x_cp * 1j * omega_e
A_mat[2, 2, :] = k_1 * 1j * omega_e + k_3
f_vec[2, :] = F_wp - rho_a * b * 1j * (k_2_AP*n_R_AP * np.exp(-1j * k * L/2) + k_2_FP*n_R_FP * np.exp(1j * k * L/2))


# ***** Fill in terms due to spatially varying pressure *****
for j in range(1, j_max + 1):

    # Computes variables dependent on j
    xi_j = Xi_j(j, rho_0, p_0, h_0, b, L, k_2_AP, k_2_FP, a_0_AP, a_0_FP, dQdp_0, x_F)  # [-] Relative damping ratio of mode j  # TODO check if I should use lcg_fan instead of x_F
    omega_j = Omega_j(j, L)
    k_4 = K_4(xi_j, h_0, omega_e, omega_j)

    # TODO: Compute A_0j(omega_e) and B_0j(omega_e) amplitudes beforehand(?)

    if j % 2 == 1:  # if j is odd
        # ***** Pitching moments due to spatially varying pressure in eq. (88) in Steen and Faltinsen (1995)
        # Heave DOF
        A_mat[1, 0, :] += 2*rho_0*b*1j * (L/j/np.pi)**2 * np.multiply(omega_e, A_3j(L, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP))

        # Pitch DOF
        A_mat[1, 1, :] += 2*rho_0*b*1j * (L/j/np.pi)**2 * np.multiply(omega_e, A_5j(j, L, omega_e, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP))

        # Uniform pressure DOF
        A_mat[1, 2, :] += 2*rho_0*b*1j * (L/j/np.pi)**2 * np.multiply(omega_e, A_0j(j, b, L, p_0, dQdp_0, x_F, k_2_AP, k_2_FP, a_0_AP, a_0_FP, k_4))

        # Wave excitation term
        f_vec[1, :] -= 2*rho_0*b*1j * (L/j/np.pi)**2 * np.multiply(omega_e, A_7j(j, k_4, L, k, omega_e, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * zeta_a

        # ***** Terms from leakage at AP in eq. (82) in Steen and Faltinsen (1995) *****
        # Heave DOF
        A_mat[2, 0, :] += rho_a*k_2_AP/2*a_0_AP * (-rho_0/p_0)*1j*np.multiply(omega_e, A_3j(L, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_g_AP, j, L)

        # Pitch DOF
        A_mat[2, 1, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, A_5j(j, L, omega_e, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_g_AP, j, L)

        # Uniform pressure DOF
        A_mat[2, 2, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, A_0j(j, b, L, p_0, dQdp_0, x_F, k_2_AP, k_2_FP, a_0_AP, a_0_FP, k_4)) * r_j(x_g_AP, j, L)

        # Wave excitation term
        f_vec[2, :] -= rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, A_7j(j, k_4, L, k, omega_e, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_g_AP, j, L) * zeta_a

        # ***** Terms from leakage at FP in eq. (82) in Steen and Faltinsen (1995) *****
        # Heave DOF
        A_mat[2, 0, :] += rho_a*k_2_FP/2*a_0_FP * (-rho_0/p_0)*1j*np.multiply(omega_e, A_3j(L, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_g_FP, j, L)

        # Pitch DOF
        A_mat[2, 1, :] += rho_a*k_2_FP/2*a_0_FP * (-rho_0/p_0)*1j*np.multiply(omega_e, A_5j(j, L, omega_e, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_g_FP, j, L)

        # Uniform pressure DOF
        A_mat[2, 2, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, A_0j(j, b, L, p_0, dQdp_0, x_F, k_2_AP, k_2_FP, a_0_AP, a_0_FP, k_4)) * r_j(x_g_FP, j, L)

        # Wave excitation term
        f_vec[2, :] -= rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, A_7j(j, k_4, L, k, omega_e, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_g_FP, j, L) * zeta_a

        # ***** Terms from fan inlet at x_F in eq. (82) in Steen and Faltinsen (1995) *****
        # Heave DOF
        A_mat[2, 0, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, A_3j(L, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_F, j, L)

        # Pitch DOF
        A_mat[2, 1, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, A_5j(j, L, omega_e, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_F, j, L)

        # Uniform pressure DOF
        A_mat[2, 2, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, A_0j(j, b, L, p_0, dQdp_0, x_F, k_2_AP, k_2_FP, a_0_AP, a_0_FP, k_4)) * r_j(x_F, j, L)

        # Wave excitation term
        f_vec[2, :] -= (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, A_7j(j, k_4, L, k, omega_e, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_F, j, L) * zeta_a


    elif j % 2 == 0:  # if j is even
        # ***** Terms from leakage at AP in eq. (82) in Steen and Faltinsen (1995) *****
        # Heave DOF
        A_mat[2, 0, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, B_3j(L, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_g_AP, j, L)

        # Pitch DOF
        A_mat[2, 1, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, B_5j(k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_g_AP, j, L)

        # Uniform pressure DOF
        A_mat[2, 2, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, B_0j(j, b, L, p_0, dQdp_0, x_F, k_2_AP, k_2_FP, a_0_AP, a_0_FP, k_4)) * r_j(x_g_AP, j, L)

        # Wave excitation term
        f_vec[2, :] -= rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, B_7j(j, k_4, L, k, omega_e, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_g_AP, j, L) * zeta_a

        # ***** Terms from leakage at FP in eq. (82) in Steen and Faltinsen (1995) *****
        # Heave DOF
        A_mat[2, 0, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, B_3j(L, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_g_FP, j, L)

        # Pitch DOF
        A_mat[2, 1, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, B_5j(k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_g_FP, j, L)

        # Uniform pressure DOF
        A_mat[2, 2, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, B_0j(j, b, L, p_0, dQdp_0, x_F, k_2_AP, k_2_FP, a_0_AP, a_0_FP, k_4)) * r_j(x_g_FP, j, L)

        # Wave excitation term
        f_vec[2, :] -= rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, B_7j(j, k_4, L, k, omega_e, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_g_FP, j, L) * zeta_a

        # ***** Terms from fan inlet at x_F in eq. (82) in Steen and Faltinsen (1995) *****
        # Heave DOF
        A_mat[2, 0, :] += (-rho_0*p_0*dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, B_3j(L, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_F, j, L)

        # Pitch DOF
        A_mat[2, 1, :] += (-rho_0*p_0*dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, B_5j(k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_F, j, L)

        # Uniform pressure DOF
        A_mat[2, 2, :] += (-rho_0*p_0*dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, B_0j(j, b, L, p_0, dQdp_0, x_F, k_2_AP, k_2_FP, a_0_AP, a_0_FP, k_4)) * r_j(x_F, j, L)

        # Wave excitation term
        f_vec[2, :] -= (-rho_0*p_0*dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, B_7j(j, k_4, L, k, omega_e, k_2_AP, k_2_FP, n_R_AP, n_R_FP)) * r_j(x_F, j, L) * zeta_a


print(A_mat[0, 0, :])

print('Hei')
