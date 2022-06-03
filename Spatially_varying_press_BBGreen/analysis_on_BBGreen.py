import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from air_cushion import air_cushion_area, read_fan_characteristics, interpolate_fan_characteristics, \
    wave_pumping_excitation_sesx
plt.rcParams['text.usetex'] = True
from veres import read_coefficients_from_veres
from utilitiesBBGreen import K_1, K_2, K_3, K_4, Xi_j, Omega_j, solve_mean_value_relation, A_0_AP, A_0_FP, A_0j, A_3j, A_5j, \
    A_7j, B_0j, B_3j, B_5j, B_7j, r_j, solve_linear_systems_of_eq, N_R, N_B, rms_leakage, Zeta_a, \
    append_spatially_varying_terms, PM_spectrum
from Spatially_varying_pressure.HyCoef import compute_hydrodynamic_coeff

"""Implementation of analysis with spatially varying pressure with the same input as in p. 35 in Steen and Faltinsen (1995). 'Cobblestone
    Oscillations of an SES with Flexible Bag Aft Seal'
"""

# ***** Define input parameters *****

# Wave parameters
H_s = 0.15  # [m] significant wave height
T_p = 1.5  # [s] wave peak period

# Physical constants
c = 343.  # [m/s] Speed of sound in air
g = 9.81  # [m/s^2] Acceleration of gravity
p_a = 101325.  # [Pa] Atmospheric pressure
# rho_0 = 1.2796  # [kg/m^3] Density of air at mean cushion pressure (p_0 + p_a)
# source: https://www.gribble.org/cycling/air_density.html at 15deg C, p_0 + p_a

rho_a = 1.225  # [kg/m^3] Density of air at atmospheric pressure
rho_w = 1025.  # [kg/m^3] Density of salt water
rho = 1000.  # [kg/m^3] Density of fresh water
gamma = 1.4  # [-] Ratio of specific heat of air at adiabatic condition

# main dimensions of BBGreen
B = 6  # [m] beam of BBGreen
Lpp = 19.2  # [m] L_pp of BBGreen
L = 18  # [m] Air cushion length
m = 25.6e3  # [kg] total mass of the vessel
r44 = 0.35 * B  # [m] radii of gyration in roll
r55 = 0.25 * Lpp  # [m] radii of gyration in pitch
r66 = 0.27 * Lpp  # [m] radii of gyration in yaw
r46 = 0  # [m] radii of gyration for inertial coupling of yaw and roll
lcg = 7.83  # [m] longitudinal center of gravity relative to AP
vcg = 1.98  # [m] vertical center of gravity relative to BL

# Properties of SES-X air cushion
l_1 = 12   # 0.0001  #     [m] length of the rectangular part of the air cushion
l_2 = 6  # 0.001  #      [m] length of the triangular part of the air cushion
b = 3.4  # [m] beam of the air cushion
p_0 = 3500  # [Pa] mean excess pressure in the air cushion
h_0 = 0.64  #2 # [m]  Height of the air cushion  # TODO: make sure this value is correct
h = p_0 / rho / g  # [m] Difference in water level inside and outside the air cushion

x_F_c = 17.7  # [m]  Distance from AP to fan inlet  # Measured in Rhino7

# Calculate initial density of air in air cushion
rho_0 = rho_a * ((p_0 + p_a) / p_a) ** (1 / gamma)  # [kg/m^3] Density of air at mean cushion pressure (p_0 + p_a)


# Compute area and longitudinal position of centroid relative to AP
A_c0, x_B_c = air_cushion_area(l_1, l_2, b)

V_c0 = A_c0 * h_0  # [m^3] Compute mean air cushion volume  # TODO: Make sure this is correct

# Read in fan characteristics
Q, P, rpm_dummy = read_fan_characteristics('C:/Users/mathi/code/repos/SES-X-response/Input files/fan characteristics/fan characteristics.csv', '1800rpm')

# Interpolates values
Q_0, dQdp_0 = interpolate_fan_characteristics(p_0, P, Q)

# ***** Read hydrodynamic coefficients for conceptual SES *****

# Read input.re7 file to get hydrodynamic coefficients

veres_formulation = 'high-speed'  # 'strip-theory'  #
if veres_formulation == 'high-speed':
    # For high-speed formulation
    path_high_speed_theory = 'C:/Users/mathi/SIMA Workspaces/Workspace_1/Task/BBGreen_uniform_pressure_model_hs_theory/DWL'
    A_temp, B_temp, C_temp, f_ex_temp, omega_0, U, beta, XMTN, ZMTN, lcg_veres, vcg_veres = read_coefficients_from_veres(path_high_speed_theory)
else:
    # For strip_theory
    path_strip_theory = 'C:/Users/mathi/SIMA Workspaces/Workspace_1/Task/BBGreen_uniform_pressure_model_strip_theory/DWL'
    A_temp, B_temp, C_temp, f_ex_temp, omega_0, U, beta, XMTN, ZMTN, lcg_veres, vcg_veres = read_coefficients_from_veres(path_strip_theory)

# Add terms from VERES
# Added mass
A_33 = A_temp[:, 2, 2]  # A_33
A_35 = A_temp[:, 2, 4]  # A_35
A_53 = A_temp[:, 4, 2]  # A_53
A_55 = A_temp[:, 4, 4]  # A_55
# Damping
B_33 = B_temp[:, 2, 2]  # B_33
B_35 = B_temp[:, 2, 4]  # B_35
B_53 = B_temp[:, 4, 2]  # B_53
B_55 = B_temp[:, 4, 4]  # B_55
# Restoring
C_33 = C_temp[2, 2]  # C_33
C_35 = C_temp[2, 4]  # C_35
C_53 = C_temp[4, 2]  # C_53
C_55 = C_temp[4, 4]  # C_55
# Adjust stiffness in pitch due to change in center of gravity relative to VERES
C_55 += m * g * (vcg_veres - vcg)  # C_55

I_55 = m * r55**2  # [kg m^2]

# Center of motion relative to BL and AP
x_prime = Lpp/2 - XMTN[0]
z_prime = ZMTN[0]

# Create mass matrix
x_G = x_prime - lcg  # [m] longitudinal position of COG relative to motion coordinate system
z_G = -(z_prime - vcg)  # [m] vertical position of COG relative to motion coordinate system
M = np.zeros([3, 3])
M[0, 0] = m  # [kg] M_33
M[0, 1] = -x_G * m  # [kg m] M_35
M[1, 0] = -x_G * m  # [kg m] M_53
M[1, 1] = r55**2 * m  # [kg m^2] I_55


# Other parameters
h_s_AP = 0.1  # [m] aft seal submergence
h_s_FP = 0.1  # [m] bow seal submergence
# Centroid of the air cushion at equilibrium relative to the motion coord. system
x_cp = x_prime - x_B_c  # [m] longitudinal position
y_cp = 0  # [m] transverse position
z_cp = -h / 2  # [m] vertical position
x_g_AP = -L / 2 - x_cp  # [m] position of leakage at AP relative to center of gravity
x_g_FP = L / 2 - x_cp  # [m] position of leakage at FP relative to center of gravity
lcg_fan = 5.6 - x_cp  # [m] Longitudinal fan position (from CG)  # TODO: Make sure this is correct

# Excitation
F_3a = f_ex_temp[2, :]
F_5a = f_ex_temp[4, :]  # TEST TO SEE WHAT HAPPENS WHEN THE SIGN IS CHANGED

k = np.power(omega_0, 2) / g  # wave number of water waves
omega_e = omega_0 + np.power(omega_0, 2) / g * U  # encounter frequencies
f_encounter = omega_e / 2 / np.pi  # [Hz] frequency of encounter
n_freq = len(omega_e)

# zeta_a = Zeta_a(omega_0, H_s, T_p)  # [m] wave amplitude dependent on encounter frequency
zeta_a = np.ones([len(omega_0)])

water_wavelength = g / 2 / np.pi * np.power(np.divide(2 * np.pi, omega_0), 2)
encounter_wavelength = g / 2 / np.pi * np.power(np.divide(2 * np.pi, omega_e), 2)

# ***** Compute wave pumping ***** # TODO: Make sure the phase is correct/adjust centroid
# Properties of SES-X shaped air cushion with origin in motion coordinate system
cushion_width = b  # [m]
x_s = x_prime  # [m]
x_f = -(l_1 - x_prime)  # [m]
y_s = b/2  # [m]
y_p = -b/2  # [m]
x_b = -(l_1 + l_2 - x_prime)  # [m]
F_wp = wave_pumping_excitation_sesx(x_f, x_s, y_s, x_b, omega_0, U, b)

# Plot for testing:
'''
plt.plot(f_encounter, np.abs(F_3a))
plt.xlabel('Encounter frequency [Hz]')
plt.show()
'''

plot_Coefficients = False

if plot_Coefficients:
    '''
    plt.plot(f_encounter, np.abs(F_3a), label='F_3a')
    plt.plot(f_encounter, A_35, label='A_35')
    plt.plot(f_encounter, B_35, label='B_35')
    plt.legend()
    plt.xlabel('Encounter frequency [Hz]')
    plt.show()


    plt.plot(encounter_wavelength, np.abs(F_3a), label='F_3a')
    plt.plot(encounter_wavelength, A_53, label='A_53')
    plt.plot(encounter_wavelength, B_53, label='B_53')
    plt.legend()
    plt.xlabel('Encounter wavelength [m]')
    plt.show()
    '''

    plt.plot(f_encounter, np.abs(F_3a), label='F_3a')
    plt.plot(f_encounter, A_35, label='A_35')
    plt.plot(f_encounter, B_35, label='B_35')
    plt.plot(f_encounter, A_53, label='A_35')
    plt.plot(f_encounter, np.ones([n_freq]) * B_53, label='B_35')
    plt.legend()
    plt.xlabel('Encounter frequency [Hz]')
    plt.show()

    plt.plot(f_encounter, np.abs(F_3a), label='F_3a')
    plt.plot(f_encounter, np.ones([n_freq]) * A_33, label='A_33')
    plt.plot(f_encounter, np.ones([n_freq]) * B_33, label='B_33')
    plt.plot(f_encounter, B_55, label='B_55')
    plt.legend()
    plt.xlabel('Encounter frequency [Hz]')
    plt.show()

    '''
    plt.plot(encounter_wavelength, np.abs(F_3a), label='F_3a')
    plt.plot(encounter_wavelength, A_33, label='A_33')
    plt.plot(encounter_wavelength, B_33, label='B_33')
    plt.legend()
    plt.xlabel('Encounter wavelength [m]')
    plt.show()
    '''

# ***** Compute constants *****
k_1 = K_1(rho_0, p_0, h_0, A_c0)  # [kg] eq. (83) in Steen and Faltinsen (1995)

k_2_AP = K_2(p_0, 1.0)  # [m/s] linearized equivalent outflow velocity constant at stern 1.0
k_2_FP = K_2(p_0, 0.61)  # [m/s] linearized equivalent outflow velocity constant at bow 0.61

k_3 = K_3(rho_0, p_0, Q_0, dQdp_0)  # air cushion flow airflow constant

# ***** Start of iteration loop to fit N_R and N_B *****
# Guess an initial value for bias and gain value of non-linear leakage
n_R_AP = 0.5
n_R_FP = 0.5
n_B_AP = 0.5
n_B_FP = 0.5

# Initialize values to contain RAOs of the current and previous step in the iteration process
eta_3a = np.zeros([1, n_freq], dtype=complex)
eta_5a = np.zeros([1, n_freq], dtype=complex)
mu_ua = np.zeros([1, n_freq], dtype=complex)
eta_3a_old = np.zeros([1, n_freq], dtype=complex)
eta_5a_old = np.zeros([1, n_freq], dtype=complex)
mu_ua_old = np.zeros([1, n_freq], dtype=complex)

epsi = 1e-6  # [-] allowed error in stopping criteria for iteration process
rel_err = -1  # [-] initialize error variable to start while-loop
abs_err = -1
counter = 0  # initialize counter in while-loop
max_iter = 100  # maximum number of iterations

eta_3_test = []

while ((rel_err > epsi) or (counter < 2)) and (counter < max_iter):

    # Solve mean value relation
    eta_3m, eta_5m, mu_um = solve_mean_value_relation(n_B_AP, n_B_FP, L, b, x_cp, A_c0, p_0, k_2_AP, k_2_FP, k_3, h_s_AP,
                                                      h_s_FP, C_33, C_55, C_35, C_53)
    print('eta_3m = ', eta_3m)
    print('eta_5m = ', eta_5m)
    print('eta_7m = ', mu_um)

    # Correct mean cushion pressure
    p_0 = (1 + mu_um) * p_0

    # Compute mean leakage area at AP and FP
    a_0_AP = np.maximum(A_0_AP(L, b, n_B_AP, eta_3m, eta_5m, h_s_AP), 0)  # , 0.5) #
    a_0_FP = np.maximum(A_0_FP(L, b, n_B_FP, eta_3m, eta_5m, h_s_FP), 0)  # , 0.5) #

    j_max = 2  # number of acoustic modes to include in the calculations

    A_mat = np.zeros([3, 3, n_freq], dtype=complex)  # initialize coefficient matrix for linear system of eq.
    f_vec = np.zeros([3, n_freq], dtype=complex)  # initialize column vector on the right hand side of the equation

    # Assign values for coefficient matrix in solving A_mat*x_vec = f_vec, with exception of terms due to spat. varying pressure
    # Heave equation, i.e. eq. (87) Steen and Faltinsen (1995)
    A_mat[0, 0, :] += -np.multiply((M[0, 0] + A_33), np.power(omega_e, 2)) + 1j * np.multiply(B_33, omega_e) + C_33
    A_mat[0, 1, :] += -np.multiply((A_35 + M[0, 1]), np.power(omega_e, 2)) + 1j * np.multiply(B_35, omega_e) + C_35
    A_mat[0, 2, :] += -p_0 * A_c0
    f_vec[0, :] += F_3a

    # Pitch equation, i.e. eq. (88) Steen and Faltinsen (1995)
    A_mat[1, 0, :] += -np.multiply(A_53 + M[1, 0], np.power(omega_e, 2)) + 1j * np.multiply(B_53, omega_e) + C_53
    A_mat[1, 1, :] += -np.multiply((M[1, 1] + A_55), np.power(omega_e, 2)) + 1j * np.multiply(B_55, omega_e) + C_55
    A_mat[1, 2, :] += p_0 * x_cp * A_c0
    f_vec[1, :] += F_5a  # F_5a  #

    # Equation of dynamic uniform pressure, i.e. eq. (82) Steen and Faltinsen (1995)
    A_mat[2, 0, :] = rho_a * b * (k_2_AP * n_R_AP + k_2_FP * n_R_FP) + rho_0 * A_c0 * 1j * omega_e
    A_mat[2, 1, :] = rho_a * b * L / 2 * (k_2_AP * n_R_AP - k_2_FP * n_R_FP) - rho_0 * A_c0 * x_cp * 1j * omega_e
    A_mat[2, 2, :] = k_1 * 1j * omega_e + k_3
    f_vec[2, :] = rho_0 * F_wp - rho_a * b * 1j * (
                k_2_AP * n_R_AP * np.exp(-1j * k * L / 2) + k_2_FP * n_R_FP * np.exp(1j * k * L / 2))

    # ***** Fill in terms due to spatially varying pressure *****
    for j in range(1, j_max + 1):

        # Computes variables dependent on j
        xi_j = Xi_j(j, rho_0, p_0, h_0, b, L, k_2_AP, k_2_FP, a_0_AP, a_0_FP, dQdp_0,
                    lcg_fan + x_cp)  # [-] Relative damping ratio of mode j  # TODO check if I should use lcg_fan instead of x_F
        omega_j = Omega_j(j, L)
        k_4 = K_4(xi_j, h_0, omega_e, omega_j)

        # Frequency dependent modal amplitudes for odd acoustic modes
        a_0j = A_0j(j, b, L, p_0, dQdp_0, lcg_fan + x_cp, k_2_AP, k_2_FP, a_0_AP, a_0_FP, k_4)
        a_3j = A_3j(L, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)
        a_5j = A_5j(j, L, omega_e, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)
        a_7j = A_7j(j, k_4, L, k, omega_e, k_2_AP, k_2_FP, n_R_AP, n_R_FP)

        # Frequency dependent modal amplitudes for even acoustic modes
        b_0j = B_0j(j, b, L, p_0, dQdp_0, lcg_fan + x_cp, k_2_AP, k_2_FP, a_0_AP, a_0_FP, k_4)
        b_3j = B_3j(L, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)
        b_5j = B_5j(k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)
        b_7j = B_7j(j, k_4, L, k, omega_e, k_2_AP, k_2_FP, n_R_AP, n_R_FP)

        '''
        # Deep copy for debugging  #TODO: Remove this
        temp_mat = A_mat.copy()
        temp_vec = f_vec.copy()

        ''''''
        append_spatially_varying_terms(temp_mat, temp_vec, omega_e, j, b, L, rho_0, p_0, dQdp_0, lcg_fan, k_2_AP, a_0_AP
                                       , x_g_AP, k_2_FP, a_0_FP, x_g_FP, zeta_a, a_0j, a_3j, a_5j, a_7j, b_0j, b_3j, b_5j, b_7j)
        '''

        if j % 2 == 1:  # if j is odd
            # ***** Pitching moments due to spatially varying pressure in eq. (88) in Steen and Faltinsen (1995)
            # Heave DOF
            A_mat[1, 0, :] += 2 * rho_0 * b * 1j * (L / j / np.pi) ** 2 * np.multiply(omega_e, a_3j)

            # Pitch DOF
            A_mat[1, 1, :] += 2 * rho_0 * b * 1j * (L / j / np.pi) ** 2 * np.multiply(omega_e, a_5j)

            # Uniform pressure DOF
            A_mat[1, 2, :] += 2 * rho_0 * b * 1j * (L / j / np.pi) ** 2 * np.multiply(omega_e, a_0j)

            # Wave excitation term
            f_vec[1, :] -= 2 * rho_0 * b * 1j * (L / j / np.pi) ** 2 * np.multiply(np.multiply(omega_e, a_7j), zeta_a)

            # ***** Terms from leakage at AP in eq. (82) in Steen and Faltinsen (1995) *****
            # Heave DOF
            A_mat[2, 0, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_3j) * r_j(
                x_g_AP + x_cp, j, L)

            # Pitch DOF
            A_mat[2, 1, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_5j) * r_j(
                x_g_AP + x_cp, j, L)

            # Uniform pressure DOF
            A_mat[2, 2, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_0j) * r_j(
                x_g_AP + x_cp, j, L)

            # Wave excitation term
            f_vec[2, :] -= rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(
                np.multiply(omega_e, a_7j) * r_j(x_g_AP + x_cp, j, L), zeta_a)

            # ***** Terms from leakage at FP in eq. (82) in Steen and Faltinsen (1995) *****
            # Heave DOF
            A_mat[2, 0, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_3j) * r_j(
                x_g_FP + x_cp, j, L)

            # Pitch DOF
            A_mat[2, 1, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_5j) * r_j(
                x_g_FP + x_cp, j, L)

            # Uniform pressure DOF
            A_mat[2, 2, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_0j) * r_j(
                x_g_FP + x_cp, j, L)

            # Wave excitation term
            f_vec[2, :] -= rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(
                np.multiply(omega_e, a_7j) * r_j(x_g_FP + x_cp, j, L), zeta_a)

            # ***** Terms from fan inlet at lcg_fan in eq. (82) in Steen and Faltinsen (1995) *****
            # Heave DOF
            A_mat[2, 0, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_3j) * r_j(
                lcg_fan + x_cp, j, L)

            # Pitch DOF
            A_mat[2, 1, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_5j) * r_j(
                lcg_fan + x_cp, j, L)

            # Uniform pressure DOF
            A_mat[2, 2, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_0j) * r_j(
                lcg_fan + x_cp, j, L)

            # Wave excitation term
            f_vec[2, :] -= (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(
                np.multiply(omega_e, a_7j) * r_j(lcg_fan + x_cp, j, L), zeta_a)

        elif j % 2 == 0:  # if j is even
            # ***** Terms from leakage at AP in eq. (82) in Steen and Faltinsen (1995) *****
            # Heave DOF
            A_mat[2, 0, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_3j) * r_j(
                x_g_AP + x_cp, j, L)

            # Pitch DOF
            A_mat[2, 1, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_5j) * r_j(
                x_g_AP + x_cp, j, L)

            # Uniform pressure DOF
            A_mat[2, 2, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_0j) * r_j(
                x_g_AP + x_cp, j, L)

            # Wave excitation term
            f_vec[2, :] -= rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(
                np.multiply(omega_e, b_7j) * r_j(x_g_AP + x_cp, j, L), zeta_a)

            # ***** Terms from leakage at FP in eq. (82) in Steen and Faltinsen (1995) *****
            # Heave DOF
            A_mat[2, 0, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_3j) * r_j(
                x_g_FP + x_cp, j, L)

            # Pitch DOF
            A_mat[2, 1, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_5j) * r_j(
                x_g_FP + x_cp, j, L)

            # Uniform pressure DOF
            A_mat[2, 2, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_0j) * r_j(
                x_g_FP + x_cp, j, L)

            # Wave excitation term
            f_vec[2, :] -= rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(
                np.multiply(omega_e, b_7j) * r_j(x_g_FP + x_cp, j, L), zeta_a)

            # ***** Terms from fan inlet at x_F in eq. (82) in Steen and Faltinsen (1995) *****
            # Heave DOF
            A_mat[2, 0, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_3j) * r_j(
                lcg_fan + x_cp, j, L)

            # Pitch DOF
            A_mat[2, 1, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_5j) * r_j(
                lcg_fan + x_cp, j, L)

            # Uniform pressure DOF
            A_mat[2, 2, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_0j) * r_j(
                lcg_fan + x_cp, j, L)

            # Wave excitation term
            f_vec[2, :] -= (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(
                np.multiply(omega_e, b_7j) * r_j(lcg_fan + x_cp, j, L), zeta_a)

        ''' # TODO: Remove this
        test_mat = temp_mat - A_mat
        test_vec = temp_vec - f_vec

        print('Compare matrices:')
        print(temp_mat == A_mat)
        print('Compare vector:')
        print(temp_vec == f_vec)
        print('hei')
        '''

        # Solves linear system of equations for each frequency
        eta_3a, eta_5a, mu_ua = solve_linear_systems_of_eq(A_mat, f_vec)

        eta_3_test.append(eta_3a[300])  # Append value to test vector
        if not counter == 0:  # compute error between current and previous if it is not the first iteration

            rel_errors = np.hstack([np.abs(np.divide(eta_3a - eta_3a_old, eta_3a_old)),
                                    np.abs(np.divide(eta_5a - eta_5a_old, eta_5a_old)),
                                    np.abs(np.divide(mu_ua - mu_ua_old, mu_ua_old))])

            rel_errors[np.isnan(rel_errors)] = 0.

            rel_err = np.amax(rel_errors)

            # gh = np.isnan(rel_err)

            # rel_err[np.isnan(rel_err)] = 0.

            abs_err = np.amax([np.abs(eta_3a - eta_3a_old), np.abs(eta_5a - eta_5a_old), np.abs(mu_ua - mu_ua_old)])

            err_eta_3a = np.abs(np.divide(eta_3a - eta_3a_old, eta_3a_old))
            i_max_err_eta_3a = np.argmax(err_eta_3a)
            err_eta_5a = np.abs(np.divide(eta_5a - eta_5a_old, eta_5a_old))
            i_max_err_eta_5a = np.argmax(err_eta_5a)
            err_mu_ua = np.abs(np.divide(mu_ua - mu_ua_old, mu_ua_old))
            i_max_err_mu_ua = np.argmax(err_mu_ua)

            print('The error largest error in heave occurs at index: ' + str(i_max_err_eta_3a))
            print('The value of eta_3a at the current step is:\t' + str(np.abs(eta_3a[i_max_err_eta_3a])))
            print('The value of eta_3a at previous step is:\t' + str(np.abs(eta_3a_old[i_max_err_eta_3a])))
            print('Giving a relative error of:\t\t\t\t\t' + str(err_eta_3a[i_max_err_eta_3a]))
            print('\nThe error largest error in pitch occurs at index: ' + str(i_max_err_eta_5a))
            print('The value of eta_5a at the current step is:\t' + str(np.abs(eta_5a[i_max_err_eta_5a])))
            print('The value of eta_5a at previous step is:\t' + str(np.abs(eta_5a_old[i_max_err_eta_5a])))
            print('Giving a relative error of:\t\t\t\t\t' + str(err_eta_5a[i_max_err_eta_5a]))
            print('\nThe error largest error in uniform pressure occurs at index: ' + str(i_max_err_mu_ua))
            print('The value of mu_ua at the current step is:\t' + str(np.abs(mu_ua[i_max_err_mu_ua])))
            print('The value of mu_ua at previous step is:\t\t' + str(np.abs(mu_ua_old[i_max_err_mu_ua])))
            print('Giving a relative error of:\t\t\t\t\t' + str(err_mu_ua[i_max_err_mu_ua]))
            print()
        # Stores current solution for comparison in next iteration
        eta_3a_old = eta_3a
        eta_5a_old = eta_5a
        mu_ua_old = mu_ua

        counter += 1  # increment counter
        print("Current iteration:", counter, 'with relative error:', rel_err)
        # Compute new value for bias and gain for the linearized leakage for the next iteration

        #Have taken into account change in sign
        b_L_AP = eta_3m - L / 2 * eta_5m - h_s_AP
        b_L_FP = eta_3m + L / 2 * eta_5m - h_s_FP

        sigma_L_AP = rms_leakage(L / 2, omega_0, eta_3a, eta_5a, H_s, T_p, zeta_a)
        sigma_L_FP = rms_leakage(-L / 2, omega_0, eta_3a, eta_5a, H_s, T_p, zeta_a)

        n_R_AP = N_R(b_L_AP, sigma_L_AP)
        n_R_FP = N_R(b_L_FP, sigma_L_FP)
        n_B_AP = N_B(b_L_AP, sigma_L_AP)
        n_B_FP = N_B(b_L_FP, sigma_L_FP)

'''
plt.plot(f_encounter, np.abs(eta_3a), 'x-', label='$\\eta_3$')
plt.plot(f_encounter, np.abs(eta_5a), 'x-',  label='$\\eta_5$')
plt.plot(f_encounter, np.abs(mu_ua), 'x-',  label='$\\mu_{ua}$')
plt.xlabel('Encounter frequency [Hz]')
plt.xlim([0, 16])
plt.legend()
plt.show()
'''

# Store results
store_results = False
if store_results:
    df_result = pd.DataFrame(
        {'f_enc [Hz]': f_encounter, 'zeta_a [m]': zeta_a, 'eta_3a [m]': np.abs(eta_3a), 'eta_5a [rad]': np.abs(eta_5a),
         'mu_ua [-]': np.abs(mu_ua)})
    df_result.to_csv('Results/results_Test1.csv')

# Compare with Steen and Faltinsen (1995)
df_uniform = pd.read_csv(
    'C:/Users/mathi/OneDrive - NTNU/Master Thesis/Spatially varying pressure/Results from Steen and Faltinsen (1995)/Rigid panel model/Uniform pressure/uniform_pressure_RAO_nice_format.csv')

# Define plot colors
color_BBGreen = '#5cb16d'
color_BBPurple = '#b15ca0'

# Divide by wave amplitude, but makes sure zeta_a is not zero
mu_ua_nondim = np.zeros([n_freq])
for i in range(n_freq):
    if zeta_a[i] > 1e-20:
        mu_ua_nondim[i] = np.abs(mu_ua[i]) / zeta_a[i]  # 1 #

save_RAOs_for_comparison = False
# plt.plot(f_encounter, np.abs(mu_ua), label='This program')

# plt.plot(f_encounter, mu_ua_nondim, label=r'\textrm{Computed}, $Hs=' + str(H_s) + '\,m$, $Tp=' + str(T_p) + '\,s$', color=color_BBGreen)
# plt.plot(df_uniform.iloc[:, 1], df_uniform.iloc[:, 2], label=r'\textrm{Steen and Faltinsen (1995)}, $Hs=0.15\,m$, $Tp=1.5\,s$', color=color_BBPurple)

plt.plot(f_encounter, mu_ua_nondim, label=r'\textrm{Computed}', color=color_BBGreen)
plt.plot(df_uniform.iloc[:, 1], df_uniform.iloc[:, 2], label=r'\textrm{Steen and O. M. Faltinsen (1995)}',
         color=color_BBPurple)
# plt.title(r'\textrm{Comparison with Steen and Faltinsen (1995)}')
plt.xlabel(r'\textrm{Encounter frequency} $[Hz]$')
plt.ylabel(r'$|\hat{\eta}_{7}|\,[-] /\zeta_a\,[m]$')
x_min = 0.
plt.xlim([x_min, 16])
plt.ylim([0, np.max(mu_ua_nondim[f_encounter > x_min])])
plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/Comparison RAOs/uniform pressure.pdf', bbox_inches='tight')
plt.show()

df_vert_acc = pd.read_csv(
    'C:/Users/mathi/OneDrive - NTNU/Master Thesis/Spatially varying pressure/Results from Steen and Faltinsen (1995)/Rigid panel model/Vert. Acc. AP/vertical_acc_AP_RAO_nice_format.csv')

vert_acc_AP_nondim = np.zeros([n_freq])
vert_motion_AP_nondim = np.zeros([n_freq])
for i in range(n_freq):
    if zeta_a[i] > 1e-20:
        vert_acc_AP_nondim[i] = np.absolute(omega_e[i] ** 2 * (eta_3a[i] + L / 2 * eta_5a[i])) / zeta_a[i]  # / 1  #
        vert_motion_AP_nondim[i] = np.absolute((eta_3a[i] + L / 2 * eta_5a[i])) / zeta_a[i]

# plt.plot(f_encounter, vert_acc_AP_nondim, label=r'\textrm{Computed}, $Hs=' + str(H_s) + '\,m$, $Tp=' + str(T_p) + '\,s$', color=color_BBGreen)
# plt.plot(df_vert_acc.iloc[:, 1], df_vert_acc.iloc[:, 2], label=r'\textrm{Steen and Faltinsen (1995)}, $Hs=0.15\,m$, $Tp=1.5\,s$', color=color_BBPurple)

plt.plot(f_encounter, vert_acc_AP_nondim, label=r'\textrm{Computed}', color=color_BBGreen)
plt.plot(df_vert_acc.iloc[:, 1], df_vert_acc.iloc[:, 2], label=r'\textrm{Steen and O. M. Faltinsen (1995)}',
         color=color_BBPurple)
plt.xlim([x_min, 16.])
plt.ylim([0, np.max(vert_acc_AP_nondim[f_encounter > x_min])])
# plt.title(r'\textrm{Comparison with Steen and Faltinsen (1995)}')
plt.xlabel(r'\textrm{Encounter frequency} $[Hz]$')
plt.ylabel(r'\textrm{Vert. acc.} $\,[m/s^2] / \zeta_a\,[m]$')
plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/Comparison RAOs/vert_acc_bow.pdf', bbox_inches='tight')
plt.show()

plt.plot(f_encounter, vert_motion_AP_nondim, label=r'\textrm{Computed}', color=color_BBGreen)
plt.xlim([x_min, 16.])
plt.ylim([0, np.max(vert_motion_AP_nondim[f_encounter > x_min])])
# plt.title(r'\textrm{Comparison with Steen and Faltinsen (1995)}')
plt.xlabel(r'\textrm{Encounter frequency} $[Hz]$')
plt.ylabel(r'\textrm{Vert. motion} $\,[m] / \zeta_a\,[m]$')
plt.legend()
plt.show()

# Plotting wave spectrum
plotWaveSpectrum = False
if plotWaveSpectrum:
    f_e_PM = np.linspace(0.1, 50, 1000)
    omega_0_PM = g / 2 / U * (np.sqrt(1 + 8 * np.pi * U / g * f_e_PM) - 1)
    plt.plot(f_e_PM, PM_spectrum(omega_0_PM, H_s, T_p), color=color_BBGreen)

    # plt.plot(f_encounter, PM_spectrum(omega_0, H_s, T_p))
    plt.xlabel(r'\textrm{Encounter frequency\,[Hz]}')
    plt.ylabel(r'$S_{\zeta}\,[m^2s]$')
    plt.xlim([0.1, 50])
    # plt.title('U = ' + str(round(U, 2)) + '[m/s], $H_s$ = ' + str(H_s) + '[m], $T_p$ = ' + str(T_p) + '[s], $\\beta=0^\\degree$')
    plt.savefig('Results/Comparison RAOs/mod_PM_spectrum.pdf', bbox_inches='tight')
    plt.show()

    # plt.plot(f_e_PM, PM_spectrum(omega_0_PM, H_s, T_p))

print('The iteration scheme converged after', counter)
print('Relative error:\t\t', rel_err)
print('Absolute error:\t\t', abs_err)

# plt.plot(np.abs(eta_3_test))
# plt.show()