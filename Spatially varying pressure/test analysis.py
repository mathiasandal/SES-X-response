import numpy as np
import matplotlib.pyplot as plt
from veres import read_re8_file, read_re7_file
from utilities import K_1, K_2, K_3, K_4, Xi_j, Omega_j, solve_mean_value_relation, A_0_AP, A_0_FP, A_0j, A_3j, A_5j, \
    A_7j, B_0j, B_3j, B_5j, B_7j, r_j, solve_linear_systems_of_eq, N_R, N_B, rms_leakage, Zeta_a, \
    append_spatially_varying_terms
"""Implementation of analysis with spatially varying pressure with the same input as in p. 35 in Steen and Faltinsen (1995). 'Cobblestone
    Oscillations of an SES with Flexible Bag Aft Seal'
"""

# ***** Define input parameters *****

# Wave parameters
H_s = 0.15  # [m] significant wave height
T_p = 1.5  # [s] wave peak period
plotWaveSpectrum = True

if plotWaveSpectrum:
    frequencies_PM = np.linspace(0, 16)
    U_PM = 50*0.514444
    omega_e_PM = 2*np.pi*frequencies_PM


# Physical constants
c = 343.  # [m/s] Speed of sound in air
g = 9.81  # [m/s^2] Acceleration of gravity
p_a = 101325.  # [Pa] Atmospheric pressure
rho_0 = 1.2796  # [kg/m^3] Density of air at mean cushion pressure (p_0 + p_a)
# source: https://www.gribble.org/cycling/air_density.html at 15deg C, p_0 + p_a

rho_a = 1.225  # [kg/m^3] Density of air at atmospheric pressure
rho_w = 1000.  # [kg/m^3] Density of fresh water

# SES main dimensions
L_oa = 35.  # [m] Length overall
L = 28.  # [m] Air cushion length
b = 8.  # [m] Air cushion beam
m = 140.  # [tons] Vessel total mass
m = m * 1e3  # [kg] Vessel total mass
h_0 = 2.0  # [m] Cushion height
U = 50.  # [knots] Velocity
U = U*0.514444  # [m/s] Velocity

# Fan characteristics
p_0 = 500.  # [mmWc] Mean cushion pressure
p_0 = rho_w*g*p_0*1e-3  # [Pa] Mean cushion pressure
Q_0 = 150.  # [m^3/s] Mean fan flow rate
dQdp_0 = -140.  # [m^2/s] Linear fan slope
dQdp_0 = dQdp_0 / rho_w / g  # [(m^3/s)/Pa] Linear fan slope

lcg_fan = 5.6  # [m] Longitudinal fan position (from CG)
x_F = 0  # [m]

# Other parameters
h_s_AP = 0.1  # [m] aft seal submergence
h_s_FP = 0.1  # [m] bow seal submergence
x_g_AP = -L/2  # [m] position of leakage at AP relative to center of gravity
x_g_FP = L/2  # [m] position of leakage at FP relative to center of gravity

# Derived parameters
A_c = L*b  # [m^2] Air cushion area  # TODO: Might want to use expression for rectangular cushion shape with triangle at the front
x_cp = 4  # [m] longitudinal centroid of air cushion relative to CoG(?)  #TODO: Make sure this is correct

# ***** Read hydrodynamic coefficients for conceptual SES *****

# Fetch hydrodynamic coefficients and loads
path = 'C:/Users/mathi/code/repos/SES-X-response/Spatially varying pressure/Input Files/Conceptual SES of 20m/Run 1 (0.1-16 [Hz])/'
path_re7 = path + 'input.re7'
path_re8 = path + 'input.re8'


# Read input.re7 file to get hydrodynamic coefficients
M, A_n, B_n, C_n, VEL_re7, HEAD_re7, FREQ_re7, XMTN_re7, ZMTN_re7, NDOF = read_re7_file(path_re7)

# Clean up matrices
A_temp = A_n[0, 0, :, :, :]
B_temp = B_n[0, 0, :, :, :]
C_temp = C_n[0, 0, :, :, :]

# Read input.re8 file to get excitation
REFORCE, IMFORCE, VEL_re8, HEAD_re8, FREQ_re8, XMTN_re8, ZMTN_re8 = read_re8_file(path_re8)
n_freq = len(FREQ_re8)

A_33 = A_temp[:, 2, 2]
B_33 = B_temp[:, 2, 2]
C_33 = C_temp[n_freq//2, 2, 2]

A_35 = A_temp[:, 2, 4]
B_35 = C_temp[:, 2, 4]
C_35 = C_temp[n_freq//2, 2, 4]

A_53 = A_temp[:, 4, 2]
B_53 = B_temp[:, 4, 2]
C_53 = B_temp[n_freq//2, 4, 2]

r_55 = 0.25 * L  # [m] radii of gyration in pitch
I_55 = m * r_55**2
A_55 = A_temp[:, 4, 4]
B_55 = B_temp[:, 4, 4]
C_55 = C_temp[n_freq//2, 4, 4]

# Excitation
# initialize array to contain complex force amplitudes

f_ex = np.zeros([6, n_freq], dtype=complex)

for i in range(6):
    for j in range(n_freq):
        f_ex[i, j] = REFORCE[i, j, 0, 0] + 1j * IMFORCE[i, j, 0, 0]

F_3a = f_ex[2, :]
F_5a = f_ex[4, :]

#omega_0 = np.linspace(1, 10, 1000)
omega_0 = FREQ_re7
k = np.power(omega_0, 2)/g  # wave number of water waves
omega_e = omega_0 + np.power(omega_0, 2)/g*U  # encounter frequencies
f_encounter = omega_e/2/np.pi  # [Hz] frequency of encounter
n_freq = len(omega_e)

# TODO: Make sure that zeta_a is treated correctly in all functions and expressions
zeta_a = Zeta_a(omega_0, H_s, T_p)  # [m] wave amplitude dependent on encounter frequency
#zeta_a = np.ones([len(omega_0)])

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

# Initialize values to contain RAOs of the current and previous step in the iteration process
eta_3a = np.zeros([1, n_freq], dtype=complex)
eta_5a = np.zeros([1, n_freq], dtype=complex)
mu_ua = np.zeros([1, n_freq], dtype=complex)
eta_3a_old = np.zeros([1, n_freq], dtype=complex)
eta_5a_old = np.zeros([1, n_freq], dtype=complex)
mu_ua_old = np.zeros([1, n_freq], dtype=complex)

epsi = 1e-6  # [-] allowed error in stopping criteria for iteration process
err = -1  # [-] initialize error variable to start while-loop
counter = 0  # initialize counter in while-loop
max_iter = 500  # maximum number of iterations

eta_3_test = []

while ((err > epsi) or (counter < 2)) and (counter < max_iter):

    # Solve mean value relation
    eta_3m, eta_5m = solve_mean_value_relation(n_B_AP, n_B_FP, L, b, x_cp, A_c, p_0, k_2_AP, k_2_FP, k_3, h_s_AP, h_s_FP, C_33, C_55, C_35, C_53)
    # Compute mean leakage area et AP and FP
    a_0_AP = A_0_AP(L, b, n_B_AP, eta_3m, eta_5m, h_s_AP)
    a_0_FP = A_0_FP(L, b, n_B_FP, eta_3m, eta_5m, h_s_FP)

    j_max = 4  # number of acoustic modes to include in the calculations

    A_mat = np.zeros([3, 3, n_freq], dtype=complex)  # initialize coefficient matrix for linear system of eq.
    f_vec = np.zeros([3, n_freq], dtype=complex)  # initialize column vector on the right hand side of the equation

    # Assign values for coefficient matrix in solving A_mat*x_vec = f_vec, with exception of terms due to spat. varying pressure
    # Heave equation, i.e. eq. (87) Steen and Faltinsen (1995)
    A_mat[0, 0, :] = -(m + A_33)*np.power(omega_e, 2) + B_33 * 1j * omega_e + C_33
    A_mat[0, 1, :] = -A_35*np.power(omega_e, 2) + B_35 * 1j * omega_e + C_35
    A_mat[0, 2, :] = -A_c * p_0
    f_vec[0, :] = F_3a

    # Pitch equation, i.e. eq. (88) Steen and Faltinsen (1995)
    A_mat[1, 0, :] = -A_53*np.power(omega_e, 2) + B_53 * 1j * omega_e + C_53
    A_mat[1, 1, :] = -(I_55+A_55)*np.power(omega_e, 2) + B_55 * 1j * omega_e + C_55
    A_mat[1, 2, :] = A_c * p_0 * x_cp
    f_vec[1, :] = F_5a

    # Equation of dynamic uniform pressure, i.e. eq. (82) Steen and Faltinsen (1995)
    A_mat[2, 0, :] = rho_a * b * (k_2_AP*n_R_AP + k_2_FP*n_R_FP) + rho_0 * A_c * 1j * omega_e
    A_mat[2, 1, :] = rho_a * b * L/2 * (k_2_AP*n_R_AP - k_2_FP*n_R_FP) - rho_0 * A_c * x_cp * 1j * omega_e
    A_mat[2, 2, :] = k_1 * 1j * omega_e + k_3
    f_vec[2, :] = F_wp - rho_a * b * 1j * (k_2_AP*n_R_AP * np.exp(-1j * k * L/2) + k_2_FP*n_R_FP * np.exp(1j * k * L/2))

    #A_mat_2 = A_mat
    #f_vec_2 = f_vec

    # ***** Fill in terms due to spatially varying pressure *****
    for j in range(1, j_max + 1):

        # Computes variables dependent on j
        xi_j = Xi_j(j, rho_0, p_0, h_0, b, L, k_2_AP, k_2_FP, a_0_AP, a_0_FP, dQdp_0, lcg_fan)  # [-] Relative damping ratio of mode j  # TODO check if I should use lcg_fan instead of x_F
        omega_j = Omega_j(j, L)
        k_4 = K_4(xi_j, h_0, omega_e, omega_j)

        # TODO: Compute A_0j(omega_e) and B_0j(omega_e) amplitudes beforehand(?)

        # Frequency dependent modal amplitudes for odd acoustic modes
        a_0j = A_0j(j, b, L, p_0, dQdp_0, lcg_fan, k_2_AP, k_2_FP, a_0_AP, a_0_FP, k_4)
        a_3j = A_3j(L, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)
        a_5j = A_5j(j, L, omega_e, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)
        a_7j = A_7j(j, k_4, L, k, omega_e, k_2_AP, k_2_FP, n_R_AP, n_R_FP)

        # Frequency dependent modal amplitudes for even acoustic modes
        b_0j = B_0j(j, b, L, p_0, dQdp_0, lcg_fan, k_2_AP, k_2_FP, a_0_AP, a_0_FP, k_4)
        b_3j = B_3j(L, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)
        b_5j = B_5j(k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP)
        b_7j = B_7j(j, k_4, L, k, omega_e, k_2_AP, k_2_FP, n_R_AP, n_R_FP)


        '''
        append_spatially_varying_terms(A_mat, f_vec, omega_e, j, b, L, rho_0, p_0, dQdp_0, lcg_fan,
                                                      k_2_AP, a_0_AP, x_g_AP, k_2_FP, a_0_FP, x_g_FP, zeta_a,
                                                      a_0j, a_3j, a_5j, a_7j, b_0j, b_3j, b_5j, b_7j)
        
        append_spatially_varying_terms(A_mat_2, f_vec_2, omega_e, j, b, L, rho_0, p_0, dQdp_0, lcg_fan,
                                                        k_2_AP, a_0_AP, x_g_AP, k_2_FP, a_0_FP, x_g_FP, zeta_a,
                                                        a_0j, a_3j, a_5j, a_7j, b_0j, b_3j, b_5j, b_7j)
        append_spatially_varying_terms(A_mat_2, f_vec_2, omega_e, j, b, L, rho_0, p_0, dQdp_0, lcg_fan,
                                                        k_2_AP, a_0_AP, x_g_AP, k_2_FP, a_0_FP, x_g_FP, zeta_a,
                                                        a_0j, a_3j, a_5j, a_7j, b_0j, b_3j, b_5j, b_7j)
        '''
        '''
        '''

        if j % 2 == 1:  # if j is odd
            # ***** Pitching moments due to spatially varying pressure in eq. (88) in Steen and Faltinsen (1995)
            # Heave DOF
            A_mat[1, 0, :] += 2*rho_0*b*1j * (L/j/np.pi)**2 * np.multiply(omega_e, a_3j)

            # Pitch DOF
            A_mat[1, 1, :] += 2*rho_0*b*1j * (L/j/np.pi)**2 * np.multiply(omega_e, a_5j)

            # Uniform pressure DOF
            A_mat[1, 2, :] += 2*rho_0*b*1j * (L/j/np.pi)**2 * np.multiply(omega_e, a_0j)

            # Wave excitation term
            f_vec[1, :] -= 2*rho_0*b*1j * (L/j/np.pi)**2 * np.multiply(np.multiply(omega_e, a_7j), zeta_a)

            # ***** Terms from leakage at AP in eq. (82) in Steen and Faltinsen (1995) *****
            # Heave DOF
            A_mat[2, 0, :] += rho_a*k_2_AP/2*a_0_AP * (-rho_0/p_0)*1j*np.multiply(omega_e, a_3j) * r_j(x_g_AP, j, L)

            # Pitch DOF
            A_mat[2, 1, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_5j) * r_j(x_g_AP, j, L)

            # Uniform pressure DOF
            A_mat[2, 2, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_0j) * r_j(x_g_AP, j, L)

            # Wave excitation term
            f_vec[2, :] -= rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(np.multiply(omega_e, a_7j) * r_j(x_g_AP, j, L), zeta_a)

            # ***** Terms from leakage at FP in eq. (82) in Steen and Faltinsen (1995) *****
            # Heave DOF
            A_mat[2, 0, :] += rho_a*k_2_FP/2*a_0_FP * (-rho_0/p_0)*1j*np.multiply(omega_e, a_3j) * r_j(x_g_FP, j, L)

            # Pitch DOF
            A_mat[2, 1, :] += rho_a*k_2_FP/2*a_0_FP * (-rho_0/p_0)*1j*np.multiply(omega_e, a_5j) * r_j(x_g_FP, j, L)

            # Uniform pressure DOF
            A_mat[2, 2, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_0j) * r_j(x_g_FP, j, L)

            # Wave excitation term
            f_vec[2, :] -= rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(np.multiply(omega_e, a_7j) * r_j(x_g_FP, j, L), zeta_a)

            # ***** Terms from fan inlet at lcg_fan in eq. (82) in Steen and Faltinsen (1995) *****
            # Heave DOF
            A_mat[2, 0, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_3j) * r_j(lcg_fan, j, L)

            # Pitch DOF
            A_mat[2, 1, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_5j) * r_j(lcg_fan, j, L)

            # Uniform pressure DOF
            A_mat[2, 2, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_0j) * r_j(lcg_fan, j, L)

            # Wave excitation term
            f_vec[2, :] -= (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(np.multiply(omega_e, a_7j) * r_j(lcg_fan, j, L), zeta_a)

        elif j % 2 == 0:  # if j is even
            # ***** Terms from leakage at AP in eq. (82) in Steen and Faltinsen (1995) *****
            # Heave DOF
            A_mat[2, 0, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_3j) * r_j(x_g_AP, j, L)

            # Pitch DOF
            A_mat[2, 1, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_5j) * r_j(x_g_AP, j, L)

            # Uniform pressure DOF
            A_mat[2, 2, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_0j) * r_j(x_g_AP, j, L)

            # Wave excitation term
            f_vec[2, :] -= rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(np.multiply(omega_e, b_7j) * r_j(x_g_AP, j, L), zeta_a)

            # ***** Terms from leakage at FP in eq. (82) in Steen and Faltinsen (1995) *****
            # Heave DOF
            A_mat[2, 0, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_3j) * r_j(x_g_FP, j, L)

            # Pitch DOF
            A_mat[2, 1, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_5j) * r_j(x_g_FP, j, L)

            # Uniform pressure DOF
            A_mat[2, 2, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_0j) * r_j(x_g_FP, j, L)

            # Wave excitation term
            f_vec[2, :] -= rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(np.multiply(omega_e, b_7j) * r_j(x_g_FP, j, L), zeta_a)

            # ***** Terms from fan inlet at x_F in eq. (82) in Steen and Faltinsen (1995) *****
            # Heave DOF
            A_mat[2, 0, :] += (-rho_0*p_0*dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_3j) * r_j(lcg_fan, j, L)

            # Pitch DOF
            A_mat[2, 1, :] += (-rho_0*p_0*dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_5j) * r_j(lcg_fan, j, L)

            # Uniform pressure DOF
            A_mat[2, 2, :] += (-rho_0*p_0*dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_0j) * r_j(lcg_fan, j, L)

            # Wave excitation term
            f_vec[2, :] -= (-rho_0*p_0*dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(np.multiply(omega_e, b_7j) * r_j(lcg_fan, j, L), zeta_a)

        '''
        print(np.size(A_mat == A_mat_2) - np.count_nonzero(A_mat == A_mat_2))
        print(np.size(f_vec == f_vec_2) - np.count_nonzero(f_vec == f_vec_2))
        '''

        # Solves linear system of equations for each frequency
        eta_3a, eta_5a, mu_ua = solve_linear_systems_of_eq(A_mat, f_vec)

        eta_3_test.append(eta_3a[4])  # Append value to test vector
        if not counter == 0:  # compute error between current and previous if it is not the first iteration
            err = np.amax([np.abs(np.divide(eta_3a - eta_3a_old, eta_3a_old)),
                           np.abs(np.divide(eta_5a - eta_5a_old, eta_5a_old)),
                           np.abs(np.divide(mu_ua - mu_ua_old, mu_ua_old))])

        # Stores current solution for comparison in next iteration
        eta_3a_old = eta_3a
        eta_5a_old = eta_5a
        mu_ua_old = mu_ua

        counter += 1  # increment counter
        print("Current iteration:", counter, 'with error:', err)
        # Compute new value for bias and gain for the linearized leakage for the next iteration

        b_L_AP = eta_3m + L/2*eta_5m - h_s_AP  # TODO: Might want to do this in seperate sub-routines
        b_L_FP = eta_3m - L/2*eta_5m - h_s_FP
        sigma_L_AP = rms_leakage(-L/2, omega_0, eta_3a, eta_5a, H_s, T_p, zeta_a)
        sigma_L_FP = rms_leakage(L/2, omega_0, eta_3a, eta_5a, H_s, T_p, zeta_a)

        n_R_AP = N_R(b_L_AP, b_L_AP)
        n_R_FP = N_R(b_L_FP, b_L_FP)
        n_B_AP = N_B(b_L_AP, b_L_AP)
        n_B_FP = N_B(b_L_FP, b_L_FP)


plt.plot(f_encounter, np.abs(eta_3a),  label='eta_3')
plt.plot(f_encounter, np.abs(eta_5a),  label='eta_5')
plt.plot(f_encounter, np.abs(mu_ua),  label='mu_ua')
plt.xlabel('Encounter frequency [Hz]')
plt.legend()
plt.show()

print("The iteration scheme converged after", counter, "iterations and the relative error was then", err)
