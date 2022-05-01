from scipy.stats import norm
import numpy as np

# Functions to calculate modal solution of Boundary value problem of the air cushion domain


def A_0j(j, b, L, p_0, dQdp_0, x_F, k_2_AP, k_2_FP, A_0_AP, A_0_FP, k_4):
    """
    Computes frequency dependent modal amplitude A_0j found in eq. (44) in Steen and Faltinsen (1995). 'Cobblestone
    Oscillations of an SES with Flexible Bag Aft Seal'
    :param j: (int)
        [-] Order of acoustic mode
    :param b: (double)
        [m] Air cushion beam
    :param L: (double)
        [m] Air cushion length
    :param p_0: (double)
        [Pa] Mean cushion pressure
    :param dQdp_0: (double)
        [m^2/s] Linear fan slope #TODO: need to check what unit to use
    :param x_F:
        [m] Longitudinal position of fan relative to geometrical center of air cushion
    :param k_2_AP: (double)
         [m/s] K_2 constant at AP
    :param k_2_FP: (double)
        [m/s] K_2 constant at FP
    :param A_0_AP: (double)
        [m^2] Mean leakage area at AP
    :param A_0_FP: (double)
         [m^2] Mean leakage area at FP
    :param k_4: (double)
        Constant K_4 found in eq. (27) (depend on j)
    :return: (double)
        Frequency dependent modal amplitude of odd mode j due to change in uniform pressure DOF
    """
    return k_4/L/b*(k_2_AP*A_0_AP - k_2_FP*A_0_FP - 2*p_0*dQdp_0*np.cos(j*np.pi/L * (x_F + L/2)))


def B_0j(j, b, L, p_0, dQdp_0, x_F, k_2_AP, k_2_FP, A_0_AP, A_0_FP, k_4):
    """
    Computes frequency dependent modal amplitude B_0j found in eq. (45) in Steen and Faltinsen (1995). 'Cobblestone
    Oscillations of an SES with Flexible Bag Aft Seal'
    :param j: (int)
        [-] Order of acoustic mode
    :param b: (double)
        [m] Air cushion beam
    :param L: (double)
        [m] Air cushion length
    :param p_0: (double)
        [Pa] Mean cushion pressure
    :param dQdp_0: (double)
        [m^2/s] Linear fan slope #TODO: need to check what unit to use
    :param x_F:
        [m] Longitudinal position of fan relative to geometrical center of air cushion
    :param k_2_AP: (double)
         [m/s] K_2 constant at AP
    :param k_2_FP: (double)
        [m/s] K_2 constant at FP
    :param A_0_AP: (double)
        [m^2] Mean leakage area at AP
    :param A_0_FP: (double)
         [m^2] Mean leakage area at FP
    :param k_4: (double)
        Constant K_4 found in eq. (27) (depend on j)
    :return: (double)
        Frequency dependent modal amplitude of even mode j due to change in uniform pressure DOF
    """
    return k_4/L/b*(k_2_AP*A_0_AP + k_2_FP*A_0_FP - 2*p_0*dQdp_0*np.cos(j*np.pi/L * (x_F + L/2)))


def A_3j(L, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP):
    """
    Computes frequency dependent modal amplitude A_3j found in eq. (48) in Steen and Faltinsen (1995). 'Cobblestone
    Oscillations of an SES with Flexible Bag Aft Seal'
    :param L: (double)
        [m] Air cushion length
    :param k_4:
        Constant K_4 found in eq. (27) (depend on j)
    :param k_2_AP:
        [m/s] K_2 constant at AP
    :param k_2_FP:
        [m/s] K_2 constant at FP
    :param n_R_AP:
        [-] gain-value of quasi-linearized variable leakage area at aft (between 0 and 1)
    :param n_R_FP:
        [-] gain-value of quasi-linearized variable leakage area at bow (between 0 and 1)
    :return: (double)
        Frequency dependent modal amplitude of odd mode j due to change heave DOF
    """
    return 2 * k_4 / L * (k_2_AP * n_R_AP - k_2_FP * n_R_FP)


def B_3j(L, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP):
    """
    Computes frequency dependent modal amplitude B_3j found in eq. (49) in Steen and Faltinsen (1995). 'Cobblestone
    Oscillations of an SES with Flexible Bag Aft Seal'
    :param L: (double)
        [m] Air cushion length
    :param k_4: (double)
        Constant K_4 found in eq. (27) (depend on j)
    :param k_2_AP: (double)
        [m/s] K_2 constant at AP
    :param k_2_FP: (double)
        [m/s] K_2 constant at FP
    :param N_R_AP: (double)
        [-] gain-value of quasi-linearized variable leakage area at aft (between 0 and 1)
    :param N_R_FP: (double)
        [-] gain-value of quasi-linearized variable leakage area at bow (between 0 and 1)
    :return: (double)
        Frequency dependent modal amplitude of even mode j due to change heave DOF
    """
    return 2 * k_4 / L * (k_2_AP * n_R_AP + k_2_FP * n_R_FP)


def A_5j(j, L, omega_e, k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP):
    """
    Computes frequency dependent modal amplitude A_5j found in eq. (50) in Steen and Faltinsen (1995). 'Cobblestone
    Oscillations of an SES with Flexible Bag Aft Seal'
    :param j: (int)
        [-] Order of acoustic mode
    :param L: (double)
        [m] Air cushion length
    :param omega_e: (double)
        [rad/s] Encounter frequency
    :param k_4: (double)
        Constant K_4 found in eq. (27) (depend on j)
    :param k_2_AP: (double)
        [m/s] K_2 constant at AP
    :param k_2_FP: (double)
        [m/s] K_2 constant at FP
    :param N_R_AP: (double)
        [-] gain-value of quasi-linearized variable leakage area at aft (between 0 and 1)
    :param N_R_FP: (double)
        [-] gain-value of quasi-linearized variable leakage area at bow (between 0 and 1)
    :return: (double)
        Frequency dependent modal amplitude of odd mode j due to change pitch DOF
    """
    return 4*L*k_4/(j*np.pi)**2 * 1j * omega_e + k_4 * (k_2_AP * n_R_AP + k_2_FP * n_R_FP)


def B_5j(k_4, k_2_AP, k_2_FP, n_R_AP, n_R_FP):
    """
    Computes frequency dependent modal amplitude B_5j found in eq. (51) in Steen and Faltinsen (1995). 'Cobblestone
    Oscillations of an SES with Flexible Bag Aft Seal'
    :param k_4: (double)
        Constant K_4 found in eq. (27) (depend on j)
    :param k_2_AP: (double)
        [m/s] K_2 constant at AP
    :param k_2_FP: (double)
        [m/s] K_2 constant at FP
    :param N_R_AP: (double)
        [-] gain-value of quasi-linearized variable leakage area at aft (between 0 and 1)
    :param N_R_FP: (double)
        [-] gain-value of quasi-linearized variable leakage area at bow (between 0 and 1)
    :return: (double)
        Frequency dependent modal amplitude of even mode j due to change pitch DOF
    """
    return k_4 * (k_2_AP * n_R_AP - k_2_FP * n_R_FP)


def A_7j(j, k_4, L, k, omega_e, k_2_AP, k_2_FP, n_R_AP, n_R_FP):
    """
    Computes frequency dependent modal amplitude A_7j found in eq. (52) in Steen and Faltinsen (1995). 'Cobblestone
    Oscillations of an SES with Flexible Bag Aft Seal'
    :param j: (int)
        [-] Order of acoustic mode
    :param k_4: (double)
        Constant K_4 found in eq. (27) (depend on j)
    :param L: (double)
        [m] Air cushion length
    :param k: (double)
        [m^-1] Wave number of water waves
    :param omega_e: (double)
        [rad/s] Encounter frequency
    :param k_2_AP: (double)
        [m/s] K_2 constant at AP
    :param k_2_FP: (double)
        [m/s] K_2 constant at FP
    :param N_R_AP: (double)
        [-] gain-value of quasi-linearized variable leakage area at aft (between 0 and 1)
    :param N_R_FP: (double)
        [-] gain-value of quasi-linearized variable leakage area at bow (between 0 and 1)
    :return: (double)
        Frequency dependent modal amplitude of odd mode j due to water waves
    """
    return -1j*4*k_4/L * k * np.cos(k*L/2) / (k**2 - (j*np.pi/L)**2) * omega_e - \
        1j*2*k_4/L * (k_2_AP*n_R_AP*np.exp(-1j*k*L/2) - k_2_FP*n_R_FP*np.exp(1j*k*L/2))


def B_7j(j, k_4, L, k, omega_e, k_2_AP, k_2_FP, n_R_AP, n_R_FP):
    """
    Computes frequency dependent modal amplitude B_7j found in eq. (53) in Steen and Faltinsen (1995). 'Cobblestone
    Oscillations of an SES with Flexible Bag Aft Seal'
    :param j: (int)
        [-] Order of acoustic mode
    :param k_4: (double)
        Constant K_4 found in eq. (27) (depend on j)
    :param L: (double)
        [m] Air cushion length
    :param k: (double)
        [m^-1] Wave number of water waves
    :param omega_e: (double)
        [rad/s] Encounter frequency
    :param k_2_AP: (double)
        [m/s] K_2 constant at AP
    :param k_2_FP: (double)
        [m/s] K_2 constant at FP
    :param N_R_AP: (double)
        [-] gain-value of quasi-linearized variable leakage area at aft (between 0 and 1)
    :param N_R_FP: (double)
        [-] gain-value of quasi-linearized variable leakage area at bow (between 0 and 1)
    :return: (double)
        Frequency dependent modal amplitude of even mode j due to water waves
    """
    return -4*k_4/L * k * np.sin(k*L/2) / (k**2 - (j*np.pi/L)**2) * omega_e - \
        1j*2*k_4/L * (k_2_AP*n_R_AP*np.exp(-1j*k*L/2) + k_2_FP*n_R_FP*np.exp(1j*k*L/2))

# Computing constants used in derivation of equations

def K_1(rho_0, p_0, h_0, A_c, p_a=105325, gamma=1.4):
    """
    Computes constant K_1 found in eq. (83) in Steen and Faltinsen (1995). 'Cobblestone Oscillations of an SES with Flexible
    Bag Aft Seal'
    :param rho_0: (double)
        [kg/m^3] Density of air at mean cushion pressure
    :param p_0: (double)
        [Pa] Mean cushion pressure
    :param h_0: (double)
        [m] Cushion height
    :param A_c: (double)
        [m^2] Air cushion area
    :param p_a: (double) default=101325
        [Pa] Atmospheric pressure
    :param gamma: (double) default=1.4
        [-] Ratio of specific heat for air, when assuming adiabatic density relation
    :return: (double)
        [kg] K_1 constant
    """
    return rho_0 * h_0 * A_c / gamma / (1 + p_a/p_0)


def K_2(p_0, c_n=0.61, rho_a=1.225):
    """
    Computes constant K_2 found in eq. (27) in Steen and Faltinsen (1995). 'Cobblestone Oscillations of an SES with
    Flexible Bag Aft Seal'
    :param p_0: (double)
        [Pa] Mean cushion pressure
    :param c_n: (double) default=0.61
        [-] Orifice constant
    :param rho_a: (double)
        [kg/m^3] Density of air at atmospheric pressure
    :return: (double)
        [m/s] K_2 constant
    """
    return c_n * np.sqrt(2 * p_0 / rho_a)


def K_3(rho_0, p_0, Q_0, dQdp_0):
    """
    Computes constant K_3 found in eq. (27) in Steen and Faltinsen (1995). 'Cobblestone Oscillations of an SES with
    Flexible Bag Aft Seal'
    :param rho_0: (double)
        [kg/m^3] Density of air at atmospheric pressure
    :param p_0: (double)
        [Pa] Mean cushion pressure
    :param Q_0: (double)
        [m^3/s] Mean fan flow rate
    :param dQdp_0: (double)
        [m^2/2] Linear fan slope
    :return: (double)
        K_3 constant
    """
    return rho_0 * Q_0 / 2 - rho_0 * p_0 * dQdp_0


def Xi_j(j, rho_0, p_0, h_0, b, L, k_2_AP, k_2_FP, A_0_AP, A_0_FP, dQdp_0, x_F, c=343):
    """
    Computes constant xi_j found in eq. (55) in Steen and Faltinsen (1995). 'Cobblestone Oscillations of an SES with
    Flexible Bag Aft Seal'
    :param j: (int)
        [-] Order of acoustic mode
    :param rho_0: (double)
        [kg/m^3] Density of air at mean cushion pressure
    :param p_0: (double)
        [Pa] Mean cushion pressure #TODO: need to check what unit to use
    :param h_0: (double)
        [m] Cushion height
    :param b: (double)
        [m] Air cushion beam
    :param L: (double)
        [m] Air cushion length
    :param k_2_AP: (double)
        [m/s] K_2 constant at AP
    :param k_2_FP: (double)
        [m/s] K_2 constant at FP
    :param A_0_AP: (double)
        [m^2] Mean leakage area at AP
    :param A_0_FP: (double)
        [m^2] Mean leakage area at FP
    :param dQdp_0: (double)
        [m^2/s] Linear fan slope #TODO: need to check what unit to use
    :param x_F: (double)
        [m] Longitudinal position of fan relative to geometrical center of air cushion
    :param c: (double) default=343
        [m/s] Speed of sound in air
    :return: (double)
        [-] Relative damping ratio of mode j
    """
    return c*rho_0/(j*np.pi*h_0*b)*((k_2_AP*A_0_FP + k_2_FP*A_0_AP)/2/p_0 - dQdp_0 * np.cos(j*np.pi/L * (x_F + L/2))**2)


def Omega_j(j, L, c=343):
    """
    Computes the eigenfrequency of the acoustic mode j
    :param j: (int)
        [-] Order of acoustic mode
    :param L: (double)
        [m] Air cushion length
    :param c: (double) default=343
        [m/s] Speed of sound in air
    :return: (double)
        [rad/s] Eigenfrequency of the acoustic mode j
    """
    return c*j*np.pi/L


def K_4(xi_j, h_0, omega_e, omega_j, c=343):
    """
    Computes constant K_4 found in eq. (54) in Steen and Faltinsen (1995). 'Cobblestone Oscillations of an SES with
    Flexible Bag Aft Seal'
    :param xi_j: (double)
        [-] Relative damping ratio of mode j
    :param h_0: (double)
        [m] Cushion height
    :param omega_e: (double)
        [rad/s] Encounter frequency
    :param omega_j: (double)
        [rad/s] Acoustic resonance of mode j
    :param c: (double) default=343
        [m/s] Speed of sound
    :return: (double)
        Constant K_4 found in eq. (27)
    """
    return c**2 / h_0 / (-omega_e**2 + 2*xi_j*omega_j*1j*omega_e + omega_j**2)

# Help functions used for calculation of the linearized variable leakage


def solve_mean_value_relation(n_B_AP, n_B_FP, L, b, x_cp, A_c, p_0, k_2, k_3, h_s_AP, h_s_FP, C_33, C_55, C_35, C_53, rho_a=1.225):
    """
    Solves linear system of equations for mean value relation in eq. (3.27) in Steen (1993) 'Cobblestone effect on SES'
    :param n_B_AP: (double)
        [-] bias of quasi-linearized leakage area aft
    :param n_B_FP: (double)
        [-] bias of quasi-linearized leakage area at the bow
    :param L: (double)
        [m] Air cushion length
    :param b: (double)
        [m] Air cushion width
    :param x_cp: (double)
        [m] Coordinate from center of gravity of air cushion center in x-direction
    :param A_c: (double)
        [m^2] Air cushion area
    :param p_0: (double)
        [Pa] Mean cushion pressure
    :param k_2: (double)
        [m/s] K_2 constant at AP
    :param k_3: (double)
        K_3 constant eq. (27)
    :param h_s_AP: (double)
        [m] aft seal submergence
    :param h_s_FP: (double)
        [m] bow seal submergence
    :param C_33: (double)
        Restoring in heave due to heave motion
    :param C_55: (double)
        Restoring in pitch due to pitch motion
    :param C_35: (double)
        Restoring in heave due to pitch motion
    :param C_53: (double)
        Restoring in pitch due to heave motion
    :param rho_a: (double)
        [kg/m^3] Density of air at atmospheric pressure
    :return: mu_3m (double), mu_5m (double)
        [m] Mean value of heave, [rad] Mean value of pitch
    """
    # Define system of equations
    a = np.array([[C_33, C_35, -A_c * p_0],
                  [C_53, C_55, -x_cp * A_c * p_0],
                  [rho_a*k_2*b*(n_B_AP + n_B_FP), rho_a*k_2*b*L/2*(n_B_AP - n_B_FP), k_3]])

    f = np.array([0, 0, rho_a*k_2*b*(n_B_AP*h_s_AP + n_B_FP*h_s_FP)])

    # Solve linear system of equations
    x = np.linalg.solve(a, f)

    mu_3m = x[0]  # [m] Mean value of heave
    mu_5m = x[1]  # [rad] Mean value of pitch
    #nu_um = x[2]  # [-] Mean value of dim. less dynamic uniform pressure

    return mu_3m, mu_5m



def N_R(b_L, sigma_L):
    """
    Computes variable N_R in eq. (3.21) in Steen (1993) 'Cobblestone effect on SES'. Uses the normal distributed
    cumulative density function from scipy.stats
    :param b_L: (double)
        Bias term of variable leakage
    :param sigma_L: (double)
        rms-value of the time dependent part of the variable leakage area
    :return: (double)
        [-] gain value of quasi-linearized variable leakage (always between 0 and 1)
    """
    return norm.cdf(b_L/sigma_L)


def N_B(b_L, sigma_L):
    """
    Computes variable N_B in eq. (3.21) in Steen (1993) 'Cobblestone effect on SES'. Uses the normal distributed
    probability and cumulative density function from scipy.stats
    :param b_L: (double)
        Bias term of variable leakage
    :param sigma_L: (double)
        rms-value of the time dependent part of the variable leakage area
    :return: (double)
        bias of the quasi-linearized leakage area (always between 0 and 1)
    """
    return norm.cdf(b_L/sigma_L) + sigma_L/b_L * norm.pdf(b_L/sigma_L)