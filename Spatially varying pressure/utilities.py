from scipy.stats import norm
from scipy import integrate
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
        [(m^3/s)/Pa] Linear fan slope
    :param x_F: (double)
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
        [(m^3/s)/Pa] Linear fan slope
    :param x_F: (double)
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
    :param n_R_AP: (double)
        [-] gain-value of quasi-linearized variable leakage area at aft (between 0 and 1)
    :param n_R_FP: (double)
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
    :param n_R_AP: (double)
        [-] gain-value of quasi-linearized variable leakage area at aft (between 0 and 1)
    :param n_R_FP: (double)
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
    :param n_R_AP: (double)
        [-] gain-value of quasi-linearized variable leakage area at aft (between 0 and 1)
    :param n_R_FP: (double)
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
    :param n_R_AP: (double)
        [-] gain-value of quasi-linearized variable leakage area at aft (between 0 and 1)
    :param n_R_FP: (double)
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
    :param n_R_AP: (double)
        [-] gain-value of quasi-linearized variable leakage area at aft (between 0 and 1)
    :param n_R_FP: (double)
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
        [(m^3/s)/Pa] Linear fan slope
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
        [Pa] Mean cushion pressure
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
        [(m^3/s)/Pa] Linear fan slope #TODO: need to check what unit to use
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


def solve_mean_value_relation(n_B_AP, n_B_FP, L, b, x_cp, A_c, p_0, k_2_AP, k_2_FP, k_3, h_s_AP, h_s_FP, C_33, C_55, C_35, C_53, rho_a=1.225):
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
    :param k_2_AP: (double)
        [m/s] K_2 constant at AP
    :param k_2_FP: (double)
        [m/s] K_2 constant at FP
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
                  [rho_a*b*(k_2_AP*n_B_AP + k_2_FP*n_B_FP), rho_a*b*L/2*(k_2_AP*n_B_AP - k_2_FP*n_B_FP), k_3]])

    f = np.array([0, 0, rho_a*b*(k_2_AP*n_B_AP*h_s_AP + k_2_FP*n_B_FP*h_s_FP)])

    # Solve linear system of equations
    x = np.linalg.solve(a, f)

    eta_3m = x[0]  # [m] Mean value of heave
    eta_5m = x[1]  # [rad] Mean value of pitch
    #mu_um = x[2]  # [-] Mean value of dim. less dynamic uniform pressure

    return eta_3m, eta_5m


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


def PM_spectrum(omega_0, H_s, T_p):  # TODO: Change to a modified wave spectrum
    """
    Computes the modified Pierson-Moskowitz spectrum at a given wave frequency and input parameters.
    Source: https://www.orcina.com/webhelp/OrcaFlex/Content/html/Waves,Wavespectra.htm#:~:text=The%20ISSC%20spectrum%20(
    also%20known,Hs%2C%20are%20data%20items.
    :param omega_0: (1xn) vector
        [rad/s] Frequency
    :param H_s: (double)
        [m] significant wave height
    :param T_p: (double)
        [s] peak wave period
    :return: (double)
        [-] modified Pierson-Moskowitz spectrum evaluated at a wave frequency of omega with the given input parameters.
    """

    return 5/16 * (2*np.pi)**2 * H_s**2 * T_p**-4 * \
           np.multiply(np.power(omega_0, -5), np.exp(-5/4*(T_p/2/np.pi)**-4 * np.power(omega_0, -4)))


def rms_leakage(x, omega_0s, eta_3_amps, eta_5_amps, H_s, T_p, zeta_a, g=9.81):
    """
    Computes integral in eq. (3.28) in Steen 'Cobblestone effect on SES' numerically using Simpsons rule.
    :param x: (double)
        Longitudinal position along the vessel
    :param omega_0s: vector (1xn)
        Vector including circular frequency of the incident water waves
    :param eta_3_amps: vector (1xn)
        Vector of heave amplitudes corresponding to the respective incident water waves
    :param eta_5_amps: vector 1xn
        Vector of pitch amplitudes corresponding to the respective incident water waves
    :param H_s: (double)
        [m] significant wave height
    :param T_p: (double)
        [s] peak wave period
    :param zeta_a: vector (1xn)
        [m] wave amplitude
    :param g: (double) default=9.81
        [m/s^2] Acceleration of gravity
    :return: vector 1xn
        Returns sigma_L in eq. (3.28) in Steen 'Cobblestone effect on SES'
    """

    k = np.divide(np.power(omega_0s, 2), g)  # computes wave number of water waves

    integrand = np.multiply(np.power(np.absolute(np.divide(eta_3_amps - x*eta_5_amps + 1j*np.multiply(zeta_a, np.exp(1j*k*x)), zeta_a)), 2), PM_spectrum(omega_0s, H_s, T_p))

    return np.sqrt(integrate.simpson(integrand, omega_0s))


def A_0_AP(L, b, n_b_AP, eta_3m, eta_5m, h_s_AP, A_c_AP=0):
    """
    Computes mean leakage area at AP from eq. (3.31) 'Cobblestone effect on SES'
    :param L: (double)
        [m] Air cushion length
    :param b: (double)
        [m] Air cushion beam
    :param n_b_AP: (double)
        [-] bias of quasi-linearized leakage area aft
    :param eta_3m: (double)
        [m] Mean value of heave
    :param eta_5m: (double)
        [rad] Mean value of pitch
    :param h_s_AP: (double)
        [m] aft seal submergence
    :param A_c_AP: (double) default=0
        [m^2] Mean leakage part independent of mean position
    :return: (double)
        Mean leakage area at AP
    """
    return b*n_b_AP*(eta_3m - h_s_AP + L/2*eta_5m) + A_c_AP


def A_0_FP(L, b, n_b_FP, eta_3m, eta_5m, h_s_FP, A_c_FP=0):
    """
    Computes mean leakage area at AP from eq. (3.32) 'Cobblestone effect on SES'
    :param L: (double)
        [m] Air cushion length
    :param b: (double)
        [m] Air cushion beam
    :param n_b_FP: (double)
        [-] bias of quasi-linearized leakage area at the front
    :param eta_3m: (double)
        [m] Mean value of heave
    :param eta_5m: (double)
        [rad] Mean value of pitch
    :param h_s_FP: (double)
        [m] aft seal submergence
    :param A_c_FP: (double) default=0
        [m^2] Mean leakage part independent of mean position
    :return: (double)
        Mean leakage area at FP
    """
    return b*n_b_FP*(eta_3m - h_s_FP - L/2*eta_5m) + A_c_FP


def r_j(x, j, L):
    """
    Computes mode shape of order j evaluated at x
    :param x: (double)
        Longitudinal position along air cushion
    :param j: (int)
        Order of mode shape
    :param L: (double)
        [m] Air cushion length
    :return: (double)
        [-] Mode shape evaluated at x
    """
    return np.cos(j*np.pi/L*(x + L/2))

'''
def append_spatially_varying_terms(A, f, omega_e, j, b, L, rho_0, p_0, dQdp_0, lcg_fan, k_2_AP, a_0_AP, x_g_AP,
                                   k_2_FP, a_0_FP, x_g_FP, zeta_a, a_0j, a_3j, a_5j, a_7j, b_0j, b_3j, b_5j, b_7j,
                                   rho_a=101325):

    if j % 2 == 1:  # if j is odd
        # ***** Pitching moments due to spatially varying pressure in eq. (88) in Steen and Faltinsen (1995)
        # Heave DOF
        A[1, 0, :] += 2 * rho_0 * b * 1j * (L / j / np.pi) ** 2 * np.multiply(omega_e, a_3j)

        # Pitch DOF
        A[1, 1, :] += 2 * rho_0 * b * 1j * (L / j / np.pi) ** 2 * np.multiply(omega_e, a_5j)

        # Uniform pressure DOF
        A[1, 2, :] += 2 * rho_0 * b * 1j * (L / j / np.pi) ** 2 * np.multiply(omega_e, a_0j)

        # Wave excitation term
        f[1, :] -= 2 * rho_0 * b * 1j * (L / j / np.pi) ** 2 * np.multiply(np.multiply(omega_e, a_7j), zeta_a)

        # ***** Terms from leakage at AP in eq. (82) in Steen and Faltinsen (1995) *****
        # Heave DOF
        A[2, 0, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_3j) * r_j(x_g_AP, j, L)

        # Pitch DOF
        A[2, 1, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_5j) * r_j(x_g_AP, j, L)

        # Uniform pressure DOF
        A[2, 2, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_0j) * r_j(x_g_AP, j, L)

        # Wave excitation term
        f[2, :] -= rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(np.multiply(omega_e, a_7j) * r_j(x_g_AP, j, L), zeta_a)

        # ***** Terms from leakage at FP in eq. (82) in Steen and Faltinsen (1995) *****
        # Heave DOF
        A[2, 0, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_3j) * r_j(x_g_FP, j, L)

        # Pitch DOF
        A[2, 1, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_5j) * r_j(x_g_FP, j, L)

        # Uniform pressure DOF
        A[2, 2, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_0j) * r_j(x_g_FP, j, L)

        # Wave excitation term
        f[2, :] -= rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(np.multiply(omega_e, a_7j) * r_j(x_g_FP, j, L), zeta_a)

        # ***** Terms from fan inlet at lcg_fan in eq. (82) in Steen and Faltinsen (1995) *****
        # Heave DOF
        A[2, 0, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_3j) * r_j(lcg_fan, j, L)

        # Pitch DOF
        A[2, 1, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_5j) * r_j(lcg_fan, j, L)

        # Uniform pressure DOF
        A[2, 2, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, a_0j) * r_j(lcg_fan, j,  L)

        # Wave excitation term
        f[2, :] -= (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(np.multiply(omega_e, a_7j) * r_j(lcg_fan, j, L), zeta_a)

    elif j % 2 == 0:  # if j is even
        # ***** Terms from leakage at AP in eq. (82) in Steen and Faltinsen (1995) *****
        # Heave DOF
        A[2, 0, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_3j) * r_j(x_g_AP, j, L)

        # Pitch DOF
        A[2, 1, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_5j) * r_j(x_g_AP, j, L)

        # Uniform pressure DOF
        A[2, 2, :] += rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_0j) * r_j(x_g_AP, j, L)

        # Wave excitation term
        f[2, :] -= rho_a * k_2_AP / 2 * a_0_AP * (-rho_0 / p_0) * 1j * np.multiply(np.multiply(omega_e, b_7j) * r_j(x_g_AP, j, L), zeta_a)

        # ***** Terms from leakage at FP in eq. (82) in Steen and Faltinsen (1995) *****
        # Heave DOF
        A[2, 0, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_3j) * r_j(x_g_FP, j, L)

        # Pitch DOF
        A[2, 1, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_5j) * r_j(x_g_FP, j, L)

        # Uniform pressure DOF
        A[2, 2, :] += rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_0j) * r_j(x_g_FP, j, L)

        # Wave excitation term
        f[2, :] -= rho_a * k_2_FP / 2 * a_0_FP * (-rho_0 / p_0) * 1j * np.multiply(np.multiply(omega_e, b_7j) * r_j(x_g_FP, j, L), zeta_a)

        # ***** Terms from fan inlet at x_F in eq. (82) in Steen and Faltinsen (1995) *****
        # Heave DOF
        A[2, 0, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_3j) * r_j(lcg_fan, j, L)

        # Pitch DOF
        A[2, 1, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_5j) * r_j(lcg_fan, j, L)

        # Uniform pressure DOF
        A[2, 2, :] += (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(omega_e, b_0j) * r_j(lcg_fan, j, L)

        # Wave excitation term
        f[2, :] -= (-rho_0 * p_0 * dQdp_0) * (-rho_0 / p_0) * 1j * np.multiply(np.multiply(omega_e, b_7j) * r_j(lcg_fan, j, L), zeta_a)

    return A, f
'''

def Zeta_a(omega_0, H_s, T_p):
    """
    Returns the wave amplitudes related to the PM spectrum with the given input parameters H_s and T_p to the selected
    frequencies.
    :param omega_0: vector 1xn
        Vector including circular frequency of the incident water waves
    :param H_s: (double)
        [m] significant wave height
    :param T_p: (double)
        [s] peak wave period
    :return: (1xn) array
        [m] Returns a vector containing wave amplitudes related to each frequency in the wave spectrum
    """
    delta_omega_0 = np.diff(omega_0)
    delta_omega_0 = np.concatenate(([delta_omega_0[0]], delta_omega_0))

    return np.sqrt(2 * np.multiply(PM_spectrum(omega_0, H_s, T_p), delta_omega_0))


def solve_linear_systems_of_eq(A_mat, f_vec):
    """
    Solves linear system of equations for several cases
    :param A_mat: (3x3xn) array
        Matrix containing
    :param f_vec: (3xn) array

    :return: (1xn) array, (1xn) array, (1xn) array
        eta_3a, eta_5a, mu_ua
    """
    n = len(A_mat[0, 0, :])  # number of frequencies

    x_vec = np.zeros([3, n], dtype=complex)  # initialize vector to contain solution of linear system of equations

    for i in range(n):  # solve linear system of equations for each frequency
        x_vec[:, i] = np.linalg.solve(A_mat[:, :, i], f_vec[:, i])

    eta_3a = x_vec[0, :]  # Heave RAO
    eta_5a = x_vec[1, :]  # Pitch RAO
    mu_ua = x_vec[2, :]  # Uniform pressure RAO

    return eta_3a, eta_5a, mu_ua
