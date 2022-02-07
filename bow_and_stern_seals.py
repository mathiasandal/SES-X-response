import numpy as np

'''
This file contains functions used to calculate restoring coefficients and wave excitation forces on the seals at the
front and aft of a conventional SES. The terms are derived in p. 22-24 of the SES Extension of the Veres manual.

'''


def c_33_finger_seal(tau, b, p_0):
    """
    Calculates the restoring coefficient in heave for a finger seal type.

    :param tau: (float)
        Angle of the finger seal in deg
    :param b: (float)
        Beam of the finger seal in m
    :param p_0: (float)
        Excess pressure in the air cushion
    :return: (float)
        Restoring coefficient in heave due to the finger seal
    """
    return p_0 * b / np.sin(np.deg2rad(tau))


def c_33_lobe_bag_seal(tau, b, p):
    """
    Calculates the restoring coefficient in heave for a lobe bag type of seal.

    :param tau: (float)
        Angle of the lobe bag seal in deg
    :param b: (float)
        Beam of the lobe bag seal in m
    :param p: (float)
        Excess pressure in the lobe bag
    :return: (float)
        Restoring coefficient in heave due to the lobe bag seal
    """
    return p * b / np.sin(np.deg2rad(tau))


def c_35_finger_seal(tau, b, p_0, x):
    """
    Calculates the coupling restoring coefficient in heave and pitch for a finger seal type.

    :param tau: (float)
        Angle of the finger seal in deg
    :param b: (float)
        Beam of the finger seal in m
    :param p_0: (float)
        Excess pressure in the air cushion in Pa
    :param x: (float)
        Longitudinal position of the finger seal relative to motion coord. system
    :return: (float)
        Coupling restoring coefficient in heave and pitch due to the finger seal
    """
    return -x * p_0 * b / np.sin(np.deg2rad(tau))


def c_35_lobe_bag_seal(tau, b, p, x):
    """
    Calculates the coupling restoring coefficient in heave and pitch for a lobe bag type of seal.

    :param tau: (float)
        Angle of the lobe bag seal in deg
    :param b: (float)
        Beam of the lobe bag seal in m
    :param p:
        Excess pressure in the lobe bag
    :param x: (float)
        Longitudinal position of the finger seal relative to motion coord. system
    :return:
        Coupling restoring coefficient in heave and pitch due to the lobe bag seal
    """
    return -x * p * b / np.sin(np.deg2rad(tau))


def c_55_finger_seal(tau, b, p_0, x):
    """
    Calculates the restoring coefficient in pitch for a finger seal type.

    :param tau: (float)
        Angle of the finger seal in deg
    :param b: (float)
        Beam of the finger seal in m
    :param p_0: (float)
        Excess pressure in the air cushion in Pa
    :param x: (float)
        Longitudinal position of the finger seal relative to motion coord. system
    :return: (float)
        Restoring coefficient in pitch due to the finger seal
    """
    return x ** 2 * p_0 * b / np.sin(np.deg2rad(tau))


def c_55_lobe_bag_seal(tau, b, p, x):
    """
    Calculates the restoring in heave for a lobe bag type of seal.

    :param tau: (float)
        Angle of the lobe bag seal in deg
    :param b: (float)
        Beam of the lobe bag seal in m
    :param p:
        Excess pressure in the lobe bag
    :param x: (float)
        Longitudinal position of the lobe bag type of seal relative to motion coord. system
    :return:
        Restoring coefficient in pitch due to the lobe bag type of seal
    """
    return x ** 2 * p * b / np.sin(np.deg2rad(tau))


def restoring_finger_at_bow_lobe_bag_at_the_stern(b, tau_b, tau_s, p_0, p_s, x_b, x_s):
    """
    Calculates restoring coefficients for a conventional SES with finger seals at the bow and a lobe bag type of seal at
    at the aft.

    :param b: (float)
        Beam of the seals in m. It is assumed that both seals have the same beam.
    :param tau_b: (float)
        Angle of the finger seal in deg
    :param tau_s: (float)
        Angle of the lobe bag seal in deg
    :param p_0: (float)
        Excess pressure in the air cushion in Pa
    :param p_s: (float)
        Excess pressure in the lobe bag
    :param x_b: (float)
        Longitudinal position of the finger seal relative to motion coord. system
    :param x_s: (float)
        Longitudinal position of the lobe bag type of seal relative to motion coord. system
    :return: (float), (float), (float)
        Restoring in heave, coupling between heave and pitch, and restoring in pitch
    """

    c_33 = c_33_finger_seal(tau_b, b, p_0) + c_33_lobe_bag_seal(tau_s, b, p_s)  # [N/m] total restoring in heave
    c_35 = c_35_finger_seal(tau_b, b, p_0, x_b) + c_35_lobe_bag_seal(tau_s, b, p_s,
                                                                     x_s)  # [N] coupling of rest. in heave and pitch
    c_55 = c_55_finger_seal(tau_b, b, p_0, x_b) + c_55_lobe_bag_seal(tau_s, b, p_s,
                                                                     x_s)  # [Nm] total restoring in pitch

    return c_33, c_35, c_55


def excitation_finger_at_bow_lobe_bag_at_the_stern(b, tau_b, tau_s, p_0, p_s, x_b, x_s, omega, beta, zeta_a=1, g=9.81):
    """
    Computes complex excitation for a SES with finger seals at the bow and lobe bag at the aft. The function is
    implemented according to the formulas presented in the Veres SES extension p. 24

    :param b: (float)
        Beam of the seals in m. It is assumed that both seals have the same beam.
    :param tau_b: (float)
        Angle of the finger seal in deg
    :param tau_s: (float)
        Angle of the lobe bag seal in deg
    :param p_0: (float)
        Excess pressure in the air cushion in Pa
    :param p_s: (float)
        Excess pressure in the lobe bag in Pa
    :param x_b: (float)
        Longitudinal position of the finger seal relative to motion coord. system
    :param x_s: (float)
        Longitudinal position of the lobe bag type of seal relative to motion coord. system
    :param omega: (float)
        frequency of encounter in rad/s
    :param beta: (float)
        Wave heading in deg. Beta = 0 means head sea
    :param zeta_a: (float)
        Wave amplitude in meter
    :param g: (float)
        Acceleration of gravity
    :return: (float), (float)
        Complex amplitude of heave excitation, complex amplitude of pitch excitation
    """

    k = omega ** 2 / g  # calculate wave number

    # excitation force amplitude in heave
    f_3 = -1j * (p_0 * b / np.sin(np.deg2rad(tau_b)) * np.exp(1j * k * x_b * np.cos(np.deg2rad(beta)))
                 + p_s * b / np.sin(np.deg2rad(tau_s)) * np.exp(1j * k * x_s * np.cos(np.deg2rad(beta)))) * zeta_a

    # excitation moment amplitude in pitch
    f_5 = 1j * (x_b * p_0 * b / np.sin(np.deg2rad(tau_b)) * np.exp(1j * k * x_b * np.cos(np.deg2rad(beta)))
                + x_s * p_s * b / np.sin(np.deg2rad(tau_s)) * np.exp(1j * k * x_s * np.cos(np.deg2rad(beta)))) * zeta_a

    return f_3, f_5


def excitation_skirts(b, tau_b, tau_s, p_0, p_s, x_b, x_s, omegas, beta, zeta_a=1, g=9.81):
    """
    Computes a list of complex excitation amplitudes for a SES with finger seals at the bow and lobe bag at the aft
    corresponding to the given list of input encounter frequencies. The function is implemented according to the
    formulas presented in the Veres SES extension p. 24

    :param b: (float)
        Beam of the seals in m. It is assumed that both seals have the same beam.
    :param tau_b: (float)
        Angle of the finger seal in deg
    :param tau_s: (float)
        Angle of the lobe bag seal in deg
    :param p_0: (float)
        Excess pressure in the air cushion in Pa
    :param p_s: (float)
        Excess pressure in the lobe bag in Pa
    :param x_b: (float)
        Longitudinal position of the finger seal relative to motion coord. system
    :param x_s: (float)
        Longitudinal position of the lobe bag type of seal relative to motion coord. system
    :param omegas: (float) (n_frequencies X 1) numpy array
        list of frequencies of encounter in rad/s
    :param beta: (float)
        Wave heading in deg. Beta = 0 means head sea
    :param zeta_a: (float)
        Wave amplitude in meter
    :param g: (float)
        Acceleration of gravity
    :return: (float) (n_frequencies X 1) numpy array, (float) (n_frequencies X 1) numpy array
        Complex amplitude of heave excitation, complex amplitude of pitch excitation
    """

    n_frequencies = len(omegas)  # number of frequencies in the input list of encounter frequencies

    # initialize lists to contain excitation force in heave an excitation moment in pitch for different frequencies
    f_3_seals = np.zeros([n_frequencies], dtype=complex)
    f_5_seals = np.zeros([n_frequencies], dtype=complex)

    for i in range(n_frequencies):
        f_3_seals[i], f_5_seals[i] = excitation_finger_at_bow_lobe_bag_at_the_stern(b, tau_b, tau_s, p_0, p_s, x_b, x_s, omegas[i], beta, zeta_a, g)

    return f_3_seals, f_5_seals
