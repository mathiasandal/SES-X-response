import numpy as np
'''
This file contains functions used to calculate restoring coefficients and wave excitation forces on the seals at the
front and aft of a conventional SES. The terms are derived in p. 22-24 of the SES Extension of the Veres manual.

'''



def c_33_finger_seal(tau, b, p_0):
    """
    Calculates the restoring in heave for a finger seal type.
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
    Calculates the restoring in heave for a lobe bag type of seal.
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


def c_55_finger_seal(tau, b, p_0, x):
    """
    Calculates the restoring in heave for a finger seal type.
    :param tau: (float)
        Angle of the finger seal in deg
    :param b: (float)
        Beam of the finger seal in m
    :param p_0: (float)
        Excess pressure in the air cushion in Pa
    :param x: (float)
        Longitudinal position of the finger seal relative to motion coord. system
    :return: (float)
        Restoring coefficient in heave due to the finger seal
    """
    return p_0 * b / np.sin(np.deg2rad(tau))


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
        Longitudinal position of the finger seal relative to motion coord. system
    :return:
        Restoring coefficient in heave due to the finger seal
    """
    return p * b / np.sin(np.deg2rad(tau))
