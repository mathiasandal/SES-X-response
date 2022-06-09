import numpy as np
import pandas as pd
from scipy.interpolate import interp1d


def air_cushion_area(l_rect, l_tri, b_c):
    """
    Computes projected area in zy-plane of a simplified air cushion and its centroid located from AP.
    It is assumed to be a rectangle with a triangle at the end.

            |________l_rect_____| |__l_tri_|
    _        ____________________
    |       |                      *
    |       |       Air cushion       *
    b       |             x              >
    |       |                         *
    _       |____________________ *

            |-----x_c----|

    Total area: b*l_1 + 0.5 * b*l_2

    :param l_rect: double [m]
        length of rectangular part of the air cushion
    :param l_tri: double [m]
        length of the triangular part of the air cushion
    :param b_c: double [m]
        beam of the air cushion
    :return:
    A_b: double
        Total area of the air cushion
    x_c: double
        Distance from centroid to AP
    """

    A_rect = b_c * l_rect  # [m^2] area of rectangle
    A_tri = 0.5 * b_c * l_tri  # [m^2] area of triangle

    # computes area of the air cushion
    S_0c = A_rect + A_tri  # [m^2] total area of air cushion

    # computes distance from AP to centroid
    x_c = (A_rect * 0.5 * l_rect + A_tri * (l_rect + l_tri / 3)) / S_0c

    return S_0c, x_c


def interpolate_fan_characteristics(p_0, p, Q):
    """
    Calculates the volume flow and the slope of the fan characteristics at a given air cushion pressure and a discrete
    array of pressures and volume flows charcterizing the fan.

    :param p_0: double
        Air cushion pressure at equilibrium
    :param p: double
        Array of pressures that the fan can operate at corresponding to a volume flow
    :param Q: double
        Array of volume flows corresponding to a air cushion pressure
    :return:
    Q_0: double
        Volume flow corresponding to the equilibrium
    dQdp_0: double
        Slope of the fan characteristics in the vicinity of the equilibrium
    """
    # TODO: Fails to interpolate the correct Q. Need to understand why.
    # Q_0 = np.interp(p_0, p, Q)  # Interpolates volume flow # This did not seem to work
    Q_0 = float(interp1d(p, Q)(p_0))  # This works at least for the first test

    dQdp_0 = 0  # initialize slope of fan characteristics

    # finds closest value to p_0 and its index in p
    p_0_closest, p_0_closest_index = find_closest_value(p, p_0)

    # Applies a central finite difference to estimate the slope near the working position (p_0, Q_0)
    if p_0_closest == p_0:
        dQdp_0 = (Q[p_0_closest_index + 1] - Q[p_0_closest_index - 1]) / (
                p[p_0_closest_index + 1] - p[p_0_closest_index - 1])
    elif p_0_closest < p_0:
        alpha = (p_0 - p[p_0_closest_index - 1]) / (p[p_0_closest_index] - p[p_0_closest_index - 1])
        dQdp_0 = (1 - alpha) * (Q[p_0_closest_index] - Q[p_0_closest_index - 2]) / (
                p[p_0_closest_index] - p[p_0_closest_index - 2]) + alpha * (
                         Q[p_0_closest_index + 1] - Q[p_0_closest_index - 1]) / (
                         p[p_0_closest_index + 1] - p[p_0_closest_index - 1])
    elif p_0_closest > p_0:
        alpha = (p_0 - p[p_0_closest_index]) / (p[p_0_closest_index + 1] - p[p_0_closest_index])
        dQdp_0 = (1 - alpha) * (Q[p_0_closest_index + 1] - Q[p_0_closest_index - 1]) / (
                p[p_0_closest_index + 1] - p[p_0_closest_index - 1]) + alpha * (
                         Q[p_0_closest_index + 2] - Q[p_0_closest_index]) / (
                         p[p_0_closest_index + 2] - p[p_0_closest_index])

    return Q_0, dQdp_0


def find_closest_value(arr, val):
    """
    Finds the closest value in a sorted array for a chosen value.

    :param arr: (1xn) numpy array
        Array which function is looping over
    :param val: double
        Value that the function will find the closest value to
    :return:
        closest_value: double
            Closest value the the specified input value
        closest_index: int
            Index of the closest value in the array
    """

    n = len(arr)
    index_closest = 0
    diff_closest = abs(100 * arr[-1] + val)

    for i in range(n):

        if abs(arr[i] - val) < diff_closest:
            index_closest = i
            diff_closest = abs(arr[i] - val)

    closest_value = arr[index_closest]

    return closest_value, index_closest


def read_fan_characteristics(filename, rpm='1800rpm'):
    """
    Reads a csv-file containing the fan characteristics and returns it for a give RPM.

    :param filename: str
        Directory of the csv file containing the fan characteristics.
    :param rpm: str
        String specifying which RPM that should be used.
    :return:
    Q: (nx1) numpy array
        Volume flow values for fan characteristic curve
    P: (nx1) numpy array
        Pressure values for fan characteristic curve
    """

    if rpm not in ['1000rpm', '1400rpm', '1600rpm', '1800rpm', '2108rpm']:
        raise TypeError

    # Reads csv-file and saves it as a pandas DataFrame
    df = pd.read_csv(filename)

    Q = df[['x']].to_numpy()  # gets all Q values from the DataFrame

    P = df[[rpm]].to_numpy()  # gets all P values from the DataFrame

    Q_cut = Q[-1]
    # Cuts the arrays where they intersect P = 0 [Pa]. Found by inspecting fan characteristics manually
    if rpm == '1000rpm':
        Q_cut = 12
    elif rpm == '1400rpm':
        Q_cut = 17
    elif rpm == '1600rpm':
        Q_cut = 19
    elif rpm == '1800rpm':
        Q_cut = 21.6
    elif rpm == '2108rpm':
        Q_cut = 25.5

    # Remove interval that is no well defined
    P = P[Q < Q_cut]
    Q = Q[Q < Q_cut]

    return Q, P, rpm


def wave_pumping_air_cushion(b, l_1, l_2, x_prime, beta, omega, g=9.81):

    # TODO: Add documentation

    k = omega**2/g  # calculate wave number

    accepted_error = 1e-9

    # Compute analytical solution to integral of spatial variation over the air cushion
    if np.abs(np.sin(np.deg2rad(beta))) < accepted_error:  # sin(beta) = 0
        I_1 = b/k/np.cos(np.deg2rad(beta)) * (np.exp(-1j * k * x_prime * np.cos(np.deg2rad(beta))) -
                                              np.exp(-1j * k * (x_prime - l_1) * np.cos(np.deg2rad(beta))))
        I_2 = 0
    elif np.abs(np.cos(np.deg2rad(beta))) < accepted_error:  # cos(beta) = 0
        I_1 = 2*np.sin(k*b/2*np.sin(np.deg2rad(beta)))/k/np.sin(np.deg2rad(beta))\
              * (x_prime*np.exp(-1j*k*x_prime*np.cos(np.deg2rad(beta))) -
                (x_prime - l_1)*np.exp(-1j*k*(x_prime - l_1)*np.cos(np.deg2rad(beta))))
        I_2 = 2*l_2/k/np.sin(np.deg2rad(beta))*(1 - np.cos(k*b/2 * np.sin(np.deg2rad(beta))))
    else:  # sin(beta) and cos(beta) is not equal zero
        I_1 = 2 * 1j * np.sin(k * b / 2 * np.sin(np.deg2rad(beta))) / (
                k ** 2 * np.sin(np.deg2rad(beta)) * np.cos(np.deg2rad(beta))) * \
                np.exp(-1j * k * x_prime * np.cos(np.deg2rad(beta)))*(1 - np.exp(1j*k*l_1*np.cos(np.deg2rad(beta))))
        I_2 = 4*l_2*np.exp(-1j*k*(x_prime - l_1)*np.cos(np.deg2rad(beta)))/(k*b**2*np.sin(np.deg2rad(beta))**2 -
              4*l_2**2*k*np.cos(np.deg2rad(beta))**2)*(b/2*np.sin(np.deg2rad(beta)) *
              (np.exp(-1j*k*l_2*np.cos(np.deg2rad(beta))) - np.cos(k*b/2*np.sin(np.deg2rad(beta)))) -
              1j*l_2*np.cos(np.deg2rad(beta))*np.sin(k*b/2*np.sin(np.deg2rad(beta))))

    f_7_hat = 1j*omega*(I_1 + I_2)

    return f_7_hat


def wave_pumping_sesx(x_f, x_s, y_s, x_b, omega_0, U, heading, zeta_a=1, g=9.81):
    """Computes wave pumping of a SES-X type of vessel assuming an undisturbed incident wave field. See section 2.6 """

    k = omega_0 ** 2 / g  # calculate wave number
    beta = np.deg2rad(heading)  # convert wave heading from degrees to radians
    omega_e = omega_0 + U/g*omega_0**2  # encounter frequency

    accepted_error = 1e-9

    if np.abs(np.sin(beta)) < accepted_error:  # case for sin(beta) close or equal zero
        # Compute first term of T_1
        T_11 = 2 * y_s / k / np.cos(beta) * (1 / (x_f - x_b) / k / np.cos(beta) + 1j) * \
              np.exp(-1j * k * x_f * np.cos(beta))
        # Compute second term of T_2
        T_12 = -2 * y_s / (x_f - x_b) / k**2 / np.cos(beta)**2 * np.exp(-1j * k * x_b * np.cos(beta))

        # Compute contribution of the triangular part of the air cushion for the wave pumping
        T_1 = T_11 + T_12

        # Compute contribution of the rectangular part of the air cushion for the wave pumping
        y_p = -y_s  # assume symmetry over xz-plane.
        T_2 = 1j * (y_s - y_p) / k / np.cos(beta) * (np.exp(-1j * k * x_s * np.cos(beta)) -
                                                     np.exp(-1j * k * x_f * np.cos(beta)))

    elif np.abs(np.cos(beta)) < accepted_error:  # case for cos(beta) close or equal zero
        # Compute contribution of the triangular part of the air cushion for the wave pumping
        T_1 = 2 * (np.cos(k * y_s * np.sin(beta)) - 1) / k**2 / (y_s / (x_f - x_b)) / np.sin(beta)**2

        # Compute contribution of the rectangular part of the air cushion for the wave pumping
        y_p = -y_s  # assume symmetry over xz-plane.
        T_2 = 1j * (x_s - x_f) / k / np.sin(beta) * (np.exp(-1j * k * y_s * np.sin(beta)) -
                                                     np.exp(-1j * k * y_p * np.sin(beta)))

    else:  # case for when sin(beta) and cos(beta) is not close or equal to zero
        # numerator of first term in T_1
        T_11_n = 2j * np.cos(beta) * np.sin(k * y_s * np.sin(beta)) * np.exp(-1j * k * x_f * np.cos(beta))
        # denominator of first term in T_1
        T_11_d = k**2 * np.sin(beta) * (np.cos(beta)**2 + (y_s / (x_f - x_b))**2 * np.sin(beta)**2)
        # numerator of second term in T_1
        T_12_n = 2 * y_s / (x_f - x_b) * (np.cos(k * y_s * np.sin(beta)) * np.exp(-1j * k * x_f * np.cos(beta)) -
                                          np.exp(-1j * k * x_b * np.cos(beta)))
        # denominator of second term in T_1
        T_12_d = k**2 * (np.cos(beta)**2 + (y_s / (x_f - x_b))**2 * np.sin(beta)**2)

        # Compute contribution of the triangular part of the air cushion for the wave pumping
        T_1 = T_11_n / T_11_d + T_12_n / T_12_d

        # Compute contribution of the rectangular part of the air cushion for the wave pumping
        y_p = -y_s  # assume symmetry over xz-plane.
        T_2 = 1/k**2/np.sin(beta)/np.cos(beta) * (- np.exp(-1j * k * (x_s * np.cos(beta) + y_s * np.sin(beta)))
                                                  + np.exp(-1j * k * (x_s * np.cos(beta) + y_p * np.sin(beta)))
                                                  + np.exp(-1j * k * (x_f * np.cos(beta) + y_s * np.sin(beta)))
                                                  - np.exp(-1j * k * (x_f * np.cos(beta) + y_p * np.sin(beta))))

    F_wp_amp = 1j * omega_e * zeta_a * (T_1 + T_2)

    return F_wp_amp


def wave_pumping_excitation_sesx(x_f, x_s, y_s, x_b, omega_0, U, heading, zeta_a=1, g=9.81):
    n = len(omega_0)  # length of list of omegas

    # initialize f_ex_7
    f_ex_7 = np.zeros([n], dtype=complex)

    # calculate one complex force amplitude per omega
    for i in range(n):
        f_ex_7[i] = wave_pumping_sesx(x_f, x_s, y_s, x_b, omega_0[i], U, heading, zeta_a, g)

    return f_ex_7
