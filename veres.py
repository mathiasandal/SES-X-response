import numpy as np
import scipy.linalg as la
import pandas as pd
from pathlib import Path
from air_cushion import find_closest_value
from Wave_response_utilities import add_row_and_column, decouple_matrix

'''Contains all functions related to things retrieved from VERES'''

def read_re1_file(filename):
    """
    Reads a *.re1 file created by VERES

    A description of a *.re1 file is found in page 86 of the ShipX Vessel Responses (VERES) User's manual.
    """
    f = open(filename, "r")  # open file for reading

    for i in range(6):  # skips first six lines
        f.readline()

    # Read in parameters
    RHOSW, GRAV = [float(i) for i in f.readline().split()]
    LPP, BREADTH, DRAUGHT = [float(i) for i in f.readline().split()]
    LCG, VCG = [float(i) for i in f.readline().split()]
    NOVEL, NOHEAD, NOFREQ, NDOF = [int(i) for i in f.readline().split()]

    # initialize more arrays
    VEL = np.zeros(NOVEL)
    HEAD = np.zeros(NOHEAD)
    FREQ = np.zeros(NOFREQ)
    SINK = np.zeros(NOVEL)
    TRIM = np.zeros(NOVEL)
    XMTN = np.zeros(NOVEL)
    ZMTN = np.zeros(NOVEL)

    # initialize lists to contain real and imaginary parts of RAO amplitudes
    RETRANS = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])
    IMTRANS = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])


    for ivel in range(NOVEL):
        VEL[ivel], SINK[ivel], TRIM[ivel], XMTN[ivel], ZMTN[ivel] = [float(i) for i in f.readline().split()]
        for ihead in range(NOHEAD):
            HEAD[ihead] = [float(i) for i in f.readline().split()][0]
            for ifreq in range(NOFREQ):
                FREQ[ifreq] = [float(i) for i in f.readline().split()][0]
                for dof in range(NDOF):
                    dummy, RETRANS[dof, ifreq, ihead, ivel], IMTRANS[dof, ifreq, ihead, ivel] = \
                        [float(i) for i in f.readline().split()]

    return RETRANS, IMTRANS, VEL, HEAD, FREQ, XMTN, ZMTN


def read_re5_file(filename):
    """
    Reads a *_ses.re5 file created by VERES

    A description of a *.re5 file is found in page 94 of the ShipX Vessel Responses (VERES) User's manual.
    """

    f = open(filename, "r")  # open file for reading

    for i in range(7):  # skips first seven lines
        f.readline()

    # Read in parameters
    RHOSW, GRAV, LSCALE = [float(i) for i in f.readline().split()]
    LPP, BREADTH, DRAUGHT = [float(i) for i in f.readline().split()]
    XCG, ZCG = [float(i) for i in f.readline().split()]
    NOVEL, NOHEAD, NOFREQ = [int(i) for i in f.readline().split()]
    NRESPS = [int(i) for i in f.readline().split()][0]

    # Initialize lists
    print(0)
    RESPID = [0] * NRESPS
    ISYM = [0] * NRESPS
    RESUNIT = [0] * NRESPS
    RESTXT1 = [0] * NRESPS
    RESTXT2 = [0] * NRESPS
    RESTXT3 = [0] * NRESPS
    RESTXT4 = [0] * NRESPS

    # initialize more arrays
    VEL = np.zeros(NOVEL)
    HEAD = np.zeros(NOHEAD)
    FREQ = np.zeros(NOFREQ)

    for ires in range(NRESPS):
        RESPID[ires], ISYM[ires], RESUNIT[ires], RESTXT1[ires], RESTXT2[ires], RESTXT3[ires], RESTXT4[ires] = \
            [str(i) for i in f.readline().split()]

    # Initialize arrays
    HYD = np.zeros([NOVEL, 6])
    TREHYD = np.zeros([NOVEL, 3])
    RREHYD = np.zeros([NOVEL, 3])

    RETRANS = np.zeros([NRESPS, NOFREQ, NOHEAD, NOVEL])
    IMTRANS = np.zeros([NRESPS, NOFREQ, NOHEAD, NOVEL])


    for ivel in range(NOVEL):
        HYD[ivel, :] = [float(i) for i in f.readline().split()]
        for ihead in range(NOHEAD):
            for ifreq in range(NOFREQ):
                VEL[ivel], HEAD[ihead], FREQ[ifreq] = [float(i) for i in f.readline().split()]
                for l in range(NRESPS):
                    RESPID[l], RETRANS[l, ifreq, ihead, ivel], IMTRANS[l, ifreq, ihead, ivel] = \
                        [float(i) for i in f.readline().split()]

    return FREQ, RETRANS[0, :, 0, 0], IMTRANS[0, :, 0, 0]


def read_re7_file(filename):
    # TODO: Fix/add documentation
    """
    Reads a *.re7 file created by VERES

    A description of a *.re7 file is found in page 114 of the ShipX Vessel Responses (VERES) User's manual.
    """

    f = open(filename, "r")  # open file for reading

    for i in range(6):  # skips first six lines
        f.readline()

    # Read in parameters
    RHOSW, GRAV = [float(i) for i in f.readline().split()]
    LPP, BREADTH, DRAUGHT = [float(i) for i in f.readline().split()]
    LCG, VCG = [float(i) for i in f.readline().split()]

    # saves line looking like "-1 FILEVER (=2)" in Veres_Manual.pdf
    FILEVER = [int(i) for i in f.readline().split()]

    # Read in number of velocities, headings, frequencies and degrees of freedom
    NOVEL, NOHEAD, NOFREQ, NDOF = [int(i) for i in f.readline().split()]

    # Initialize vectors to contain velocities, headings and frequencies
    VEL = np.zeros(NOVEL)  # [m/s] Velocity of the vessel
    HEAD = np.zeros(NOHEAD)  # [deg] Wave heading
    FREQ = np.zeros(NOFREQ)  # [rad/s] Wave frequency

    # Initialize values to be stores for each velocity
    SINK = np.zeros(NOVEL)  # [m] Sink
    TRIM = np.zeros(NOVEL)  # [deg] Trim
    XMTN = np.zeros(NOVEL)  # [m] X-pos. of the motion coordinate system (relative to Lpp/2)
    ZMTN = np.zeros(NOVEL)  # [m] Z-pos. of the mo􀆟on coordinate system (relative to BL)

    # Read in mass matrix
    VMAS = np.zeros([6, 6])  # Initialize mass matrix. Has the size (6x6) regardless if there is an air cushion or not.

    for j in range(6):  # Read and converts each line from string to floats
        VMAS[j, :] = [float(i) for i in f.readline().split()]

    # Have one (NDOF x NDOF) matrix for each combinations of velocity, heading and wave frequency
    # "Additional" matrices has the size (6x6) regardless if there is an air cushion or not.
    ADDMAS = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Hydrodynamic added mass
    ADDADDMAS = np.zeros([NOVEL, NOHEAD, NOFREQ, 6, 6])  # Additional added mass
    DAMP = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Hydrodynamic damping
    ADDDAMP = np.zeros([NOVEL, NOHEAD, NOFREQ, 6, 6])  # Additional damping
    REST = np.zeros([NOVEL, NOHEAD, NOFREQ, NDOF, NDOF])  # Restoring
    ADDREST = np.zeros([NOVEL, NOHEAD, NOFREQ, 6, 6])  # Additional restoring
    # "Additional" matrices has the size (6x6) regardless if there is an air cushion or not.

    # Initialize viscous roll damping
    VISCDL = np.zeros([NOVEL, NOHEAD, NOFREQ])  # Viscous roll damping, linear part
    VISCDN = np.zeros([NOVEL, NOHEAD, NOFREQ])  # Viscous roll damping, nonlinear part
    VISCDNL = np.zeros([NOVEL, NOHEAD, NOFREQ])  # Viscous roll damping, linearized

    for i in range(NOVEL):
        VEL[i], SINK[i], TRIM[i], XMTN[i], ZMTN[i] = [float(m) for m in f.readline().split()]
        for j in range(NOHEAD):
            HEAD[j] = float(f.readline())  # Should only contain one element
            for k in range(NOFREQ):
                FREQ[k] = float(f.readline())  # Should only contain one element

                for m in range(NDOF):  # NDOF = 7 if the air cushion is on and NDOF = 7 if the air cushion is off.
                    ADDMAS[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(6):
                    ADDADDMAS[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(NDOF):  # NDOF = 7 if the air cushion is on and NDOF = 7 if the air cushion is off.
                    DAMP[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(6):
                    ADDDAMP[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(NDOF):  # NDOF = 7 if the air cushion is on and NDOF = 7 if the air cushion is off.
                    REST[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                for m in range(6):
                    ADDREST[i, j, k, m, :] = [float(i) for i in f.readline().split()]

                # Read in viscous roll damping
                VISCDL[i, j, k], VISCDN[i, j, k], VISCDNL[i, j, k] = [float(i) for i in f.readline().split()]
    return VMAS, ADDMAS, DAMP, REST, VEL, HEAD, FREQ, XMTN, ZMTN, NDOF


# TODO: Create seperate function that works on 2D strip theory and change the name of this function to something with
#  high-speed formulation
def read_re8_file(filename):
    # TODO: Fix/add documentation
    """
    Reads a *.re8 file created by VERES

    A description of a *.re8 file is found in page 116 of the ShipX Vessel Responses (VERES) User's manual.

    Only works for High-Speed formulation in VERES
    """

    f = open(filename, 'r')

    for i in range(6):  # skips first six lines
        f.readline()

    # Read in parameters
    RHOSW, GRAV = [float(i) for i in f.readline().split()]
    LPP, BREADTH, DRAUGHT = [float(i) for i in f.readline().split()]
    LCG, VCG = [float(i) for i in f.readline().split()]

    # Read in number of velocities, headings, frequencies and degrees of freedom
    NOVEL, NOHEAD, NOFREQ, NDOF = [int(i) for i in f.readline().split()]

    # Initialize vectors to contain velocities, headings and frequencies
    VEL = np.zeros(NOVEL)  # [m/s] Velocity of the vessel
    HEAD = np.zeros(NOHEAD)  # [deg] Wave heading
    FREQ = np.zeros(NOFREQ)  # [rad/s] Wave frequency

    # Initialize values to be stores for each velocity
    SINK = np.zeros(NOVEL)  # [m] Sink
    TRIM = np.zeros(NOVEL)  # [deg] Trim
    XMTN = np.zeros(NOVEL)  # [m] X-pos. of the motion coordinate system (relative to Lpp/2)
    ZMTN = np.zeros(NOVEL)  # [m] Z-pos. of the mo􀆟on coordinate system (relative to BL)

    # Initialize Force components
    # For High-Speed formulation only the total excitation force is given
    REFORCE = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])  # Real part of excitation load
    IMFORCE = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])  # Imaginary part of excitation load

    # For conventional strip theory the Froude-Krylov and diffraction loads are also stored separately
    REFORCE_FK = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])  # Real part of Froude-Krylov part of excitation load
    IMFORCE_FK = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])  # Imaginary part of Froude-Krylov part of excitation load
    REFORCE_D1 = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])  # Real part of diffraction load
    IMFORCE_D1 = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])  # Imaginary part of diffraction load
    REFORCE_D2 = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])  # Real part of diffraction load without speed terms
    IMFORCE_D2 = np.zeros([NDOF, NOFREQ, NOHEAD, NOVEL])  # Imaginary part of diffraction load without speed terms

    for i in range(NOVEL):
        VEL[i], SINK[i], TRIM[i], XMTN[i], ZMTN[i] = [float(m) for m in f.readline().split()]
        for j in range(NOHEAD):
            HEAD[j] = float(f.readline())  # Should only contain one element
            for k in range(NOFREQ):
                FREQ[k] = float(f.readline())  # Should only contain one element

                temp = [float(m) for m in f.readline().split()][1:]

                if len(temp) == 2:
                    REFORCE[0, k, j, i], IMFORCE[0, k, j, i] = temp
                elif len(temp) == 8:
                    REFORCE[0, k, j, i], IMFORCE[0, k, j, i], REFORCE_FK[0, k, j, i], IMFORCE_FK[0, k, j, i], \
                    REFORCE_D1[0, k, j, i], IMFORCE_D1[0, k, j, i], REFORCE_D2[0, k, j, i], IMFORCE_D2[0, k, j, i] = \
                        temp

                for m in range(1, NDOF):
                    if len(temp) == 2:
                        REFORCE[m, k, j, i], IMFORCE[m, k, j, i] = [float(m) for m in f.readline().split()][1:]
                    elif len(temp) == 8:
                        REFORCE[m, k, j, i], IMFORCE[m, k, j, i], REFORCE_FK[m, k, j, i], IMFORCE_FK[m, k, j, i], \
                        REFORCE_D1[m, k, j, i], IMFORCE_D1[m, k, j, i], REFORCE_D2[m, k, j, i], IMFORCE_D2[m, k, j, i] = \
                            [float(m) for m in f.readline().split()][1:]

    return REFORCE, IMFORCE, VEL, HEAD, FREQ, XMTN, ZMTN


def iterate_natural_frequencies(wave_frequencies, velocity, heading, added_mass, mass, restoring, g=9.81, tolerance=1e-5):

    """
    Iterates to the correct natural undamped natural frequencies of a system with n degrees of freedom for a SES-X vessel.

    :param wave_frequencies:
    :param velocity:
    :param heading:
    :param added_mass: (NVEL x NHEAD x NFREQ x 6 x 6) Matrix
        Added mass of from the hull calculated using VERES.
    :param mass: (7x7) or (3x3) matrix
        if (7x7)
            contains mass properties of all 7dof, i.e. surge, sway, heave, roll, pitch, yaw and uniform cushion pressure
        if (3x3)
            contains mass properties of three degrees of freedom. Heave, pitch and uniform cushion pressure.
    :param restoring: (3x3) matrix
        if (3x3)
            contains restoring properties of three degrees of freedom. Heave, pitch and uniform cushion pressure.
    :param g: double
        accelerations of gravity
    :param tolerance:

    :return:
    """

    n = len(wave_frequencies)
    m = len(restoring)

    nat_frequencies = np.ones([m], dtype=complex)
    eigen_modes = np.zeros([m, m], dtype=complex)

    # calculates encounter frequency corresponding to each wave frequency and the vessel velocity and wave heading
    encounter_frequencies = wave_frequencies + velocity / g * np.cos(np.deg2rad(heading)) * np.power(wave_frequencies, 2)

    for i in range(m):
        counter = 0  # initialize
        omega_real = 0  # initialize
        index_frequency_upper = int(np.ceil(n/2))
        index_frequency_lower = index_frequency_upper - 1
        err = -1.

        # Get correct size of

        while (counter <= 100 and err >= tolerance) or err == -1:
            # Create matrices for the current encounter frequency
            '''
            M = decouple_matrix(mass, [2, 4, 6])
            A = decouple_matrix(add_row_and_column(added_mass[0, 0, index_frequency, :, :]), [2, 4, 6])
            C = decouple_matrix(restoring, [2, 4, 6])
            '''
            M = mass

            if m == 3:
                if err == -1.:
                    A = decouple_matrix(add_row_and_column(added_mass[0, 0, index_frequency_upper, :, :]), [2, 4, 6])
                else:
                    A_l = decouple_matrix(add_row_and_column(added_mass[0, 0, index_frequency_lower, :, :]), [2, 4, 6])
                    A_u = decouple_matrix(add_row_and_column(added_mass[0, 0, index_frequency_upper, :, :]), [2, 4, 6])
                    A = interpolate_matrices(omega_real, encounter_frequencies[index_frequency_lower], encounter_frequencies[index_frequency_upper], A_l, A_u)

            else:
                if err == -1.:
                    A = add_row_and_column(added_mass[0, 0, index_frequency_upper, :, :])
                else:
                    A_l = add_row_and_column(added_mass[0, 0, index_frequency_lower, :, :])
                    A_u = add_row_and_column(added_mass[0, 0, index_frequency_upper, :, :])
                    A = interpolate_matrices(omega_real, encounter_frequencies[index_frequency_lower], encounter_frequencies[index_frequency_upper], A_l, A_u)

            C = restoring

            nat_freq_temp, eigen_modes_temp = la.eig(C, M + A)

            '''
            if nat_freq_temp[i].real <= 0:
                raise ValueError
            '''

            omega_real = np.sqrt(nat_freq_temp[i].real**2 + nat_freq_temp[i].imag**2)

            # Finds next guess
            dummy_omega, index_frequency = find_closest_value(encounter_frequencies, omega_real)

            if dummy_omega < omega_real:
                index_frequency_lower = index_frequency
                index_frequency_upper = index_frequency + 1
            elif dummy_omega > omega_real:
                index_frequency_lower = index_frequency - 1
                index_frequency_upper = index_frequency

            # Computes relative error
            err = abs((np.sqrt(nat_frequencies[i].real**2 + nat_frequencies[i].imag**2) - omega_real) /
                      np.sqrt(nat_frequencies[i].real**2 + nat_frequencies[i].imag**2))

            # err = abs((encounter_frequencies[index_frequency] - omega_real)/encounter_frequencies[index_frequency])

            # Assigns new natural frequency to array
            nat_frequencies[i] = nat_freq_temp[i]
            eigen_modes[i, :] = eigen_modes_temp[i, :]

            if omega_real == float('inf'):  # infinite frequency
                break  # Breaks if the natural frequency has gone to infinity
            elif index_frequency_lower < 0:
                break  # Breaks if natural frequency goes to zero

            counter += 1  # increment counter

    return nat_frequencies, eigen_modes, encounter_frequencies


def read_veres_input(path):
    # TODO: Add documentation

    # Reads *.re7 file
    VMAS, ADDMAS, DAMP, REST, VEL_re7, HEAD_re7, FREQ_re7, XMTN_re7, ZMTN_re7 = read_re7_file(path + '//input.re7')

    # Reads *.re8 file
    REFORCE, IMFORCE, VEL_re8, HEAD_re8, FREQ_re8, XMTN_re8, ZMTN_re8 = read_re8_file(path + '//input.re8')

    # Store the necessary data
    A_h = ADDMAS
    B_h = DAMP
    C_h = REST
    F_ex_real = REFORCE
    F_ex_im = IMFORCE

    # Input handling. Checks that the hydrodynamic coefficients from *.re7 have the same velocity, heading and
    # frequencies as the excitation forces in *.re8 file.
    if VEL_re7.all() == VEL_re8.all() and HEAD_re7.all() == HEAD_re8.all() and FREQ_re7.all() == FREQ_re8.all() and \
            XMTN_re7.all() == XMTN_re8.all() and ZMTN_re7.all() == ZMTN_re8.all():
        VEL = VEL_re7
        HEAD = HEAD_re7
        FREQ = FREQ_re7
        XMTN = XMTN_re7
        ZMTN = ZMTN_re7
    else:
        raise ValueError

    return A_h, B_h, C_h, F_ex_real, F_ex_im, VEL, HEAD, FREQ, XMTN, ZMTN


def read_group_of_re7_input(filepath):

    if not Path(filepath + '/Run ' + str(1) + '/input.re7').is_file():
        return

    M, A_temp, B_temp, C_temp, VEL, HEAD, FREQ_temp, XMTN, ZMTN, NDOF = \
        read_re7_file(filepath + '/Run ' + str(1) + '/input.re7')

    # Clean up matrices
    A = A_temp[0, 0, :, :, :]
    B = B_temp[0, 0, :, :, :]
    C = C_temp[0, 0, :, :, :]
    omega_0 = FREQ_temp

    counter = 2

    while Path(filepath + '/Run ' + str(counter) + '/input.re7').is_file():

        M, A_temp, B_temp, C_temp, VEL, HEAD, FREQ_temp, XMTN, ZMTN, NDOF = \
            read_re7_file(filepath + '/Run ' + str(counter) + '/input.re7')

        # Append hydrodynamic coefficients and frequencies
        A = np.vstack((A_temp[0, 0, :, :, :].copy(), A))
        B = np.vstack((B_temp[0, 0, :, :, :].copy(), B))
        C = np.vstack((C_temp[0, 0, :, :, :].copy(), C))
        omega_0 = np.append(FREQ_temp.copy(), omega_0)

        counter += 1  # Increment for next iteration

    return M, A, B, C, VEL, HEAD, omega_0, XMTN, ZMTN, NDOF


def read_group_of_re8_input(filepath):
    if not Path(filepath + '/Run ' + str(1) + '/input.re8').is_file():
        return

    REFORCE, IMFORCE, VEL, HEAD, FREQ_temp, XMTN, ZMTN = \
        read_re8_file(filepath + '/Run ' + str(1) + '/input.re8')

    f_ex = np.zeros([6, len(FREQ_temp)], dtype=complex)
    omega_0 = FREQ_temp

    for i in range(6):
        for j in range(len(FREQ_temp)):
            f_ex[i, j] = REFORCE[i, j, 0, 0] + 1j * IMFORCE[i, j, 0, 0]

    counter = 2

    while Path(filepath + '/Run ' + str(counter) + '/input.re8').is_file():
        REFORCE, IMFORCE, VEL, HEAD, FREQ_temp, XMTN, ZMTN = \
            read_re8_file(filepath + '/Run ' + str(counter) + '/input.re8')

        f_ex_temp = np.zeros([6, len(FREQ_temp)], dtype=complex)

        for i in range(6):
            for j in range(len(FREQ_temp)):
                f_ex_temp[i, j] = REFORCE[i, j, 0, 0] + 1j * IMFORCE[i, j, 0, 0]

        f_ex = np.hstack((f_ex_temp.copy(), f_ex))
        omega_0 = np.append(FREQ_temp.copy(), omega_0)

        counter += 1  # Increment for next iteration

    return f_ex, VEL, HEAD, omega_0, XMTN, ZMTN

def interpolate_matrices(omega, omega_lower, omega_upper, mat_lower, mat_upper):
    """
    Interpolates two (nxn) numpy arrays at two different frequencies.

    :param omega: double
        Frequency of the evaluated frequency
    :param omega_lower: double
        Frequency at the lower end for the interpolation
    :param omega_upper: double
        Frequency at the upper end for the interpolation
    :param mat_lower: (nxn) matrix
        Matrix evaluated at the lower end frequency
    :param mat_upper: (nxn) matrix
        Matrix evaluated at the upper end frequency
    :return:
        mat: (nxn) matrix
            interpolated matrix evaluated at omega.
    """

    # input handling
    if len(mat_upper) != len(mat_lower):
        raise ValueError
    elif omega_lower > omega_upper:
        raise ValueError
    elif omega < omega_lower or omega > omega_upper:
        raise ValueError

    # interpolate
    mat = mat_lower + (omega - omega_lower) / (omega_upper - omega_lower) * (mat_upper - mat_lower)

    return mat


def compute_RAOs(velocity, heading, wave_frequencies, M, A, B_c, B_h, C, F_ex_real, F_ex_imag, f_ex_7,
                 force_combination=1, g=9.81):
    # TODO: Add documentation

    n = len(wave_frequencies)

    encounter_frequencies = wave_frequencies + velocity / g * np.cos(np.deg2rad(heading)) * np.power(wave_frequencies,
                                                                                                     2)

    rao = np.zeros([7, n], dtype=complex)

    # force_combination = 2  # 1: Hydro, 2: Wave pumping, 3: both

    for i in range(n):
        A_temp = add_row_and_column(A[0, 0, i, :, :])
        B_temp = add_row_and_column(B_h[0, 0, i, :, :]) + B_c
        # adds a empty column to model no excitation in the air cushion

        if force_combination == 1:
            F_ex_real_temp = np.r_[F_ex_real[:, i, 0, 0], np.array([0])]
            F_ex_imag_temp = np.r_[F_ex_imag[:, i, 0, 0], np.array([0])]
        elif force_combination == 2:
            F_ex_real_temp = np.zeros([7])
            F_ex_imag_temp = np.zeros([7])
            F_ex_real_temp[6] = f_ex_7[i].real
            F_ex_imag_temp[6] = f_ex_7[i].imag
        elif force_combination == 3:
            F_ex_real_temp = np.r_[F_ex_real[:, i, 0, 0], np.array([f_ex_7[i].real])]
            F_ex_imag_temp = np.r_[F_ex_imag[:, i, 0, 0], np.array([f_ex_7[i].imag])]

        # Solves linear system of equations
        rao[:, i] = np.linalg.solve(
            -encounter_frequencies[i] ** 2 * (M + A_temp) + 1j * encounter_frequencies[i] * B_temp + C,
            F_ex_real_temp + F_ex_imag_temp * 1j)

    return encounter_frequencies, rao


def print_natfrequencies_and_eigenmodes(nat_frequencies, eigen_modes, type="radians_per_sec"):
    decimal_precision = 2
    keep_complex = True

    if type == 'periods':
        for i in range(len(nat_frequencies)):
            if nat_frequencies[i] == 0 or nat_frequencies[i] < 1e-6:
                nat_frequencies[i] = float('inf')
            elif nat_frequencies[i] == float('inf'):
                nat_frequencies[i] = 0
            else:
                nat_frequencies[i] = 2 * np.pi / nat_frequencies[i]
    elif type == 'Hz':
        for i in range(len(nat_frequencies)):
            nat_frequencies[i] = nat_frequencies[i] / 2 / np.pi

    # create list containing eigenmodes DOFs
    vertical_list = []
    horizontal_list = []
    for i in range(len(nat_frequencies)):
        vertical_list.append('eta_' + str(i + 1))

        if type == 'periods':
            horizontal_list.append(str(round(nat_frequencies[i], decimal_precision)) + ' [s]')
        elif type == 'Hz':
            horizontal_list.append(str(round(nat_frequencies[i], decimal_precision)) + ' [Hz]')
        else:
            horizontal_list.append(str(round(nat_frequencies[i], decimal_precision)) + ' [rad/s]')

    df = pd.DataFrame(np.around((eigen_modes), decimal_precision), vertical_list, horizontal_list)

    print(df)

    return df


def uncoupled_natural_frequency(encounter_frequency, M, C, A, tol=10 - 8):
    n = len(encounter_frequency)

    omega_temp = encounter_frequency[n // 2]

    A_temp = np.interp(omega_temp, encounter_frequency, A)

    nat_frequency = np.sqrt(C / (M + A_temp))

    err = (nat_frequency - omega_temp) / omega_temp  # computes relative error

    counter = 1  # initialize counter

    while err < tol and counter < 100:
        omega_temp = nat_frequency
        A_temp = np.interp(omega_temp, encounter_frequency, A)

        nat_frequency = np.sqrt(C / (M + A_temp))

        err = (nat_frequency - omega_temp) / omega_temp  # computes relative error

        counter += 1

    return nat_frequency


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    path_re5 = "Input files/Conceptual SES/with air cushion/input_ses.re5"
    path_re7 = "Input files/Conceptual SES/with air cushion/input.re7"
    path_re1 = "Input files/Conceptual SES/with air cushion/input.re1"

    # Define some parameters
    p_0 = 3500  # [Pa]
    rho = 1025  # [kg/m^3]
    g = 9.81  # [m/s^2]
    L = 20  # [m]
    zeta_a = 1  # [m]

    frequency_nond, real_trans, imag_trans = read_re5_file(path_re5)

    VMAS, ADDMAS, DAMP, REST, VEL, HEAD, FREQ, XMTN, ZMTN, NDOF = read_re7_file(path_re7)

    trans_nond = np.sqrt(np.power(real_trans, 2) + np.power(imag_trans, 2))
    trans = trans_nond * rho * g * zeta_a / p_0
    frequency = frequency_nond * np.sqrt(g / L)

    plt.plot(frequency, trans, '-x')
    plt.xlabel('encounter frequency [rad/s]')
    plt.ylabel('$\\eta_{7}$')
    plt.title('RAO in uniform pressure')
    plt.show()

    # Read RAOs
    RETRANS, IMTRANS, VEL_1, HEAD_1, FREQ_1, XMTN_1, ZMTN_1 = read_re1_file(path_re1)

    rao_re = RETRANS[:, :, 0, 0]
    rao_im = IMTRANS[:, :, 0, 0]
    rao = rao_re + 1j * rao_im

    print('hi')
