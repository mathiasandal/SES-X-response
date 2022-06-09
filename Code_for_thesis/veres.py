import numpy as np
from pathlib import Path

'''Contains all functions related to things retrieved from VERES'''


def read_re7_file(filename):
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
    return VMAS, ADDMAS, DAMP, REST, VEL, HEAD, FREQ, XMTN, ZMTN, LCG, VCG, NDOF


def read_re8_file(filename):
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
                        REFORCE_D1[m, k, j, i], IMFORCE_D1[m, k, j, i], REFORCE_D2[m, k, j, i], IMFORCE_D2[m, k, j, i]=\
                            [float(m) for m in f.readline().split()][1:]

    return REFORCE, IMFORCE, VEL, HEAD, FREQ, XMTN, ZMTN, LCG, VCG


def read_group_of_re7_input(filepath):
    """Reads input of multiple *.re7 files for retrieving values from several runs in VERES"""

    if not Path(filepath + '/Run ' + str(1) + '/input.re7').is_file():
        return

    M, A_temp, B_temp, C_temp, VEL, HEAD, FREQ_temp, XMTN, ZMTN, lcg_veres, vcg_veres, NDOF = \
        read_re7_file(filepath + '/Run ' + str(1) + '/input.re7')

    # Clean up matrices
    A = A_temp[0, 0, :, :, :]
    B = B_temp[0, 0, :, :, :]
    C = C_temp[0, 0, :, :, :]
    omega_0 = FREQ_temp

    counter = 2

    while Path(filepath + '/Run ' + str(counter) + '/input.re7').is_file():

        M, A_temp, B_temp, C_temp, VEL, HEAD, FREQ_temp, XMTN, ZMTN, lcg_veres, vcg_veres, NDOF = \
            read_re7_file(filepath + '/Run ' + str(counter) + '/input.re7')

        # Append hydrodynamic coefficients and frequencies
        A = np.vstack((A_temp[0, 0, :, :, :].copy(), A))
        B = np.vstack((B_temp[0, 0, :, :, :].copy(), B))
        C = np.vstack((C_temp[0, 0, :, :, :].copy(), C))
        omega_0 = np.append(FREQ_temp.copy(), omega_0)

        counter += 1  # Increment for next iteration

    return M, A, B, C, VEL, HEAD, omega_0, XMTN, ZMTN, lcg_veres, vcg_veres, NDOF


def read_group_of_re8_input(filepath):
    """Reads input of multiple *.re7 files for retrieving values from several runs in VERES"""
    if not Path(filepath + '/Run ' + str(1) + '/input.re8').is_file():
        return

    REFORCE, IMFORCE, VEL, HEAD, FREQ_temp, XMTN, ZMTN, lcg_veres, vcg_veres = \
        read_re8_file(filepath + '/Run ' + str(1) + '/input.re8')

    f_ex = np.zeros([6, len(FREQ_temp)], dtype=complex)
    omega_0 = FREQ_temp

    for i in range(6):
        for j in range(len(FREQ_temp)):
            f_ex[i, j] = REFORCE[i, j, 0, 0] + 1j * IMFORCE[i, j, 0, 0]

    counter = 2

    while Path(filepath + '/Run ' + str(counter) + '/input.re8').is_file():
        REFORCE, IMFORCE, VEL, HEAD, FREQ_temp, XMTN, ZMTN, lcg_veres, vcg_veres = \
            read_re8_file(filepath + '/Run ' + str(counter) + '/input.re8')

        f_ex_temp = np.zeros([6, len(FREQ_temp)], dtype=complex)

        for i in range(6):
            for j in range(len(FREQ_temp)):
                f_ex_temp[i, j] = REFORCE[i, j, 0, 0] + 1j * IMFORCE[i, j, 0, 0]

        f_ex = np.hstack((f_ex_temp.copy(), f_ex))
        omega_0 = np.append(FREQ_temp.copy(), omega_0)

        counter += 1  # Increment for next iteration

    return f_ex, VEL, HEAD, omega_0, XMTN, ZMTN, lcg_veres, vcg_veres


def read_coefficients_from_veres(path):
    """Read hydrodynamic coefficients and loads from several rund in VERES"""

    # Read hydrodynamic coefficients
    M, A, B, C, U, beta, omega_0, XMTN, ZMTN, lcg_veres, vcg_veres, NDOF = read_group_of_re7_input(path)

    # Read hydrodynamic loads
    f_ex, U, beta, omega_0, XMTN, ZMTN, lcg_veres, vcg_veres = read_group_of_re8_input(path)

    return A, B, C[len(omega_0)//2, :, :], f_ex, omega_0, U, beta, XMTN, ZMTN, lcg_veres, vcg_veres
