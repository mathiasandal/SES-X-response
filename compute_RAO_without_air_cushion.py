from veres import read_re8_file, read_re7_file, iterate_natural_frequencies, print_natfrequencies_and_eigenmodes
from air_cushion import wave_pumping_rect
from Wave_response_utilities import solve_eq_motion_steady_state, add_row_and_column
import numpy as np
import pandas as pd
import scipy.linalg as la
import matplotlib.pyplot as plt


# Fetch hydrodynamic coefficients and loads
path = 'Input files/Conceptual SES/without air cushion/'
path_re7 = path + 'input.re7'
path_re8 = path + 'input.re8'


# Read input.re7 file to get hydrodynamic coefficients
M, A_n, B_n, C_n, VEL_re7, HEAD_re7, FREQ_re7, XMTN_re7, ZMTN_re7, NDOF = read_re7_file(path_re7)

# Clean up matrices
A = A_n[0, 0, :, :, :]
B = B_n[0, 0, :, :, :]
C = C_n[0, 0, :, :, :]
omegas = FREQ_re7

# Read input.re8 file to get excitation
REFORCE, IMFORCE, VEL_re8, HEAD_re8, FREQ_re8, XMTN_re8, ZMTN_re8 = read_re8_file(path_re8)

n_frequencies = len(FREQ_re7)

# initialize array to contain complex force amplitudes
raos = np.zeros([6, n_frequencies], dtype=complex)

# initialize array to contain complex force amplitudes
f_ex = np.zeros([6, n_frequencies], dtype=complex)

for i in range(6):
    for j in range(n_frequencies):
        f_ex[i, j] = REFORCE[i, j, 0, 0] + 1j * IMFORCE[i, j, 0, 0]


# Compute and store complex response amplitudes

for i in range(n_frequencies):
    raos[:, i] = la.solve(-omegas[i]**2 * (M + A[i, :, :]) + 1j * omegas[i] * B[i, :, :] + C[i, :, :], f_ex[:, i])

raos_magnitude = np.transpose(np.abs(raos))  # store the magnitude of the complex response amplitude
raos_phase = np.transpose(np.rad2deg(np.angle(raos)))  # store the phase of the complex response amplitude

k = np.power(omegas, 2)/9.81

# Organize RAO data in a similar manner as in the Excel Veres postprocessor
dat = np.array([raos_magnitude[:, 0], raos_phase[:, 0], raos_magnitude[:, 1], raos_phase[:, 1],
                raos_magnitude[:, 2], raos_phase[:, 2], np.divide(raos_magnitude[:, 3], k), raos_phase[:, 3],
                np.divide(raos_magnitude[:, 4], k), raos_phase[:, 4], np.divide(raos_magnitude[:, 5], k), raos_phase[:, 5]]).transpose()

df = pd.DataFrame(dat, index=omegas, columns=['eta_1 mag', 'eta_1 phase', 'eta_2 mag', 'eta_2 phase',
                                              'eta_3 mag', 'eta_3 phase', 'eta_4 mag', 'eta_4 phase',
                                              'eta_5 mag', 'eta_5 phase', 'eta_6 mag', 'eta_6 phase'])

# Plot RAOs

plot_raos = True

if plot_raos:
    plt.plot(omegas, raos_magnitude[:, 2], '-x')
    plt.xlabel('encounter frequency [rad/s]')
    plt.ylabel('$\\eta_{3}/\\zeta_a$')
    plt.title('RAO in heave')
    plt.show()

    plt.plot(omegas, np.multiply(np.power(omegas, 2), raos_magnitude[:, 2]), '-x')
    plt.xlabel('encounter frequency [rad/s]')
    plt.ylabel('$|\\omega^2\\eta_{3}|/\\zeta_a[s^{-2}]$')
    plt.title('RAO for acceleration in heave')
    plt.show()
