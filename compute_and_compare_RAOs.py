from veres import read_re1_file, read_re8_file, read_re7_file, read_re5_file, iterate_natural_frequencies, print_natfrequencies_and_eigenmodes
from air_cushion import wave_pumping_rect
from Wave_response_utilities import solve_eq_motion_steady_state, add_row_and_column
import numpy as np
import pandas as pd
import scipy.linalg as la
import matplotlib.pyplot as plt


# Fetch hydrodynamic coefficients and loads
path = 'Input files/Conceptual SES/with air cushion/'
path_re7 = path + 'input.re7'
path_re8 = path + 'input.re8'
path_re1 = path + 'input.re1'

# Read input.re7 file to get hydrodynamic coefficients
M, A_n, B_n, C_n, VEL_re7, HEAD_re7, FREQ_re7, XMTN_re7, ZMTN_re7, NDOF = read_re7_file(path_re7)

# Clean up matrices
M = add_row_and_column(M)
A = A_n[0, 0, :, :, :]
B = B_n[0, 0, :, :, :]
C = C_n[0, 0, :, :, :]
omegas = FREQ_re7

# Read input.re8 file to get excitation
REFORCE, IMFORCE, VEL_re8, HEAD_re8, FREQ_re8, XMTN_re8, ZMTN_re8 = read_re8_file(path_re8)

n_frequencies = len(FREQ_re7)

# initialize array to contain complex force amplitudes
raos = np.zeros([7, n_frequencies], dtype=complex)

# initialize array to contain complex force amplitudes
f_ex = np.zeros([7, n_frequencies], dtype=complex)

# Properties
L = 20  # [m] length of the vessel
b = 7  # [m] beam of air cushion

# Compute wave pumping excitation
x_s = L/2
x_f = -L/2
y_s = b/2
y_p = -b/2

f_ex_7 = wave_pumping_rect(x_f, x_s, y_p, y_s, omegas, HEAD_re7[0])

for i in range(6):
    for j in range(n_frequencies):
        f_ex[i, j] = REFORCE[i, j, 0, 0] + 1j * IMFORCE[i, j, 0, 0]

f_ex[6, :] = 1j * f_ex_7  # add the wave pumping to the force vectors

# Compute and store complex response amplitudes

for i in range(n_frequencies):
    raos[:, i] = la.solve(-omegas[i]**2 * (M + A[i, :, :]) + 1j * omegas[i] * B[i, :, :] + C[i, :, :], f_ex[:, i])
    # raos[:, i] = solve_eq_motion_steady_state(M + A[i, :, :], B[i, :, :], C[i, :, :], f_ex[:, i], omegas[i])

raos_magnitude = np.transpose(np.abs(raos))  # store the magnitude of the complex response amplitude
raos_phase = np.transpose(np.rad2deg(np.angle(raos)))  # store the phase of the complex response amplitude
raos_python = np.transpose(raos)

# Organize RAO data in a similar manner as in the Excel Veres postprocessor
dat = np.array([raos_magnitude[:, 0], raos_phase[:, 0], raos_magnitude[:, 1], raos_phase[:, 1],
                raos_magnitude[:, 2], raos_phase[:, 2], raos_magnitude[:, 3], raos_phase[:, 3],
                raos_magnitude[:, 4], raos_phase[:, 4], raos_magnitude[:, 5], raos_phase[:, 5],
                raos_magnitude[:, 6], raos_phase[:, 6]]).transpose()

df = pd.DataFrame(dat, index=omegas, columns=['eta_1 mag', 'eta_1 phase', 'eta_2 mag', 'eta_2 phase',
                                              'eta_3 mag', 'eta_3 phase', 'eta_4 mag', 'eta_4 phase',
                                              'eta_5 mag', 'eta_5 phase', 'eta_6 mag', 'eta_6 phase',
                                              'eta_7 mag', 'eta_7 phase'])

# Compute natural frequencies and eigenmodes
nat_frequencies, eigen_modes = la.eig(C[40, :, :], M + A[40, :, :])

df_nat = print_natfrequencies_and_eigenmodes(np.power(nat_frequencies, 0.5), eigen_modes)

k = np.power(omegas, 2) / 9.81
wavelength = np.divide(2*np.pi, k)

# Read RAO computed in VERES


# Compute and compare heave accelerations
# Read RAOs
RETRANS, IMTRANS, VEL_1, HEAD_1, FREQ_1, XMTN_1, ZMTN_1 = read_re1_file(path_re1)

rao_veres_re = RETRANS[:, :, 0, 0]
rao_veres_im = IMTRANS[:, :, 0, 0]
rao_veres = np.transpose(rao_veres_re + 1j * rao_veres_im)

# Compute RAO for uniform pressure calculated in Veres
# Plot and compare uniform pressure RAO
path_re5 = "Input files/Conceptual SES/with air cushion/input_ses.re5"
path_re7 = "Input files/Conceptual SES/with air cushion/input.re7"

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

a = np.array([real_trans + 1j * imag_trans]).transpose()
rao_veres = np.hstack((rao_veres, a * rho * g * zeta_a / p_0))


# Plot RAOs

plot_raos = True

if plot_raos:
    plt.plot(omegas, raos_magnitude[:, 2], '-x', label='python')
    plt.plot(omegas, np.abs(rao_veres[:, 2]), '-x', label='veres')
    plt.xlabel('encounter frequency [rad/s]')
    plt.ylabel('$\\eta_{3}/\\zeta_a$')
    plt.title('RAO in heave')
    plt.show()

    plt.plot(omegas, np.multiply(np.power(omegas, 2), raos_magnitude[:, 2]), '-x', label='python')
    plt.plot(omegas, np.abs(np.multiply(np.power(omegas, 2), rao_veres[:, 2])), '-x', label='veres')
    plt.xlabel('encounter frequency [rad/s]')
    plt.ylabel('$|\\omega^2\\eta_{3}|/\\zeta_a[s^{-2}]$')
    plt.title('RAO for acceleration in heave')
    plt.show()

    plot_type = 2
    if plot_type == 1:
        plt.plot(frequency, trans, '-rx', label='Veres')
        plt.plot(omegas, raos_magnitude[:, 6], fillstyle='none', marker='o', linestyle='-', label='Python')
        plt.xlabel('encounter frequency [rad/s]')
    elif plot_type == 2:
        plt.plot(wavelength, trans, '-rx', label='Veres')
        plt.plot(wavelength, raos_magnitude[:, 6], fillstyle='none', marker='o', linestyle='-', label='Python')
        plt.xlabel('$\\lambda[m]$')

    plt.ylabel('$\\eta_{7}$')
    plt.title('RAO in uniform pressure')
    plt.legend()
    plt.show()

    plot_dof = 2

    # plot RAO
    if plot_type == 1:
        plt.plot(omegas, raos_magnitude[:, plot_dof], '-x', label='python')
        plt.plot(omegas, np.abs(rao_veres[:, plot_dof]), fillstyle='none', marker='o', linestyle='-', label='veres')
        plt.xlabel('encounter frequency [rad/s]')
    elif plot_type == 2:
        plt.plot(wavelength, raos_magnitude[:, plot_dof], '-x', label='python')
        plt.plot(wavelength, np.abs(rao_veres[:, plot_dof]), fillstyle='none', marker='o', linestyle='-', label='veres')
        plt.xlabel('$\\lambda[m]$')

    plt.ylabel('$\\eta_{' + str(plot_dof + 1) + '}/\\zeta_a$')
    plt.legend()
    plt.show()

    # plot phase
    if plot_type == 1:
        plt.plot(omegas, np.angle(raos_python[:, plot_dof]) * 180, '-x', label='python')
        plt.plot(omegas, np.angle(rao_veres[:, plot_dof]) * 180, fillstyle='none', marker='o', linestyle='-', label='veres')
        plt.xlabel('encounter frequency [rad/s]')
    elif plot_type == 2:
        plt.plot(wavelength, np.rad2deg(np.angle(raos_python[:, plot_dof])), '-x', label='python')
        plt.plot(wavelength, np.rad2deg(np.angle(rao_veres[:, plot_dof])), fillstyle='none', marker='o', linestyle='-', label='veres')
        plt.xlabel('$\\lambda[m]$')

    plt.ylabel('Phase shift, ' + '$\\eta_{' + str(plot_dof + 1) + '}/\\zeta_a$')
    plt.legend()
    plt.show()


comparison_magnitude = np.divide(np.abs(np.abs(rao_veres) - np.abs(raos_python)), np.abs(rao_veres))
comparison_phase = np.divide(np.abs(np.angle(rao_veres) - np.angle(raos_python)), np.abs(np.angle(rao_veres)))


print('1')
