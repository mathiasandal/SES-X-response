from veres import read_re7_file, read_re8_file
from air_cushion import wave_pumping_excitation_sesx, air_cushion_area, read_fan_characteristics, \
    interpolate_fan_characteristics, wave_pumping_rect
from mass_matrix import create_mass_matrix
from Wave_response_utilities import add_row_and_column
import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

# ***** Physical parameters *****
rho = 1025  # [kg/m^3] density of sea water
p_a = 101325  # [Pa] atmospheric pressure
g = 9.81  # [m/s^2] acceleration of gravity
gamma = 1.4  # [-] Ratio of specific heat for air, when assuming adiabatic density relation

# ***** Input parameters for BBGreen *****

# main dimensions of BBGreen
B = 6  # [m] beam of BBGreen
Lpp = 19.2  # [m] L_pp of BBGreen
total_mass = 25.6e3  # [kg] total mass of the vessel
r44 = 0.35 * B  # [m] radii of gyration in roll
r55 = 0.25 * Lpp  # [m] radii of gyration in pitch
r66 = 0.27 * Lpp  # [m] radii of gyration in yaw
r46 = 0  # [m] radii of gyration for inertial coupling of yaw and roll
lcg = 7.83  # [m] longitudinal center of gravity relative to AP
vcg = 1.98  # [m] vertical center of gravity relative to BL

# Properties of SES-X air cushion
l_1 = 12   # 0.0001  #     [m] length of the rectangular part of the air cushion
l_2 = 6  # 0.001  #      [m] length of the triangular part of the air cushion
b = 3.4  # [m] beam of the air cushion
p_0 = 3500  # [Pa] mean excess pressure in the air cushion
h_b = 0.64  #2 # [m]  Height of the air cushion  # TODO: make sure this value is correct
h = p_0 / rho / g  # [m] Difference in water level inside and outside the air cushion

# Compute area and longitudinal position of centroid relative to AP
A_b, x_B_c = air_cushion_area(l_1, l_2, b)

V_c0 = A_b * h_b  # [m^3] Compute mean air cushion volume  # TODO: Make sure this is correct

# Read in fan characteristics
Q, P, rpm_dummy = read_fan_characteristics('Input files/fan characteristics/fan characteristics.csv', '1800rpm')

# Interpolates values
Q_0, dQdp_0 = interpolate_fan_characteristics(p_0, P, Q)

# Fetch hydrodynamic coefficients and loads
path = 'Input files/BBGreen/22kn,0.17draft,0.33trim,z-coord at waterline/'
path_re7 = path + 'input.re7'
path_re8 = path + 'input.re8'

VMAS_veres, A_h, B_h, C_h, VEL_re7, HEAD_re7, FREQ_re7, XMTN_re7, ZMTN_re_7, NDOF = read_re7_file(path_re7)

REFORCE, IMFORCE, VEL_re8, HEAD_re8, FREQ_re8, XMTN_re8, ZMTN_re8 = read_re8_file(path_re8)

# Center of motion relative to BL and AP
x_prime = Lpp/2 - XMTN_re7[0]
z_prime = ZMTN_re_7[0]


# Properties of SES-X shaped air cushion with origin in motion coordinate system
cushion_width = b  # [m]
x_s = x_prime  # [m]
x_f = -(l_1 - x_prime)  # [m]
y_s = b/2  # [m]
y_p = -b/2  # [m]
x_b = -(l_1 + l_2 - x_prime)  # [m]

beta = HEAD_re7[0]  # [deg] wave heading
U = VEL_re7[0]  # [m/s] forward speed
omega_0 = FREQ_re7  # [rad/s] wave frequency

omega_e = omega_0 + U / g * np.multiply(np.power(omega_0, 2), np.cos(np.rad2deg(beta)))   # [rad/s] encounter frequency
n_frequencies = len(FREQ_re7)  # number of wave frequencies

# Initialize numpy array to contain complex excitation amplitudes
f_ex_1 = np.zeros([6, n_frequencies], dtype=complex)

# fill the array with Froude Krylov and diffraction loads
for i in range(n_frequencies):
    for j in range(6):
        f_ex_1[j, i] = REFORCE[j, i, 0, 0] + 1j * IMFORCE[j, i, 0, 0]

# compute wave pumping
f_7_test = 1j*wave_pumping_rect(x_b, x_s, y_p, y_s, omega_e, beta)
#f_7_test = np.ones([n_frequencies])*30
f_7_wp = wave_pumping_excitation_sesx(x_f, x_s, y_s, x_b, omega_e, beta)

f_7 = f_7_wp

# plot wavepumping
encounter_frequency_Hz = omega_e / 2 / np.pi  # compute encounter frequency in Hz
k_encounter = np.power(omega_e, 2) / g  # encounter wave number
encounter_wavelength = np.divide(2*np.pi, k_encounter)


fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(encounter_frequency_Hz, np.abs(f_7_wp), label='SES-X shape')
ax1.plot(encounter_frequency_Hz, np.abs(f_7_test), label='Test')
fig.suptitle(
    'Wave pumping excitation for SES-X air cushion. $l_1 = %.0f$m, ' % l_1 + '$l_2 = %.0f$m, ' % l_2 + '$b = %0.1f$m, ' % b + '$\\beta = %0.1f^{\circ}$' % beta)
ax1.set_xlabel("encounter frequency [Hz]")
ax1.set_ylabel('$F_{wp} [m^3/s]$')
ax1.set_xlim([0, 12])
ax1.legend()
ax2.plot(encounter_frequency_Hz, np.angle(f_7_wp), label='SES-X shape')
ax2.plot(encounter_frequency_Hz, np.angle(f_7_test), label='Test')
ax2.set_xlabel("encounter frequency [Hz]")
ax2.set_ylabel('Phase shift')
ax2.set_xlim([0, 12])
ax2.legend()
plt.show()


f_ex = np.vstack([f_ex_1, f_7])  # append the wave pumping as a seventh excitation
f_ex_test = np.vstack([f_ex_1, f_7_test])
# TODO: Notice signs of these values. Should be fixed in content from project thesis
x_G = x_prime - lcg  # [m] longitudinal position of COG relative to the motion coordinate system
z_G = -(z_prime - vcg)  # [m] vertical position of COG relative to the motion coordinate system

# Creates mass matrix
M = create_mass_matrix(total_mass, r44, r55, r66, r46, x_G, z_G)
M = add_row_and_column(M)  # Add empty row and column for the seventh DOF

# Centroid of the air cushion at equilibrium relative to the motion coord. system
x_c = x_prime - x_B_c  # [m] longitudinal position
y_c = 0  # [m] transverse position
z_c = -h / 2  # [m] vertical position

# Compute restoring terms due to the presence of the air cushion
C_44_c = rho * g * h * A_b * z_c  # [Nm] restoring coefficient in roll due to the air cushion
C_55_c = rho * g * h * A_b * z_c  # [Nm] restoring coefficient in pitch due to the air cushion
C_57_c = -rho * g * h * A_b * x_c  # [Nm] coupling term in pitch due to change in the air cushion pressure
C_37_c = -rho * g * h * A_b  # [N] coupling term in heave due to change in the air cushion pressure
C_77_c = 0.5 * Q_0 - p_0 * dQdp_0  # [m^3/s] equivalent restoring term in mass continuity eq. for air inside air cushion

# Compute damping terms du to the presence of the air cushion
B_73_c = A_b  # [m] equivalent damping coefficient in uniform pressure DOF due to velocity in heave.
B_75_c = -A_b * x_c  # [m^2] equivalent damping coefficient in uniform pressure DOF due to velocity in pitch.
B_77_c = p_0 * V_c0 / gamma / (p_0 + p_a)  # [m^3] equivalent damping coefficient in uniform pressure DOF.

# Remove dimensions in arrays connected to velocity and heading
A_h = A_h[0, 0, :, :, :]
B_h = B_h[0, 0, :, :, :]
C_h = C_h[0, 0, :, :, :]

# Initialize added mass-, damping- and restoring matrices for 7DOF
A = np.zeros([n_frequencies, 7, 7])
B = np.zeros([n_frequencies, 7, 7])
C = np.zeros([n_frequencies, 7, 7])

# Add empty row and columns in order to add the seventh DOF
for i in range(n_frequencies):
    A[i] = add_row_and_column(A_h[i])
    B[i] = add_row_and_column(B_h[i])
    C[i] = add_row_and_column(C_h[i])

# Add cushion terms to restoring- and damping matrices
# Restoring
C[:, 3, 3] += C_44_c
C[:, 4, 4] += C_55_c
C[:, 4, 6] += C_57_c
C[:, 2, 6] += C_37_c
C[:, 6, 6] += C_77_c
# Damping
B[:, 6, 2] += B_73_c
B[:, 6, 4] += B_75_c
B[:, 6, 6] += B_77_c

# Correct stiffness matrix for change of center of gravity
vcg_veres = 0.5  # [m] vertical center of gravity used in Veres  # TODO: Get this automatically
C[:, 3, 3] -= total_mass * g * (vcg_veres - vcg)
C[:, 4, 4] -= total_mass * g * (vcg_veres - vcg)

# Compute RAOs

# initialize array to contain complex force amplitudes
raos = np.zeros([7, n_frequencies], dtype=complex)
raos_test = np.zeros([7, n_frequencies], dtype=complex)

# Compute and store complex response amplitudes

for i in range(n_frequencies):
    raos[:, i] = la.solve(-omega_e[i]**2 * (M + A[i]) + 1j * omega_e[i] * B[i] + C[i], f_ex[:, i])
    raos_test[:, i] = la.solve(-omega_e[i] ** 2 * (M + A[i]) + 1j * omega_e[i] * B[i] + C[i], f_ex_test[:, i])

raos_magnitude = np.transpose(np.abs(raos))  # store the magnitude of the complex response amplitude
raos_phase = np.transpose(np.rad2deg(np.angle(raos)))  # store the phase of the complex response amplitude
raos_magnitude_test = np.transpose(np.abs(raos_test))  # store the magnitude of the complex response amplitude
raos_phase_test = np.transpose(np.rad2deg(np.angle(raos_test)))  # store the phase of the complex response amplitude


# Plot RAOs

plot_raos = True

max_encounter_frequency = 20  # Hz


if plot_raos:
    plt.plot(encounter_frequency_Hz, raos_magnitude[:, 2], '-x', label='SES-X shape')
    plt.plot(encounter_frequency_Hz, raos_magnitude_test[:, 2], '-x', label='Test')
    plt.xlabel('encounter frequency [Hz]')
    plt.ylabel('$\\eta_{3}/\\zeta_a$')
    plt.title('RAO in heave at 22kn and head sea')
    plt.xlim([0, max_encounter_frequency])
    plt.legend()
    plt.show()

    plt.plot(encounter_frequency_Hz, np.multiply(np.power(omega_e, 2), raos_magnitude[:, 2]), '-x', label='SES-X shape')
    plt.plot(encounter_frequency_Hz, np.multiply(np.power(omega_e, 2), raos_magnitude_test[:, 2]), '-x', label='Test')
    plt.xlabel('encounter frequency [Hz]')
    plt.ylabel('$|\\omega^2\\eta_{3}|/\\zeta_a[s^{-2}]$')
    plt.title('RAO for acceleration in heave')
    plt.xlim([0, max_encounter_frequency])
    plt.ylim([0, np.max(np.multiply(np.power(omega_e, 2), raos_magnitude_test[:, 2])[encounter_frequency_Hz < max_encounter_frequency])])
    plt.legend()
    plt.show()

    plt.plot(encounter_frequency_Hz, raos_magnitude[:, 6], fillstyle='none', marker='o', linestyle='-', label='SES-X shape')
    plt.plot(encounter_frequency_Hz, raos_magnitude_test[:, 6], fillstyle='none', marker='o', linestyle='-', label='Test')
    plt.xlabel('encounter frequency [Hz]')
    plt.title('RAO in uniform pressure at 22kn and head sea')
    plt.ylabel('$\\eta_{7} [-]$')
    plt.xlim([0, max_encounter_frequency])
    plt.legend()
    plt.show()

    plt.plot(encounter_frequency_Hz, np.divide(raos_magnitude[:, 4], k_encounter), '-x', label='SES-X shape')
    plt.plot(encounter_frequency_Hz, np.divide(raos_magnitude_test[:, 4], k_encounter), '-x', label='Test')
    plt.xlabel('encounter frequency [Hz]')
    plt.ylabel('$\\eta_{5}/k\\zeta_a$')
    plt.title('RAO in pitch at 22kn and head sea')
    plt.xlim([0, max_encounter_frequency])
    plt.legend()
    plt.show()

    plt.plot(encounter_wavelength, raos_magnitude[:, 2], '-x', label='SES-X shape')
    plt.plot(encounter_wavelength, raos_magnitude_test[:, 2], '-x', label='Test')
    plt.xlabel('encounter wavelength [m]')
    plt.ylabel('$\\eta_{3}/\\zeta_a$')
    plt.title('RAO in heave at 22kn and head sea')
    plt.xlim([0, max_encounter_frequency])
    plt.legend()
    plt.show()

    plt.plot(encounter_wavelength, np.multiply(np.power(omega_e, 2), raos_magnitude[:, 2]), '-x', label='python')
    plt.plot(encounter_wavelength, np.multiply(np.power(omega_e, 2), raos_magnitude_test[:, 2]), '-x', label='Test')
    plt.xlabel('encounter wavelength [m]')
    plt.ylabel('$|\\omega^2\\eta_{3}|/\\zeta_a[s^{-2}]$')
    plt.title('RAO for acceleration in heave')
    plt.xlim([0, max_encounter_frequency])
    plt.ylim([0, np.max(
        np.multiply(np.power(omega_e, 2), raos_magnitude[:, 2])[encounter_frequency_Hz < max_encounter_frequency])])
    plt.legend()
    plt.show()

    plt.plot(encounter_wavelength, raos_magnitude[:, 6], fillstyle='none', marker='o', linestyle='-',
             label='SES-X shape')
    plt.plot(encounter_wavelength, raos_magnitude_test[:, 6], fillstyle='none', marker='o', linestyle='-',
             label='Test')
    plt.xlabel('encounter wavelength [m]')
    plt.title('RAO in uniform pressure at 22kn and head sea')
    plt.ylabel('$\\eta_{7} [-]$')
    plt.xlim([0, max_encounter_frequency])
    plt.legend()
    plt.show()

    plt.plot(encounter_wavelength, np.divide(raos_magnitude[:, 4], k_encounter), '-x', label='SES-X shape')
    plt.plot(encounter_wavelength, np.divide(raos_magnitude_test[:, 4], k_encounter), '-x', label='Test')
    plt.xlabel('encounter wavelength [m]')
    plt.ylabel('$\\eta_{5}/k\\zeta_a$')
    plt.title('RAO in pitch at 22kn and head sea')
    plt.xlim([0, max_encounter_frequency])
    plt.show()
