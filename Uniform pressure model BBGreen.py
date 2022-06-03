from veres import read_coefficients_from_veres
from air_cushion import air_cushion_area, read_fan_characteristics, interpolate_fan_characteristics, \
    wave_pumping_excitation_sesx
import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

# ***** Physical parameters *****
rho = 1025  # [kg/m^3] density of sea water
p_a = 101325  # [Pa] atmospheric pressure
g = 9.81  # [m/s^2] acceleration of gravity
gamma = 1.4  # [-] Ratio of specific heat for air, when assuming adiabatic density relation

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
h_0 = 0.64  #2 # [m]  Height of the air cushion  # TODO: make sure this value is correct
h = p_0 / rho / g  # [m] Difference in water level inside and outside the air cushion

# Compute area and longitudinal position of centroid relative to AP
A_c0, x_B_c = air_cushion_area(l_1, l_2, b)

V_c0 = A_c0 * h_0  # [m^3] Compute mean air cushion volume  # TODO: Make sure this is correct

# Read in fan characteristics
Q, P, rpm_dummy = read_fan_characteristics('Input files/fan characteristics/fan characteristics.csv', '1800rpm')

# Interpolates values
Q_0, dQdp_0 = interpolate_fan_characteristics(p_0, P, Q)

# Read hydrodynamic coefficients and loads from VERES

veres_formulation = 'high-speed'  # 'strip-theory'  #
if veres_formulation == 'high-speed':
    # For high-speed formulation
    path_high_speed_theory = 'C:/Users/mathi/SIMA Workspaces/Workspace_1/Task/BBGreen_uniform_pressure_model_hs_theory/DWL'
    A_temp, B_temp, C_temp, f_ex_temp, omega_0, U, beta, XMTN, ZMTN, lcg_veres, vcg_veres = read_coefficients_from_veres(path_high_speed_theory)
else:
    # For strip_theory
    path_strip_theory = 'C:/Users/mathi/SIMA Workspaces/Workspace_1/Task/BBGreen_uniform_pressure_model_strip_theory/DWL'
    A_temp, B_temp, C_temp, f_ex_temp, omega_0, U, beta, XMTN, ZMTN, lcg_veres, vcg_veres = read_coefficients_from_veres(path_strip_theory)

k = np.power(omega_0, 2) / g
omega_e = omega_0 + U / g * np.multiply(np.power(omega_0, 2), np.cos(np.rad2deg(beta)))  # [rad/s] encounter freq.
f_enc = omega_e / 2 / np.pi  # [Hz] encounter frequency
k_enc = np.power(omega_e, 2) / g

# Initialize vectors for hydrodynamic coeffs. and loads
A = np.zeros([len(omega_0), 3, 3])
B = np.zeros([len(omega_0), 3, 3])
C = np.zeros([3, 3])
f_ex = np.zeros([3, len(omega_0)], dtype=complex)

# Add terms from VERES
# Added mass
A[:, 0, 0] = A_temp[:, 2, 2]  # A_33
A[:, 0, 1] = A_temp[:, 2, 4]  # A_35
A[:, 1, 0] = A_temp[:, 4, 2]  # A_53
A[:, 1, 1] = A_temp[:, 4, 4]  # A_55
# Damping
B[:, 0, 0] = B_temp[:, 2, 2]  # B_33
B[:, 0, 1] = B_temp[:, 2, 4]  # B_35
B[:, 1, 0] = B_temp[:, 4, 2]  # B_53
B[:, 1, 1] = B_temp[:, 4, 4]  # B_55
# Restoring
C[0, 0] = C_temp[2, 2]  # C_33
C[0, 1] = C_temp[2, 4]  # C_35
C[1, 0] = C_temp[4, 2]  # C_53
C[1, 1] = C_temp[4, 4]  # C_55
# Excitation(diffraction and Froude-Kriloff)
f_ex[0, :] = f_ex_temp[2, :]  # f_ex_3
f_ex[1, :] = f_ex_temp[4, :]  # f_ex_5

# Center of motion relative to BL and AP
x_prime = Lpp/2 - XMTN[0]
z_prime = ZMTN[0]

# Create mass matrix
x_G = x_prime - lcg  # [m] longitudinal position of COG relative to motion coordinate system
z_G = -(z_prime - vcg)  # [m] vertical position of COG relative to motion coordinate system
M = np.zeros([3, 3])
M[0, 0] = total_mass  # [kg] M_33
M[0, 1] = -x_G * total_mass  # [kg m] M_35
M[1, 0] = -x_G * total_mass  # [kg m] M_53
M[1, 1] = r55**2 * total_mass  # [kg m^2] I_55

# Centroid of the air cushion at equilibrium relative to the motion coord. system
x_cp = x_prime - x_B_c  # [m] longitudinal position
y_cp = 0  # [m] transverse position
z_cp = -h / 2  # [m] vertical position

# Properties of SES-X shaped air cushion with origin in motion coordinate system
cushion_width = b  # [m]
x_s = x_prime  # [m]
x_f = -(l_1 - x_prime)  # [m]
y_s = b/2  # [m]
y_p = -b/2  # [m]
x_b = -(l_1 + l_2 - x_prime)  # [m]

# Append wave pumping to excitation vector
f_ex[2, :] += wave_pumping_excitation_sesx(x_f, x_s, y_s, x_b, omega_0, U, b)  # f_ex_wp

# Append damping terms due to air cushion
B[2, 0] += A_c0  # B_73
B[2, 1] += x_cp * A_c0  # B_75
B[2, 2] += p_0 * V_c0 / gamma / (p_a + p_0)  # B_77

# Append restoring terms due to air cushion
C[1, 1] += p_0 * z_cp * A_c0  # C_55
C[0, 2] += -p_0 * A_c0  # C_37
C[1, 2] += p_0 * x_cp * A_c0  # C_57
C[2, 2] += 0.5*Q_0 - p_0 * dQdp_0  # C_77

# Adjust stiffness in pitch due to change in center of gravity relative to VERES
C[1, 1] += total_mass * g * (vcg_veres - vcg)  # C_55

# Initialize array to contain RAOs
raos = np.zeros([3, len(omega_0)], dtype=complex)

for i in range(len(omega_0)):
    raos[:, i] = np.linalg.solve(-omega_e[i]**2*(A[i] + M) + 1j * omega_e[i] * B[i] + C, f_ex[:, i])

# Iterate natural frequencies
nat_freq, eigen_modes = la.eig(C, M + A[999])

f_03 = np.sqrt(nat_freq[0])/2/np.pi
f_05 = np.sqrt(nat_freq[1])/2/np.pi

# Heave and pitch resonance
def it_nat_freq(M, A_vec, C, omega_e, g=9.81):

    A_temp = A_vec[len(omega_e) // 2]

    omega_nat_old = np.sqrt(C / (M + A_temp))

    err = -1
    epsi = 1e-6
    counter = 1

    while err > epsi or counter == 1 or counter > 100:

        A_temp = np.interp(omega_nat_old, omega_e, A_vec)

        omega_nat = np.sqrt(C / (M + A_temp))

        err = np.abs((omega_nat - omega_nat_old)/omega_nat_old)

        omega_nat_old = omega_nat

        counter += 1

    return omega_nat

# heave - iterated
omega_03_BBGreen_it = it_nat_freq(M[0, 0], A[:, 0, 0], C[0, 0], omega_e)
f_03_it = omega_03_BBGreen_it / 2 / np.pi

# pitch - iterated
omega_05_it = it_nat_freq(M[1, 1], A[:, 1, 1], C[1, 1], omega_e)
f_05_it = omega_05_it / 2 / np.pi

# Plotting RAOs

# Define plot colors
color_BBGreen = '#5cb16d'
color_BBPurple = '#b15ca0'

plt.plot(f_enc, np.absolute(A[:, 0, 0]), color=color_BBGreen)
plt.xlabel(r'\textrm{Encounter frequency} $[Hz]$')
plt.show()

# heave

plt.plot(f_enc, np.absolute(raos[0]), color=color_BBGreen)
plt.xlabel(r'\textrm{Encounter frequency} $[Hz]$')
plt.ylabel(r'$|\hat{\eta}_3|\,/\,\zeta_a\,[-]$')
plt.show()

plt.plot(f_enc, np.divide(np.absolute(raos[1]), k), color=color_BBGreen)
plt.xlabel(r'\textrm{Encounter frequency} $[Hz]$')
plt.ylabel(r'$|\hat{\eta}_5|\,/\,k \zeta_a\,[-]$')
plt.show()

plt.plot(f_enc, np.absolute(raos[2]), color=color_BBGreen)
plt.xlabel(r'\textrm{Encounter frequency} $[Hz]$')
plt.ylabel(r'$|\hat{\eta}_7|\,[-]$')
plt.show()

plt.plot(f_enc, np.absolute(np.multiply(raos[0], np.power(omega_e, 2))), color=color_BBGreen)
plt.xlabel(r'\textrm{Encounter frequency} $[Hz]$')
plt.ylabel(r'$|\omega_e^2\hat{\eta}_3|\,/\,\zeta_a\,[-]$')
plt.show()



