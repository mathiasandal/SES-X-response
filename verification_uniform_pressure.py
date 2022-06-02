"""
In this script RAOs are plotted for a conceptual SES with and without an air cushion
"""
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
from veres import read_group_of_re1_input
from air_cushion import wave_pumping_rect, wave_pumping_steen

save_plots = False

# Define plot colors
color_BBGreen = '#5cb16d'
color_BBPurple = '#b15ca0'

filepath_no_air = 'C:/Users/mathi/SIMA Workspaces/Workspace_1/Task/conceptual_SES_20m_no_air_cushion/DWL'
filepath_w_air = 'C:/Users/mathi/SIMA Workspaces/Workspace_1/Task/conceptual_SES_20m_w_air_cushion/DWL'
raos_no_air, VEL1, HEAD1, omega_01, XMTN1, ZMTN1 = read_group_of_re1_input(filepath_no_air)
raos_w_air, VEL2, HEAD2, omega_02, XMTN2, ZMTN2 = read_group_of_re1_input(filepath_w_air)

g = 9.81
beta = HEAD1[0]
frequencies = omega_01
frequency_Hz = frequencies / 2 / np.pi  # compute encounter frequency in Hz
k = np.power(frequencies, 2) / g  # encounter wave number
wavelength = np.divide(2*np.pi, k)


plt.plot(omega_01/2/np.pi, np.absolute(raos_no_air[2, :]), label='no air cushion', color=color_BBPurple, fillstyle='none', marker='o', linestyle='-')
plt.plot(omega_01/2/np.pi, np.absolute(raos_w_air[2, :]), label='with air cushion', color=color_BBGreen, marker='x', linestyle='-')
plt.legend()
plt.xlim([0, 3])
plt.ylabel(r'$|\hat{\eta}_3|\,/\,\zeta_a\,[-]$')
plt.xlabel(r'$\textrm{Encounter frequency}\,[Hz]$')
if save_plots:
    plt.savefig('Results/Verification VERES/heave against frequency.pdf', bbox_inches='tight')
plt.show()

plt.plot(omega_01/2/np.pi, np.divide(np.absolute(raos_no_air[4, :]), k), label='no air cushion', color=color_BBPurple, fillstyle='none', marker='o', linestyle='-')
plt.plot(omega_01/2/np.pi, np.divide(np.absolute(raos_w_air[4, :]), k), label='with air cushion', color=color_BBGreen, marker='x', linestyle='-')
plt.legend()
plt.xlim([0, 3])
plt.ylabel(r'$|\hat{\eta}_5|\,/\, k \zeta_a\,[-]$')
plt.xlabel(r'$\textrm{Encounter frequency}\,[Hz]$')
plt.show()

plt.plot(omega_01/2/np.pi, np.absolute(np.multiply(np.power(omega_01, 2), raos_no_air[2, :])), label='no air cushion', color=color_BBPurple, fillstyle='none', marker='o', linestyle='-')
plt.plot(omega_01/2/np.pi, np.absolute(np.multiply(np.power(omega_01, 2), raos_w_air[2, :])), label='with air cushion', color=color_BBGreen, marker='x', linestyle='-')
plt.legend()
plt.xlim([0, 7])
plt.ylabel(r'$\textrm{Vert. acc.}\,/\,\zeta_a\,[s^{-2}]$')
plt.xlabel(r'$\textrm{Encounter frequency}\,[Hz]$')
if save_plots:
    plt.savefig('Results/Verification VERES/vert acc against frequency.pdf', bbox_inches='tight')
plt.show()

plt.plot(omega_01/2/np.pi, np.divide(np.absolute(np.multiply(np.power(omega_01, 2), raos_no_air[2, :])), k), label='no air cushion', color=color_BBPurple, fillstyle='none', marker='o', linestyle='-')
plt.plot(omega_01/2/np.pi, np.divide(np.absolute(np.multiply(np.power(omega_01, 2), raos_w_air[2, :])), k), label='with air cushion', color=color_BBGreen, marker='x', linestyle='-')
plt.legend()
plt.xlim([0, 7])
plt.ylabel(r'$\textrm{Pitch acc.}\,/\,k\zeta_a\,[s^{-2}]$')
plt.xlabel(r'$\textrm{Encounter frequency}\,[Hz]$')
plt.show()

plt.plot(wavelength, np.absolute(raos_no_air[2, :]), label='no air cushion', color=color_BBPurple, fillstyle='none', marker='o', linestyle='-')
plt.plot(wavelength, np.absolute(raos_w_air[2, :]), label='with air cushion', color=color_BBGreen, marker='x', linestyle='-')
plt.legend()
plt.xlim([0, 40])
plt.ylabel(r'$|\hat{\eta}_3|\,/\,\zeta_a\,[-]$')
plt.xlabel(r'$\textrm{Wavelength}\,[m]$')
if save_plots:
    plt.savefig('Results/Verification VERES/heave against wavelength.pdf', bbox_inches='tight')
plt.show()


plt.plot(wavelength, np.absolute(np.multiply(np.power(omega_01, 2), raos_no_air[2, :])), label='no air cushion', color=color_BBPurple, fillstyle='none', marker='o', linestyle='-')
plt.plot(wavelength, np.absolute(np.multiply(np.power(omega_01, 2), raos_w_air[2, :])), label='with air cushion', color=color_BBGreen, marker='x', linestyle='-')
plt.legend()
plt.xlim([0, 40])
plt.ylabel(r'$\textrm{Vert. acc.}\,/\,\zeta_a\,[s^{-2}]$')
plt.xlabel(r'$\textrm{Wavelength}\,[m]$')
if save_plots:
    plt.savefig('Results/Verification VERES/vert acc against wavelength.pdf', bbox_inches='tight')
plt.show()


# Compute and plot wave pumping
# properties
L = 20  # [m] length of the vessel
b = 7  # [m] beam of air cushion
p_s = 0  # [Pa] Membrane seal pressure
p_0 = 3500  # [Pa] excess pressure in the air cushion

# Compute wave pumping excitation
x_s = L/2
x_f = -L/2
y_s = b/2
y_p = -b/2

omegas = np.linspace(0.1*2*np.pi, 2*np.pi*7, 10000)
f_ex_7 = wave_pumping_rect(x_f, x_s, y_p, y_s, omegas, HEAD1)
#f_ex_7_2 = wave_pumping_steen(L, b, omegas)
plt.xlabel(r'$\textrm{Encounter frequency}\,[Hz]$')
plt.ylabel(r'$|\hat{F}_{WP}|\,[m^3/s]$')
plt.plot(omegas/2/np.pi, np.absolute(f_ex_7), color=color_BBGreen)
if save_plots:
    plt.savefig('Results/Verification VERES/F_wp against frequency.pdf', bbox_inches='tight')
plt.show()


frequencies = omegas
k = np.power(frequencies, 2) / g  # encounter wave number
wavelength = np.divide(2*np.pi, k)
plt.xlabel(r'$\textrm{Wavelength}\,[m]$')
plt.ylabel(r'$|\hat{F}_{WP}|\,[m^3/s]$')
plt.xlim([0, 40])
plt.plot(wavelength, np.absolute(f_ex_7), color=color_BBGreen)
if save_plots:
    plt.savefig('Results/Verification VERES/F_wp against wavelength.pdf', bbox_inches='tight')
plt.show()



