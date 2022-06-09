import pandas as pd
import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

U = 22  # [kn]
U = 0.5144*U  # [m/s]
g = 9.81
L = 18

x_size = 3.6


# Read RAOs from csv-files
df_1 = pd.read_csv('Results/RAOs/RAO_with_wave_pumping.csv')
df_2 = pd.read_csv('Results/RAOs/RAO_without_wave_pumping.csv')
df_3 = pd.read_csv('Results/RAOs/RAO_without_wave_pumping_and_hull_excitation.csv')

omega_0 = df_1['omega_0 [rad/s]'].astype(float)
k = np.power(omega_0, 2) / g  # wave number of water waves
omega_e = omega_0 + np.power(omega_0, 2) / g * U  # encounter frequencies
f_encounter = omega_e / 2 / np.pi  # [Hz] frequency of encounter
wavelength = np.divide(2*np.pi, k)
n_freq = len(omega_e)
b = np.power(omega_0, 2)

eta_3a_1 = df_1['eta_3a [m]'].astype(complex)
eta_3a_2 = df_2['eta_3a [m]'].astype(complex)
eta_3a_3 = df_3['eta_3a [m]'].astype(complex)

eta_5a_1 = df_1['eta_5a [rad]'].astype(complex)
eta_5a_2 = df_2['eta_5a [rad]'].astype(complex)
eta_5a_3 = df_3['eta_5a [rad]'].astype(complex)

eta_7a_1 = df_1['mu_ua [-]'].astype(complex)
eta_7a_2 = df_2['mu_ua [-]'].astype(complex)
eta_7a_3 = df_3['mu_ua [-]'].astype(complex)

# Define plot colors
color_BBGreen = '#5cb16d'
color_BBPurple = '#b15ca0'


save_RAOs_for_comparison = False

f_encounter_new_range = np.linspace(0.1001, 16.-0.01, 150)
f = interpolate.interp1d(f_encounter, np.absolute(eta_7a_3))
eta_7a_3_new_range = f(f_encounter_new_range)

plt.plot(f_encounter, np.absolute(eta_7a_1), label=r'\textrm{All excitations}', color=color_BBGreen)
plt.plot(f_encounter, np.absolute(eta_7a_2), label=r'\textrm{No wave pumping}', color=color_BBPurple)
plt.plot(f_encounter_new_range, eta_7a_3_new_range, label=r'\textrm{No wave pumping and hull excitation}', markersize=x_size, fillstyle='none', color='k', linestyle='', marker='x')
plt.xlabel(r'\textrm{Encounter frequency} $[Hz]$')
plt.ylabel(r'$|\hat{\eta}_{7}|\,/\,\zeta_a\,[m^{-1}]$')
x_min = 0.
plt.xlim([x_min, 16])
plt.ylim([0, np.max(np.abs(eta_7a_1[f_encounter > x_min]))])
plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/RAOs/uniform pressure.pdf', bbox_inches='tight')
plt.show()

# against wavelength
wavelength_new_range = np.linspace(min(wavelength), 25.-0.01, 150)
f = interpolate.interp1d(wavelength, np.absolute(eta_7a_3))
eta_7a_3_new_range = f(wavelength_new_range)

plt.plot(wavelength, np.absolute(eta_7a_1), label=r'\textrm{All excitations}', color=color_BBGreen)
plt.plot(wavelength, np.absolute(eta_7a_2), label=r'\textrm{No wave pumping}', color=color_BBPurple)
plt.plot(wavelength_new_range, eta_7a_3_new_range, label=r'\textrm{No wave pumping and hull excitation}', markersize=x_size, fillstyle='none', color='k', linestyle='', marker='x')
plt.xlabel(r'\textrm{Wavelength} $[m]$')
plt.ylabel(r'$|\hat{\eta}_{7}|\,/\,\zeta_a\,[m^{-1}]$')
x_min = 0.
plt.xlim([x_min, 25.])
plt.ylim([0, np.max(np.abs(eta_7a_1[wavelength > x_min]))])
plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/RAOs/uniform pressure against wavelength.pdf', bbox_inches='tight')
plt.show()


df_vert_acc = pd.read_csv(
    'C:/Users/mathi/OneDrive - NTNU/Master Thesis/Spatially varying pressure/Results from Steen and Faltinsen (1995)/Rigid panel model/Vert. Acc. AP/vertical_acc_AP_RAO_nice_format.csv')

vert_acc_FP_1 = np.absolute(omega_e ** 2 * (eta_3a_1 + L / 2 * eta_5a_1))
vert_acc_FP_2 = np.absolute(omega_e ** 2 * (eta_3a_2 + L / 2 * eta_5a_2))
vert_acc_FP_3 = np.absolute(omega_e ** 2 * (eta_3a_3 + L / 2 * eta_5a_3))

vert_motion_FP_1 = np.absolute((eta_3a_1 + L / 2 * eta_5a_1))
vert_motion_FP_2 = np.absolute((eta_3a_2 + L / 2 * eta_5a_2))
vert_motion_FP_3 = np.absolute((eta_3a_3 + L / 2 * eta_5a_3))

f_encounter_new_range = np.linspace(0.1001, 16.-0.01, 150)
f = interpolate.interp1d(f_encounter, np.absolute(vert_acc_FP_3))
vert_acc_FP_3_new_range = f(f_encounter_new_range)

plt.plot(f_encounter, np.absolute(vert_acc_FP_1), label=r'\textrm{All excitations}', color=color_BBGreen)
plt.plot(f_encounter, np.absolute(vert_acc_FP_2), label=r'\textrm{No wave pumping}', color=color_BBPurple)
plt.plot(f_encounter_new_range, vert_acc_FP_3_new_range, label=r'\textrm{No wave pumping and hull excitation}', markersize=x_size, fillstyle='none', color='k', linestyle='', marker='x')
plt.xlim([x_min, 16.])
#plt.ylim([0, np.max(np.abs(vert_acc_FP_1[f_encounter > x_min]))])
plt.gca().set_ylim(bottom=0)
plt.xlabel(r'\textrm{Encounter frequency} $[Hz]$')
plt.ylabel(r'\textrm{Vert. acc. at the bow} $\,/\,\zeta_a\,[s^{-2}]$')
plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/RAOs/vert acc bow against encounter frequency.pdf', bbox_inches='tight')
plt.show()

# against wavelength
wavelength_new_range = np.linspace(min(wavelength), 25.-0.01, 150)
f = interpolate.interp1d(wavelength, vert_acc_FP_3)
vert_acc_FP_3_new_range = f(wavelength_new_range)

plt.plot(wavelength, vert_acc_FP_1, label=r'\textrm{All excitations}', color=color_BBGreen)
plt.plot(wavelength, vert_acc_FP_2, label=r'\textrm{No wave pumping}', color=color_BBPurple)
plt.plot(wavelength_new_range, vert_acc_FP_3_new_range, label=r'\textrm{No wave pumping and hull excitation}', markersize=x_size, fillstyle='none', color='k', linestyle='', marker='x')
plt.xlim([x_min, 25.])
#plt.ylim([0, np.max(np.abs(vert_acc_FP_1[f_encounter > x_min]))])
plt.gca().set_ylim(bottom=0)
plt.xlabel(r'\textrm{Wavelength} $[m]$')
plt.ylabel(r'\textrm{Vert. acc. at the bow} $\,/\,\zeta_a\,[s^{-2}]$')
plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/RAOs/vert acc bow against wavelength.pdf', bbox_inches='tight')
plt.show()


f_encounter_new_range = np.linspace(0.1001, 16.-0.01, 100)
f = interpolate.interp1d(f_encounter, np.absolute(vert_motion_FP_3))
vert_motion_FP_3_new_range = f(f_encounter_new_range)

plt.plot(f_encounter, np.absolute(vert_motion_FP_1), label=r'\textrm{All excitations}', color=color_BBGreen)
plt.plot(f_encounter, np.absolute(vert_motion_FP_2), label=r'\textrm{No wave pumping}', color=color_BBPurple)
plt.plot(f_encounter_new_range, vert_motion_FP_3_new_range, label=r'\textrm{No wave pumping and hull excitation}', markersize=x_size, fillstyle='none', color='k', linestyle='', marker='x')
plt.xlim([x_min, 16.])
#plt.ylim([0, np.max(vert_motion_FP_1[f_encounter > x_min])])
plt.gca().set_ylim(bottom=0)
plt.xlabel(r'\textrm{Encounter frequency} $[Hz]$')
plt.ylabel(r'\textrm{Vert. motion at the bow} $\,/\,\zeta_a\,[-]$')
plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/RAOs/vert motion bow against encounter frequency.pdf', bbox_inches='tight')
plt.show()

# against wavelength
wavelength_new_range = np.linspace(min(wavelength), 25.-0.01, 150)
f = interpolate.interp1d(wavelength, vert_motion_FP_3)
vert_motion_FP_3_new_range = f(wavelength_new_range)

plt.plot(wavelength, vert_motion_FP_1, label=r'\textrm{All excitations}', color=color_BBGreen)
plt.plot(wavelength, vert_motion_FP_2, label=r'\textrm{No wave pumping}', color=color_BBPurple)
plt.plot(wavelength_new_range, vert_motion_FP_3_new_range, label=r'\textrm{No wave pumping and hull excitation}', markersize=x_size, fillstyle='none', color='k', linestyle='', marker='x')
plt.xlim([x_min, 25.])
#plt.ylim([0, np.max(np.abs(vert_acc_FP_1[f_encounter > x_min]))])
plt.gca().set_ylim(bottom=0)
plt.xlabel(r'\textrm{Wavelength} $[m]$')
plt.ylabel(r'\textrm{Vert. motion at the bow} $\,/\,\zeta_a\,[s^{-2}]$')
plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/RAOs/vert motion bow against wavelength.pdf', bbox_inches='tight')
plt.show()

# Heave RAO

f_encounter_new_range = np.linspace(0.1001, 16.-0.01, 100)
f = interpolate.interp1d(f_encounter, np.absolute(eta_3a_3))
eta_3a_3_new_range = f(f_encounter_new_range)

plt.plot(f_encounter, np.absolute(eta_3a_1), label=r'\textrm{All excitations}', color=color_BBGreen)
plt.plot(f_encounter, np.absolute(eta_3a_2), label=r'\textrm{No wave pumping}', color=color_BBPurple)
plt.plot(f_encounter_new_range, eta_3a_3_new_range, label=r'\textrm{No wave pumping and hull excitation}', markersize=x_size, fillstyle='none', color='k', linestyle='', marker='x')
plt.xlim([x_min, 16.])
plt.gca().set_ylim(bottom=0)
plt.xlabel(r'\textrm{Encounter frequency} $[Hz]$')
plt.ylabel(r'$|\hat{\eta}_3|\,/\,\zeta_a\,[-]$')
plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/RAOs/heave against encounter frequency.pdf', bbox_inches='tight')
plt.show()

# Against wavelength

# against wavelength
wavelength_new_range = np.linspace(min(wavelength), 25.-0.01, 150)
f = interpolate.interp1d(wavelength, np.absolute(eta_3a_3))
eta_3a_3_new_range = f(wavelength_new_range)

plt.plot(wavelength, np.absolute(eta_3a_1), label=r'\textrm{All excitations}', color=color_BBGreen)
plt.plot(wavelength, np.absolute(eta_3a_2), label=r'\textrm{No wave pumping}', color=color_BBPurple)
plt.plot(wavelength_new_range, eta_3a_3_new_range, label=r'\textrm{No wave pumping and hull excitation}', markersize=x_size, fillstyle='none', color='k', linestyle='', marker='x')
plt.xlim([x_min, 25.])
plt.gca().set_ylim(bottom=0)
plt.xlabel(r'\textrm{Wavelength} $[m]$')
plt.ylabel(r'$|\hat{\eta}_3|\,/\,\zeta_a\,[-]$')
plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/RAOs/heave against wavelength.pdf', bbox_inches='tight')
plt.show()

'''
plt.plot(water_wavelength, np.absolute(eta_3a), label=r'\textrm{Computed}', color=color_BBGreen)
plt.xlim([x_min, 50.])
plt.xlabel(r'$\textrm{Wavelength}\,[m]$')
plt.ylabel(r'$\hat{\eta}_3\,/\,\zeta_a\,[-]$')
#plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/RAOs/no WP heave against wavelength.pdf', bbox_inches='tight')
plt.show()
'''

# Pitch RAO
f_encounter_new_range = np.linspace(0.1001, 16.-0.01, 100)
f = interpolate.interp1d(f_encounter, np.absolute(eta_5a_3))
eta_5a_3_new_range = f(f_encounter_new_range)

f = interpolate.interp1d(f_encounter, k)
k_new_range = f(f_encounter_new_range)

plt.plot(f_encounter, np.absolute(eta_5a_1), label=r'\textrm{All excitations}', color=color_BBGreen)
plt.plot(f_encounter, np.absolute(eta_5a_2), label=r'\textrm{No wave pumping}', color=color_BBPurple)
plt.plot(f_encounter_new_range, eta_5a_3_new_range, label=r'\textrm{No wave pumping and hull excitation}', markersize=x_size, fillstyle='none', color='k', linestyle='', marker='x')
plt.xlabel(r'\textrm{Encounter frequency} $[Hz]$')
plt.ylabel(r'$|\hat{\eta}_5|\,\,/\,\zeta_a\,[rad/m]$')
plt.xlim([0.0, 16.])
plt.gca().set_ylim(bottom=0.0)
#plt.ylim([0, np.max(np.divide(np.absolute(eta_5a_2), k)[f_encounter > 0.0])])
plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/RAOs/pitch against encounter frequency.pdf', bbox_inches='tight')
plt.show()

# against wavelength
wavelength_new_range = np.linspace(min(wavelength), 25.-0.01, 150)
f = interpolate.interp1d(wavelength, np.absolute(eta_5a_3))
eta_5a_3_new_range = f(wavelength_new_range)

f = interpolate.interp1d(wavelength, k)
k_new_range = f(wavelength_new_range)

plt.plot(wavelength, np.absolute(eta_5a_1), label=r'\textrm{All excitations}', color=color_BBGreen)
plt.plot(wavelength, np.absolute(eta_5a_2), label=r'\textrm{No wave pumping}', color=color_BBPurple)
plt.plot(wavelength_new_range, eta_5a_3_new_range, label=r'\textrm{No wave pumping and hull excitation}', markersize=x_size, fillstyle='none', color='k', linestyle='', marker='x')
plt.xlim([x_min, 25.])
plt.gca().set_ylim(bottom=0)
#plt.ylim([0, np.max(vert_motion_FP_nondim[f_encounter > x_min])])
plt.xlabel(r'\textrm{Wavelength} $[m]$')
plt.ylabel(r'$|\hat{\eta}_5|\,\,/\,\zeta_a\,[rad/m]$')
plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/RAOs/pitch against wavelength.pdf', bbox_inches='tight')
plt.show()


'''
plt.plot(water_wavelength, np.absolute(np.divide(eta_5a, k)), label=r'\textrm{Computed}', color=color_BBGreen)
plt.xlim([x_min, 50.])
plt.xlabel(r'$\textrm{Wavelength}\,[m]$')
plt.ylabel(r'$\hat{\eta}_5\,\,/\,k\zeta_a\,[-]$')
#plt.legend()
if save_RAOs_for_comparison:
    plt.savefig('Results/RAOs/no WP pitch against wavelength.pdf', bbox_inches='tight')
plt.show()

print('The iteration scheme converged after', counter)
print('Relative error:\t\t', rel_err)
print('Absolute error:\t\t', abs_err)
'''