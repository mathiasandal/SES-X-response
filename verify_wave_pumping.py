import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
from air_cushion import wave_pumping_rect, wave_pumping_excitation, wave_pumping_excitation_sesx, wave_pumping_steen

g = 9.81  # [m/s^2] acceleration of gravity

heading = 0  # [deg] wave heading
zeta_a = 1  # [m] wave amplitude
U = 50.  # [kn]
U = U * 0.51444  # [m/s]

# Properties of SES-X air cushion
l_1 = 12   # 0.0001  #     [m] length of the rectangular part of the air cushion
l_2 = 6  # 0.001  #      [m] length of the triangular part of the air cushion
b = 3.4  # [m] beam of the air cushion
# Location of motion coordinate system relative to intersection of AP, CL and BL
x_prime = 4  # [m]
z_prime = 3  # [m]

beta = heading  # [deg] wave heading relative to positive x-direction

# compute effective lengths
b_eff = (b * l_1 + 0.5 * b * l_2) / (l_1 + l_2)
l_eff = l_1 + 0.5 * l_2

# Properties of rectangular cushion
cushion_width = b  # [m]
x_s = l_1/2  # [m]
x_f = -l_1/2  # [m]
y_s = b/2  # [m]
y_p = -b/2  # [m]
x_b = -l_2 - l_1/2  # [m]

U = 22.  # [kn]
U = U * 0.51444  # [m/s]
omega_0 = np.linspace(0.5, 200, 1000000)
k = np.power(omega_0, 2) / g
omega_e = omega_0 + U / g * np.power(omega_0, 2)
f_enc = omega_e / 2 / np.pi
wavelength_enc = g / 2 / np.pi * (1 / f_enc)**2

wavelength = np.divide(2 * np.pi, k)

f_wp_rect = wave_pumping_steen(l_eff, b, omega_0, U)
f_wp_sesx = wave_pumping_excitation_sesx(x_f, x_s, y_s, x_b, omega_0, U, beta)

# Define plot colors
color_BBGreen = '#5cb16d'
color_BBPurple = '#b15ca0'

save_plots = True

plt.plot(f_enc, np.absolute(f_wp_rect), label=r'Rectangle', color=color_BBPurple)
plt.plot(f_enc, np.absolute(f_wp_sesx), label=r'BBGreen', color=color_BBGreen)
plt.xlabel(r'$\textrm{Encounter frequency}\,[Hz]$')
plt.ylabel(r'$|\hat{F}_{\textit{WP}}|\,/\,\zeta_a\,[m^2/s]$')
plt.xlim([0., 16])
plt.gca().set_ylim(bottom=0)
plt.legend()
if save_plots:
    plt.savefig('Results/Results and discussion/Wave pumping BBGreen/WP against encounter frequency.pdf', bbox_inches='tight')
plt.show()

plt.plot(wavelength, np.absolute(f_wp_rect), label='Rectangle', color=color_BBPurple)
plt.plot(wavelength, np.absolute(f_wp_sesx), label='BBGreen', color=color_BBGreen)
plt.xlabel(r'$\textrm{Wavelength}\,[m]$')
plt.ylabel(r'$|\hat{F}_{\textit{WP}}|\,/\,\zeta_a\,[m^2/s]$')
plt.xlim([0., 25])
plt.gca().set_ylim(bottom=0)
plt.legend()
if save_plots:
    plt.savefig('Results/Results and discussion/Wave pumping BBGreen/WP against wavelength.pdf', bbox_inches='tight')
plt.show()

"""
fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(wavelength, np.abs(f_wp_sesx), label='SES-X shape')
ax1.plot(wavelength, np.abs(f_wp_rect), label='Rectangle')
fig.suptitle(
    'Wave pumping excitation for SES-X air cushion. $l_1 = %.0f$m, ' % l_1 + '$l_2 = %.0f$m, ' % l_2 + '$b = %0.1f$m, ' % b + '$\\beta = %0.1f^{\circ}$' % beta)
ax1.set_xlabel("$\\lambda [m]$")
ax1.set_ylabel('$F_{wp} [m^3/s]$')
ax1.set_xlim([0, 12])
ax1.legend()

ax2.plot(wavelength, np.angle(f_wp_sesx), label='SES-X shape')
ax2.plot(wavelength, np.angle(f_wp_rect), label='Rectangle')
#ax2.title(
#    'Wave pumping excitation for SES-X air cushion. $l_1 = %.0f$m, ' % l_1 + '$l_2 = %.0f$m, ' % l_2 + '$b = %0.1f$m, ' % b + '$\\beta = %0.1f^{\circ}$' % beta)
'''
# Get fractions of pi as y label
y_pi = np.angle(f_wp_sesx)/np.pi
unit = 0.25
y_tick = np.arange(-0.5, 0.5+unit, unit)

y_label = [r"$-\frac{\pi}{2}$", r"$-\frac{\pi}{4}$", r"$0$", r"$+\frac{\pi}{4}$",   r"$+\frac{\pi}{2}$"]
ax2.set_yticks(y_tick*np.pi)
ax2.set_yticklabels(y_label, fontsize=12)
'''
ax2.set_xlabel("$\\lambda [m]$")
ax2.set_ylabel('Phase shift')
ax2.set_xlim([0, 12])
ax2.legend()
plt.show()


fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(frequencies_Hz, np.abs(f_wp_sesx), label='SES-X shape')
#ax1.plot(frequencies_Hz, np.abs(f_wp_rect), label='Rectangle')
fig.suptitle(
    'Wave pumping excitation for SES-X air cushion. $l_1 = %.0f$m, ' % l_1 + '$l_2 = %.0f$m, ' % l_2 + '$b = %0.1f$m, ' % b + '$\\beta = %0.1f^{\circ}$' % beta)
ax1.set_xlabel("encounter frequency [Hz]")
ax1.set_ylabel('$F_{wp} [m^3/s]$')
ax1.set_xlim([0, 12])
ax1.legend()

ax2.plot(frequencies_Hz, np.angle(f_wp_sesx), label='SES-X shape')
#ax2.plot(frequencies_Hz, np.angle(f_wp_rect), label='Rectangle')
ax2.set_xlabel("encounter frequency [Hz]")
ax2.set_ylabel('Phase shift')
ax2.set_xlim([0, 12])
ax2.legend()
plt.show()
"""