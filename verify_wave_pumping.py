import numpy as np
import matplotlib.pyplot as plt
from air_cushion import wave_pumping_rect, wave_pumping_excitation, wave_pumping_excitation_sesx


g = 9.81  # [m/s^2] acceleration of gravity

heading = 90  # [deg] wave heading
zeta_a = 1  # [m] wave amplitude

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

frequencies = np.linspace(0.5, 17, 10000)
k = np.power(frequencies, 2) / g
wavelength = np.divide(2 * np.pi, k)


f_wp_rect = 1j*wave_pumping_rect(x_b, x_s, y_p, y_s, frequencies, beta)
f_wp_sesx = wave_pumping_excitation_sesx(x_f, x_s, y_s, x_b, frequencies, beta)

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
