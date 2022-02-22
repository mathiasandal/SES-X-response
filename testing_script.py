import numpy as np
import matplotlib.pyplot as plt
from air_cushion import wave_pumping_rectangle, wave_pumping_rect, wave_pumping_excitation

g = 9.81
omega = 1.083308

beta = 0  # [deg] wave heading
zeta_a = 1  # [m] wave amplitude

# Properties of SES-X air cushion
l_1 = 12  # [m] length of the rectangular part of the air cushion
l_2 = 6  # [m] length of the triangular part of the air cushion
b = 3.4  # [m] beam of the air cushion
# Location of motion coordinate system relative to intersection of AP, CL and BL
x_prime = 8  # [m]
z_prime = 3  # [m]

# compute effective lengths
b_eff = (b * l_1 + 0.5 * b * l_2) / (l_1 + l_2)
l_eff = l_1 + 0.5 * l_2
print('Effective beam: ' + str(b_eff) + '[m]')
print('Effective length: ' + str(l_eff) + '[m]')

# Properties of rectangular cushion
L = 18  # [m]
cushion_width = 3.4  # [m]
x_s = L/2  # [m]
x_f = -L/2  # [m]
y_s = cushion_width/2  # [m]
y_p = -cushion_width/2  # [m]

frequencies = np.linspace(0.5, 17, 10000)
k = np.power(frequencies, 2) / g
wavelength = np.divide(2 * np.pi, k)


f_wp_rect = wave_pumping_rect(x_f, x_s, y_p, y_s, frequencies, beta)
f_wp_sesx = wave_pumping_excitation(b, l_1, l_2, x_prime, beta, frequencies)


plot_type = 3  # determines what the x-axis should be
plot_rect = True
# Plot for rectangular air cushion
if plot_rect:
    if plot_type == 1:
        plt.plot(frequencies, np.abs(f_wp_rect))
        plt.xlabel('Encounter frequency $[rad/s]$')
    elif plot_type == 2:
        plt.plot(k, np.abs(f_wp_rect))
        plt.xlabel('wave number $[m^{-1}]$')
    elif plot_type == 3:
        plt.plot(wavelength, np.abs(f_wp_rect))
        plt.xlabel("$\\lambda [m]$")
        plt.ylim([0, np.amax(np.abs(f_wp_rect))])
        plt.xlim([1, 60])

    plt.ylabel('$F_{wp} [m^3/s]$')
    plt.title('Wave pumping excitation for ' + str(x_s - x_f) + '[m]$\\times$' + str(y_s - y_p) + '[m] air cushion. $\\beta=$' + str(beta) + '$^{\circ}$')
    plt.show()

# Plot for SES-X air cushion
plot_sesx = True

if plot_sesx:
    if plot_type == 1:
        plt.plot(frequencies, np.abs(f_wp_sesx))
        plt.xlabel('Encounter frequency $[rad/s]$')
    elif plot_type == 2:
        plt.plot(k, np.abs(f_wp_rect))
        plt.xlabel('wave number $[m^{-1}]$')
    elif plot_type == 3:
        plt.plot(wavelength, np.abs(f_wp_sesx))
        plt.xlabel("$\\lambda [m]$")
        plt.ylim([0, np.amax(np.abs(f_wp_sesx))])
        plt.xlim([0, 25])

    plt.ylabel('$F_{wp} [m^3/s]$')
    plt.title('Wave pumping excitation for SES-X air cushion. $l_1 = %.0f$m, ' % l_1 + '$l_2 = %.0f$m, ' % l_2 + '$b = %0.1f$m, ' % b + '$\\beta = %0.1f^{\circ}$' % beta)
    plt.show()

    plt.plot(wavelength, np.abs(f_wp_sesx), label='SES-X shape')
    plt.plot(wavelength, np.abs(f_wp_rect), label='Rectangle')
    plt.title('Wave pumping excitation for SES-X air cushion. $l_1 = %.0f$m, ' % l_1 + '$l_2 = %.0f$m, ' % l_2 + '$b = %0.1f$m, ' % b + '$\\beta = %0.1f^{\circ}$' % beta)
    plt.xlabel("$\\lambda [m]$")
    plt.xlim([0, 40])
    plt.legend()
    plt.show()

"""
f_wp_fasit = wave_pumping_rectangle(x_f, x_s, y_p, y_s, omega, beta, zeta_a, g)

a = 1j * omega * zeta_a * (y_s - y_p) / k / np.cos(beta)

theta = k * x_s * np.cos(beta)
f_wp_1 = a * (np.exp(-1j * theta) - np.exp(1j * theta))
f_wp_2 = -2j * a * np.sin(theta)
"""

