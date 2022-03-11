import numpy as np
import matplotlib.pyplot as plt
from air_cushion import wave_pumping_rect
from basic_units import radians, degrees, cos

g = 9.81  # [m/s^2] acceleration of gravity

heading = 0  # [deg] wave heading
zeta_a = 1  # [m] wave amplitude

# Incoming wave properties
frequencies = np.linspace(0.01, 200, 100000)  # [rad/s] wave frequencies
k = np.power(frequencies, 2) / g  # [m^-1] wave number
wavelength = np.divide(2 * np.pi, k)  # [m] wavelength
frequencies_Hz = frequencies / 2 / np.pi  # [Hz] wave frequencies
beta = heading  # [deg] wave heading relative to positive x-direction

# Properties
L = 20  # [m] length of the vessel
b_c = 7  # [m] beam of air cushion

# Compute wave pumping excitation
x_s = L/2
x_f = -L/2
y_s = b_c/2
y_p = -b_c/2

f_wp = 1j * wave_pumping_rect(x_f, x_s, y_p, y_s, frequencies, beta)

x_min = 0
x_max = 24

# Define plot colors
color_BBGreen = '#5cb16d'
color_BBPurple = '#b15ca0'
plt.rc('font', size=15)  # sets size

# Plot amplitude against wavelength
plot_wavelength = True
if plot_wavelength:
    plt.plot(wavelength, np.abs(f_wp), label='Conceptual SES', color=color_BBGreen)
    plt.title('$\\beta = %0.1f^{\\circ}$, ' % beta + '$\\zeta_a = %.0f$ $[m]$' % zeta_a)
    plt.xlabel('$\\lambda [m]$')
    plt.ylabel('$\\left|\\hat{F}_{7}\\right|$ $[m^3/s]$')
    plt.xlim([x_min, x_max])
    plt.ylim([np.min(np.abs(f_wp)[(wavelength > x_min) & (wavelength < x_max)]),
              np.max(np.abs(f_wp)[(wavelength > x_min) & (wavelength < x_max)]) + 1])

    save_plot = True
    if save_plot:
        plot_path = 'Results/Conceptual SES/wave pumping/'
        plot_name = 'magnitude,wavelength,beta=%0.1fdeg, ' % beta + 'zeta_a=%.0fm' % zeta_a + '.pdf'
        plt.tight_layout()
        plt.savefig(plot_path + plot_name)
    plt.show()

# Plot amplitude against frequency
plot_wavelength = False
if plot_wavelength:
    plt.plot(frequencies, np.abs(f_wp), label='Conceptual SES', color=color_BBGreen)
    plt.title('$\\beta = %0.1f^{\\circ}$, ' % beta + '$\\zeta_a = %.0f$ $[m]$' % zeta_a)
    plt.xlabel('wave frequency [rad/s]')
    plt.ylabel('$\\left|\\hat{F}_{7}\\right|$ $[m^3/s]$')
    plt.xlim([x_min, x_max])
    plt.ylim([np.min(np.abs(f_wp)[(frequencies > x_min) & (frequencies < x_max)]),
              np.max(np.abs(f_wp)[(frequencies > x_min) & (frequencies < x_max)]) + 1])

    save_plot = False
    if save_plot:
        plot_path = 'Results/Conceptual SES/wave pumping/'
        plot_name = 'magnitude,frequency,beta=%0.1fdeg, ' % beta + 'zeta_a=%.0fm' % zeta_a + '.pdf'
        plt.tight_layout()
        plt.savefig(plot_path + plot_name)
    plt.show()


plt.plot(wavelength, np.angle(f_wp), yunits=degrees, label='Conceptual SES', color=color_BBPurple)
plt.title('$\\beta = %0.1f^{\\circ}$, ' % beta + '$\\zeta_a = %.0f$ $[m]$' % zeta_a)
plt.xlabel('$\\lambda [m]$')
plt.ylabel('$\\arg\\left(\\hat{F}_{7}\\right)$ $[rad]$')
plt.xlim([x_min, x_max])
#plt.show()

'''
fig, (ax1, ax2) = plt.subplots(2)
ax1.plot(wavelength, np.abs(f_wp), label='Conceptual SES')
fig.suptitle('$L = %.0f$ $[m]$, ' % L + '$b_c = %.0f$ $[m]$, ' % b_c + '$\\beta = %0.1f^{\\circ}$, ' % beta + '$\\zeta_a = %.0f$ $[m]$' % zeta_a)
ax1.set_xlabel('$\\lambda [m]$')
ax1.set_ylabel('$\\left|\\hat{F}_{7}\\right|$ $[m^3/s]$')
ax1.set_xlim([x_min, x_max])

ax2.plot(wavelength, np.angle(f_wp), label='Conceptual SES')
ax2.set_xlabel('$\\lambda [m]$')
ax2.set_ylabel('$\\arg\\left(\\hat{F}_{7}\\right)$ $[rad]$')
ax2.set_xlim([x_min, x_max])
plt.show()
'''
