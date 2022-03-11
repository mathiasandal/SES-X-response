"""
In this script RAOs are plotted for a conceptual SES with and without an air cushion
"""
import numpy as np
import matplotlib.pyplot as plt
from veres import read_re1_file

g = 9.81  # [m^2/s]

# paths
path_air_cushion = 'Input files/Conceptual SES/with air cushion no skirts/input.re1'
path_no_air_cushion = 'Input files/Conceptual SES/without air cushion/input.re1'

RETRANS_a, IMTRANS_a, VEL_a, HEAD_a, FREQ_a, XMTN_a, ZMTN_a = read_re1_file(path_air_cushion)

RETRANS_na, IMTRANS_na, VEL_na, HEAD_na, FREQ_na, XMTN_na, ZMTN_na = read_re1_file(path_no_air_cushion)

# initialize array to contain complex force amplitudes
beta = HEAD_na[0]
frequencies = FREQ_na
frequency_Hz = frequencies / 2 / np.pi  # compute encounter frequency in Hz
k = np.power(frequencies, 2) / g  # encounter wave number
wavelength = np.divide(2*np.pi, k)
n_frequencies = len(FREQ_a)


raos_air = np.transpose(RETRANS_a[:, :, 0, 0] + 1j * IMTRANS_a[:, :, 0, 0])
raos_no_air = np.transpose(RETRANS_na[:, :, 0, 0] + 1j * IMTRANS_na[:, :, 0, 0])

# ***** Plot RAOs *****

x_min = 0
x_max = 50

# Define plot colors
color_BBGreen = '#5cb16d'
color_BBPurple = '#b15ca0'
plt.rc('font', size=15)  # sets size

# Plot amplitude against wavelength
plot_wavelength = True
if plot_wavelength:
    # Heave response
    plt.plot(wavelength, np.abs(raos_air[:, 2]), fillstyle='none', marker='o', linestyle='-', label='Air cushion', color=color_BBGreen)
    plt.plot(wavelength, np.abs(raos_no_air[:, 2]), marker='x', linestyle='-', label='No air cushion', color=color_BBPurple)
    plt.title('Heave RAO, ' + '$\\beta = %0.1f^{\\circ}$' % beta)
    plt.xlabel('$\\lambda [m]$')
    plt.ylabel('$\\left|\\eta_{3}\\right|/\\zeta_a$ $[-]$')
    plt.legend(loc='best')
    #plt.xlim([x_min, x_max])
    #plt.ylim([np.min(np.abs(raos_air[:, 2])[(wavelength > x_min) & (wavelength < x_max)]),
    #          np.max(np.abs(raos_air[:, 2])[(wavelength > x_min) & (wavelength < x_max)]) + 1])

    save_plot = True
    if save_plot:
        plot_path = 'Results/Conceptual SES/RAOs/'
        plot_name = 'heave,wavelength,beta=%0.1fdeg' % beta + '.pdf'
        plt.tight_layout()
        plt.savefig(plot_path + plot_name)
    plt.show()

    # Heave acceleration
    plt.plot(wavelength, np.multiply(np.power(frequencies, 2), np.abs(raos_air[:, 2])), fillstyle='none', marker='o', linestyle='-', label='Air cushion',
             color=color_BBGreen)
    plt.plot(wavelength, np.multiply(np.power(frequencies, 2), np.abs(raos_no_air[:, 2])), marker='x', linestyle='-', label='No air cushion',
             color=color_BBPurple)
    plt.title('Heave RAO acceleration, ' + '$\\beta = %0.1f^{\\circ}$' % beta)
    plt.xlabel('$\\lambda [m]$')
    plt.ylabel('$\\left|\\omega^2 \\eta_{3}\\right|/\\zeta_a$ $[rad^2/s^2]$')
    plt.legend(loc='best')
    #plt.xlim([x_min, x_max])
    #plt.ylim([np.min(np.abs(raos_air[:, 2])[(wavelength > x_min) & (wavelength < x_max)]),
    #          np.max(np.abs(raos_air[:, 2])[(wavelength > x_min) & (wavelength < x_max)]) + 1])

    save_plot = True
    if save_plot:
        plot_path = 'Results/Conceptual SES/RAOs/'
        plot_name = 'heave_acc,wavelength,beta=%0.1fdeg' % beta + '.pdf'
        plt.tight_layout()
        plt.savefig(plot_path + plot_name)
    plt.show()


# Plot amplitude against frequency
plot_frequencies = True
if plot_wavelength:
    # Heave response
    plt.plot(frequencies, np.abs(raos_air[:, 2]), fillstyle='none', marker='o', linestyle='-', label='Air cushion', color=color_BBGreen)
    plt.plot(frequencies, np.abs(raos_no_air[:, 2]), marker='x', linestyle='-', label='No air cushion', color=color_BBPurple)
    plt.title('Heave RAO, ' + '$\\beta = %0.1f^{\\circ}$' % beta)
    plt.xlabel('wave frequencies $[rad/s]$')
    plt.ylabel('$\\left|\\eta_{3}\\right|/\\zeta_a$ $[-]$')
    plt.legend(loc='best')
    #plt.xlim([x_min, x_max])
    # plt.ylim([np.min(np.abs(raos_air[:, 2])[(frequencies > x_min) & (frequencies < x_max)]),
    #          np.max(np.abs(raos_air[:, 2])[(frequencies > x_min) & (frequencies < x_max)]) + 1])

    save_plot = True
    if save_plot:
        plot_path = 'Results/Conceptual SES/RAOs/'
        plot_name = 'heave,frequencies,beta=%0.1fdeg' % beta + '.pdf'
        plt.tight_layout()
        plt.savefig(plot_path + plot_name)
    plt.show()

    # heave acceleration
    plt.plot(frequencies, np.multiply(np.power(frequencies, 2), np.abs(raos_air[:, 2])), fillstyle='none', marker='o', linestyle='-', label='Air cushion', color=color_BBGreen)
    plt.plot(frequencies, np.multiply(np.power(frequencies, 2), np.abs(raos_no_air[:, 2])), marker='x', linestyle='-', label='No air cushion', color=color_BBPurple)
    plt.title('Heave RAO acceleration, ' + '$\\beta = %0.1f^{\\circ}$' % beta)
    plt.xlabel('wave frequencies $[rad/s]$')
    plt.ylabel('$\\left|\\omega^2 \\eta_{3}\\right|/\\zeta_a$ $[rad^2/s^2]$')
    plt.legend(loc='best')
    #plt.xlim([x_min, x_max])
    #plt.ylim([np.min(np.abs(raos_air[:, 2])[(frequencies > x_min) & (frequencies < x_max)]),
    #          np.max(np.abs(raos_air[:, 2])[(frequencies > x_min) & (frequencies < x_max)]) + 1])

    save_plot = True
    if save_plot:
        plot_path = 'Results/Conceptual SES/RAOs/'
        plot_name = 'heave_acc,frequencies,beta=%0.1fdeg' % beta + '.pdf'
        plt.tight_layout()
        plt.savefig(plot_path + plot_name)
    plt.show()
