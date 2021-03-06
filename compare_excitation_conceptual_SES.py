import numpy as np
import pandas as pd
from veres import read_re8_file
from bow_and_stern_seals import excitation_skirts
from air_cushion import wave_pumping_rect, wave_pumping_rectangle

"""
This script compares excitation forces from Veres on a simplified SES geometry with rectangular side hulls. 

The center of gravity, and the motion coordinate system is located at the centroid of the water plane area at the 
mean free surface.
"""

# ------- Bow and stern seals ------
# The vessel has a finger seal at the bow and a lobe bag type of seal at the stern

# properties
L = 20  # [m] length of the vessel
b = 7  # [m] beam of air cushion
b_seals = 7  # [m] Beam of the seals
tau_b = 60  # [deg] angle of the bow finger seal
tau_s = 30  # [deg] angle of the stern lobe bag seal
p_s = 0  # [Pa] Membrane seal pressure
x_lobe_bag_seal = L / 2  # [m] Longitudinal position of the lobe bag seal at the stern relative to motion coord. system
x_finger_seal = -L / 2  # [m] Longitudinal position of the finger seal at the bow relative to motion coord. system
p_0 = 3500  # [Pa] excess pressure in the air cushion

# ------- Comparison with Veres ------

path_re8_air_cushion = 'Input files/Conceptual SES/with air cushion no skirts/input.re8'
path_re8_no_air_cushion = 'Input files/Conceptual SES/without air cushion/input.re8'

REFORCE_na, IMFORCE_na, VEL_na, HEAD_na, FREQ_na, XMTN_na, ZMTN_na = read_re8_file(path_re8_no_air_cushion)
REFORCE_a, IMFORCE_a, VEL_a, HEAD_a, FREQ_a, XMTN_a, ZMTN_a = read_re8_file(path_re8_air_cushion)


# compute excitation loads on skirts

# excitation_skirts(b, tau_b, tau_s, p_0, p_s, x_b, x_s, omegas, beta, zeta_a=1, g=9.81):
f_3_skirts, f_5_skirts = excitation_skirts(b_seals, tau_b, tau_s, p_0, p_s, x_finger_seal, x_lobe_bag_seal, FREQ_na, HEAD_na[0])


# Append excitation loads due to skirts into the excitation terms without any excitation loads due to the skirts.

n_frequencies = len(FREQ_na)

# Initialize arrays to hold the force components
f_ex_re_without_air = np.zeros([n_frequencies, 6])
f_ex_im_without_air = np.zeros([n_frequencies, 6])
f_ex_re_with_air = np.zeros([n_frequencies, 6])
f_ex_im_with_air = np.zeros([n_frequencies, 6])

# Copy values from input files with and without an air cushion present
for i in range(n_frequencies):
    for j in range(6):
        f_ex_re_without_air[i, j] = REFORCE_na[j, i, 0, 0]
        f_ex_im_without_air[i, j] = IMFORCE_na[j, i, 0, 0]
        f_ex_re_with_air[i, j] = REFORCE_a[j, i, 0, 0]
        f_ex_im_with_air[i, j] = IMFORCE_a[j, i, 0, 0]

f_ex_with_air = f_ex_re_with_air + 1j * f_ex_im_with_air
f_ex_without_air = f_ex_re_without_air + 1j * f_ex_im_without_air

# Append excitation from skirts

f_3_python = np.zeros([n_frequencies, 1], dtype=complex)
f_5_python = np.zeros([n_frequencies, 1], dtype=complex)

for i in range(n_frequencies):
    f_3_python[i] = f_ex_re_without_air[i, 2] + 1j * f_ex_im_without_air[i, 2] + f_3_skirts[i, 0]

    f_5_python[i] = f_ex_re_without_air[i, 4] + 1j * f_ex_im_without_air[i, 4] + f_5_skirts[i, 0]

    '''
    f_ex_re_without_air[i, 2] += f_3_skirts[i].real
    f_ex_im_without_air[i, 2] += f_3_skirts[i].imag
    
    f_ex_re_without_air[i, 4] += f_5_skirts[i].real
    f_ex_im_without_air[i, 4] += f_5_skirts[i].imag
    '''

print('Values from Python calculation:')
print('Absolute value: ', np.abs(f_3_skirts[-1][0]))
print('Phase: ', np.rad2deg(np.angle(f_3_skirts[-1][0])))
print()
print('Values from Veres calculation:')
print('Absolute value: ', np.abs(f_ex_with_air[-1][2]))
print('Phase', np.rad2deg(np.angle(f_ex_with_air[-1][2])))
print()
print('Relative error: ', 100*np.abs((np.abs(f_3_skirts[-1][0]) - np.abs(f_ex_with_air[-1][2])) / np.abs(f_ex_with_air[-1][2])), '%')

print(f_3_python[:, 0])
print()
print(f_5_python[:, 0])
print()
print(f_ex_with_air[:, 2])
print()
print(f_ex_with_air[:, 4])


dat = np.array([f_ex_with_air[:, 2], f_3_python[:, 0], f_ex_with_air[:, 4], f_5_python[:, 0]]).transpose()

print(dat)

df = pd.DataFrame(dat, index=FREQ_a, columns=['f_3_veres', 'f_3_python', 'f_5_veres', 'f_5_python'])

# Compute relative errors between Veres and Python calculations.
df['f_3_rel_err'] = np.abs((df['f_3_veres'] - df['f_3_python']) / df['f_3_veres'])
df['f_5_rel_err'] = np.abs((df['f_5_veres'] - df['f_5_python']) / df['f_5_veres'])


# Compute wave pumping excitation
x_s = L/2
x_f = -L/2
y_s = b/2
y_p = -b/2

f_ex_7 = wave_pumping_rect(x_f, x_s, y_p, y_s, FREQ_a, HEAD_a[0])

a = f_ex_with_air.transpose()

f_ex = np.vstack([f_ex_with_air.transpose(), f_ex_7])

# Expressions for comparisons
zeta_a = 1
g = 9.81
beta = HEAD_na[0]
omegas = FREQ_a
k = np.power(omegas, 2) / g
theta = k * x_s * np.cos(beta)
a = 1j * zeta_a * (y_s - y_p) / np.cos(beta) * np.divide(omegas, k)
f_wp_comparison = np.multiply(a, (np.exp(-1j * theta) - np.exp(1j * theta)))

comp = f_wp_comparison - f_ex_7

print(df)
