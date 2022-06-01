import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
from veres import read_group_of_re7_input

g = 9.81

M_BBGreen = 26*1e3
L_pp_BBGreen = 19.2
r_55_BBGreen = 0.25 * L_pp_BBGreen
I_55_BBGreen = M_BBGreen * r_55_BBGreen**2

M_35m = 140*1e3
I_55_35m = 2860000.

M_20m = 10.3  # 60.2*1e3 #
L_pp_20m = 20
r_55_20m = 1.
I_55_20m = M_20m * r_55_20m**2


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


# read restoring and added mass coefficients
# BBGreen at 22kn
path_BBGreen = 'C:/Users/mathi/SIMA Workspaces/Workspace_1/Task/BBGreen/DWL'
M_temp, A_BBGreen, B_BBGreen, C_BBGreen, VEL_BBGreen, HEAD_BBGreen, FREQ_BBGreen, XMTN_BBGreen, ZMTN_BBGreen, NDOF_BBGreen = read_group_of_re7_input(path_BBGreen)
omega_enc_BBGreen = FREQ_BBGreen + VEL_BBGreen / g * np.power(FREQ_BBGreen, 2)
f_enc_BBGreen = omega_enc_BBGreen / 2 / np.pi

# Conceptual 35m at 50kn
path_35m = 'C:/Users/mathi/SIMA Workspaces/Workspace_1/Task/conceptual_35m_fine_contour_strip_theory/DWL'
M_temp, A_35m, B_35m, C_35m, VEL_35m, HEAD_35m, FREQ_35m, XMTN_35m, ZMTN_35m, NDOF_35m = read_group_of_re7_input(path_35m)
omega_enc_35m = FREQ_35m + VEL_35m / g * np.power(FREQ_35m, 2)
f_enc_35m = omega_enc_35m / 2 / np.pi

# Conceptual 20m at 0kn
path_20m = 'C:/Users/mathi/SIMA Workspaces/Workspace_1/Task/conceptual_SES_20m_no_air_cushion/DWL'
M_temp, A_20m, B_20m, C_20m, VEL_20m, HEAD_20m, FREQ_20m, XMTN_20m, ZMTN_20m, NDOF_20m = read_group_of_re7_input(path_20m)
omega_enc_20m = FREQ_20m + VEL_20m / g * np.power(FREQ_20m, 2)
f_enc_20m = omega_enc_20m / 2 / np.pi

# compute natural frequencies in heave and pitch

# BBGreen
# heave - inf added mass
omega_03_BBGreen = np.sqrt(C_BBGreen[len(f_enc_BBGreen)//2, 2, 2] / (M_BBGreen + A_BBGreen[-1, 2, 2]))
f_03_BBGreen = omega_03_BBGreen / 2 / np.pi

# heave - iterated
omega_03_BBGreen_it = it_nat_freq(M_BBGreen, A_BBGreen[:, 2, 2], C_BBGreen[len(f_enc_BBGreen)//2, 2, 2], omega_enc_BBGreen)
f_03_BBGreen_it = omega_03_BBGreen_it / 2 / np.pi

# pitch - inf added mass
omega_05_BBGreen = np.sqrt(C_BBGreen[len(f_enc_BBGreen)//2, 4, 4] / (I_55_BBGreen + A_BBGreen[-1, 4, 4]))
f_05_BBGreen = omega_05_BBGreen / 2 / np.pi

# pitch - iterated
omega_05_BBGreen_it = it_nat_freq(I_55_BBGreen, A_BBGreen[:, 4, 4], C_BBGreen[len(f_enc_BBGreen)//2, 4, 4], omega_enc_BBGreen)
f_05_BBGreen_it = omega_05_BBGreen_it / 2 / np.pi

# Conceptual 35m
# heave - inf added mass
omega_03_35m = np.sqrt(C_35m[len(f_enc_35m)//2, 2, 2] / (M_35m + A_35m[-1, 2, 2]))
f_03_35m = omega_03_35m / 2 / np.pi

# heave - iterated
omega_03_35m_it = it_nat_freq(M_35m, A_35m[:, 2, 2], C_35m[len(f_enc_35m)//2, 2, 2], omega_enc_35m)
f_03_35m_it = omega_03_35m_it / 2 / np.pi

# pitch - inf added mass
omega_05_35m = np.sqrt(C_35m[len(f_enc_35m)//2, 4, 4] / (I_55_35m + A_35m[-1, 4, 4]))
f_05_35m = omega_05_35m / 2 / np.pi

# pitch - iterated
omega_05_35m_it = it_nat_freq(I_55_35m, A_35m[:, 4, 4], C_35m[len(f_enc_35m)//2, 4, 4], omega_enc_35m)
f_05_35m_it = omega_05_35m_it / 2 / np.pi

# Conceptual 20m
# heave - inf added mass
omega_03_20m = np.sqrt(C_20m[len(f_enc_20m)//2, 2, 2] / (M_20m + A_20m[-1, 2, 2]))
f_03_20m = omega_03_20m / 2 / np.pi

# heave - iterated
omega_03_20m_it = it_nat_freq(M_20m, A_20m[:, 2, 2], C_20m[len(f_enc_20m)//2, 2, 2], omega_enc_20m)
f_03_20m_it = omega_03_20m_it / 2 / np.pi

# pitch - inf added mass
omega_05_20m = np.sqrt(C_20m[len(f_enc_20m)//2, 4, 4] / (I_55_20m + A_20m[-1, 4, 4]))
f_05_20m = omega_05_20m / 2 / np.pi

# pitch - iterated
omega_05_20m_it = it_nat_freq(I_55_20m, A_20m[:, 4, 4], C_20m[len(f_enc_20m)//2, 4, 4], omega_enc_20m)
f_05_20m_it = omega_05_20m_it / 2 / np.pi

# Plot A_33

# Define plot colors
color_BBGreen = '#5cb16d'
color_BBPurple = '#b15ca0'

plt.plot(f_enc_BBGreen, A_BBGreen[:, 2, 2], label='BBGreen', color=color_BBGreen)
plt.plot(f_enc_20m, A_20m[:, 2, 2], label='Conceptual 20m', color=color_BBPurple)
plt.plot(f_enc_35m, A_35m[:, 2, 2], label='Conceptual 35m', color='k')
plt.xlabel(r'$\textrm{Encounter frequency}\,[Hz]$')
plt.ylabel(r'$A_{33}$')
plt.legend()
plt.show()

print('BBGreen')
print('Infinite added mass')
print('f_03 =', str(round(f_03_BBGreen, 3)), '[Hz]')
print('f_05 =', str(round(f_05_BBGreen, 3)), '[Hz]')
print('Iterated')
print('f_03 =', str(round(f_03_BBGreen_it, 3)), '[Hz]')
print('f_05 =', str(round(f_05_BBGreen_it, 3)), '[Hz]')

print('Conceptual 35m')
print('Infinite added mass')
print('f_03 =', str(round(f_03_35m, 3)), '[Hz]')
print('f_05 =', str(round(f_05_35m, 3)), '[Hz]')
print('Iterated')
print('f_03 =', str(round(f_03_35m_it, 3)), '[Hz]')
print('f_05 =', str(round(f_05_35m_it, 3)), '[Hz]')

print('Conceptual 20m')
print('Infinite added mass')
print('f_03 =', str(round(f_03_20m, 3)), '[Hz]')
print('f_05 =', str(round(f_05_20m, 3)), '[Hz]')
print('Iterated')
print('f_03 =', str(round(f_03_20m_it, 3)), '[Hz]')
print('f_05 =', str(round(f_05_20m_it, 3)), '[Hz]')

