import numpy as np
import matplotlib.pyplot as plt
from veres import read_group_of_re7_input, read_group_of_re8_input
from Spatially_varying_pressure.HyCoef import compute_hydrodynamic_coeff

# Get result from ShipX
filepath = 'C:/Users/mathi/SIMA Workspaces/Workspace_1/Task/BBGreen/DWL'
M, A_temp, B_temp, C_temp, VEL, HEAD, FREQ_re7, XMTN_re7, ZMTN_re7, NDOF = read_group_of_re7_input(filepath)

A_33_shipx = A_temp[:, 2, 2]
B_33_shipx = B_temp[:, 2, 2]

A_35_shipx = A_temp[:, 2, 4]
B_35_shipx = B_temp[:, 2, 4]

A_53_shipx = A_temp[:, 4, 2]
B_53_shipx = B_temp[:, 4, 2]

A_55_shipx = A_temp[:, 4, 4]
B_55_shipx = B_temp[:, 4, 4]

# Formulas of Salvesen, Tuck and Faltinsen (1970)
L_c = 35.  # [m]
Bs = 0.5  # [m] Beam of the side hulls  # TODO: Check what this variable means in Salvesen, Tuck and Faltinsen (1970)

g = 9.81  # [m/s^2]  Acceleration of gravity
U = 50.  # [knots]  Vessel velocity
U = U*0.514444  # [m/s]  Vessel velocity

rho_w = 1025.  # [kg/m^3]  Density of water

omega_0 = FREQ_re7
k = np.power(omega_0, 2)/g  # wave number of water waves
omega_e = omega_0 + np.power(omega_0, 2)/g*U  # encounter frequencies
f_encounter = omega_e/2/np.pi  # [Hz] frequency of encounter
n_freq = len(omega_e)

A_33, A_35, A_53, A_55, B_33, B_35, B_53, B_55 = compute_hydrodynamic_coeff(L_c, Bs, U, omega_0)

# Simplified formula from Sea Loads
Beam = 0.5  # [m]
Draught = 0.7  # [m]

Beam_draught_ratio = Beam/2/Draught  # [-]

A_33_mass_ratio = 0.75*0.95/1.65  # [m]

A_33_2D_simplified = A_33_mass_ratio * rho_w * Beam * Draught
A_33_simplified = 2*A_33_2D_simplified * L_c


# Plot results
plt.plot(f_encounter, A_33_shipx, label='ShipX')
#plt.plot(f_encounter, np.ones([n_freq])*A_33, label='Salvesen, Tuck and Faltinsen (1970)')
#plt.plot(f_encounter, np.ones([n_freq])*A_33_simplified, label='Sea Loads p. 52')
plt.xlabel('Encounter frequency [Hz]')
plt.ylabel('$A_{33}$ [kg]')
plt.title('')
plt.legend()
plt.show()

plt.plot(f_encounter, A_55_shipx, label='ShipX')
#plt.plot(f_encounter, A_55, label='Salvesen, Tuck and Faltinsen (1970)')
plt.xlabel('Encounter frequency [Hz]')
plt.ylabel('$A_{55}$ [kg m^2]')
plt.title('')
plt.legend()
plt.show()

plt.plot(f_encounter, -A_35_shipx, label='ShipX')
#plt.plot(f_encounter, A_35, label='Salvesen, Tuck and Faltinsen (1970)')
plt.xlabel('Encounter frequency [Hz]')
plt.ylabel('$A_{35}$ [kg m]')
plt.title('')
plt.legend()
plt.show()

plt.plot(f_encounter, -A_53_shipx, label='ShipX')
#plt.plot(f_encounter, A_53, label='Salvesen, Tuck and Faltinsen (1970)')
plt.xlabel('Encounter frequency [Hz]')
plt.ylabel('$A_{53}$ [kg m]')
plt.title('')
plt.legend()
plt.show()

plt.plot(f_encounter, B_33_shipx, label='ShipX')
#plt.plot(f_encounter, np.ones([n_freq])*B_33, label='Salvesen, Tuck and Faltinsen (1970)')
plt.xlabel('Encounter frequency [Hz]')
plt.ylabel('$B_{33}$ [kg / s]')
plt.title('')
plt.legend()
plt.show()

plt.plot(f_encounter, -B_35_shipx, label='ShipX')
#plt.plot(f_encounter, B_35, label='Salvesen, Tuck and Faltinsen (1970)')
plt.xlabel('Encounter frequency [Hz]')
plt.ylabel('$B_{35}$ [kg m / s]')
plt.title('')
plt.legend()
plt.show()

plt.plot(f_encounter, -B_53_shipx, label='ShipX')
#plt.plot(f_encounter, np.ones([n_freq])*B_53, label='Salvesen, Tuck and Faltinsen (1970)')
plt.xlabel('Encounter frequency [Hz]')
plt.ylabel('$B_{53}$ [kg m / s]')
plt.title('')
plt.legend()
plt.show()

plt.plot(f_encounter, B_55_shipx, label='ShipX')
#plt.plot(f_encounter, B_55, label='Salvesen, Tuck and Faltinsen (1970)')
plt.xlabel('Encounter frequency [Hz]')
plt.ylabel('$B_{55}$ [kg m/s]')
plt.title('')
plt.legend()
plt.show()
