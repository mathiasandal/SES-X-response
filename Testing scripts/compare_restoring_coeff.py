import numpy as np
import matplotlib.pyplot as plt
from veres import read_group_of_re7_input
from Spatially_varying_pressure.HyCoef import compute_hydrostatic_coeff


# Get result from ShipX
filepath = 'C:/Users/mathi/SIMA Workspaces/Workspace_1/Task/conceptual_35m_fine_contour_strip_theory/DWL'
M, A_temp, B_temp, C_temp, VEL, HEAD, FREQ_re7, XMTN_re7, ZMTN_re7, NDOF = read_group_of_re7_input(filepath)

n_freq = len(FREQ_re7)

C_33_shipx = C_temp[n_freq//2, 2, 2]

C_35_shipx = C_temp[n_freq//2, 2, 4]

C_53_shipx = C_temp[n_freq//2, 4, 2]

C_55_shipx = C_temp[n_freq//2, 4, 4]

# Formulas of Salvesen, Tuck and Faltinsen (1970)
L_c = 35.  # [m]
L_oa = 35.  # [m]
Bs = 0.5  # [m] Beam of the side hulls  # TODO: Check what this variable means in Salvesen, Tuck and Faltinsen (1970)

m = 140*1e3  # [kg]
g = 9.81  # [m/s^2]  Acceleration of gravity
U = 50.  # [knots]  Vessel velocity
U = U*0.514444  # [m/s]  Vessel velocity

rho_w = 1025.  # [kg/m^3]  Density of water

C_33, C_35, C_53, C_55 = compute_hydrostatic_coeff(m, L_oa, L_c, Bs)

print(C_33)

print('hei')