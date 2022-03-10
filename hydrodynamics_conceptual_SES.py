from air_cushion import read_fan_characteristics, interpolate_fan_characteristics
from veres import read_re7_file
from bow_and_stern_seals import restoring_finger_at_bow_lobe_bag_at_the_stern
from Wave_response_utilities import add_row_and_column
import numpy as np

"""
This script compares calculations from Veres on a simplified SES geometry with rectangular side hulls. 

The center of gravity, and the motion coordinate system is located at the centroid of the water plane area at the 
mean free surface.
"""

# ------ Parameters -------
# Main parameters of conceptual SES with rectangular side hulls
L = 20  # [m] length of the vessel
b_G = 0.5  # [m] beam of side hulls (gondolas)
b_c = 7  # [m] beam of the air cushion/length between the side hulls
h_hull = 2  # [m] distance from bottom of side hull to roof of air cushion
d = 0.5  # [m] draught with excess pressure in the air cushion
p_0 = 3500  # [Pa] excess pressure in the air cushion
B = b_c + 2 * b_G  # [m] total beam of the vessel

# ------- Physical parameters -------
rho = 1025  # [kg/m^3] density of sea water
p_a = 101325  # [Pa] atmospheric pressure
g = 9.81  # [m/s^2] acceleration of gravity
gamma = 1.4  # [-] Ratio of specific heat for air, when assuming adiabatic density relation

# Some derived parameters
A_wp = 2 * b_G * L  # [m^2] water plane area of the side hulls
I_x = 2 * (b_G ** 3 * L / 12 + (b_c / 2 + b_G / 2) ** 2 * b_G * L)  # [m^4] second moment of area of A_wp about x-axis.
I_y = 2 * (L ** 3 * b_G / 12)  # [m^4] second moment of area of A_wp about y-axis
V_0 = A_wp * d  # [m^3] volume displaced by the side hulls at equilibrium

# ------- Air cushion properties -------
h = p_0 / rho / g  # [m] Difference in water level inside and outside the air cushion
h_b1 = h_hull + h - d  # [m] height of the air cushion after the water level inside the cushion i suppressed.
h_b = h_hull - d  # [m] height of the air cushion before the water level inside the cushion i suppressed.
A_b = L * b_c  # [m] Area of the air cushion
x_B_c = 0  # [m] longitudinal position of the air cushion centroid relative to the L_pp/2 in front of AP
V_c0 = A_b * h_b  # [m^3] Volume inside air cushion at equilibrium

# Read in fan characteristics
Q, P, rpm_dummy = read_fan_characteristics('Input files/fan characteristics/fan characteristics.csv', '1800rpm')

# Interpolates values
Q_0, dQdp_0 = interpolate_fan_characteristics(p_0, P, Q)

# ------- Calculate mass of the vessel -------
M = V_0 * rho + p_0 * A_b / g  # [kg] total mass of the vessel

# ------- Coordinate systems -------
# Motion coordinate system located relative to the intersection of AP, CL and BL.
x_m = L / 2  # [m] longitudinal position
y_m = 0  # [m] transverse position
z_m = d  # [m] vertical position

# Center of gravity (CoG) located at (x_G, y_G, z_G)
x_G = 0  # [m] longitudinal position of the center of gravity relative to the motion coord. system
y_G = 0  # [m] transverse position of the center of gravity relative to the motion coord. system
z_G = 0  # [m] vertical position of the center of gravity relative to the motion coord. system

# Centroid of the air cushion at equilibrium relative to the motion coord. system
x_c = 0  # [m] longitudinal position
y_c = 0  # [m] transverse position
z_c = -h / 2  # [m] vertical position

# Centroid of the submerged hull volume relative to the motion coord. system
x_B = 0  # [m] longitudinal position
y_B = 0  # [m] transverse position
z_B = -d / 2  # [m] vertical position

# ------- Bow and stern seals ------
# The vessel has a finger seal at the bow and a lobe bag type of seal at the stern

# properties
b_seals = b_c  # [m] Beam of the seals
tau_b = 89  # [deg] angle of the bow finger seal
tau_s = 89  # [deg] angle of the stern lobe bag seal
p_s = 0  # [Pa] Membrane seal pressure
x_lobe_bag_seal = L / 2  # [m] Longitudinal position of the lobe bag seal at the stern relative to motion coord. system
x_finger_seal = -L / 2  # [m] Longitudinal position of the finger seal at the bow relative to motion coord. system

# Computes restoring coefficients
C_33_seal, C_35_seal, C_55_seal = restoring_finger_at_bow_lobe_bag_at_the_stern(b_seals, tau_b, tau_s, 0, p_s,
                                                                                x_finger_seal, x_lobe_bag_seal)

# ------- Hydrodynamic coefficients -------

# Restoring
C_33_h = rho * g * A_wp  # [N/m] hydrostatic restoring coefficient in heave
C_44_h = rho * g * V_0 * z_B + rho * g * I_x - M * g * z_G  # [Nm] hydrostatic restoring coefficient in roll
C_55_h = rho * g * V_0 * z_B + rho * g * I_y - M * g * z_G  # [Nm] hydrostatic restoring coefficient in pitch
C_53_h = -rho * g * A_wp * x_B  # [N] hydrostatic restoring coefficient in pitch due to heave motion
C_35_h = C_53_h  # [m] hydrostatic restoring in heave due to pitch motion

C_44_c = rho * g * h * A_b * z_c  # [Nm] restoring coefficient in roll due to the air cushion
C_55_c = rho * g * h * A_b * z_c  # [Nm] restoring coefficient in pitch due to the air cushion
C_57_c = -rho * g * h * A_b * x_c  # [Nm] coupling term in pitch due to change in the air cushion pressure
C_37_c = -rho * g * h * A_b  # [N] coupling term in heave due to change in the air cushion pressure

C_77_c = 0.5 * Q_0 - p_0 * dQdp_0  # [m^3/s] equivalent restoring term in mass continuity eq. for air inside air cushion

C_33 = C_33_h + C_33_seal  # [N/m] total restoring coefficient in heave
C_44 = C_44_h + C_44_c  # [Nm] total restoring coefficient in roll
C_55 = C_55_h + C_55_c + C_55_seal  # [Nm] total restoring coefficient in pitch
C_35 = C_35_seal + C_35_h  # [N] total restoring coefficient in heave due to pitch motion
C_53 = C_35_seal + C_53_h  # [N] Total restoring coefficient in pitch due to heave motion

# Damping
B_73_c = A_b  # [m] equivalent damping coefficient in uniform pressure DOF due to velocity in heave.
B_75_c = -A_b * x_c  # [m^2] equivalent damping coefficient in uniform pressure DOF due to velocity in pitch.
B_77_c = p_0 * V_c0 / gamma / (p_0 + p_a)  # [m^3] equivalent damping coefficient in uniform pressure DOF.

# ------- Establish stiffness and damping matrix calculated in python -------
C_python = np.array([[0,     0,     0,    0,    0,   0,      0],
                     [0,     0,     0,    0,    0,   0,      0],
                     [0,     0,  C_33,    0, C_35,   0, C_37_c],
                     [0,     0,     0, C_44,    0,   0,      0],
                     [0,     0,  C_53,    0, C_55,   0, C_57_c],
                     [0,     0,     0,    0,    0,   0,      0],
                     [0,     0,     0,    0,    0,   0, C_77_c]])

B_python = np.array([[0,     0,      0,    0,      0,   0,      0],
                     [0,     0,      0,    0,      0,   0,      0],
                     [0,     0,      0,    0,      0,   0,      0],
                     [0,     0,      0,    0,      0,   0,      0],
                     [0,     0,      0,    0,      0,   0,      0],
                     [0,     0,      0,    0,      0,   0,      0],
                     [0,     0, B_73_c,    0, B_75_c,   0, B_77_c]])


# ------- Comparison with Veres ------

path_re7_air_cushion = 'Input files/Conceptual SES/with air cushion no skirts/input.re7'
path_re7_no_air_cushion = 'Input files/Conceptual SES/without air cushion/input.re7'

# Get Veres results with air cushion turned on
M_veres_air_cushion, A_veres_air_cushion, B_veres_air_cushion, C_veres_air_cushion, VEL_a, HEAD_a, FREQ_a, XMTN_a, \
    ZMTN_a, NDOF_a = read_re7_file(path_re7_air_cushion)

# Get Veres results with air cushion turned off
M_veres_no_air_cushion, A_veres_no_air_cushion, B_veres_no_air_cushion, C_veres_no_air_cushion, VEL_na, HEAD_na, \
    FREQ_na, XMTN_na, ZMTN_na, NDOF_na = read_re7_file(path_re7_no_air_cushion)

# Combine Veres without air cushion with air cushion terms

# initialize stiffness matrix to be added air cushion terms and adds a row and a column such that it is a 7x7 matrix
C_veres_plus_cushion = add_row_and_column(C_veres_no_air_cushion[0][0][0])

# Add air cushion terms
C_veres_plus_cushion[2, 2] += C_33_seal  # heave
C_veres_plus_cushion[2, 4] += C_35_seal  # heave due to pitch
C_veres_plus_cushion[4, 2] += C_35_seal  # pitch due to heave

C_veres_plus_cushion[3, 3] += C_44_c  # roll
C_veres_plus_cushion[4, 4] += C_55_c + C_55_seal  # pitch
C_veres_plus_cushion[2, 6] += C_37_c  # heave due to uniform pressure
C_veres_plus_cushion[4, 6] += C_57_c  # pitch due to uniform pressure
C_veres_plus_cushion[6, 6] += C_77_c  # Uniform pressure DOF

# computes relative errors
C_rel_error = np.zeros([7, 7])  # initialize array

C_rel_error[2, 2] = np.abs((C_veres_plus_cushion[2, 2] - C_veres_air_cushion[0][0][0][2][2])/C_veres_air_cushion[0][0][0][2][2])  # heave
C_rel_error[2, 4] = np.abs((C_veres_plus_cushion[2, 4] - C_veres_air_cushion[0][0][0][2][4])/C_veres_air_cushion[0][0][0][2][4])  # heave due to pitch
C_rel_error[4, 2] = np.abs((C_veres_plus_cushion[4, 2] - C_veres_air_cushion[0][0][0][4][2])/C_veres_air_cushion[0][0][0][4][2])  # pitch due to heave
C_rel_error[3, 3] = np.abs((C_veres_plus_cushion[3, 3] - C_veres_air_cushion[0][0][0][3][3])/C_veres_air_cushion[0][0][0][3][3])  # roll
C_rel_error[4, 4] = np.abs((C_veres_plus_cushion[4, 4] - C_veres_air_cushion[0][0][0][4][4])/C_veres_air_cushion[0][0][0][4][4])  # pitch
C_rel_error[2, 6] = np.abs((C_veres_plus_cushion[2, 6] - C_veres_air_cushion[0][0][0][2][6])/C_veres_air_cushion[0][0][0][2][6])  # heave due to uniform pressure
C_rel_error[4, 6] = np.abs((C_veres_plus_cushion[4, 6] - C_veres_air_cushion[0][0][0][4][6])/C_veres_air_cushion[0][0][0][4][6])  # pitch due to uniform pressure
C_rel_error[6, 6] = np.abs((C_veres_plus_cushion[6, 6] - C_veres_air_cushion[0][0][0][6][6])/C_veres_air_cushion[0][0][0][6][6])  # uniform pressure DOF

print('hi')
