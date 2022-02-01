from air_cushion import read_fan_characteristics, interpolate_fan_characteristics
from veres import read_re7_file
from bow_and_stern_seals import restoring_finger_at_bow_lobe_bag_at_the_stern

"""
This script compares calculations from Veres on a simplified SES geometry with rectangular side hulls. 

The center of gravity, and the motion coordinate system is located at the centroid of the water plane area at the 
mean free surface.
"""

# ------ Parameters -------
# Main parameters of conceptual SES with rectangular side hulls
L = 20  # [m] length of the vessel
B = 10  # [m] total beam of the vessel
b_G = 0.5  # [m] beam of side hulls (gondolas)
b_c = 7  # [m] beam of the air cushion/beam between the side hulls
h_hull = 2  # [m] distance from bottom of side hull to roof of air cushion
d = 0.5  # [m] draught with excess pressure in the air cushion
p_0 = 3500  # [Pa] excess pressure in the air cushion

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
# The vessel has a finger seal at the bow and a lobe bag type of seal at the back

# properties
b_seals = b_c  # [m] Beam of the seals
tau_b = 60  # [deg] angle of the bow finger seal
tau_s = 30  # [deg] angle of the stern lobe bag seal
p_s = 0  # [Pa] Membrane seal pressure
x_lobe_bag_seal = -L / 2  # [m] Longitudinal position of the lobe bag seal relative to motion coord. system
x_finger_seal = L / 2  # [m] Longitudinal position of the finger seal relative to motion coord. system

# Computes restoring coefficients
c_33_seal, c_35_seal, c_55_seal = restoring_finger_at_bow_lobe_bag_at_the_stern(b_seals, tau_b, tau_s, p_0, p_s,
                                                                                x_finger_seal, x_lobe_bag_seal)

# ------- Hydrodynamic coefficients -------

# Restoring
C_33_h = rho * g * A_wp  # [N/m] hydrostatic restoring coefficient in heave
C_44_h = rho * g * V_0 * z_B + rho * g * I_x - M * g * z_G  # [Nm] hydrostatic restoring coefficient in roll
C_55_h = rho * g * V_0 * z_B + rho * g * I_y - M * g * z_G  # [Nm] hydrostatic restoring coefficient in pitch

C_44_c = rho * g * h * A_b * z_c  # [Nm] restoring coefficient in roll due to the air cushion
C_55_c = rho * g * h * A_b * z_c  # [Nm] restoring coefficient in pitch due to the air cushion
C_57_c = -rho * g * h * A_b * x_B  # [Nm] coupling term in pitch due to change in the air cushion pressure
C_37_c = -rho * g * h * A_b  # [N] coupling term in heave due to change in the air cushion pressure

C_77_c = 0.5 * Q_0 - p_0 * dQdp_0  # [m^3/s] equivalent restoring term in mass continuity eq. for air inside air cushion

C_44 = C_44_h + C_44_c  # [Nm] total restoring coefficient
C_55 = C_55_h + C_55_c  # [Nm] total restoring coefficient

# Damping
B_73_c = A_b  # [m] equivalent damping coefficient in uniform pressure DOF due to velocity in heave.
B_75_c = -A_b * x_B  # [m^2] equivalent damping coefficient in uniform pressure DOF due to velocity in pitch.
B_77_c = p_0 * V_c0 / gamma / (p_0 + p_a)  # [m^3] equivalent damping coefficient in uniform pressure DOF.

# ------- Comparison with Veres ------

path_re7 = 'Input files/Conceptual SES/with air cushion/input.re7'

VMAS, ADDMAS, DAMP, REST, VEL, HEAD, FREQ, XMTN, ZMTN, NDOF = read_re7_file(path_re7)

print('hi')
