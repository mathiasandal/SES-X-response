from veres import read_re8_file

"""
This script compares excitation forces from Veres on a simplified SES geometry with rectangular side hulls. 

The center of gravity, and the motion coordinate system is located at the centroid of the water plane area at the 
mean free surface.
"""

# ------- Bow and stern seals ------
# The vessel has a finger seal at the bow and a lobe bag type of seal at the stern

# properties
L = 20  # [m] length of the vessel
b_seals = 7  # [m] Beam of the seals
tau_b = 60  # [deg] angle of the bow finger seal
tau_s = 30  # [deg] angle of the stern lobe bag seal
p_s = 0  # [Pa] Membrane seal pressure
x_lobe_bag_seal = L / 2  # [m] Longitudinal position of the lobe bag seal at the stern relative to motion coord. system
x_finger_seal = -L / 2  # [m] Longitudinal position of the finger seal at the bow relative to motion coord. system


# ------- Comparison with Veres ------

path_re8_air_cushion = 'Input files/Conceptual SES/with air cushion/input.re8'
path_re8_no_air_cushion = 'Input files/Conceptual SES/without air cushion/input.re8'

REFORCE_na, IMFORCE_na, VEL_na, HEAD_na, FREQ_na, XMTN_na, ZMTN_na = read_re8_file(path_re8_no_air_cushion)
REFORCE_a, IMFORCE_a, VEL_a, HEAD_a, FREQ_a, XMTN_a, ZMTN_a = read_re8_file(path_re8_air_cushion)

print('hi')



