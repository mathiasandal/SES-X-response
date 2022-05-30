import numpy as np

U = 0.  # [knots] Vessel velocity
U = U * 0.514444  # [m/s] Vessel velocity
g = 9.81  # [m/s^2] acceleration of gravity
number_of_runs = 10


def f_e_to_T_0(f_e, U, g=9.81):

    if U == 0:
        return 1/f_e
    else:
        return (4*np.pi*U/g)/(np.sqrt(1 + 8*np.pi*U/g*f_e) - 1)


f_e_min = 0.1  # [Hz] minimum encounter frequency in range
f_e_max = 7.  # [Hz] maximum encounter frequency in range
f_e = np.linspace(f_e_min, f_e_max, 100*number_of_runs)

T_0_1 = np.sort(f_e_to_T_0(f_e, U))

T_0_max = f_e_to_T_0(f_e_min, U)
T_0_min = f_e_to_T_0(f_e_max, U)

T_0_2 = np.linspace(T_0_min, T_0_max, 100)

for i in range(len(T_0_1)):
    print(T_0_1[i])
