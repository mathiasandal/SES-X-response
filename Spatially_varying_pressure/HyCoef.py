"""
Implement the added mass and damping coefficients in the same manner as Sverre Steen
in the HyCoef.for.
"""
import numpy as np
import matplotlib.pyplot as plt

def compute_hydrodynamic_coeff(L_c, Bs, U, omega_0, rho_w=1025., g=9.81):
    # Calculate encounter frequency
    omega_e = omega_0 + np.power(omega_0, 2) * U / g  # [rad/s] encounter frequency

    # Calculate two-dimensional hydrodynamic coefficients

    # 2D added mass of a rectangular cross-section
    A_33_2D = 0.3 * rho_w * Bs**2*2.  # [kg/m]  # TODO: Check if this is correct with the literature
    # 2D damping term of a rectangular cross-section
    B_33_2D = 2 * rho_w * Bs**2 * omega_e * (1 - np.tanh(4*omega_e**np.sqrt(Bs/2/g) - np.pi))

    include_end_terms = True
    if include_end_terms:
        A_33_0 = A_33_2D * L_c
        A_33 = A_33_0 - np.divide(U * B_33_2D, np.power(omega_e, 2))

        B_33_0 = B_33_2D * L_c
        B_33 = B_33_0 + U * A_33_2D

        A_35 = np.divide(-U * B_33_0, np.power(omega_e, 2)) - np.divide(U*L_c*B_33_2D, 2*np.power(omega_e, 2)) - \
               np.divide(U**2, A_33_2D, np.power(omega_e, 2))

        B_35 = U * A_33_0 + U*L_c*A_33_2D/2. - np.divide(U**2*B_33_2D, np.power(omega_e, 2))

        A_53 = np.divide(U * B_33_0, np.power(omega_e, 2)) - np.divide(U*L_c*B_33_2D, 2*np.power(omega_e, 2))

        B_53 = -U * A_33_0 + U*L_c*A_33_2D/2.

        A_55 = A_33_2D * L_c**3/12. + np.divide(U**2 * A_33_0, np.power(omega_e, 2)) - \
               np.divide(U*L_c**2*B_33_2D, np.power(2*omega_e, 2)) + np.divide(U**2*L_c*A_33_2D, 2*np.power(omega_e, 2))

        B_55 = B_33_2D*L_c**3/12. + np.divide(U**2 * B_33_0, np.power(omega_e, 2)) + U*L_c**2*A_33_2D/4. \
               - np.divide(U**2 * L_c * B_33_2D, 2*np.multiply(omega_e, 2))

    else:
        A_33_0 = A_33_2D * L_c
        A_33 = A_33_0

        B_33_0 = B_33_2D * L_c
        B_33 = B_33_0

        A_35 = np.divide(-U * B_33_0, np.power(omega_e, 2))

        B_35 = U * A_33_0

        A_53 = np.divide(U * B_33_0, np.power(omega_e, 2))

        B_53 = -U * A_33_0

        A_55 = A_33_2D * L_c ** 3 / 12. + np.divide(U ** 2 * A_33_0, np.power(omega_e, 2))

        B_55 = B_33_2D * L_c ** 3 / 12. + np.divide(U ** 2 * B_33_0, np.power(omega_e, 2))

    # TODO: Make sure if I need to multiply by 2 to include both side hulls
    return 2*A_33, 2*A_35, 2*A_53, 2*A_55, 2*B_33, 2*B_35, 2*B_53, 2*B_55


def compute_hydrostatic_coeff(m, L_oa, L_c, b_s, kb, zb, zg, rho_w=1025., g=9.81):
    A_w = 2*L_oa*b_s
    C_33 = rho_w * g * A_w
    C_35 = 0.0
    C_53 = 0.0
    C_55 = 0.0  # TODO: Need to implement this correctly
    C_55 = m * g * ((1 - kb) * zb - zg) + rho_w * g * (L_c ** 3 * b_s / 6.)

    return C_33, C_35, C_53, C_55


if __name__ == "__main__":
    L_c = 35.  # [m]
    Bs = 0.5  # [m] Beam of the side hulls  # TODO: Check what this variable means in Salvesen, Tuck and Faltinsen (1970)

    g = 9.81  # [m/s^2]  Acceleration of gravity
    U = 50.  # [knots]  Vessel velocity
    U = U*0.514444  # [m/s]  Vessel velocity

    rho_w = 1025.  # [kg/m^3]  Density of water

    omega_0 = np.linspace(0.1, 4)  # [rad/s] water wave frequency
    omega_e = omega_0 + np.power(omega_0, 2)*U/g  # [rad/s] encounter frequency
    f_encounter = omega_e/2/np.pi  # [Hz] encounter frequency

    # Calculate two-dimensional hydrodynamic coefficients

    '''
    # 2D added mass of a rectangular cross-section
    A_33_2D = 0.3 * rho_w * Bs**2*2.  # [kg/m]  # TODO: Check if this is correct with the literature
    B_33_2D = 2 * rho_w * Bs**2 * omega_e * (1 - np.tanh(4*omega_e**np.sqrt(Bs/2/g) - np.pi))

    A_33_0 = A_33_2D * L_c
    A_33 = A_33_0

    B_33_0 = B_33_2D * L_c
    B_33 = B_33_0

    A_35 = np.divide(-U * B_33_0, np.power(omega_e, 2))

    B_35 = U * A_33_0

    A_53 = np.divide(U * B_33_0, np.power(omega_e, 2))

    B_53 = -U * A_33_0

    A_55 = A_33_2D * L_c**3/12. + np.divide(U**2 * A_33_0, np.power(omega_e, 2))  # TODO:  Check what terms that should be included

    B_55 = B_33_2D*L_c**3/12. + np.divide(U**2 * B_33_0, np.power(omega_e, 2))
    '''

    A_33, A_35, A_53, A_55, B_33, B_35, B_53, B_55 = compute_hydrodynamic_coeff(L_c, Bs, U, omega_0)

    plt.plot(f_encounter, A_55)
    plt.xlabel('Encounter frequency [Hz]')
    plt.ylabel('A_55')
    plt.show()

    print('Hei')
