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

    # TODO: Make sure if I need to multiply by 2 to include both side hulls
    return A_33, A_35, A_53, A_55, B_33, B_35, B_53, B_55


if __name__ == "__main__":
    L_c = 28.  # [m]
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