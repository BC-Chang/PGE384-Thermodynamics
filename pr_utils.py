import numpy as np

def da_dT(Tc, Pc, T, alpha, kappa):
    """
    da/dT
    :param Tc: Critical temperature
    :param Pc: Critical pressure
    :param T: Temperature
    :param alpha:
    :return: dadt
    """
    R = 8.314
    dadT = -0.45724 * R**2 * Tc**2 / Pc * kappa * np.sqrt(alpha / (T * Tc))

    return dadT

def enthalpy_departure(T, Z, dadT, a, b, B):
    """
    Calculate enthalpy departure from equation 6.4-29 in Sandler 5th Edition
    :param T: Temperature (K)
    :param Z: Critical compressibility factor
    :param dadT: dadT
    :param a:
    :param b:
    :return: enthalpy departure Hd
    """
    R = 8.314
    h_dep = R * T * (Z - 1) + (T * dadT - a) / (2 * np.sqrt(2) * b) * \
           (np.log(Z + (1 + np.sqrt(2)) * B) - np.log(Z + (1 - np.sqrt(2)) * B))

    return h_dep

def entropy_departure(Z, dadT, b, B):
    """
    Calculate entropy departure from equation 6.4-30 in Sandler 5th Edition
    :param Z: Critical pressure
    :param dadt: dadt
    :param b:
    :param B:
    :return:
    """

    R = 8.314
    s_dep = R*np.log(Z-B) + dadT/(2*np.sqrt(2)*b) * (np.log(Z+(1+np.sqrt(2))*B) - np.log(Z+(1-np.sqrt(2))*B))

    return s_dep


def fugacity_coefficient(Z, A, B):
    """
    Calculate fugacity coefficient using equations 7.4-14(a,b) from Sandler 5th edition
    :param Z: Compressibility given by root (either liquid or vapor compressibility)
    :param A: Dimensionless attraction
    :param B: Dimensionless covolume
    :return: Fugacity coefficient
    """
    # First term of eq 7.4-14
    fc1 = Z - 1

    # Second term of eq. 7.4-14
    fc2 = np.log(Z - B)

    # Third term of eq. 7.4-14
    fc3 = A / (2*np.sqrt(2)*B) * (np.log(Z + (1 + np.sqrt(2))*B) - np.log(Z + (1 - np.sqrt(2))*B))

    # Put them all together
    fc = fc1 - fc2 - fc3

    return fc
