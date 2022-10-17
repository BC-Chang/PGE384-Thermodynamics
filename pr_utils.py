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

