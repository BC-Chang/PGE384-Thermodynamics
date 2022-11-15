import numpy as np

def pr_eos(T, mol_v, a, b):
    """
    :param T: Temperature
    :param mol_v: Molar volume
    :param a: dimensional attraction parameter
    :param b: dimensional covolume parameter
    :return: pressure at the given temperature and molar volume using PR EoS
    """

    return 8.3144598 * T / (mol_v - b) - a / (mol_v * (mol_v + b) + b * (mol_v - b))


def da_dT(Tc, Pc, T, alpha, kappa, R = 8.3144598):
    """
    da/dT
    :param Tc: Critical temperature
    :param Pc: Critical pressure
    :param T: Temperature
    :param alpha:
    :return: dadt
    """

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

def fugacity_coefficient_phi(Z, P, T, a, b, R=8.3144598):
    """
    Calculate fugacity coefficient (phi) by rearranging equations 7.4-14(a,b). This avoids using logs
    :param Z: Compressibility factor
    :param P: Pressure
    :param T: Temperature
    :param a: dimensional attraction parameter
    :param b: dimensional covolume parameter
    :param R: Universal Gas constant
    :return: Fugacity coefficient f/P (no log)
    """
    # Dimensionless covolume parameter
    B = b * P / (R * T)

    return (np.exp(Z - 1) /
            (Z - B)) / \
           (np.exp(a / (2 * np.sqrt(2) * b * R * T)) ** np.log((Z + (1 + np.sqrt(2)) * B) / (Z + (1 - np.sqrt(2)) * B)))

def fugacity_coefficient_multicomponent(Z, a_mix, b_mix, P, T, a_ij, b_i, y_i, R= 8.3144598):
    """
    Fugacity coefficient of a species with multiple components, assuming Peng-Robinson equation of state
    :param Z: Compressibility (root of cubic equation)
    :param a_mix: attraction of the mixture
    :param b_mix: covolume of the mixture
    :param P: System pressure
    :param T: System temperature
    :param a_ij: binary attraction parameter (aij)
    :param b_i: species covolume
    :param y_i: Species molar volumes
    :param R: Universal Gas constant
    :return: Fugacity coefficient of the component (ln(phi))
    """
    term_1 = -np.log(Z - b_mix * P / (R * T))
    term_2 = b_i / b_mix * (Z - 1)
    term_3a = a_mix / (2 * np.sqrt(2) * R * T * b_mix)
    my_sum = np.sum(y_i * a_ij, axis=1)

    term_3b = (2 * my_sum / a_mix) - (b_i / b_mix)
    term_3c = np.log(Z + (1 + np.sqrt(2)) * b_mix * P / (R * T)) - np.log(Z + (1 - np.sqrt(2)) * b_mix * P / (R * T))

    ln_phi = term_1 + term_2 - term_3a * term_3b * term_3c

    f = np.exp(ln_phi) * P * y_i

    return ln_phi, f
