# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 14:24:46 2022

@author: bchan
"""

import numpy as np
import matplotlib.pyplot as plt


def cubic_eos(P, T,
              Pc: np.float32 = 4.600E6, Tc: np.float32 = 190.6, w: np.float32 = 0.008,
              eos: str = 'PR'):
    """
    Get coefficients for the specified equation of state. All variables in SI
    :param P: Pressure (Pa)
    :param V: Molar volume (m3 / mol)
    :param T: Temperature (K)
    :param Pc: Pressure at critical point (Pa). Default = 4.6 MPa (Methane)
    :param Tc: Temperature at critical point (K). Default = 190.6 K (Methane)
    :param w: Acentric factor. Default = 0.008 (Methane)
    :param eos: Desired equation of state. Default = PR (Peng-Robinson)
    :return: A dictionary of eos parameters, including cubic equation coefficients, dimensional and dimensionless attraction and covolume,
    and alpha and kappa parameters for PR
    """

    eos_dict = {'vdw': _van_der_waals_cubic,
                'PR': _peng_robinson_cubic,
                'RK': _redlich_kwong_cubic,
                'S': _soave_cubic}

    ab_dict = {'vdw': _van_der_waals_ab,
               'PR': _peng_robinson_ab}

    a, b, alph, kapp = ab_dict[eos](Pc, Tc, T, w)

    alpha, beta, gamma, A, B = eos_dict[eos](P, T, a, b)

    eos_params = {'alpha': alpha, 'beta': beta, 'gamma': gamma,
                'a': a, 'b': b, 'A': A, 'B': B,
                'alph': alph, 'kapp': kapp}

    return eos_params


def _get_compressibility_factor(P, V, T, R=8.3144598):
    """
    Get compression factor, Z, with P, V, and T in SI
    :param P: Pressure (Pa)
    :param V: Molar volume (m3 / mol)
    :param T: Temperature (K)
    :return: Compression factor (Z)
    """

    return P * V / (R * T)


def _van_der_waals_cubic(P, T, a, b):
    """
    Coefficients for the cubic form of Redlich-Kwong equation of state
    :param P: Pressure (Pa)
    :param T: Temperature (K)
    :param a:
    :param b:
    :return:
    """

    A, B = _get_A_B(P, T, a, b)

    alpha = -1 - B
    beta = A
    gamma = -1 * A * B

    return alpha, beta, gamma


def _peng_robinson_cubic(P, T, a, b):
    """
    Coefficients for the cubic form of van der Waals equation of state
    :param P: Pressure (Pa)
    :param T: Temperature (K)
    :param a:
    :param b:
    :return:
    """

    A, B = _get_A_B(P, T, a, b)

    alpha = -1 + B
    beta = A - 3* B**2 - 2 * B
    gamma = -1 * A * B + B**2 + B**3

    return alpha, beta, gamma, A, B


def _redlich_kwong_cubic(P, T, a, b):
    """
    Coefficients for the cubic form of van der Waals equation of state
    :param P: Pressure (Pa)
    :param T: Temperature (K)
    :param a:
    :param b:
    :return:
    """

    A, B = _get_A_B(P, T, a, b, is_RK=True)

    alpha = -1
    beta = A - B - B**2
    gamma = -1 * A * B

    return alpha, beta, gamma, A, B


def _soave_cubic(P, T, a, b):
    """
    Coefficients for the cubic form of van der Waals equation of state
    :param P: Pressure (Pa)
    :param T: Temperature (K)
    :param a:
    :param b:
    :return:
    """

    A, B = _get_A_B(P, T, a, b)

    alpha = -1
    beta = A - B - B ** 2
    gamma = -1 * A * B

    return alpha, beta, gamma, A, B


def _get_A_B(P, T, a, b, is_RK: bool = False):
    """
    Get the intermediate A & B parameters for cubic forms of equations of state
    :param P: Pressure (Pa)
    :param T: Temperature (K)
    :param a:
    :param b:
    :param is_RK: True if EoS is Redlich-Kwong, False otherwise
    :return: A and B
    """
    R = 8.3144598
    B = b * P / (R * T)
    if is_RK:
        A = a * P / (R**2 * T**2.5)

    else:
        A = a * P / (R * T) ** 2

    return A, B


def _van_der_waals_ab(Pc, Tc, *args, **kwargs):
    """
    Calculate a & b for van der Waals equation of state
    :param Pc:
    :param Tc:
    :return: a , b
    """
    R = 8.3144598  # Universal Gas Constant
    a = 27 * R**2 * Tc**2 * 0.015625 / Pc
    b = R * Tc * 0.125 / Pc
    alpha = 0.
    kappa = 0.
    return a, b, alpha, kappa


def _peng_robinson_ab(Pc, Tc, T, w):
    """
    Calculate a & b for Peng-Robinson equation of state
    :param Pc: Critical Pressure (Pa)
    :param Tc: Critical Temperature (K)
    :param w: Acentric factor
    :return:
    """
    R = 8.3144598
    if w <= 0.491:
        kappa = 0.37464 + 1.54226 * w - 0.26992 * w**2
    else:
        kappa = 0.379642 + 1.48503 * w - 0.164423 * w**2 + 0.016666*w**3

    alpha = (1 + kappa * (1 - np.sqrt(T / Tc)))
    a = 0.45724 * R**2 * Tc**2 * alpha**2 / Pc
    b = 0.07780 * R * Tc / Pc

    return a, b, alpha**2, kappa


def wilson_correlation(input_dict, i):
    """
    :param input_dict: Dictionary with input values from YML file
    :param i: Index of component
    :return: Vapor pressure as calculated by Wilson's correlation
    """
    return input_dict["Pc"][i] * np.exp(5.373*(1 + input_dict["w"][i])*(1 - input_dict["Tc"][i]/input_dict["T"]))

def get_molar_volume(z, T, P):
    """
    Calculate molar volume of a given phase
    :param z: Compressibility Factor
    :param input_dict: Dictionary of inputs containing pressure and temperature
    :return: Molar Volume
    """
    return z * 8.3144598 * T / P


if __name__ == "__main__":
    pass