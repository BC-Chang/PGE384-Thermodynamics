import numpy as np
import warnings
from eos import cubic_eos, get_molar_volume
from pr_utils import fugacity_coefficient
from time import perf_counter_ns

def _get_delta(a: np.float32, b: np.float32, c: np.float32, d: np.float32) -> list:
    """
    Find Delta to return combination of real and imaginary roots
    :param a: x^3 coefficient
    :param b: x^2 coefficient
    :param c: x coefficient
    :param d: units coefficient
    :return: A list of length 3 corresponding to D, E, and Delta from the appendix of Wilcezek-Vera et al. (2015)
    """
    D = (b/3)**3 - b*c/6 + d*0.5
    E = c/3 - (b/3)**2
    delta = D**2 + E**3
    return D, E, delta

def _delta_0(b, D, *args, **kwargs):
    """
    Find roots when delta = 0. There are 3 real roots with at least two equal roots
    :param b: x^2 coefficient
    :param D: D from appendix of Wilcezek-Vera et al. (2015)
    :return: List of 3 real roots with at least 2 equal roots
    """
    x1 = 2*(-D)**(1/3) - (b/3)
    x2 = -(-D)**(1/3) - (b/3)
    x3 = x2.copy()
    
    return x1, x2, x3

def _delta_positive(b, D, delta, *args, **kwargs):
    """
    Find roots when delta > 0. There is 1 real root and 2 complex roots
    :param b: x^2 coefficient
    :param D: D from appendix of Wilcezek-Vera et al. (2015)
    :return: List of 1 real root and 2 complex roots
    """
    F = ((-D) + delta**0.5)**(1/3)
    G = ((-D) - delta**0.5)
    if G < 0:
        G = -(-G)**(1./3)
    else:
        G = G**(1./3)
    
    x1 = F + G - (b/3)
    x2 = -(0.5*(F + G) + b/3) + 3**(0.5)*0.5j*(F-G)#*(-1)**0.5
    x3 = -(0.5*(F + G) + b/3) - 3**(0.5)*0.5j*(F-G)#*(-1)**0.5
    
    return x1, x2, x3

def _delta_negative(b, D, E, *args, **kwargs):
    """
    Find roots when delta < 0. There are 3 real and unequal roots
    :param b: x^2 coefficient
    :param D: D from appendix of Wilcezek-Vera et al. (2015)
    :param E:
    :param args:
    :param kwargs:
    :return:
    """
    theta = np.arccos(-D/(-E**3)**0.5)
    
    x1 = 2*(-E)**(0.5)*np.cos(theta / 3) - b/3
    x2 = 2*(-E)**(0.5)*np.cos(theta / 3 + 2/3*np.pi) - b/3
    x3 = 2*(-E)**(0.5)*np.cos(theta / 3 + 4/3*np.pi) - b/3
    
    return x1, x2, x3

def _get_real_roots(*args):
    """
    Check the roots to sort whether they correspond to vapor or liquid phase
    :param x1: root 1
    :param x2: root 2
    :param x3: root 3
    :return: Root corresponding to liquid phase, root corresponding to vapor phase
    """
    # Figure out which roots are real numbers and ignore complex roots
    real_roots = [x for x in args if isinstance(x, float)]

    # Make sure there is at least one real root:
    assert len(real_roots) >= 1, "No real roots, check inputs."

    # Root corresponding to liquid phase = minimum of real roots
    xl = np.amin(real_roots)

    # Root corresponding to vapor phase = maximum of real roots
    xv = np.amax(real_roots)

    if len(real_roots) < 2:
        warnings.warn("Warning.... Only 1 real root.")

    return xl, xv




def solve_cardanos(a: np.float32, b: np.float32, c: np.float32, d: np.float32) -> tuple:
    """
    Find cubic roots from polynomial of form a*x^3 + b*x^2 + c*x + d using Cardano's equation
    :param a: x^3 coefficient
    :param b: x^2 coefficient
    :param c: x coefficient
    :param d: units coefficient
    :return: A list of length 3 corresponding to cubic roots
    """
    assert a != 0, "a cannot be 0. Use quadratic formula if it is"

    D, E, delta = _get_delta(a, b, c, d)

    
    if delta == 0:
        return _delta_0(b, D)
    elif delta > 0:
        return _delta_positive(b, D, delta)
    elif delta < 0:
        return _delta_negative(b, D, E)
    else:
        raise Exception("Delta is not a real number :(")
        
    # R = (9*a*b*c - 27*a**2*d - 2*b**3) / (54 * a**3)
    # Q = (3*a*c - b**2) / (9*a**2)
    # T = (R - (Q**3 + R**2) ** (1/2)) ** (1/3)
    # S =  (R + (Q**3 + R**2)**(1/2)) ** (1/3)

    # # Calculate 3 roots
    # x1 = S + T - b / (3*a)
    # x2 = - (S + T) * 0.5 - b/(3*a) + (-3)**(1/2) * 0.5 * (S - T)
    # x3 = - (S + T) * 0.5 - b/(3*a) - (-3)**(1/2) * 0.5 * (S - T)

    # return x1, x2, x3

def get_vapor_pressure(input_dict):
    """
    Function to iteratively calculate the vapor pressure at a given temperature
    :param input_dict: Dictionary of input values. Must contain P, T, eos, Pc, Tc, and w keys
    :return: Input dictionary with updated vapor pressure value and compressibility factors of each root
    """

    # Set initial error to be very large
    err = 1.E9
    # Initialize counter
    i = 0
    # Start timer
    tic = perf_counter_ns()

    # Calculating the vapor pressure
    while err > input_dict["eps"]:
        # Solve Cardano's method and process to get compressibility factors of liquid and vapor pressure.
        roots, eos_params = get_roots(input_dict["P"], input_dict["T"], input_dict)
        zl, zv = roots

        # Calculate the fugacity coefficients
        # Liquid fugacity coefficient
        fc_l = fugacity_coefficient(zl, eos_params["A"], eos_params["B"])

        # Vapor phase fugacity coefficient
        fc_v = fugacity_coefficient(zv, eos_params["A"], eos_params["B"])

        # Calculate error using fugacity coefficients
        err = abs(fc_l - fc_v)

        # Update pressure if error is greater than convergence criterion
        if err > input_dict["eps"]:
            input_dict["P"] = input_dict["P"] * np.exp(fc_l) / np.exp(fc_v)

        if i == input_dict["maxiter"]:
            print("Maximum number of iterations reached")
            break

        i += 1

    toc = perf_counter_ns()

    # Assign compressibility factors to input dictionary
    input_dict["zl"] = zl
    input_dict["zv"] = zv

    # Assign fugacity coefficients to input_dict
    input_dict["fc_l"] = fc_l
    input_dict["fc_v"] = fc_v

    return input_dict, eos_params

def get_roots(p, T, input_dict):
    eos_params = cubic_eos(P=p, T=T, eos=input_dict['eos'],
                                        Pc=input_dict["Pc"], Tc=input_dict["Tc"], w=input_dict["w"])
    x1, x2, x3 = solve_cardanos(1, eos_params["alpha"], eos_params["beta"], eos_params["gamma"])
    roots = _get_real_roots(x1, x2, x3)

    return roots, eos_params


def pt_flash_nonaqueous(p, T, p_vap, input_dict):
    roots, eos_params = get_roots(p, T, input_dict)
    if p < p_vap:
        desired_root = max(roots)
    else:
        desired_root = min(roots)

    return get_molar_volume(desired_root, T, p), eos_params, desired_root



if __name__ == '__main__':

    from eos import cubic_eos
    import matplotlib.pyplot as plt
    from io import read_input
    
    alpha, beta, gamma, Z = cubic_eos(P=2.E6, V=1.293E-3, T=310.93, eos='vdw')
    V = np.linspace(-0, 1.6, 1000)
    P = V**3 + alpha*V**2 + beta*V + gamma
    plt.plot(V, P, '-')

    print(alpha, beta, gamma)
    x1, x2, x3 = solve_cardanos(1, alpha, beta, gamma)
    print(x1, x2, x3)

    plt.show()