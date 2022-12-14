import numpy as np
import warnings
from eos import cubic_eos, get_molar_volume

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
    F = ((-D) + delta**0.5)#**(1./3)
    if F < 0:
        F = -(-F)**(1./3)
    else:
        F = F**(1./3)

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