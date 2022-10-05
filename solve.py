import numpy as np
from eos import cubic_eos
import matplotlib.pyplot as plt

def solve_cardanos(a: np.float32, b: np.float32, c: np.float32, d: np.float32) -> list:
    """
    Find cubic roots from polynomial of form a*x^3 + b*x^2 + c*x + d using Cardano's equation
    :param a: x^3 coefficient
    :param b: x^2 coefficient
    :param c: x coefficient
    :param d: units coefficient
    :return: A list of length 3 corresponding to cubic roots
    """
    assert a != 0, "a cannot be 0. Use quadratic formula if it is"
    R = (9*a*b*c - 27*a**2*d - 2*b**3) / (54 * a**3)
    Q = (3*a*c - b**2) / (9*a**2)
    T = (R - (Q**3 + R**2) ** (1/2)) ** (1/3)
    S =  (R + (Q**3 + R**2)**(1/2)) ** (1/3)

    # Calculate 3 roots
    x1 = S + T - b / (3*a)
    x2 = - (S + T) * 0.5 - b/(3*a) + (-3)**(1/2) * 0.5 * (S - T)
    x3 = - (S + T) * 0.5 - b/(3*a) - (-3)**(1/2) * 0.5 * (S - T)

    return x1, x2, x3

if __name__ == '__main__':
    alpha, beta, gamma, Z = cubic_eos(P=2.E6, V=1.293E-3, T=310.93, eos='vdw')
    V = np.linspace(-100, 100, 1000)
    P = V**3 + alpha*V**2 + beta*V + gamma
    plt.plot(V, P, 'o')

    print(alpha, beta, gamma)
    x1, x2, x3 = solve_cardanos(1, alpha, beta, gamma)
    print(x1, x2, x3)

    plt.show()