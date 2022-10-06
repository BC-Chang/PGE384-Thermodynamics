import numpy as np
from eos import cubic_eos
from solve import solve_cardanos
import matplotlib.pyplot as plt
from scipy.optimize import root

def main():

    global R
    R = 8.314  # Universal gas constant (J/K-mol)

    alpha, beta, gamma, a, b = cubic_eos(P=2.E6, T=310.93, eos='vdw')
    print(f"{alpha = :.5f}\n{beta = :.5f}\n{gamma = :.5f}")

    x1, x2, x3 = solve_cardanos(1, alpha, beta, gamma)
    # print(x1, x2, x3)


    Z = np.linspace(0, 1.6, 1000)
    P = Z**3 + alpha*Z**2 + beta*Z + gamma
    plt.plot(Z, P, '-')
    plt.grid(True)
    plt.show()




if __name__ == '__main__':
    main()
