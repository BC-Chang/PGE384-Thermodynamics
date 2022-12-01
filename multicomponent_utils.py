import numpy as np
from eos import _peng_robinson_ab, _peng_robinson_cubic
from itertools import product
from solve import solve_cardanos, _get_real_roots
from pr_utils import fugacity_coefficient_multicomponent

def raoults_law(pvap, xi):
    """
    Calculate vapor pressure of a component assuming Raoults law
    :param pvap: Equilibrium vapor pressure of component i assuming pure component
    :param xi: Mole fraction of component i
    :return: Partial pressure of component i assuming Raoults law
    """

    return pvap * xi

def get_bubble_point(pvap, z):
    """
    Calculate bubble point of a component assuming Raoults law
    :param pvap: Equilibrium vapor pressure of component i assuming pure component
    :param z: Mole fraction of component i
    :return: Bubble point assuming Raoults law
    """

    return np.sum(pvap * z)


def get_dew_point(pvap, z):
    """
    Calculate dew point of a component assuming Raoults law
    :param pvap: Equilibrium vapor pressure of component i assuming pure component
    :param z: Mole fraction of component i
    :return: Dew point assuming Raoults law
    """

    return 1 / np.sum(z / pvap)

def get_phase_compositions(beta_v, Ki, Zi):
    """
    Calculate phase compositions of a component using root of Raoults law
    :param beta_v: Root of Raoults law
    :param Ki: Binary interactions assuming Raoults law
    :param Zi: List of compressibility (roots of cubic) for each pure component. [Zl, Zv]
    :return: Liquid and vapor phase compositions
    """

    xi = Zi / (1 - (1 - Ki) * beta_v)  # np.zeros_like(Zl)
    yi = Ki * xi  # np.zeros_like(Zl)

    # for i in range(len(Zl)):
    #     xi[i] = Zl[i] / (1  - (1 - Ki[i]) * beta_v)
    #     yi[i] = Ki[i] * xi[i]

    return xi, yi

def get_purecomponent_a_b(Pc: np.ndarray, Tc: np.ndarray, T: float, w: np.ndarray, Nc: int):
    """
    Get attraction and covolume of each component assuming pure components
    :param input_dict: Input dictionary
    :return: Lists of attraction and covolume parameters for each component
    """
    # Calculate vapor pressure at the given temperature assuming pure fluid
    pure_fluid_a_b = [_peng_robinson_ab(Pc[i], Tc[i],
                                        T, w[i]) for i in range(Nc)]
    a_ii = np.array([pf_a_b[0] for pf_a_b in pure_fluid_a_b])
    b_ii = np.array([pf_a_b[1] for pf_a_b in pure_fluid_a_b])

    return a_ii, b_ii

def get_multicomponent_A_B(a_ii, b_ii, zi, K_ij, P, T, Nc):
    # Create a matrix for calculating a_mix
    a_mix = 0
    a_matrix = np.eye(Nc) * a_ii
    a_ij = np.eye(Nc) * a_ii
    for i, j in product(range(Nc), range(Nc)):
        if i != j:
            a_ij[i, j] = np.sqrt(a_matrix[i, i] * a_matrix[j, j]) * (1 - K_ij[i, j])
        # a_ij[i, j] = a_matrix[i, j]
        a_mix += a_ij[i, j] * zi[i] * zi[j]

    # Get attraction and covolume parameters of mixture
    # a_mix = np.sum(a_matrix)
    b_mix = np.sum(zi * b_ii)

    alpha, beta, gamma, A_mix, B_mix = _peng_robinson_cubic(P, T, a_mix, b_mix)

    return alpha, beta, gamma, a_mix, b_mix, A_mix, B_mix, a_ij

def get_phase_fugacity(a_ii, b_ii, zi, P, input_dict):
    alpha, beta, gamma, a_mix, b_mix, A_mix, B_mix, a_ij = get_multicomponent_A_B(a_ii, b_ii, zi, input_dict['K_ij'],
                                                                            P, input_dict['T'], input_dict['Nc'])


    x1, x2, x3 = solve_cardanos(1, alpha, beta, gamma)
    roots = _get_real_roots(x1, x2, x3)
    roots = np.unique(roots)
    # Get fugacity coefficient
    fc_a, f_a = fugacity_coefficient_multicomponent(max(roots), a_mix, b_mix, P, input_dict['T'],
                                                a_ij, b_ii, zi)
    fc_b, f_b = fugacity_coefficient_multicomponent(min(roots), a_mix, b_mix, P, input_dict['T'],
                                                a_ij, b_ii, zi)

    fc, f, root = choose_root(zi, [fc_a, fc_b], [f_a, f_b], roots)

    return fc, f, root

def choose_root(zi, fc, f, roots):
    """
    Choose the root that results in the lower Gibbs Free Energy
    :param zi: Mole fraction of component i
    :param fc: List of fugacity coefficients
    :param f: List of fugacities
    :return: Fugacity coefficient, fugacity, and root to use
    """
    fc_a, fc_b = fc
    f_a, f_b = f
    G = np.sum(zi * (fc_a - fc_b))
    if G > 0:
        fc = fc_b
        f = f_b
        root = min(roots)
    else:
        fc = fc_a
        f = f_a
        root = max(roots)

    return fc, f, root

def check_TBD_sign(Xi):
    """
    Check TBD sign for stability. 0 or Positive TBD means the phase is assumed to be stable. Negative TBD means the phase is unstable.
    :param Xi:
    :return: True if assumed stable, False if unstable
    """
    Xi_sum = np.sum(Xi)
    if Xi_sum > 1.0:
        return False
    else:
        return True

def get_total_composition(xi, yi, beta_v):
    """
    Get the total composition (z_i) from phase compositions (x_i and y_i)
    :param xi: Composition of liquid phase
    :param yi: Composition of vapor phase
    :param beta_v: Phase composition (root of RR equation)
    :return: Total composition
    """

    beta = np.array([1-beta_v, beta_v])
    x_ij = np.array([xi, yi])
    zi = np.sum(beta[:, None] * x_ij, axis=0)

    return zi

