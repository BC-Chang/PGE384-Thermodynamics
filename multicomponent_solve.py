import numpy as np
from multicomponent_utils import get_phase_compositions, get_purecomponent_a_b, get_phase_fugacity, check_TBD_sign

def newton_raphson(f, df, lower_bound, upper_bound, x0=None, tol=1e-6, maxiter=1000, **fkwargs):
    """
    Find a root using Augmented Newton-Raphson iteration.
    :param f: function to find root of
    :param df: derivative of f
    :param lower_bound: lower bound of x
    :param upper_bound: upper bound of x
    :param x0: initial guess
    :param tol: convergence tolerance
    :param maxiter: maximum number of iterations
    :param fwargs: additional arguments to be passed to lambda functions
    :return x where f(x) = 0
    """

    # Set initial error to be very large
    err = 1E9
    counter = 0

    if x0 is None:
        x_old = (lower_bound + upper_bound)/2
    else:
        x_old = x0

    while err > tol:
        x_new = x_old - f(**fkwargs, x=x_old) / df(**fkwargs, x=x_old)

        # Under relaxation conditions
        if x_new < lower_bound:
            x_new = 0.5 * (x_old + lower_bound)
        if x_new > upper_bound:
            x_new = 0.5 * (x_old + upper_bound)

        # Compute error
        err = np.abs(f(**fkwargs, x=x_new) - f(**fkwargs, x=x_old))

        # Check if iteration is larger than max iterations
        if counter == maxiter:
            print("Maximum number of iterations reached")
            break

        # Update guess
        x_old = x_new

        counter += 1
    return x_old

def get_rachford_rice():
    """
    Lambda function for Rachford-Rice and its derivative
    :return: List of lambda functions
    """
    f = lambda K, z, x: np.sum(((1 - K) * z) / (1 - (1 - K) * x))
    f_deriv = lambda K, z, x: np.sum((z * (1 - K) ** 2) / (1 + x * (K - 1)) ** 2)

    return f, f_deriv

def rachford_rice_root(K, zi, input_dict):
    """
    Find the root of the Rachford-Rice function.
    :param K: List of Binary interactions between components
    :param input_dict: dictionary of input arguments
    :return: root of the Rachford-Rice function found using Newton Raphson's method
    """

    # f and derivative of f. Here, x = beta_V
    f, f_deriv = get_rachford_rice()

    lower_bound = 1 / (1 - np.amax(K))
    upper_bound = 1 / (1 - np.amin(K))
    root = newton_raphson(f, f_deriv, lower_bound, upper_bound, tol=input_dict['eps'], maxiter=input_dict['maxiter'],
                          K=K, z=zi)

    return root

def two_phase_flash(input_dict, Ki=None, P=None, T=None, zi=None):
    """
    Perform two phase PT flash calculation
    :param input_dict: Dictionary of input parameters
    :param Ki: Initial guess for Ki parameter. Default: None uses Wilson's correlation
    :return: Final Ki of each component
    """

    if P is None:
        P = input_dict['P']
    if T is None:
        T = input_dict['T']
    if zi is None:
        zi = input_dict['zi']
    if Ki is None:
        # print("Initial Ki set using Wilson's correlation")
        Pri = input_dict["Pc"] / P
        Tri = input_dict["Tc"] / T
        Ki = Pri * np.exp((5.373 * (1 + input_dict['w']) * (1 - Tri)))

    # Set initial error to be very large
    err = 1E9
    count = 0
    while err > input_dict["eps"] and count <= input_dict["maxiter"]:
        RR_root = rachford_rice_root(Ki, zi, input_dict)
        xi, yi = get_phase_compositions(RR_root, Ki, zi)
        Ki = yi / xi

        a_ii, b_ii = get_purecomponent_a_b(input_dict['Pc'], input_dict['Tc'], T, input_dict['w'],
                                           input_dict['Nc'])

        # Get attraction and covolume parameters of liquid phase
        phi_l, f_l, _ = get_phase_fugacity(a_ii, b_ii, xi, P, input_dict)
        phi_v, f_v, _ = get_phase_fugacity(a_ii, b_ii, yi, P, input_dict)

        err = np.max(np.log(xi) + phi_l - np.log(yi) - phi_v)

        Ki = np.exp(phi_l - phi_v)

        # if err < input_dict['eps']:
        #     print(f"Flash calculation converged in {count} iterations")

        if count >= input_dict['maxiter']:
            print(f"Flash calculation did not converge after {input_dict['maxiter']} iterations")
            flash_params = None

        count += 1

    # Populate flash parameters dictionary
    flash_params = {'a_ii': a_ii, 'b_i': b_ii,
                    'xi': xi, 'yi': yi, 'beta_v': RR_root,
                    'phi_l': phi_l, 'phi_v': phi_v, 'f_l': f_l, 'f_v': f_v}

    return Ki, flash_params

def single_phase_stability(Xi_guess, a_ii, b_ii, fc_z, zi, P, input_dict: dict) -> bool:
    """
    Perform single-phase stability analysis
    :param Xi: Iterable with initial guesses for Xi parameter
    :param a_ii: Attraction parameter assuming pure component
    :param b_ii: Covolume parameter assuming pure component
    :param fc_z: Fugacity coefficient of phase z
    :param input_dict: Dictionary of input parameters
    :return: Boolean with True indicating assumed stability
    """

    # Initialize stability for each phase contained in Xi guesses
    if Xi_guess.ndim == 1:
        Xi_guess = Xi_guess[None, :]
    stability = np.zeros(Xi_guess.shape[0], dtype=bool)

    # Loop through each Xi guess individually
    for i, Xi in enumerate(Xi_guess):
        count = 0
        err = 1E9
        while err > input_dict['eps'] and count < input_dict['maxiter']:
            # Get mole fractions
            xi = Xi / np.sum(Xi)
            # Calculate compressibility factor and fugacity coefficients of x with root selection
            fc_x, f_x, root_x = get_phase_fugacity(a_ii, b_ii, xi, P, input_dict)
            # Calculate error
            err = np.max(np.log(Xi) + fc_x - np.log(zi) - fc_z)
            # If converged, check TBD sign for stability
            if err < input_dict['eps']:
                stability[i] = check_TBD_sign(Xi)
            else:
                Xi = zi * np.exp(fc_z) / np.exp(fc_x)

            count += 1
            if count == input_dict['maxiter']:
                print(f'Stability analysis did not converge after {count} iterations.')

    return stability.all()

def run_flash_main(P, T, zi, input_dict):
    """
    Run program outlined in Problem 1 of Project 2. Perform single phase stability analysis.
    If the given phase is unstable, perform 2 phase flash calculation
    :return: Either 2 phase flash results or None
    """

    # Calculate compressibility factors and fugacity coefficients for phase z
    a_ii, b_i = get_purecomponent_a_b(input_dict['Pc'], input_dict['Tc'], input_dict['T'], input_dict['w'],
                                      input_dict['Nc'])
    fc, f, root = get_phase_fugacity(a_ii, b_i, zi, P, input_dict)

    # Initialize Ki with wilson's correlation
    Pri = input_dict["Pc"] / P
    Tri = input_dict["Tc"] / T
    Ki = Pri * np.exp((5.373 * (1 + input_dict['w']) * (1 - Tri)))

    # Initial guesses for Xi assuming emerging phase is vapor-like, then liquid-like
    Xi = np.array([Ki * zi, zi / Ki])
    # print("Initial guesses for Xi: ", Xi)
    stable_single_phase = single_phase_stability(Xi, a_ii, b_i, fc, zi, P, input_dict)

    if stable_single_phase:
        return None, None
        # TODO perform single phase flash

    else:
        Ki, flash_params = two_phase_flash(input_dict, P=P, zi=zi)
        return Ki, flash_params




