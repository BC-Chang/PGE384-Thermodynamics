import numpy as np

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

def rachford_rice_root(K, input_dict):
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
                          K=K, z=input_dict["zi"])

    return root

