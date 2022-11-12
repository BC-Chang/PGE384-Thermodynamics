from eos import cubic_eos, get_molar_volume
from time import perf_counter_ns
from solve import solve_cardanos, _get_real_roots
import numpy as np
from pr_utils import fugacity_coefficient

def pt_flash_nonaqueous(p, T, p_vap, input_dict):
    roots, eos_params = get_roots(p, T, input_dict)
    if p < p_vap:
        desired_root = max(roots)
    else:
        desired_root = min(roots)

    return get_molar_volume(desired_root, T, p), eos_params, desired_root

def get_vapor_pressure(input_dict):
    """
    Function to iteratively calculate the vapor pressure at a given temperature
    :param input_dict: Dictionary of input values. Must contain P, T, eos, Pc, Tc, and w keys
    :return: Input dictionary with updated vapor pressure value and compressibility factors of each root
    """

    # Start timer
    tic = perf_counter_ns()

    input_dict["zl"] = np.zeros(input_dict["Nc"])
    input_dict["zv"] = np.zeros(input_dict["Nc"])

    # Assign fugacity coefficients to input_dict
    input_dict["fc_l"] = np.zeros(input_dict["Nc"])
    input_dict["fc_v"] = np.zeros(input_dict["Nc"])
    eos_params = {}

    for j in range(input_dict["Nc"]):
        # Set initial error to be very large
        err = 1.E9
        # Initialize counter
        i = 0
        # Calculating the vapor pressure
        while err > input_dict["eps"]:
            # Solve Cardano's method and process to get compressibility factors of liquid and vapor pressure.
            roots, eos_params[f"component_{j}"] = get_roots(input_dict["Pvap"][j], input_dict["T"],
                                                            input_dict["Pc"][j], input_dict["Tc"][j], input_dict["w"][j],
                                                            input_dict)
            zl, zv = roots

            # Calculate the fugacity coefficients
            # Liquid fugacity coefficient
            fc_l = fugacity_coefficient(zl, eos_params[f"component_{j}"]["A"], eos_params[f"component_{j}"]["B"])

            # Vapor phase fugacity coefficient
            fc_v = fugacity_coefficient(zv, eos_params[f"component_{j}"]["A"], eos_params[f"component_{j}"]["B"])

            # Calculate error using fugacity coefficients
            err = abs(fc_l - fc_v)

            # Update pressure if error is greater than convergence criterion
            if err > input_dict["eps"]:
                input_dict["Pvap"][j] = input_dict["Pvap"][j] * np.exp(fc_l) / np.exp(fc_v)

            if i == input_dict["maxiter"]:
                print("Maximum number of iterations reached")
                break

            i += 1

        toc = perf_counter_ns()

        # Assign compressibility factors to input dictionary
        input_dict["zl"][j] = zl
        input_dict["zv"][j] = zv

        # Assign fugacity coefficients to input_dict
        input_dict["fc_l"][j] = fc_l
        input_dict["fc_v"][j] = fc_v

    return input_dict, eos_params

def get_roots(p, T, Pc, Tc, w, input_dict):
    eos_params = cubic_eos(P=p, T=T, eos=input_dict['eos'], Pc=Pc, Tc=Tc, w=w)
    x1, x2, x3 = solve_cardanos(1, eos_params["alpha"], eos_params["beta"], eos_params["gamma"])
    roots = _get_real_roots(x1, x2, x3)

    return roots, eos_params