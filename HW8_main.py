import numpy as np
import pandas as pd
import os
from eos import cubic_eos, get_molar_volume, _get_compressibility_factor, _peng_robinson_ab, _peng_robinson_cubic
from singlecomponent_utils import get_vapor_pressure
import matplotlib.pyplot as plt
from io_utils import read_input, pout, redirect_stdout, close_output
import os
from multicomponent_utils import raoults_law, get_bubble_point, get_dew_point
from multicomponent_solve import rachford_rice_root, get_rachford_rice
from pr_utils import fugacity_coefficient_multicomponent
from itertools import product
from solve import solve_cardanos, _get_real_roots


def p2_main():
    # Create a HW output directory if one does not already exist
    output_path = "HW8_Output/Problem_2"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Redirect output to the output directory
    redirect_stdout(f"{output_path}/output_file.txt")

    # Read in input file
    input_dict = read_input(filename="Input_Files/hw8_input_file_p2.yml")

    print("*"*50)
    print(f"Initial pressure guess using Wilson's correlation:")
    print(input_dict['Pvap'])
    print("K values using Wilson's correlation")
    K_wilson = input_dict['Pvap']/input_dict['P']
    print(K_wilson)

    # Calculate vapor pressure at the given temperature assuming pure fluid
    input_dict, eos_params = get_vapor_pressure(input_dict)
    print(f"Vapor Pressure of components [1, 2, 3] at T = {input_dict['T']}K:")
    print(input_dict['Pvap'])


    print("K values using Raoult's Law")
    K_raoults = input_dict['Pvap']/input_dict['P']
    print(K_raoults)

    f, _ = get_rachford_rice()

    beta_vals = np.linspace(1/(1-np.amax(K_raoults))+1E-6, 1/(1-np.amin(K_raoults))-1E-6, 1000)
    f_vals = np.zeros_like(beta_vals)
    for i, b_val in enumerate(beta_vals):
        f_vals[i] = f(K_raoults, input_dict["zi"], b_val)

    # Rachford-Rice Internal loop using Newton Iteration
    RR_root = rachford_rice_root(K_raoults, input_dict)

    bubble_pt = get_bubble_point(input_dict["Pvap"], input_dict["zi"])
    dew_pt = get_dew_point(input_dict["Pvap"], input_dict["zi"])

    print(f"Beta_V: ", RR_root)
    print(f"Lower bound of Beta_V: {1 / (1 - max(K_raoults)) : .3f}")
    print(f"Upper bound of Beta_V: {1 / (1 - min(K_raoults)) : .3f}")
    print(f"Dew point:  {dew_pt:.3f} Pa")
    print(f"Bubble point: {bubble_pt:.3f} Pa")

    plt.style.use('seaborn')
    plt.figure()
    plt.plot(beta_vals[2:-1], f_vals[2:-1], '-b', alpha=0.8, label=r'$f(\beta_V)$')
    # plt.plot([beta_vals[1], beta_vals[-2]], [0,0], '--k')
    plt.plot(RR_root, 0, '*y', label=r'$\beta_V$')
    plt.ylim([-40, 20])
    plt.xlabel(r'$\beta_V$')
    plt.ylabel(r'$f(\beta_V)$')
    plt.legend()
    plt.savefig(os.path.join(output_path, 'RR.png'))
    plt.show()



    # f_tmp = f(K_raoults, input_dict["xi"], 0.5)
    # f_deriv_tmp = f_deriv(K_raoults, input_dict["xi"], 0.5)

    # Close the output file
    close_output(f"{output_path}/output_file.txt")
    return


# Code for problem 3
def p3_main():
    R = 8.3144598

    # Create a HW output directory if one does not already exist
    output_path = "HW8_Output/Problem_3"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Redirect output to the output directory
    redirect_stdout(f"{output_path}/output_file.txt")

    # Get a & b for each component assuming pure fluid
    input_dict = read_input(filename="Input_Files/hw8_input_file_p3.yml")
    print(f"Initial pressure guess using Wilson's correlation:")
    print(input_dict['Pvap'])

    # Create binary interaction matrix
    K_matrix = np.zeros((input_dict['Nc'], input_dict['Nc']))
    K_matrix[0, 1:] = input_dict['Kij']['Component 0']
    K_matrix = K_matrix + K_matrix.T

    # Calculate vapor pressure at the given temperature assuming pure fluid
    pure_fluid_a_b = [_peng_robinson_ab(input_dict['Pc'][i], input_dict['Tc'][i],
                                        input_dict['T'], input_dict['w'][i]) for i in range(input_dict["Nc"])]
    a_ii = np.array([pf_a_b[0] for pf_a_b in pure_fluid_a_b])
    b_ii = np.array([pf_a_b[1] for pf_a_b in pure_fluid_a_b])

    # Create a matrix for calculating a_mix
    a_mix = 0
    a_matrix = np.eye(input_dict['Nc']) * a_ii
    a_ij = np.eye(input_dict['Nc']) * a_ii
    for i, j in product(range(input_dict["Nc"]), range(input_dict["Nc"])):
        if i != j:
            a_ij[i, j] = np.sqrt(a_matrix[i, i] * a_matrix[j, j]) * (1 - K_matrix[i, j])
        # a_ij[i, j] = a_matrix[i, j]
        a_mix += a_ij[i, j] * input_dict['zi'][i] * input_dict['zi'][j]

    # Get attraction and covolume parameters of mixture
    # a_mix = np.sum(a_matrix)
    b_mix = np.sum(input_dict["zi"] * b_ii)

    alpha, beta, gamma, A_mix, B_mix = _peng_robinson_cubic(input_dict["P"], input_dict["T"], a_mix, b_mix)

    x1, x2, x3 = solve_cardanos(1, alpha, beta, gamma)
    roots = _get_real_roots(x1, x2, x3)

    roots = np.unique(roots)

    # Get fugacity coefficient
    fc, f = fugacity_coefficient_multicomponent(roots[0], a_mix, b_mix, input_dict['P'], input_dict['T'],
                                                a_ij, b_ii, input_dict['zi'])
    print('Fugacity coefficient')
    print(fc)

    # Fugacity
    print('Fugacity')
    print(f * 1.451e-4, "psia")

    # Close the output file
    close_output(f"{output_path}/output_file.txt")

    return


if __name__ == "__main__":
    # Uncomment the following line to run problem 2
    p2_main()

    # Uncomment the following line to run problem 3
    p3_main()






