import numpy as np
import pandas as pd
import os
from eos import cubic_eos, get_molar_volume, _get_compressibility_factor, _peng_robinson_ab, _peng_robinson_cubic
from singlecomponent_utils import get_vapor_pressure
import matplotlib.pyplot as plt
from io_utils import read_input, pout, redirect_stdout, close_output
import os
from multicomponent_utils import get_bubble_point, get_dew_point, get_phase_compositions
from multicomponent_solve import rachford_rice_root, get_rachford_rice
from pr_utils import fugacity_coefficient_multicomponent
from itertools import product
from solve import solve_cardanos, _get_real_roots
from unit_conversions import Unit_Converter


def p2_main():
    # Create a HW output directory if one does not already exist
    output_path = "HW8_Output/Problem_2"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Initialize a unit converter object
    unit_converter = Unit_Converter()

    # Redirect output to the output directory
    # redirect_stdout(f"{output_path}/output_file.txt")

    # Read in input file
    input_dict = read_input(filename="Input_Files/hw8_input_file_p2.yml")

    print("*"*50)
    print("Homework 8 Problem 2: ")
    print("\nProblem 2b:")
    print("\tInitial vapor pressure guesses using Wilson's correlation:")
    print(f"\t\t{unit_converter.Pa_to_psi(input_dict['Pvap'])} psi")
    print("\tK values using Wilson's correlation (Problem 2b)")
    K_wilson = input_dict['Pvap']/input_dict['P']
    print("\t\t", K_wilson)

    # Calculate vapor pressure at the given temperature assuming pure fluid
    input_dict, eos_params = get_vapor_pressure(input_dict)
    print("\nProblem 2a:")
    print(f"\tVapor Pressure of components [1, 2, 3] at T = {unit_converter.K_to_F(input_dict['T']):.2f}F")
    print(f"\t\t{unit_converter.Pa_to_psi(input_dict['Pvap'])} psi")


    print("\tK values using Raoult's Law")
    K_raoults = input_dict['Pvap']/input_dict['P']
    print("\t\t", K_raoults)



    # Rachford-Rice Internal loop using Newton Iteration
    RR_root = rachford_rice_root(K_raoults, input_dict)
    print('RR_root', RR_root)

    xi, yi = get_phase_compositions(RR_root, K_raoults, input_dict['zi'])
    print("\nProblem 2c:")
    print('\tPhase compositions:')
    print(f"\t\tLiquid Phase Compositions (x_i) of components [1, 2, 3]: \n\t\t\t{xi}")
    print(f"\t\tVapor Phase Compositions (y_i) of components [1, 2, 3]: \n\t\t\t{yi}")

    # Get Dew and Bubble points assuming Raoults law
    bubble_pt = get_bubble_point(input_dict["Pvap"], input_dict["zi"])
    dew_pt = get_dew_point(input_dict["Pvap"], input_dict["zi"])

    print("\tPoles of RR:")
    print(f"\t\tMole fraction of vapor phase (beta_v): ", RR_root)
    print(f"\t\tLower bound of Beta_V (Pole 1): {1 / (1 - max(K_raoults))}")
    print(f"\t\tUpper bound of Beta_V (Pole 2): {1 / (1 - min(K_raoults))}")
    print("\nProblem 2d:")
    print(f"\tDew point:  {unit_converter.Pa_to_psi(dew_pt)} psi")
    print(f"\tBubble point: {unit_converter.Pa_to_psi(bubble_pt)} psi")



    # Create Rachford-Rice Array for Plotting
    # Get lambda function for Rachford-Rice
    f, _ = get_rachford_rice()

    beta_vals = np.linspace(1 / (1 - np.amax(K_raoults)) + 1E-6, 1 / (1 - np.amin(K_raoults)) - 1E-6, 1000)
    f_vals = np.zeros_like(beta_vals)
    for i, b_val in enumerate(beta_vals):
        f_vals[i] = f(K_raoults, input_dict["zi"], b_val)
    # Plot Rachford-Rice for values of T between lower bound and upper bound
    plt.style.use('seaborn')
    plt.figure()
    plt.plot(beta_vals[2:-1], f_vals[2:-1], '-b', alpha=0.8, label=r'$f(\beta_V)$')
    plt.plot(RR_root, 0, '*y', label=r'$\beta_V$')
    plt.ylim([-40, 20])
    plt.xlabel(r'$\beta_V$')
    plt.ylabel(r'$f(\beta_V)$')
    plt.legend()
    plt.savefig(os.path.join(output_path, 'Problem2C_RR.png'))
    plt.show()

    # Close the output file
    # close_output(f"{output_path}/output_file.txt")

    return


# Code for problem 3
def p3_main():
    R = 8.3144598

    # Create a HW output directory if one does not already exist
    output_path = "HW8_Output/Problem_3"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    unit_converter = Unit_Converter()

    # Redirect output to the output directory
    redirect_stdout(f"{output_path}/output_file.txt")

    # Get a & b for each component assuming pure fluid
    input_dict = read_input(filename="Input_Files/hw8_input_file_p3.yml")
    print("*"*50)
    print("Problem 3: ")

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
    print(f'\tFugacity coefficients of components [1, 2, 3, 4]: \n\t\t{fc}')

    # Fugacity
    print('\tFugacity of components [1, 2, 3, 4]: ')
    print("\t\t", unit_converter.Pa_to_psi(f), "psia")

    # Close the output file
    close_output(f"{output_path}/output_file.txt")

    return


if __name__ == "__main__":
    # Uncomment the following line to run problem 2
    p2_main()

    # Uncomment the following line to run problem 3
    p3_main()






