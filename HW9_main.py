import numpy as np
import pandas as pd
import os

from IPython.core.pylabtools import figsize

from eos import cubic_eos, get_molar_volume, _get_compressibility_factor, _peng_robinson_ab, _peng_robinson_cubic
from singlecomponent_utils import get_vapor_pressure
import matplotlib.pyplot as plt
from io_utils import read_input, pout, redirect_stdout, close_output
import os
from multicomponent_utils import get_bubble_point, get_dew_point, get_phase_compositions
from multicomponent_solve import rachford_rice_root, get_rachford_rice
from pr_utils import fugacity_coefficient_multicomponent, fugacity_coefficient_phi
from itertools import product, combinations
from solve import solve_cardanos, _get_real_roots
from unit_conversions import Unit_Converter
from scipy.signal import argrelextrema
from scipy.optimize import fsolve, least_squares, minimize


if __name__ == "__main__":

    # Create a HW output directory if one does not already exist
    output_path = "HW9_Output/Problem_1"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Initialize a unit converter object
    unit_converter = Unit_Converter()

    # Read in input file
    input_dict = read_input(filename="Input_Files/hw9_input_file.yml")

    input_dict['P'] = unit_converter.bar_to_Pa(input_dict['P'])
    input_dict['Pc'] = unit_converter.bar_to_Pa(input_dict['Pc'])

    # Create binary interaction matrix
    K_matrix = np.zeros((input_dict['Nc'], input_dict['Nc']))
    K_matrix[0, 1:] = input_dict['Kij']['Component 0']
    K_matrix = K_matrix + K_matrix.T

    # Calculate vapor pressure at the given temperature assuming pure fluid
    pure_fluid_a_b = [_peng_robinson_ab(input_dict['Pc'][i], input_dict['Tc'][i],
                                        input_dict['T'], input_dict['w'][i]) for i in range(input_dict["Nc"])]
    a_ii = np.array([pf_a_b[0] for pf_a_b in pure_fluid_a_b])
    b_ii = np.array([pf_a_b[1] for pf_a_b in pure_fluid_a_b])

    # Create a range of mole fractions of component 1
    n_pts = 10000
    input_dict["zi"][0] = np.linspace(0+1E-9, 1 - 1E-9, n_pts)
    input_dict["zi"][1] = np.linspace(1-1E-9, 0+1E-9, n_pts)

    # Get multicomponent fugacity coefficients
    fc_multicomponent = np.zeros_like(input_dict["zi"])
    fc_purecomponent = np.zeros_like(a_ii)

    for i, (a, b) in enumerate(zip(a_ii, b_ii)):
        alpha, beta, gamma, A, B = _peng_robinson_cubic(input_dict["P"], input_dict["T"], a, b)

        x1, x2, x3 = solve_cardanos(1, alpha, beta, gamma)
        roots = _get_real_roots(x1, x2, x3)

        roots = np.unique(roots)


        fc_purecomponent[i] = fugacity_coefficient_phi(roots, input_dict["P"], input_dict["T"], a, b)


    for count in range(len(input_dict["zi"][0])):
        zi_tmp_list = [input_dict["zi"][nc][count] for nc in range(len(input_dict["zi"]))]
        # Create a matrix for calculating a_mix
        a_mix = 0
        a_matrix = np.eye(input_dict['Nc']) * a_ii
        a_ij = np.eye(input_dict['Nc']) * a_ii

        for i, j in product(range(input_dict["Nc"]), range(input_dict["Nc"])):
            if i != j:
                a_ij[i, j] = np.sqrt(a_matrix[i, i] * a_matrix[j, j]) * (1 - K_matrix[i, j])
            a_mix += a_ij[i, j] * input_dict['zi'][i][count] * input_dict['zi'][j][count]


        # Get attraction and covolume parameters of mixture
        b_mix = np.sum(zi_tmp_list * b_ii)

        alpha, beta, gamma, A_mix, B_mix = _peng_robinson_cubic(input_dict["P"], input_dict["T"], a_mix, b_mix)

        x1, x2, x3 = solve_cardanos(1, alpha, beta, gamma)
        roots = _get_real_roots(x1, x2, x3)

        roots = np.unique(roots)
        # Get fugacity coefficient
        fc_tmp, _ = fugacity_coefficient_multicomponent(roots[0], a_mix, b_mix, input_dict['P'], input_dict['T'],
                                                    a_ij, b_ii, zi_tmp_list)
        fc_multicomponent[0][count], fc_multicomponent[1][count] = np.exp(fc_tmp) #* zi_tmp_list
    # print(f'\tFugacity coefficients of components [1, 2, 3, 4]: \n\t\t{fc_multicomponent}')

    # np.sum(input_dict['zi'][0] * np.log(fc_purecomponent)
    # print(input_dict['zi'] * np.log(input_dict['zi']  * fc_multicomponent))
    delG = np.sum(input_dict['zi'] * np.log(input_dict['zi'] * fc_multicomponent), axis=0) - np.sum(input_dict['zi'] * np.log(fc_purecomponent[:, None]), axis=0)

    plt.figure(dpi=400)
    plt.style.use('seaborn')
    plt.plot(input_dict['zi'][0], delG, '-k')
    plt.xlabel('$C_1$')
    plt.ylabel('$\Delta_{mix} G_{molar} / RT$')
    plt.ylim([-0.1, 0.01])

    d_delG = np.gradient(delG, input_dict['zi'][0])
    # Indices of local minima
    delG_minima_id = argrelextrema(delG, np.less)[0].tolist()

    # Count number of lobes
    n_lobes = len(delG_minima_id)
    print(f'Number of lobes = {n_lobes}')

    # Local minima
    C1_minima = input_dict['zi'][0][delG_minima_id]
    delG_minima = delG[delG_minima_id]

    # Plot local minima
    for minima in C1_minima:
        plt.plot([minima, minima], [-1, 1], '--k')
    # plt.show()
    #
    #
    # # Take derivative of of del G with respect to C1
    # d_delG = np.gradient(delG, input_dict['zi'][0])
    # plt.figure()
    # plt.plot(input_dict['zi'][0], d_delG)
    # plt.show()

    # id_1 = abs(input_dict['zi'][0] - 0.12).argmin()
    # id_2 = abs(input_dict['zi'][0] - 0.9).argmin()
    # tangent_slopes, tangent_intercepts = np.polyfit([input_dict['zi'][0][id_1], input_dict['zi'][0][id_2]],
    #                                                 [delG[id_1], delG[id_2]], 1)
    # plt.plot(input_dict['zi'][0], input_dict['zi'][0] * tangent_slopes + tangent_intercepts)
    # plt.show()





    # Find common tangent lines
    # Define an objective function to minimize. Here, we minimize the difference between the derivative of delta G
    def F(x):
        x1, x2 = x[0], x[1]
        id_1 = int(x1)
        id_2 = int(x2)
        E1 = abs(d_delG[id_1] - d_delG[id_2])

        return E1

    # Find the tangent lines between the different lobes
    tangent_slopes = np.zeros(n_lobes)
    tangent_intercepts = np.zeros(n_lobes)
    for i, (min_a, min_b) in enumerate(combinations(delG_minima_id, 2)):
        # Set upper and lower bounds
        lower_bound = min_a + 0.03*n_pts
        upper_bound = min_b + 0.03*n_pts
        # Minimize the objective function. Set initial guess to be at the local minima.
        res = minimize(F, x0=[min_a, min_b], bounds=[(lower_bound, None), (None, upper_bound)], method='SLSQP')
        # Extract the indices corresponding to the solution to the minimization
        id_1 = int(res.x[0])
        id_2 = int(res.x[1])
        # Fit a line to the guessed points
        tangent_slopes[i], tangent_intercepts[i] = np.polyfit([input_dict['zi'][0][id_1], input_dict['zi'][0][id_2]],
                                                             [delG[id_1], delG[id_2]], 1)
        # Plot the approximate common tangents
        plt.plot(input_dict['zi'][0], input_dict['zi'][0] * tangent_slopes[i] + tangent_intercepts[i], '--',
                 linewidth=1.5, label=f'Tangent {i+1}')

    plt.legend()
    plt.show()

    # Problem 1c:
    # Dimensionless tangent plane distance
    # Find approximate indices of reference points
    zi_id = [abs(input_dict['zi'][nc] - zi_ref).argmin() for nc, zi_ref in enumerate(input_dict['zi_ref'])]
    # Get fugacity coefficients of components at reference compositions
    fc_ref = np.array([fc_multicomponent[nc, zi_idx] for nc, zi_idx in enumerate(zi_id)])

    Dr = np.sum(input_dict['zi'] * np.log(input_dict['zi'] * fc_multicomponent), axis=0) - np.sum(input_dict['zi'] * np.log(input_dict['zi_ref']*fc_ref)[:, None], axis=0)

    plt.figure(dpi=400)
    plt.plot(input_dict['zi'][0], Dr, label='D/RT')
    plt.plot([input_dict['zi_ref'][0], input_dict['zi_ref'][0]], [1, -1], '--k', label='C1 Composition')
    plt.plot([-1, 1.1], [0, 0], '-k')
    plt.ylim([-0.05, 0.35])
    plt.xlim([0, 1])
    plt.xlabel('$C_1$')
    plt.ylabel('$D/RT$')
    plt.legend()
    plt.show()

    # Problem 1d:
    # Dimensionless tangent plane distance
    # Find approximate indices of reference points







