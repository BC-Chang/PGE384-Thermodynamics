import numpy as np
import os
from eos import cubic_eos, get_molar_volume, _get_compressibility_factor, _peng_robinson_ab, _peng_robinson_cubic
from singlecomponent_utils import get_vapor_pressure
import matplotlib.pyplot as plt
from io_utils import read_input, pout, redirect_stdout, close_output
import os
from multicomponent_utils import get_bubble_point, get_dew_point, get_phase_compositions, get_purecomponent_a_b, \
    get_multicomponent_A_B, get_phase_fugacity, check_TBD_sign, get_total_composition
from multicomponent_solve import two_phase_flash, single_phase_stability, run_flash_main
from pr_utils import fugacity_coefficient_multicomponent, fugacity_coefficient_phi
from solve import solve_cardanos, _get_real_roots
from unit_conversions import Unit_Converter
from scipy.optimize import minimize
from itertools import product, combinations
from scipy.signal import argrelextrema



def initialize(output_path="Project2_Output/Problem2", input_filepath="Input_Files/project2_input_file_2.yml",
                         write_out=False):
    # Create a HW output directory if one does not already exist
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Initialize a unit converter object
    unit_converter = Unit_Converter()

    # Redirect output to the output directory
    if write_out:
        redirect_stdout(f"{output_path}/output_file.txt")

    # Read in input file
    input_dict = read_input(filename=input_filepath)

    # Convert bar to Pa
    input_dict['P'] = unit_converter.bar_to_Pa(input_dict['P'])
    input_dict['Pc'] = unit_converter.bar_to_Pa(input_dict['Pc'])

    return output_path, unit_converter, input_dict


# Two-Phase PT Flash
def p2a_main(output_path, unit_converter, input_dict):
    # Number of points to test
    n_pts = 100

    # Create a range of pressures to test
    P = np.linspace(7.91, 60, n_pts)
    P = unit_converter.bar_to_Pa(P)

    xi_flash = np.zeros_like(P)
    yi_flash = np.zeros_like(P)

    print("*" * 50)
    for i, Pi in enumerate(P):
        Ki, flash_params = run_flash_main(Pi, input_dict['T'], input_dict['zi'], input_dict)
        if flash_params is None:
            flash_params = {'xi': [0, 0], 'yi': [0, 0]}
        xi_flash[i] = flash_params['xi'][0]
        yi_flash[i] = flash_params['yi'][0]

    P_plot = unit_converter.Pa_to_bar(P)

    fig = plt.figure()
    plt.plot(xi_flash[xi_flash > 0], P_plot[xi_flash > 0], '-b', label='Flash calculation bubble point')
    plt.plot(yi_flash[yi_flash > 0], P_plot[yi_flash > 0], '-r', label='Flash calculation dew point')
    plt.plot([xi_flash[n_pts//4::n_pts//4], yi_flash[n_pts//4::n_pts//4]], [P_plot[n_pts//4::n_pts//4], P_plot[n_pts//4::n_pts//4]], 'o--k')
    plt.xlabel('Mole fraction of Component 1')
    plt.ylabel('Pressure (bar)')
    plt.ylim([0, 75])
    plt.xlim([0, 1.1])
    plt.legend()
    return fig


def p2b_main(output_path, unit_converter, input_dict):

    fig = p2a_main(output_path, unit_converter, input_dict)

    joyce_P = np.array([7.91, 14.80, 25.14, 35.49, 45.83, 56.17, 62.38, 64.51])
    joyce_xA = 1 - np.array([0.765, 0.626, 0.449, 0.330, 0.240, 0.160, 0.116, 0.0702])
    joyce_yA = 1 - np.array([0.00466, 0.00347, 0.00364, 0.00493, 0.00812, 0.0157, 0.0306, 0.0702])
    plt.plot(joyce_xA[:-1], joyce_P[:-1], '^b', label='Joyce et al. (2000) liquid compositions')
    plt.plot(joyce_yA[:-1], joyce_P[:-1], '^r', label='Joyce et al. (2000) vapor compositions')
    plt.plot(joyce_xA[-1], joyce_P[-1], '*y', markersize=13, label='Joyce et al. (2000) critical point')
    plt.legend()

    return fig


def p2c_main(output_path, unit_converter, input_dict):
    fig = p2b_main(output_path, unit_converter, input_dict)

    return fig


def p2d_main(input_filepath, fig=None, **plot_kwargs):
    output_path, unit_converter, input_dict = initialize(input_filepath=input_filepath)
    # Calculate vapor pressure at the given temperature assuming pure fluid
    pure_fluid_a_b = [_peng_robinson_ab(input_dict['Pc'][i], input_dict['Tc'][i],
                                        input_dict['T'], input_dict['w'][i]) for i in range(input_dict["Nc"])]
    a_ii = np.array([pf_a_b[0] for pf_a_b in pure_fluid_a_b])
    b_ii = np.array([pf_a_b[1] for pf_a_b in pure_fluid_a_b])

    # Create a range of mole fractions of component 1
    n_pts = 10000
    input_dict["zi"][0] = np.linspace(0 + 1E-9, 1 - 1E-9, n_pts)
    input_dict["zi"][1] = np.linspace(1 - 1E-9, 0 + 1E-9, n_pts)

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
                a_ij[i, j] = np.sqrt(a_matrix[i, i] * a_matrix[j, j]) * (1 - input_dict['K_ij'][i, j])
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
        fc_multicomponent[0][count], fc_multicomponent[1][count] = np.exp(fc_tmp)  # * zi_tmp_list
    # print(f'\tFugacity coefficients of components [1, 2, 3, 4]: \n\t\t{fc_multicomponent}')

    # np.sum(input_dict['zi'][0] * np.log(fc_purecomponent)
    # print(input_dict['zi'] * np.log(input_dict['zi']  * fc_multicomponent))
    delG = np.sum(input_dict['zi'] * np.log(input_dict['zi'] * fc_multicomponent), axis=0) - np.sum(
        input_dict['zi'] * np.log(fc_purecomponent[:, None]), axis=0)

    if fig is None:
        fig = plt.figure()
    plt.plot(input_dict['zi'][0], delG, **plot_kwargs)
    plt.xlabel('$C_1$')
    plt.ylabel('$\Delta_{mix} G_{molar} / RT$')
    # plt.ylim([-0.1, 0.01])
    plt.legend()

    # Get equilibrium positions
    # Put in arbitrary values for overall composition
    Ki, flash_params = run_flash_main(input_dict['P'], input_dict['T'], [0.8, 0.2], input_dict)
    z_id_a = np.abs(flash_params['xi'][0] - input_dict['zi'][0]).argmin()
    z_id_b = np.abs(flash_params['yi'][0] - input_dict['zi'][0]).argmin()
    print(delG[z_id_a])
    plt.plot(flash_params['xi'][0], delG[z_id_a], 'o')
    plt.plot(flash_params['yi'][0], delG[z_id_b], 'o')


    return fig



if __name__ == "__main__":

    plt.style.use('seaborn')
    # # Problem 2a
    # output_path, unit_converter, input_dict = initialize()
    # p2a_main(output_path, unit_converter, input_dict)
    # plt.show()
    #
    # # Problem 2b
    # output_path, unit_converter, input_dict = initialize()
    # p2b_main(output_path, unit_converter, input_dict)
    # plt.show()
    #
    # # Problem 2c
    # output_path, unit_converter, input_dict = initialize(input_filepath='Input_Files/project2_input_file_2c.yml')
    # p2c_main(output_path, unit_converter, input_dict)
    # plt.show()

    # Problem 2d
    fig = p2d_main(input_filepath='Input_Files/project2_input_file_2.yml', color='black', label='BIP=0')
    fig = p2d_main(input_filepath='Input_Files/project2_input_file_2c.yml', fig=fig, color=[0, 153/255, 51/255], label='BIP=0.09')

    plt.show()

