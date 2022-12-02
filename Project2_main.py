import numpy as np
import os

from astropy.wcs.docstrings import mix

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
import ternary
from ternary_plot_utils import initialize_ternary_diagram, compositions_to_coords, get_axis_intercept


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
def p2a_main(output_path, unit_converter, input_dict, parse_P_max=56.4):
    # Number of points to test
    n_pts = 100

    # Create a range of pressures to test
    P = np.linspace(7.91, parse_P_max, n_pts)
    P = unit_converter.bar_to_Pa(P)

    xi_flash = np.zeros_like(P)
    yi_flash = np.zeros_like(P)

    print("*" * 50)
    for i, Pi in enumerate(P):
        Ki, flash_params = two_phase_flash(input_dict, P=Pi, zi=input_dict['zi'])#Pi, input_dict['T'], input_dict['zi'], input_dict)
        if flash_params is None:
            flash_params = {'xi': [0, 0], 'yi': [0, 0]}
        xi_flash[i] = flash_params['xi'][0]
        yi_flash[i] = flash_params['yi'][0]

    P_plot = unit_converter.Pa_to_bar(P)

    fig = plt.figure()
    plt.plot(xi_flash[xi_flash > 0], P_plot[xi_flash > 0], '-b', label='Flash calculation bubble point')
    plt.plot(yi_flash[yi_flash > 0], P_plot[yi_flash > 0], '-r', label='Flash calculation dew point')
    plt.plot([xi_flash[n_pts//4::n_pts//4], yi_flash[n_pts//4::n_pts//4]], [P_plot[n_pts//4::n_pts//4], P_plot[n_pts//4::n_pts//4]], 'o--k')
    plt.plot(xi_flash[xi_flash > 0][-1], P_plot[xi_flash > 0][-1], '*m', markersize=10, label='Flash Critical point')
    plt.xlabel('Mole fraction of Component 1')
    plt.ylabel('Pressure (bar)')
    plt.ylim([0, 75])
    plt.xlim([0, 1.1])
    plt.legend()
    return fig


def p2b_main(output_path, unit_converter, input_dict, parse_P_max=56.4):

    fig = p2a_main(output_path, unit_converter, input_dict, parse_P_max)

    joyce_P = np.array([7.91, 14.80, 25.14, 35.49, 45.83, 56.17, 62.38, 64.51])
    joyce_xA = 1 - np.array([0.765, 0.626, 0.449, 0.330, 0.240, 0.160, 0.116, 0.0702])
    joyce_yA = 1 - np.array([0.00466, 0.00347, 0.00364, 0.00493, 0.00812, 0.0157, 0.0306, 0.0702])
    plt.plot(joyce_xA[:-1], joyce_P[:-1], '^b', label='Joyce et al. (2000) liquid compositions')
    plt.plot(joyce_yA[:-1], joyce_P[:-1], '^r', label='Joyce et al. (2000) vapor compositions')
    plt.plot(joyce_xA[-1], joyce_P[-1], '*y', markersize=13, label='Joyce et al. (2000) critical point')
    plt.legend()

    return fig


def p2c_main(output_path, unit_converter, input_dict):
    fig = p2b_main(output_path, unit_converter, input_dict, parse_P_max=63.5)

    return fig


def p2d_main(input_filepath, fig=None, *pt_args, **curve_kwargs):
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
    plt.plot(input_dict['zi'][0], delG, **curve_kwargs)
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
    plt.plot(flash_params['xi'][0], delG[z_id_a], 'o', color=curve_kwargs['color'])
    plt.plot(flash_params['yi'][0], delG[z_id_b], 'o', color=curve_kwargs['color'])

    return fig

def p3a_main(output_path, unit_converter, input_dict):
    # Problem 3
    # Create a HW output directory if one does not already exist
    output_path = "Project2_Output/Problem3"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Initialize a unit converter object
    unit_converter = Unit_Converter()

    # Read in input file
    input_dict = read_input(filename='Input_Files/project2_input_file_3.yml')
    # Preprocess system pressure and critical pressures
    input_dict['P'] = unit_converter.psi_to_Pa(input_dict['P'])
    input_dict['Pc'][0] = unit_converter.psi_to_Pa(input_dict['Pc'][0])
    input_dict['T'] = unit_converter.F_to_K(input_dict['T'])
    input_dict['Pc'][1:] = unit_converter.bar_to_Pa(input_dict['Pc'][1:])
    input_dict['Tc'][0] = unit_converter.F_to_K(input_dict['Tc'][0])

    _, tax = initialize_ternary_diagram()

    zC4 = np.linspace(0., 0.2, 100)
    # x_coords = np.zeros((len(zC4), 3))
    # y_coords = np.zeros((len(zCO), 3))
    # zi = np.array([[0.75, 0.1, 0.150], [0.7, 0.05, 0.25], [0.8, 0.1, 0.1], [0.850, 0.05, 0.1], [0.9, 0.05, 0.05]])
    for i, z in enumerate(zC4):
        # calculate c4 composition
        zCO = (1 - z) / 1.142  # / 1.5
        zi = [zCO, z, 1 - z - zCO]

        Ki, flash_params = two_phase_flash(input_dict, zi=zi)
        if flash_params is None:
            continue
        x_coords = compositions_to_coords(flash_params['xi'])
        y_coords = compositions_to_coords(flash_params['yi'])
        tax.scatter([x_coords, y_coords], color='black', marker='.')
    tax.scatter([x_coords], color='magenta', marker='*')
    print('Estimated Critical Point:')
    print(f'\t\t{[x_coords[1], x_coords[0], x_coords[2]]}')

    # Get some tie lines
    zi = np.array([[0.75, 0.15, 0.1], [0.8, 0.1, 0.1], [0.85, 0.1, 0.05], [0.87, 0.07, 0.06], [0.92, 0.03, 0.05]])
    line_styles = ['-', '--', '-.', ':', '-']
    colors = ['black', 'red', 'green', 'blue', 'magenta']
    for i, z in enumerate(zi):
        Ki, flash_params = run_flash_main(input_dict['P'], input_dict['T'], z, input_dict)
        print(Ki, flash_params)
        if flash_params is None:
            continue
        x_coords = compositions_to_coords(flash_params['xi'])
        y_coords = compositions_to_coords(flash_params['yi'])
        tax.line(x_coords, y_coords, color=colors[i], linestyle=line_styles[i],
                 label=f'Tie line {z[0]:.2f}/{z[1]:.2f}/{z[2]:.2f}')

    tax.legend()
    tax.show()

    return tax

def p3b_main(output_path, unit_converter, input_dict):
    # Problem 3
    # Create a HW output directory if one does not already exist
    output_path = "Project2_Output/Problem3"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Initialize a unit converter object
    unit_converter = Unit_Converter()

    # Read in input file
    input_dict = read_input(filename='Input_Files/project2_input_file_3.yml')
    # Preprocess system pressure and critical pressures
    input_dict['P'] = unit_converter.psi_to_Pa(input_dict['P'])
    input_dict['Pc'][0] = unit_converter.psi_to_Pa(input_dict['Pc'][0])
    input_dict['T'] = unit_converter.F_to_K(input_dict['T'])
    input_dict['Pc'][1:] = unit_converter.bar_to_Pa(input_dict['Pc'][1:])
    input_dict['Tc'][0] = unit_converter.F_to_K(input_dict['Tc'][0])

    tax = p3a_main(output_path, unit_converter, input_dict)

    # Problem 3b
    Ki, flash_params = run_flash_main(input_dict['P'], input_dict['T'], input_dict['zi'], input_dict)
    x_coords = compositions_to_coords(flash_params['xi'])
    y_coords = compositions_to_coords(flash_params['yi'])
    composition_a, composition_b = get_axis_intercept(flash_params['xi'], flash_params['yi'])
    right_axis_intercept = compositions_to_coords(composition_a)
    bottom_axis_intercept = compositions_to_coords(composition_b)

    tax.line(right_axis_intercept, bottom_axis_intercept, color='black', linestyle='-',
             label=f'Tie line {input_dict["zi"][0]:.2f}/{input_dict["zi"][1]:.2f}/{input_dict["zi"][2]:.2f}')

    tax.legend()
    tax.show()

    # Calculate vapor pressure at the given temperature assuming pure fluid
    pure_fluid_a_b = [_peng_robinson_ab(input_dict['Pc'][i], input_dict['Tc'][i],
                                        input_dict['T'], input_dict['w'][i]) for i in range(input_dict["Nc"])]
    a_ii = np.array([pf_a_b[0] for pf_a_b in pure_fluid_a_b])
    b_ii = np.array([pf_a_b[1] for pf_a_b in pure_fluid_a_b])

    # Initialize mixing ratio
    mixing_ratio = np.linspace(0, 1, 100000)
    xi = np.array([mixing_ratio * composition_a[0],
                   composition_b[1] + (composition_a[1] - composition_b[1]) * mixing_ratio,
                   composition_b[2] * (1 - mixing_ratio)])
    fc_multicomponent = np.zeros_like(xi)  # ((mixing_ratio.shape[0], 3))
    GFE = np.zeros_like(mixing_ratio)

    for count in range(len(xi[0])):
        xi_tmp_list = [xi[nc][count] for nc in range(len(xi))]
        # Create a matrix for calculating a_mix
        a_mix = 0
        a_matrix = np.eye(input_dict['Nc']) * a_ii
        a_ij = np.eye(input_dict['Nc']) * a_ii

        for i, j in product(range(input_dict["Nc"]), range(input_dict["Nc"])):
            if i != j:
                a_ij[i, j] = np.sqrt(a_matrix[i, i] * a_matrix[j, j]) * (1 - input_dict['K_ij'][i, j])
            a_mix += a_ij[i, j] * xi[i][count] * xi[j][count]

        # Get attraction and covolume parameters of mixture
        b_mix = np.sum(xi_tmp_list * b_ii)

        alpha, beta, gamma, A_mix, B_mix = _peng_robinson_cubic(input_dict["P"], input_dict["T"], a_mix, b_mix)

        x1, x2, x3 = solve_cardanos(1, alpha, beta, gamma)
        roots = _get_real_roots(x1, x2, x3)

        roots = np.unique(roots)
        # Get fugacity coefficient
        fc_tmp, _ = fugacity_coefficient_multicomponent(roots[0], a_mix, b_mix, input_dict['P'], input_dict['T'],
                                                        a_ij, b_ii, xi_tmp_list)
        fc_multicomponent[:, count] = np.exp(fc_tmp)

    GFE = np.sum(xi * np.log(xi * fc_multicomponent), axis=0)

    # Calculate gradient of GFE
    dGFE = np.gradient(GFE, mixing_ratio)

    fig, ax1 = plt.subplots(dpi=400)
    ax1.set_xlabel(r'Mixing Ratio ($\beta$)')
    ax1.set_ylabel(r'$G\underbar_R / RT$')
    ax1.plot(mixing_ratio, GFE, '-b', label='Gr/RT')
    ax1.grid(False)

    ax2 = ax1.twinx()
    ax2.set_ylabel(r'$dG\underbar_R / RT$')
    ax2.plot(mixing_ratio, dGFE, '-.r', label='dGr/RT')
    ax2.grid(False)

    # plt.legend()
    # plt.grid('off')

    # Get gradient of liquid and vapor phase compositions
    liquid_phase_beta = flash_params['xi'][0] / composition_a[0]
    vapor_phase_beta = flash_params['yi'][0] / composition_a[0]
    ax2.plot([liquid_phase_beta, liquid_phase_beta], [0, 10], '-k', label='Equilibrium liquid phase')
    ax2.plot([vapor_phase_beta, vapor_phase_beta], [0, 10], '--k', label='Equilibrium vapor phase')
    print(liquid_phase_beta, vapor_phase_beta)

    # Find index of beta values and return dGFE value at these indices
    x_beta_ind = np.abs(liquid_phase_beta - mixing_ratio).argmin()
    y_beta_ind = np.abs(vapor_phase_beta - mixing_ratio).argmin()

    print(f"Liquid phase Gibbs gradient = {dGFE[x_beta_ind]}")
    print(f"Vapor phase Gibbs gradient = {dGFE[y_beta_ind]}")

    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, loc=0)

    return fig


if __name__ == "__main__":

    plt.style.use('seaborn')
    # # Problem 2a
    # output_path, unit_converter, input_dict = initialize()
    # p2a_main(output_path, unit_converter, input_dict)
    # plt.show()
    # #
    # # Problem 2b
    # output_path, unit_converter, input_dict = initialize()
    # p2b_main(output_path, unit_converter, input_dict)
    # plt.show()
    #
    # # Problem 2c
    # output_path, unit_converter, input_dict = initialize(input_filepath='Input_Files/project2_input_file_2c.yml')
    # p2c_main(output_path, unit_converter, input_dict)
    # plt.show()
    #
    # # Problem 2d
    # fig = p2d_main('Input_Files/project2_input_file_2.yml', color='black', label='BIP=0')
    # fig = p2d_main('Input_Files/project2_input_file_2c.yml', fig, color=[0, 153/255, 51/255], label='BIP=0.09')
    # plt.show()
    #
    # Problem 3a
    # p3a_main(output_path, unit_converter, input_dict)
    # tax.show()
    #
    # Problem 3b
    # p3b_main(output_path, unit_converter, input_dict)
    # tax.show()
    # plt.show()
