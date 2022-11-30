import numpy as np
import os
from eos import cubic_eos, get_molar_volume, _get_compressibility_factor, _peng_robinson_ab, _peng_robinson_cubic
from singlecomponent_utils import get_vapor_pressure
import matplotlib.pyplot as plt
from io_utils import read_input, pout, redirect_stdout, close_output
import os
from multicomponent_utils import get_bubble_point, get_dew_point, get_phase_compositions, get_purecomponent_a_b, \
    get_multicomponent_A_B, get_phase_fugacity, check_TBD_sign, get_total_composition
from multicomponent_solve import two_phase_flash, single_phase_stability
from pr_utils import fugacity_coefficient_multicomponent
from itertools import product
from solve import solve_cardanos, _get_real_roots
from unit_conversions import Unit_Converter

# Two-Phase PT Flash
def p2_main():
    # Create a HW output directory if one does not already exist
    output_path = "Project2_Output/Problem_2"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Initialize a unit converter object
    unit_converter = Unit_Converter()

    # Redirect output to the output directory
    # redirect_stdout(f"{output_path}/output_file.txt")

    # Read in input file
    input_dict = read_input(filename="Input_Files/project2_input_file_2.yml")
    n_pts = 100

    # Create a range of pressures to test
    # P = np.linspace(5, 75, n_pts)

    # Convert bar to Pa
    input_dict['P'] = unit_converter.bar_to_Pa(input_dict['P'])
    input_dict['Pc'] = unit_converter.bar_to_Pa(input_dict['Pc'])
    # P = unit_converter.bar_to_Pa(P)
    P = np.linspace(7.91, 54, n_pts) #input_dict['P'].copy()
    P = unit_converter.bar_to_Pa(P)

    xi_flash = np.zeros_like(P)
    yi_flash = np.zeros_like(P)
    print("*" * 50)
    for i, Pi in enumerate(P):
        input_dict['P'] = Pi
        # input_dict['P'] = P[i]
        # input_dict['zi'] = zi[:, i]
        # Calculate compressibility factors and fugacity coefficients for phase z
        a_ii, b_i = get_purecomponent_a_b(input_dict['Pc'], input_dict['Tc'], input_dict['T'], input_dict['w'], input_dict['Nc'])
        fc, f, root = get_phase_fugacity(a_ii, b_i, input_dict['zi'], input_dict)

        # Initialize Ki with wilson's correlation
        Pri = input_dict["Pc"] / input_dict["P"]
        Tri = input_dict["Tc"] / input_dict["T"]
        Ki = Pri * np.exp((5.373 * (1 + input_dict['w']) * (1 - Tri)))

        # Initial guesses for Xi assuming emerging phase is vapor-like, then liquid-like
        Xi = np.array([Ki*input_dict['zi'], input_dict['zi']/Ki])
        # print("Initial guesses for Xi: ", Xi)
        stable_single_phase = single_phase_stability(Xi, a_ii, b_i, fc, input_dict['zi'], input_dict)

        if stable_single_phase:
            print(P[i], "Single phase assumed stable.")
            #TODO perform single phase flash

        else:
            Ki, flash_params = two_phase_flash(input_dict, zi=input_dict['zi'])
            xi_flash[i] = flash_params['xi'][0]
            yi_flash[i] = flash_params['yi'][0]

    P_plot = unit_converter.Pa_to_bar(P)
    plt.style.use('seaborn')
    plt.figure()
    plt.plot(xi_flash, P_plot, 'o-b', label='Bubble point curve')
    plt.plot(yi_flash, P_plot, 'o-r', label='Dew point curve')
    plt.plot([xi_flash[n_pts//4::n_pts//4], yi_flash[n_pts//4::n_pts//4]], [P_plot[n_pts//4::n_pts//4], P_plot[n_pts//4::n_pts//4]], 'o--k')
    plt.xlabel('Mole fraction of Component 1')
    plt.ylabel('Pressure (bar)')
    plt.ylim([0, 75])
    plt.xlim([0, 1.1])
    # plt.plot(np.concatenate((zi_flash, zi_flash[::-1])), np.concatenate((unit_converter.Pa_to_bar(bp_flash), unit_converter.Pa_to_bar(dew_flash)[::-1])), 'ob')
    # plt.plot(zi_flash[2:-1], unit_converter.Pa_to_bar(dew_flash)[2:-1], 'or')
    plt.legend()
    plt.show()

def p3_main():
    # Create a HW output directory if one does not already exist
    output_path = "Project2_Output/Problem_2"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Initialize a unit converter object
    unit_converter = Unit_Converter()

    # Redirect output to the output directory
    # redirect_stdout(f"{output_path}/output_file.txt")

    # Read in input file
    input_dict = read_input(filename="Input_Files/project2_input_file_2.yml")
    n_pts = 1000

    # Create a range of pressures to test
    P = np.linspace(5, 75, n_pts)

    # Convert bar to Pa
    input_dict['P'] = unit_converter.bar_to_Pa(input_dict['P'])
    input_dict['Pc'] = unit_converter.bar_to_Pa(input_dict['Pc'])
    P = unit_converter.bar_to_Pa(P)
    # P = input_dict['P'].copy()

    xi_flash = np.zeros_like(P)
    yi_flash = np.zeros_like(P)
    print("*" * 50)
    for i, Pi in enumerate(P):
        input_dict['P'] = Pi
        # Calculate compressibility factors and fugacity coefficients for phase z
        a_ii, b_i = get_purecomponent_a_b(input_dict['Pc'], input_dict['Tc'], input_dict['T'], input_dict['w'], input_dict['Nc'])
        fc, f, root = get_phase_fugacity(a_ii, b_i, input_dict['zi'], input_dict)

        # Initialize Ki with wilson's correlation
        Pri = input_dict["Pc"] / input_dict["P"]
        Tri = input_dict["Tc"] / input_dict["T"]
        Ki = Pri * np.exp((5.373 * (1 + input_dict['w']) * (1 - Tri)))

        # Initial guesses for Xi assuming emerging phase is vapor-like, then liquid-like
        Xi = np.array([Ki*input_dict['zi'], input_dict['zi']/Ki])
        # print("Initial guesses for Xi: ", Xi)
        stable_single_phase = single_phase_stability(Xi, a_ii, b_i, fc, input_dict['zi'], input_dict)

        if stable_single_phase:
            print(P[i], "Single phase assumed stable.")
            #TODO perform single phase flash

        else:
            Ki, flash_params = two_phase_flash(input_dict, zi=input_dict['zi'])
            xi_flash[i] = flash_params['xi'][0]
            yi_flash[i] = flash_params['yi'][0]

    P_plot = unit_converter.Pa_to_bar(P)
    plt.style.use('seaborn')
    fig = plt.figure()
    plt.plot(xi_flash, P_plot, '-b', label='Bubble point curve')
    plt.plot(yi_flash, P_plot, '-r', label='Dew point curve')
    plt.plot([xi_flash[n_pts//4::n_pts//4], yi_flash[n_pts//4::n_pts//4]], [P_plot[n_pts//4::n_pts//4], P_plot[n_pts//4::n_pts//4]], 'o--k')
    plt.xlabel('Mole fraction of Component 1')
    plt.ylabel('Pressure (bar)')
    plt.ylim([0, 75])
    plt.xlim([0, 1.1])


    joyce_P = np.array([7.91, 14.80, 25.14, 35.49, 45.83, 56.17, 62.38, 64.51])
    joyce_xA = 1 - np.array([0.765, 0.626, 0.449, 0.330, 0.240, 0.160, 0.116, 0.0702])
    joyce_yA = 1 - np.array([0.00466, 0.00347, 0.00364, 0.00493, 0.00812, 0.0157, 0.0306, 0.0702])
    plt.plot(joyce_xA[:-1], joyce_P[:-1], '^b', label='Joyce et al. (2000) liquid compositions')
    plt.plot(joyce_yA[:-1], joyce_P[:-1], '^r', label='Joyce et al. (2000) vapor compositions')
    plt.plot(joyce_xA[-1], joyce_P[-1], '*y', markersize=13, label='Joyce et al. (2000) critical point')
    plt.legend()
    plt.show()

def p4_main():
        # Create a HW output directory if one does not already exist
        output_path = "Project2_Output/Problem_2"
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        # Initialize a unit converter object
        unit_converter = Unit_Converter()

        # Redirect output to the output directory
        # redirect_stdout(f"{output_path}/output_file.txt")

        # Read in input file
        input_dict = read_input(filename="Input_Files/project2_input_file_2c.yml")
        n_pts = 100

        # Create a range of pressures to test
        P = np.linspace(5, 75, n_pts)

        # Convert bar to Pa
        input_dict['P'] = unit_converter.bar_to_Pa(input_dict['P'])
        input_dict['Pc'] = unit_converter.bar_to_Pa(input_dict['Pc'])
        P = unit_converter.bar_to_Pa(P)
        # P = input_dict['P'].copy()

        xi_flash = np.zeros_like(P)
        yi_flash = np.zeros_like(P)
        print("*" * 50)
        for i, Pi in enumerate(P):
            input_dict['P'] = Pi
            # Calculate compressibility factors and fugacity coefficients for phase z
            a_ii, b_i = get_purecomponent_a_b(input_dict['Pc'], input_dict['Tc'], input_dict['T'], input_dict['w'],
                                              input_dict['Nc'])
            fc, f, root = get_phase_fugacity(a_ii, b_i, input_dict['zi'], input_dict)

            # Initialize Ki with wilson's correlation
            Pri = input_dict["Pc"] / input_dict["P"]
            Tri = input_dict["Tc"] / input_dict["T"]
            Ki = Pri * np.exp((5.373 * (1 + input_dict['w']) * (1 - Tri)))

            # Initial guesses for Xi assuming emerging phase is vapor-like, then liquid-like
            Xi = np.array([Ki * input_dict['zi'], input_dict['zi'] / Ki])
            # print("Initial guesses for Xi: ", Xi)
            stable_single_phase = single_phase_stability(Xi, a_ii, b_i, fc, input_dict['zi'], input_dict)

            if stable_single_phase:
                print(P[i], "Single phase assumed stable.")
                # TODO perform single phase flash

            else:
                Ki, flash_params = two_phase_flash(input_dict, zi=input_dict['zi'])
                xi_flash[i] = flash_params['xi'][0]
                yi_flash[i] = flash_params['yi'][0]

        P_plot = unit_converter.Pa_to_bar(P)
        plt.style.use('seaborn')
        fig = plt.figure()
        plt.plot(xi_flash[xi_flash > 0], P_plot[xi_flash > 0], '-b', label='Bubble point curve')
        plt.plot(yi_flash[yi_flash > 0], P_plot[yi_flash > 0], '-r', label='Dew point curve')
        # plt.plot([xi_flash[n_pts // 4::n_pts // 4], yi_flash[n_pts // 4::n_pts // 4]],
        #          [P_plot[n_pts // 4::n_pts // 4], P_plot[n_pts // 4::n_pts // 4]], 'o--k')
        plt.xlabel('Mole fraction of Component 1')
        plt.ylabel('Pressure (bar)')
        plt.ylim([0, 75])
        plt.xlim([0, 1.1])

        joyce_P = np.array([7.91, 14.80, 25.14, 35.49, 45.83, 56.17, 62.38, 64.51])
        joyce_xA = 1 - np.array([0.765, 0.626, 0.449, 0.330, 0.240, 0.160, 0.116, 0.0702])
        joyce_yA = 1 - np.array([0.00466, 0.00347, 0.00364, 0.00493, 0.00812, 0.0157, 0.0306, 0.0702])
        plt.plot(joyce_xA[:-1], joyce_P[:-1], '^b', label='Joyce et al. (2000) liquid compositions')
        plt.plot(joyce_yA[:-1], joyce_P[:-1], '^r', label='Joyce et al. (2000) vapor compositions')
        plt.plot(joyce_xA[-1], joyce_P[-1], '*y', markersize=13, label='Joyce et al. (2000) critical point')
        plt.legend()
        plt.show()


if __name__ == "__main__":
    # p2_main()
    p3_main()
    # p4_main()