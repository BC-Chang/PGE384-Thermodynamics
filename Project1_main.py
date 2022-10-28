import numpy as np
from eos import cubic_eos, get_molar_volume, _get_compressibility_factor
from solve import get_vapor_pressure, solve_cardanos, pt_flash_nonaqueous, get_roots
import matplotlib.pyplot as plt
from io_utils import read_input, pout
import os
from pr_utils import fugacity_coefficient
import pandas as pd

def main_1():
    # Create a HW output directory if one does not already exist
    output_path = "Project1_Output/Problem1"
    if not os.path.exists(output_path):
        os.makedirs(output_path)


    # Read in input file
    input_dict = read_input(filename="Input_Files/project1_input_file_1.yml")
    print(f"Initial pressure guess using Wilson's correlation = {input_dict['P']: 0.3f} Pa")

    # Calculate vapor pressure at the given temperature
    input_dict, _ = get_vapor_pressure(input_dict)

    # Calculate molar volumes of each phases
    mol_v_l = get_molar_volume(input_dict["zl"], input_dict["T"], input_dict["P"])
    mol_v_v = get_molar_volume(input_dict["zv"], input_dict["T"], input_dict["P"])
    print(f"Molar volume of liquid phase = {mol_v_l: 0.3e} m3/mol")
    print(f"Molar volume of vapor phase = {mol_v_v: 0.3e} m3/mol")

    # Create an array of pressure values to test
    pressures = np.linspace(0.5, 2.0, 20) * 1.E6
    pressures_refined = np.linspace(input_dict["P"] - 1E4, input_dict["P"] + 1E4, 10)
    pressures = np.append(pressures, pressures_refined)
    pressures = np.sort(pressures)

    # Initialize molar volume array
    mol_volumes = np.zeros(len(pressures))

    for i, p in enumerate(pressures):
        mol_volumes[i], _ = pt_flash_nonaqueous(p, input_dict["T"], input_dict["P"], input_dict)

    plt.figure(dpi=300)
    plt.plot(mol_volumes, pressures, '-ob', alpha=0.8, label="T = 313.15K")
    plt.xlabel('Molar Volume [m3/mol]')
    plt.ylabel('Pressure [Pa]')
    plt.legend()
    plt.show()

    # Calculate molar volume
    input_dict["Vc"] = get_molar_volume(0.307, input_dict["Tc"], input_dict["Pc"])
    # #%% =================================================================
    # #Read PV data from excel files
    # df_b = pd.read_excel(f"{output_path}/pv_isotherm.xlsx")
    # df_c = pd.read_excel(f"{output_path}/pv_isotherm_1.xlsx")
    #
    #
    # plt.figure(dpi=400)
    # plt.plot(df_b['mol_volumes_t313'], df_b['P'], '-ob', alpha=0.8, label="T = 313.15K")
    # plt.plot(df_c['mol_volumes_t343'], df_c['P'], '-or', alpha=0.8, label="T = 343.15K")
    # plt.plot(input_dict["Vc"], input_dict["Pc"], '*y', label="Critical Point")
    #
    # plt.xlabel('Molar Volume [m3/mol]')
    # plt.ylabel('Pressure [Pa]')
    # plt.legend()
    # plt.show()

    return

# # Problem 2:
# def main_2():
#     # Create a HW output directory if one does not already exist
#     output_path = "Project1_Output/Problem2"
#     if not os.path.exists(output_path):
#         os.makedirs(output_path)
#
#     # Read in input file
#     input_dict = read_input(filename="Input_Files/project1_input_file_2.yml")
#
#     # Calculate vapor pressure at the given temperature
#     input_dict = get_vapor_pressure(input_dict)
#
#     # Calculate PV Isotherm
#     pressures = np.linspace(0.5, 2.0, 20) * 1.E6
#     pressures_refined = np.linspace(input_dict["P"] - 1E4, input_dict["P"] + 1E4, 10)
#     pressures = np.append(pressures, pressures_refined)
#     pressures = np.sort(pressures)


if __name__ == '__main__':
    # Create a HW output directory if one does not already exist
    output_path = "Project1_Output/Problem2"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Read in input file
    input_dict = read_input(filename="Input_Files/project1_input_file_2.yml")

    # Calculate vapor pressure at the given temperature
    input_dict, _ = get_vapor_pressure(input_dict)

    # Calculate PV Isotherm
    pressures = np.linspace(2.5, 4, 1000) * 1.E6
    # pressures_refined = np.linspace(input_dict["P"] - 1E4, input_dict["P"] + 1E4, 10)
    # pressures = np.append(pressures, pressures_refined)
    # pressures = np.sort(pressures)

    # Initialize molar volume array
    mol_volumes = np.zeros(len(pressures))
    Del_molar_G = np.zeros(len(pressures))
    # pref = 100000  # Reference pressure
    # root = np.linspace(0.07, 0.8, len(pressures))
    # mol_volumes = np.linspace(1E-4, 7E-4, 50)
    # Z = _get_compressibility_factor(input_dict["P"], mol_volumes, input_dict["T"], R=8.3144598)
    # _, eos_params = get_roots(input_dict["P"], input_dict["T"], input_dict)
    # A = eos_params["a"]
    # B = eos_params["b"]
    # Del_molar_G = fugacity_coefficient(Z, eos_params["A"], eos_params["B"]) + np.log(input_dict["P"]) - np.log(pref)
    #
    for i, p in enumerate(pressures):
        roots, eos_params = get_roots(p, input_dict["T"], input_dict)
        mol_volumes[i], _, _ = pt_flash_nonaqueous(p, input_dict["T"], input_dict["P"], input_dict)
        # Calculate fugacity coefficients
        # Calculate the fugacity coefficients
        Del_molar_G[i] = fugacity_coefficient(min(roots), eos_params["A"], eos_params["B"]) + np.log(p) - np.log(input_dict["pref"])


    # plt.figure(dpi=300)
    # plt.plot(mol_volumes, pressures, '-ob', alpha=0.8, label=f"T = {input_dict['T']}")
    # plt.xlabel('Molar Volume [m3/mol]')
    # plt.ylabel('Pressure [Pa]')
    # plt.legend()

    plt.figure(dpi=300)
    plt.plot(mol_volumes, Del_molar_G, 'or', alpha=0.8)
    plt.xlabel('Molar Volume [m3/mol]')
    plt.ylabel('$\Delta G_{molar} / RT$')
    # plt.legend()
    plt.show()




