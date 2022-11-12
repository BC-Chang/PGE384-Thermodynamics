import numpy as np
import pandas as pd
import os
from eos import cubic_eos, get_molar_volume, _get_compressibility_factor
from solve import get_vapor_pressure, solve_cardanos, pt_flash_nonaqueous, get_roots
import matplotlib.pyplot as plt
from io_utils import read_input, pout, redirect_stdout, close_output
import os
from pr_utils import fugacity_coefficient_phi, pr_eos


def main_1():
    # Create a HW output directory if one does not already exist
    output_path = "Project1_Output/Problem1"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Redirect output to the output directory
    redirect_stdout(f"{output_path}/output_file.txt")

    # Read in input file
    input_dict = read_input(filename="Input_Files/project1_input_file_1.yml")

    print("*"*50)
    print(f"Working Fluid: Propane, T = {input_dict['T']}K")
    print(f"Initial pressure guess using Wilson's correlation:")
    print(f"\t{input_dict['P']: 0.3f} Pa")

    # Calculate vapor pressure at the given temperature
    input_dict, _ = get_vapor_pressure(input_dict)
    print(f"Vapor Pressure of Propane at {input_dict['T']}K:")
    print(f"\t{input_dict['P']: 0.3f} Pa")

    # Calculate molar volumes of each phases
    mol_v_l = get_molar_volume(input_dict["zl"], input_dict["T"], input_dict["P"])
    mol_v_v = get_molar_volume(input_dict["zv"], input_dict["T"], input_dict["P"])
    print(f"Molar volume of liquid phase")
    print(f"\t{mol_v_l: 0.3e} m3/mol")
    print(f"Molar volume of vapor phase")
    print(f"\t{mol_v_v: 0.3e} m3/mol")

    # Calculate PV isotherm using PT flash calculations
    # Create an array of pressure values to test
    pressures = np.linspace(input_dict["P"]*0.5, input_dict["P"]*1.15, 30) #* 1.E6
    pressures_refined = np.linspace(input_dict["P"] - 1E4, input_dict["P"] + 1E4, 10)
    pressures = np.append(pressures, pressures_refined)
    pressures = np.sort(pressures)

    # Initialize molar volume array
    mol_volumes = np.zeros(len(pressures))

    for i, p in enumerate(pressures):
        mol_volumes[i], _, _ = pt_flash_nonaqueous(p, input_dict["T"], input_dict["P"], input_dict)

    plt.figure(dpi=300)
    plt.plot(mol_volumes, pressures, '-ob', alpha=0.8, label=f"T = {input_dict['T']}K")
    plt.xlabel('Molar Volume [m3/mol]')
    plt.ylabel('Pressure [Pa]')
    plt.title(f"PV Isotherm at T = {input_dict['T']}K")

    # Calculate critical molar volume
    input_dict["Vc"] = get_molar_volume(0.307, input_dict["Tc"], input_dict["Pc"])
    plt.plot(input_dict["Vc"], input_dict["Pc"], '*y', label=f'Critical Point')

    plt.legend()
    plt.savefig(f'{output_path}/PV_isotherm_{input_dict["T"]}K.png')

    # Write PV data to Excel file
    df = pd.DataFrame({"Molar Volume": mol_volumes, "Pressure": pressures, "Temperature": input_dict["T"]})
    df.to_excel(f'{output_path}/PV_isotherm_{input_dict["T"]}K.xlsx')

    # Close the output file
    print("*"*50)
    close_output(f"{output_path}/output_file.txt")

    """
    =================================================================
    Plot both isotherms together, the below code reads in any Excel
    file in the output directory and plots the PV Isotherm
    =================================================================
    """
    isotherms = os.listdir(f"{output_path}/")
    isotherms = [iso for iso in isotherms if "xlsx" in iso]
    plt.figure(dpi=400)
    for iso in isotherms:
        df_tmp = pd.read_excel(f"{output_path}/{iso}")
        plt.plot(df_tmp["Molar Volume"], df_tmp["Pressure"], label=f"T = {df_tmp['Temperature'][0]} K")

    # Also plot the critical point
    plt.plot(input_dict["Vc"], input_dict["Pc"], '*y', label="Critical Point")
    plt.xlabel('Molar Volume [m3/mol]')
    plt.ylabel('Pressure [Pa]')
    plt.legend()
    plt.savefig(f'{output_path}/Combined_PV_isotherm.png')
    # plt.show()

    plt.show()
    return

# Problem 2:
def main_2():
    # Create a HW output directory if one does not already exist
    output_path = "Project1_Output/Problem2"
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    # Redirect stdout to the output file
    redirect_stdout(f"{output_path}/output_file.txt")

    # Read in input file
    input_dict = read_input(filename="Input_Files/project1_input_file_2.yml")

    """
    =======================================================================
    Problem 2b: Plot \Delta G/RT for 30F Isotherm.
    Identify stable, unstable, and metastable regions
    Delete regions that are not stable to get single-phase regions
    =======================================================================
    """

    # Calculate vapor pressure at the given temperature
    input_dict, _ = get_vapor_pressure(input_dict)
    print("*" * 50)
    print(f"Working Fluid: CO2, T = {input_dict['T']}K")
    print("Vapor pressure: ")
    print(f"\t{input_dict['P']} Pa")

    # Get liquid and vapor molar volumes at the given vapor pressure
    mol_v_l = get_molar_volume(input_dict["zl"], input_dict["T"], input_dict["P"])
    mol_v_v = get_molar_volume(input_dict["zv"], input_dict["T"], input_dict["P"])
    print("Equilibrium molar volume of liquid phase:")
    print(f"\t{mol_v_l} m3/mol")
    print("Equilibrium molar volume of vapor phase:")
    print(f"\t{mol_v_v} m3/mol")

    # Calculate pressure at specified molar volumes
    mol_v = np.arange(mol_v_l - 0.3 * mol_v_l, mol_v_v + 0.3 * mol_v_v, 1E-8)
    # Get dimensional attraction and covolume parameters
    eos_params = cubic_eos(1, input_dict["T"], input_dict["Pc"], input_dict["Tc"], input_dict["w"])
    # Calculate pressures using PR EoS
    P = pr_eos(T=input_dict["T"], mol_v=mol_v, a=eos_params["a"], b=eos_params["b"])
    Z = _get_compressibility_factor(P, mol_v, input_dict["T"])

    # Calculate fugacity coefficient
    phi = fugacity_coefficient_phi(Z, P, input_dict["T"], eos_params["a"], eos_params["b"])

    # Calculate delta molar Gibbs Free Energy
    del_G = np.log(phi * P / input_dict["pref"])

    # Plot isotherm with Delta molar Gibbs free energy
    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Molar Volume (m3/mol)')
    ax1.set_ylabel("Pressure (Pa)", color=color)
    ax1.plot(mol_v, P, color=color, label="Pressure")
    ax1.plot([mol_v[0], mol_v[-1]], [input_dict["P"], input_dict["P"]], "--k", label="Vapor Pressure")
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('$\Delta G_{molar} / RT$', color=color)
    ax2.plot(mol_v, del_G, color=color)
    ax2.tick_params(axis='y', labelcolor=color)

    fig.tight_layout()

    # Find stable, metastable and unstable regions of \Delta G curve
    # First deriviative of del G term
    del_G_deriv = np.diff(del_G)
    # Find where derivative changes sign
    zero_del_G_deriv = np.argwhere(np.diff(np.sign(del_G_deriv)) != 0)

    ax2.plot(mol_v[zero_del_G_deriv], del_G[zero_del_G_deriv], '^y',
             markersize=6, label="Minima and Maxima")

    ax2.legend()
    plt.title('PV Isotherm and Molar Gibbs Free Energy')
    fig.savefig(f"{output_path}/PV_Molar_Gibbs_{input_dict['T']}K.png", dpi=300)
    plt.show()


    # Plot stable, metastable and unstable regions
    mol_v_v_id = np.amax(np.where(np.isclose(mol_v, mol_v_v, atol=1e-7)))
    mol_v_l_id = np.amin(np.where(np.isclose(mol_v, mol_v_l, atol=1e-7)))

    stable_id_l = np.arange(0, mol_v_l_id + 1)
    stable_id_v = np.arange(mol_v_v_id, len(mol_v))
    metastable_id_1 = np.arange(mol_v_l_id + 1, zero_del_G_deriv[0])
    metastable_id_2 = np.arange(zero_del_G_deriv[-1], mol_v_v_id)
    unstable_id = np.arange(zero_del_G_deriv[0], zero_del_G_deriv[-1])

    plt.figure()
    plt.plot(mol_v[stable_id_l], del_G[stable_id_l], '-g')
    plt.plot(mol_v[stable_id_v], del_G[stable_id_v], '-g', label="Stable region")

    plt.plot(mol_v[metastable_id_1], del_G[metastable_id_1], '-y')
    plt.plot(mol_v[metastable_id_2], del_G[metastable_id_2], '-y', label="Metastable")

    plt.plot(mol_v[unstable_id], del_G[unstable_id], '-r', label="Unstable region")
    plt.legend()
    plt.xlabel('Molar Volume [m3/mol]')
    plt.ylabel('$\Delta G_{molar} / RT$')
    plt.savefig(f'{output_path}/stability_{input_dict["T"]}K.png')
    plt.show()

    # Plot only stable regions
    plt.figure()
    stable_id = np.append(stable_id_l, stable_id_v)
    plt.plot(mol_v, del_G, '-r', label="Total curve")
    plt.plot(mol_v[stable_id], del_G[stable_id], '-g', label="Stable curve")

    plt.xlabel('Molar Volume [m3/mol]')
    plt.ylabel('$\Delta G_{molar} / RT$')
    plt.title('Stable $\Delta G_{molar} / RT$ curve')
    plt.legend()
    plt.savefig(f'{output_path}/stable_regions_{input_dict["T"]}K.png')


    """
    =======================================================================
    Problem 2c: Determine \Delta G/RT at equilibrium molar volumes.
    =======================================================================
    """
    print("Delta G/RT at molar volume of liquid at equilibrium: ")
    print(f"\t{del_G[stable_id_l[-1]]}")

    print("Delta G/RT at molar volume of vapor at equilibrium: ")
    print(f"\t{del_G[stable_id_v[0]]}")

    print("Molar Gibbs free energy of liquid and vapor phases are approximately equal at equilibrium.")


    """
    =======================================================================
    Problem 2d: Is CO2 a liquid or vapor at 25 bars and 30F?
    =======================================================================
    """
    p = 2.5E6  # Pressure in Pa

    # Calculate vapor pressure at the given temperature
    input_dict, _ = get_vapor_pressure(input_dict)
    print(f"Z_l = \n\t{input_dict['zl']}")
    print(f"Z_v = \n\t{input_dict['zv']}")

    mol_volumes, eos_params, desired_root = pt_flash_nonaqueous(p, input_dict["T"], input_dict["P"], input_dict)
    print(f"Equilibrium molar volume: \n\t{mol_volumes} m3/mol")

    # Calculate fugacity coefficient of liquid phase
    phi_l = fugacity_coefficient_phi(input_dict["zl"], p, input_dict["T"], eos_params["a"], eos_params["b"])
    # Calculate molar Gibbs Free Energy
    del_G_l = np.log(phi_l * p / input_dict["pref"])

    # Calculate fugacity coefficient of vapor phase
    phi_v = fugacity_coefficient_phi(input_dict["zv"], p, input_dict["T"], eos_params["a"], eos_params["b"])
    # Calculate delta molar Gibbs Free Energy
    del_G_v = np.log(phi_v * p / input_dict["pref"])
    print(f"Molar Gibbs Free Energy of Liquid Phase: \n\t{del_G_l}")
    print(f"Molar Gibbs Free Energy of Vapor Phase:  \n\t{del_G_v}")

    if del_G_v < del_G_l:
        phase = "vapor"
    elif del_G_l < del_G_v:
        phase = "liquid"
    else:
        phase = "unknown"

    print(f"CO2 is a {phase} phase at {p} Pa and {input_dict['T']}K")
    print(f"Compressibility factor of {phase} phase results in lower Gibbs Free Energy")

    plt.show()

    close_output(f"{output_path}/output_file.txt")

if __name__ == '__main__':
    # Uncomment the following line to run problem 1
    main_1()

    # Uncomment the following line to run problem 2
    # main_2()
