import numpy as np
from eos import cubic_eos, get_molar_volume
from solve import get_vapor_pressure, solve_cardanos
import matplotlib.pyplot as plt
from io_utils import read_input, pout
import os
import pandas as pd

# Create a HW output directory if one does not already exist
output_path = "HW5_Output"
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Outline for vapor pressure calculation for single component system
# 1. Set Tc, Pc, and w
# 2. Set T and guess a value for P
# 3. Calculate a & b (A & B)
# 4. Solve PR EoS for Z_L and Z_V
# 5. Calculate fugacity coefficients
# 6. Check error in fugacity coefficients
# 7. Update P and loop back to 4

# Read in input file
input_dict = read_input(filename="Input_Files/project1_input_file.yml")
# input_dict["T"] = 70 + 273.15
print(f"Initial pressure guess using Wilson's correlation = {input_dict['P']: 0.3f} Pa")

# Calculate vapor pressure at the given temperature
input_dict = get_vapor_pressure(input_dict)

# Calculate molar volumes of each phases
mol_v_l = get_molar_volume(input_dict["zl"], input_dict["T"], input_dict["P"])
mol_v_v = get_molar_volume(input_dict["zv"], input_dict["T"], input_dict["P"])
print(f"Molar volume of liquid phase = {mol_v_l: 0.3e} m3/mol")
print(f"Molar volume of vapor phase = {mol_v_v: 0.3e} m3/mol")

# Create an array of pressure values
pressures = np.linspace(100, 1E3, 50)
mol_volumes = np.empty_like(pressures)
for i, p in enumerate(pressures):
    # Calculate molar volume at that pressure
    # Calculate roots using Cardano's method
    alpha, beta, gamma, _, _ = cubic_eos(P=p, T=input_dict["T"], eos=input_dict['eos'],
                                        Pc=input_dict["Pc"], Tc=input_dict["Tc"], w=input_dict["w"])
    x1, x2, x3 = solve_cardanos(1, alpha, beta, gamma)
    print(x1, x2, x3)

    mol_volumes[i] = get_molar_volume((x1, x2, x3), input_dict["T"], p)

plt.figure(dpi=300)
plt.plot(mol_volumes, pressures, 'ob', alpha=0.8)
plt.xlabel('Molar Volume [m3/mol]')
plt.ylabel('Pressure [Pa]')

# plt.ylim(0.0, 10.0)
plt.show()

