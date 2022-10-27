import numpy as np
from eos import cubic_eos, get_molar_volume
from solve import get_vapor_pressure, solve_cardanos, pt_flash_nonaqueous
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

# Create an array of pressure values to test
pressures = np.linspace(1.2, 1.5, 20) * 1.E6
pressures_refined = np.linspace(input_dict["P"] - 1E4, input_dict["P"] + 1E4, 10)
pressures = np.append(pressures, pressures_refined)
pressures = np.sort(pressures)

# Initialize molar volume array
mol_volumes = np.zeros(len(pressures))

for i, p in enumerate(pressures):

    mol_volumes[i] = pt_flash_nonaqueous(p, input_dict["T"], input_dict["P"], input_dict)

plt.figure(dpi=300)
plt.plot(mol_volumes, pressures, '-ob', alpha=0.8, label="T = 313.15K")
plt.xlabel('Molar Volume [m3/mol]')
plt.ylabel('Pressure [Pa]')

# Change Temperature to 70C
input_dict["T"] = 70. + 273.15
mol_volumes = np.zeros(len(pressures))
for i, p in enumerate(pressures):
    mol_volumes[i] = pt_flash_nonaqueous(p, input_dict["T"], input_dict["P"], input_dict)
# # plt.figure(dpi=300)
plt.plot(mol_volumes, pressures, '-or', alpha=0.8, label='T = 343.15K')
# plt.xlabel('Molar Volume [m3/mol]')
# plt.ylabel('Pressure [Pa]')
# plt.ylim(0.0, 10.0)
plt.legend()
plt.show()

