import numpy as np
import os
from eos import cubic_eos, get_molar_volume, _get_compressibility_factor, _peng_robinson_ab, _peng_robinson_cubic
from singlecomponent_utils import get_vapor_pressure
import matplotlib.pyplot as plt
from io_utils import read_input, pout, redirect_stdout, close_output
import os
from multicomponent_utils import get_bubble_point, get_dew_point, get_phase_compositions, get_purecomponent_a_b, \
    get_multicomponent_A_B, get_phase_fugacity, check_TBD_sign
from multicomponent_solve import two_phase_flash, single_phase_stability
from pr_utils import fugacity_coefficient_multicomponent
from itertools import product
from solve import solve_cardanos, _get_real_roots
from unit_conversions import Unit_Converter

# Two-Phase PT Flash

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
# Convert bar to Pa
input_dict['P'] = unit_converter.bar_to_Pa(input_dict['P'])
input_dict['Pc'] = unit_converter.bar_to_Pa(input_dict['Pc'])

print("*" * 50)
print("Project 2 Problem 1: ")
# Calculate compressibility factors and fugacity coefficients for phase z
a_ii, b_i = get_purecomponent_a_b(input_dict['Pc'], input_dict['Tc'], input_dict['T'], input_dict['w'], input_dict['Nc'])
fc, f, root = get_phase_fugacity(a_ii, b_i, input_dict['zi'], input_dict)

# Initialize Ki with wilson's correlation
Pri = input_dict["Pc"] / input_dict["P"]
Tri = input_dict["Tc"] / input_dict["T"]
Ki = Pri * np.exp((5.373 * (1 + input_dict['w']) * (1 - Tri)))

# Initial guesses for Xi assuming emerging phase is vapor-like, then liquid-like
Xi = np.array([Ki*input_dict["zi"], input_dict["zi"]/Ki])
print("Initial guesses for Xi: ", Xi)
stable_single_phase = single_phase_stability(Xi, a_ii, b_i, fc, input_dict)

if stable_single_phase:
    print("Single phase assumed stable. Perform single phase flash")
    #TODO perform single phase flash
else:
    Ki, flash_params = two_phase_flash(input_dict)


# Close the output file
# close_output(f"{output_path}/output_file.txt")