import numpy as np
import pandas as pd
import os
from eos import cubic_eos, get_molar_volume, _get_compressibility_factor
from singlecomponent_utils import get_vapor_pressure
import matplotlib.pyplot as plt
from io_utils import read_input, pout, redirect_stdout, close_output
import os
from multicomponent_utils import raoults_law
from pr_utils import fugacity_coefficient_phi, pr_eos

# Create a HW output directory if one does not already exist
output_path = "HW8_Output/"
if not os.path.exists(output_path):
    os.makedirs(output_path)

# Redirect output to the output directory
# redirect_stdout(f"{output_path}/output_file.txt")

# Read in input file
input_dict = read_input(filename="Input_Files/hw8_input_file.yml")

print("*"*50)
print(f"Initial pressure guess using Wilson's correlation:")
print(input_dict['Pvap'])
print("K values using Wilson's correlation")
K_wilson = input_dict['Pvap']/input_dict['P']
print(K_wilson)

# Calculate vapor pressure at the given temperature assuming pure fluid
input_dict, eos_params = get_vapor_pressure(input_dict)
print(f"Vapor Pressure of components [1, 2, 3] at T = {input_dict['T']}K:")
print(input_dict['Pvap'])


print("K values using Raoult's Law")
K_raoults = input_dict['Pvap']/input_dict['P']
print(K_raoults)

# Close the output file
# close_output(f"{output_path}/output_file.txt")
