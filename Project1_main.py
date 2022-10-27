import numpy as np
from eos import cubic_eos
from solve import solve_cardanos, check_roots
import matplotlib.pyplot as plt
from io_utils import read_input, pout
from pr_utils import fugacity_coefficient
import os
import pandas as pd
from time import perf_counter_ns

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
print(f"Initial pressure guess using Wilson's correlation = {input_dict['P']: 0.3f} Pa")

# Set initial error to be very large
err = 1.E9
i = 0

tic = perf_counter_ns()
# Calculating the vapor pressure
while err > input_dict["eps"]:
    # Find the coefficients of the specified cubic EoS at given pressure and temperature
    alpha, beta, gamma, A, B = cubic_eos(P=input_dict["P"], T=input_dict["T"], eos=input_dict['eos'],
                                   Pc=input_dict["Pc"], Tc=input_dict["Tc"], w=input_dict["w"])

    # Calculate roots using Cardano's method
    x1, x2, x3 = solve_cardanos(1, alpha, beta, gamma)


    # Get the roots corresponding to liquid and vapor phase, respectively.
    zl, zv = check_roots(x1, x2, x3)
    # print(zl, zv)

    # Calculate the fugacity coefficients
    # Liquid fugacity coefficient
    fc_l = fugacity_coefficient(zl, A, B)

    # Vapor phase fugacity coefficient
    fc_v = fugacity_coefficient(zv, A, B)

    # Calculate error using fugacity coefficients
    err = abs(fc_l - fc_v)

    # Update pressure if error is greater than convergence criterion
    if err > input_dict["eps"]:
        input_dict["P"] = input_dict["P"] * np.exp(fc_l) / np.exp(fc_v)

    if i == input_dict["maxiter"]:
        print("Maximum number of iterations reached")
        break

    i += 1

toc = perf_counter_ns()
print(f"Final vapor pressure = {input_dict['P'] :.3f} Pa")
print(f"Elapsed Time = {(toc - tic)*1E-9 :.5f} s")
