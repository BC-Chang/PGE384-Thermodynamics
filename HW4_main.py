import numpy as np
from eos import cubic_eos
from solve import solve_cardanos
import matplotlib.pyplot as plt
from io_utils import read_input, pout
import os
import pandas as pd


def main():
    # Create a HW4 output directory if one does not already exist
    if not os.path.exists("HW4_Output"):
        os.makedirs("HW4_Output")
    
    
    # Read in input file
    input_dict = read_input()
    
    # Find the coefficients of the specified cubic EoS at given pressure and temperature
    alpha, beta, gamma, a, b = cubic_eos(P=input_dict["P"], T=input_dict["T"], eos=input_dict['eos'])
    
    # Print the cubic coefficients to output file
    pout("-"*53)
    pout(f"Homework 4 Output")
    pout(f"P = {input_dict['P']} Pa")
    pout(f"T = {input_dict['T']} K")
    pout(f"EoS = {input_dict['eos']}")
    pout("-"*53)
    pout(f"Coefficients: ")
    pout(f"\t{alpha = :.5f}")
    pout(f"\t{beta = :.5f}")
    pout(f"\t{gamma = :.5f}")
    pout("-"*53)
    
    # Calculate roots using Cardano's method
    x1, x2, x3 = solve_cardanos(1, alpha, beta, gamma)
    
    # Print roots to the output file
    pout(f"Roots: ")
    pout(f"\t{x1 = :.6f}")
    pout(f"\t{x2 = :.6f}")
    pout(f"\t{x3 = :.6f}")
    pout("-"*53)
    pout(f"Compressibility Factor = {x1:.6f}")
    
    # Calculate molar volume
    molar_V = x1*8.314*input_dict['T'] / input_dict['P']
    
    # Print molar volume to output file
    pout(f"Molar Volume = {molar_V:.6f} m3/mol")
    pout("-"*53)
    
    # Create array of cubic EoS vs. Z
    Z = np.linspace(0, 1.6, 1000)
    P = Z**3 + alpha*Z**2 + beta*Z + gamma
    
    # Write compressibility factor and P to an Excel File
    df = pd.DataFrame(np.array([Z, P]).T, columns=['Z', 'Cubic EoS'])
    df.to_excel("HW4_Output/Graphical_Solution.xlsx")
    
    plt.figure(dpi=400)
    plt.plot(Z, P, 'b-')
    plt.plot([x1, x1], [-.25, 0], 'r--')  # Plot red line for Z root
    plt.plot([0, 1.6], [0, 0], 'k-')  # Plot y = 0 axis
    plt.grid(True)
    plt.xlabel('Compressibility Factor (Z)')
    plt.ylabel('Cubic EoS')
    plt.xlim([0, 1.6])
    plt.ylim([-0.2, 1.5])
    plt.savefig(f"HW4_Output/cubic_eos_plot_{input_dict['P']*10**-6}MPa_{input_dict['T']}K.png")

    return




if __name__ == '__main__':
    main()
