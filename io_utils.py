import yaml
import sys
import datetime
from eos import wilson_correlation
import numpy as np

# Read in yaml file
def read_input(filename: str='input_file.yml',) -> dict:
    """
    Read in yaml file
    :param filename: A string with path to input file
    :return: Dictionary with input parameters
    """
    with open(filename, 'r') as f:
        input_dict = yaml.load(f, Loader=yaml.FullLoader)

    # Number of components and check that critical values are present for all components
    input_dict["Nc"] = len(input_dict["Pc"])
    assert len(input_dict["Tc"]) == input_dict["Nc"]
    assert len(input_dict["w"]) == input_dict["Nc"]

    for i, component_pvap in enumerate(input_dict["Pvap"]):
        if component_pvap == "None":
            input_dict["Pvap"][i] = wilson_correlation(input_dict, i)

    # Construct K_ij matrix
    if input_dict["Nc"] > 1:
        input_dict["K_ij"] = construct_Kij(input_dict)

    # Convert list to arrays
    for key in ["Pc", "Tc", "w", "Pvap"]:
        if not isinstance(input_dict[key], np.ndarray):
            input_dict[key] = np.array(input_dict[key])

    return input_dict

def construct_Kij(input_dict):
    # Construct binary interaction matrix
    K_ij = np.zeros((input_dict['Nc'], input_dict['Nc']))
    Kij_nonzero_components = None
    if 'Kij' in input_dict:
        Kij_nonzero_components = [*input_dict['Kij'].keys()]

    if Kij_nonzero_components is not None:
        for key in Kij_nonzero_components:
            print(key)
            row = int(key.rsplit(" ")[-1])
            K_ij[row, row+1:] = input_dict['Kij'][key]
    K_ij = K_ij + K_ij.T

    return K_ij


# Write to output file
def pout(msg: str, filename: str='HW5_Output/outfile.txt') -> None:
    """
    Write to output file
    :param msg: Message to be logged
    :param filename: A string with path to output file
    :return: None
    """
    with open(filename, 'a') as f:
        f.write(datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S') + '\t\t' + msg + '\n')
    
    return

# Redirect stdout to the output file
def redirect_stdout(filename: str) -> None:
    """
    Redirect stdout to the output file
    :param filename: A string with path to output file
    :return: None
    """
    sys.stdout = open(filename, 'a')

def close_output(filename: str) -> None:
    """
    Close the output file
    :param filename: A string with path to output file
    :return: None
    """
    sys.stdout.close()
    sys.stdout = sys.__stdout__