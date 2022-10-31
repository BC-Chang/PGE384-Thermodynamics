import yaml
import sys
import datetime
from eos import wilson_correlation

# Read in yaml file
def read_input(filename: str='input_file.yml',) -> dict:
    """
    Read in yaml file
    :param filename: A string with path to input file
    :return: Dictionary with input parameters
    """
    with open(filename, 'r') as f:
        input_dict = yaml.load(f, Loader=yaml.FullLoader)

    if input_dict["P"] == "None":
        input_dict["P"] = wilson_correlation(input_dict)

    return input_dict




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