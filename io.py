import yaml
import sys


# Read in yaml file
def read_yaml(filename: str='input_file.yml',) -> dict:
    """
    Read in yaml file
    :param filename: A string with path to input file
    :return: Dictionary with input parameters
    """
    with open(filename, 'r') as f:
        return yaml.load(f, Loader=yaml.FullLoader)


# Write to output file
def pout(msg: str, filename: str='outfile.txt') -> None:
    """
    Write to output file
    :param msg: Message to be logged
    :param filename: A string with path to output file
    :return: None
    """
    with open(filename, 'a') as f:
        f.write(msg + '\n')
        return
