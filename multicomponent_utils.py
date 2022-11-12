import numpy as np

def raoults_law(pvap, xi):
    """
    Calculate vapor pressure of a component assuming Raoults law
    :param pvap: Equilibrium vapor pressure of component i assuming pure component
    :param xi: Mole fraction of component i
    :return: Partial pressure of component i assuming Raoults law
    """

    return pvap * xi
