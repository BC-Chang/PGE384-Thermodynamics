import numpy as np

def raoults_law(pvap, xi):
    """
    Calculate vapor pressure of a component assuming Raoults law
    :param pvap: Equilibrium vapor pressure of component i assuming pure component
    :param xi: Mole fraction of component i
    :return: Partial pressure of component i assuming Raoults law
    """

    return pvap * xi

def get_bubble_point(pvap, z):
    """
    Calculate bubble point of a component assuming Raoults law
    :param pvap: Equilibrium vapor pressure of component i assuming pure component
    :param z: Mole fraction of component i
    :return: Bubble point assuming Raoults law
    """

    return np.sum(pvap * z)


def get_dew_point(pvap, z):
    """
    Calculate dew point of a component assuming Raoults law
    :param pvap: Equilibrium vapor pressure of component i assuming pure component
    :param z: Mole fraction of component i
    :return: Dew point assuming Raoults law
    """

    return 1 / np.sum(z / pvap)
