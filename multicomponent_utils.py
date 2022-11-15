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

def get_phase_compositions(beta_v, Ki, Zv):
    """
    Calculate phase compositions of a component using root of Raoults law
    :param beta_v: Root of Raoults law
    :param Ki: Binary interactions assuming Raoults law
    :param Zi: List of compressibility (roots of cubic) for each pure component. [Zl, Zv]
    :return: Liquid and vapor phase compositions
    """

    xi = Zv / (1 - (1 - Ki) * beta_v)  # np.zeros_like(Zl)
    yi = Ki * xi  # np.zeros_like(Zl)

    # for i in range(len(Zl)):
    #     xi[i] = Zl[i] / (1  - (1 - Ki[i]) * beta_v)
    #     yi[i] = Ki[i] * xi[i]

    return xi, yi
