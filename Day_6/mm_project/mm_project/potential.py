"""
Functions for potential energies.
"""

import numpy as np

def lennard_jones_potential(rij2):
    """
    Calculate the Lennard Jones pairwise potential between two particles based on a separation distance.
    
    Parameters
    ----------
    rij2 : float
        The squared distance between two particles.
    
    Returns
    -------
    float
        The Lennard Jones interaction energy for two particles.
    """

    sig_by_r6 = np.power(1 / rij2, 3)
    sig_by_r12 = np.power(sig_by_r6, 2)
    return 4.0 * (sig_by_r12 - sig_by_r6)
