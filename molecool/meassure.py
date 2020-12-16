"""
This module is for functions which perform meassurments.
"""

import numpy as np
from .atom_data import atomic_weights 

def calculate_angle(rA, rB, rC, degrees=False):
    # Calculate the angle between three points. Answer is given in radians by default, but can be given in degrees
    # by setting degrees=True
    AB = rB - rA
    BC = rB - rC
    theta=np.arccos(np.dot(AB, BC)/(np.linalg.norm(AB)*np.linalg.norm(BC)))

    if degrees:
        return np.degrees(theta)
    else:
        return theta

def calculate_distance(rA, rB):
    # This function calculates the distance between two points given as numpy arrays.
    d=(rA-rB)
    dist=np.linalg.norm(d)
    return dist

def calculate_molecular_mass(symbols):
    """Calculate the mass of a molecule.

    Parameters
    ----------
    symbols : list
        A list of elements.

    Returns
    -------
    mass : float
        The mass of the molecule
    """
    molecular_mass =  sum(atomic_weights[atom] for atom in symbols)

    return molecular_mass

def calculate_center_of_mass(symbols, coordinates):
   """Calculate the center of mass of a molecule.
   
   The center of mass is weighted by each atom's weight.
   
   Parameters
   ----------
   symbols : list
       A list of elements for the molecule
   coordinates : np.ndarray
       The coordinates of the molecule.
   
   Returns
   -------
   center_of_mass: np.ndarray
       The center of mass of the molecule.

   Notes
   -----
   The center of mass is calculated with the formula
   
   .. math:: \\vec{R}=\\frac{1}{M} \\sum_{i=1}^{n} m_{i}\\vec{r_{}i}
   
   """


   molecular_mass =  sum(atomic_weights[atom] for atom in symbols) 
   
   center_of_mass = 1/molecular_mass*sum(atomic_weights[symbols[i]] * coordinates[i] for i in range(len(symbols)))
   
   return center_of_mass