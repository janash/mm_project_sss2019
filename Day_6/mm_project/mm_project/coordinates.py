"""
Functions for building ininital coordinates
"""

import numpy as np
import math

def generate_initial_coordinates(method, **kwargs):
    """
    Generate initial coordinates for MC system.

    Parameters
    ----------
    method : str
        Method for creating initial configuration. Valid options are 'random' or 'file'. 'random' will place particles in the simulation box randomly, depending on the num_particles argument, while 'file' will read coordinates from an xyz file.

    Returns
    -------
    coordinates : numpy array
        Array of coordinates (x, y, z)
    """
    
    if method == "random":
        coord_method = _random
    elif method == "file":
        coord_method = _from_file
    else:
        raise ValueError("Method not found.")
    
    return coord_method(**kwargs)

def _random(**kwargs):

    # Make sure kwargs are what we want
    try:
        num_particles = kwargs['num_particles']
        box_length = kwargs['box_length']
    except:
        raise KeyError('Incorrect arguments for generage initial coordinates with method=random!')

    if num_particles is None:
        raise ValueError('generate_initial_coordinates - "random" particle placement chosen, please input the number of particles using the num_particles argument.')
    if box_length is None:
        raise ValueError('generate_initial_coordinates - "random" particle placement chosen, please input the box length using the box_length argument.')
    # Randomly placing particles in a box
    coordinates = (0.5 - np.random.rand(num_particles, 3)) * box_length
    
    return coordinates, box_length

def _from_file(**kwargs):
    
     # Make sure kwargs are what we want
    try:
        fname = kwargs['fname']
    except:
        raise KeyError('Incorrect arguments for generage initial coordinates with method=file!')

    
    try:
        # Reading a reference configuration from NIST
        coordinates = np.loadtxt(fname, skiprows=2, usecols=(1,2,3))
        with open(fname) as f:
            f.readline()
            box_length = float(f.readline().split()[0])
    except ValueError:
        if fname is None:
            raise ValueError("generate_initial_coordinates: Method set to 'file', but no filepath given. Please specify an input file")
        else:
            raise ValueError
    except OSError:
        raise OSError(F'File {fname} not found.')
    except Exception as e:
        print(e)
        raise
    
    return coordinates, box_length

class Box:
    def __init__(self, box_length, particles=None):
        self.box_length = box_length
        self.particles = particles # if None, no particles in MCSystem
    
    @property
    def number_particles(self):
        """Calculate the number of particles"""
        if self.particles is None:
            return 0
        else:
            return len(self.particles)
    
    def minimum_image_distance(self, r_i, r_j):
        """
        Calculate the distance between two points using the minimum image convention.

        The minimum image convention is for calculating distances in boxes with periodic boundary conditions.

        Parameters
        ----------
        r_i : numpy array
            The coordinates of first point np.array(x,y,z) 
        r_j : numpy array
            The coordinates of second point np.array(x,y,z)
        
        Returns
        -------
        rij2 : float
            The square distance between two particles.
        """

        rij = r_i - r_j
        rij = rij - self.box_length * np.round(rij / self.box_length)
        rij2 = np.dot(rij, rij)
        return rij2
    
    @property
    def volume(self):
        """Calculate the volume of the box based on the box length"""
        return self.box_length ** 3
    
    def wrap(self):
        """Wrap points which are outside of box into box"""
        # TODO
        pass