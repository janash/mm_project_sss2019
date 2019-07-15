import os
import numpy as np
import mm_cpp

import math

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

class MCState:
    def __init__(self, box, cutoff, max_displacement, reduced_temperature):
        self.box = box
        self.cutoff = cutoff
        self.beta = 1 / reduced_temperature
        self.max_displacement = max_displacement
        
        self.total_pair_energy = self.calculate_total_pair_energy()
        self.tail_correction = self.calculate_tail_correction()

        self.total_energy = self.total_pair_energy + self.tail_correction
            
    def calculate_total_pair_energy(self):
        """
        Calculate the total potential energy of the system using a pairwise Lennard Jones potential.

        Parameters
        ----------
        coordinates: numpy array
            The coordinates of all particles in the system. 
        box_length: float
            Length of simulation box
        
        Returns
        -------
        e_total : float
            The total pairwise potential energy of the system.
        """

        e_total = mm_cpp.calculate_total_pair_energy(self.box.box_length, self.box.particles, self.cutoff)
        return e_total

    
    def calculate_tail_correction(self):
        """
        Calculate the long range tail correction for a system.

        Parameters
        ----------
        box_length : float   
            The length of the simulation box
        cutoff : float
            The simulation cutoff. Beyond this distance, interactions are not calculated.
        number_particles:
            The number of particles in the simulation
        
        Returns
        -------
        e_correction : float
            The long distance tail correction.
        """

        volume = self.box.volume
        sig_by_cutoff3 = np.power(1.0 / self.cutoff, 3)
        sig_by_cutoff9 = np.power(sig_by_cutoff3, 3)
        e_correction = sig_by_cutoff9 - 3.0 * sig_by_cutoff3
        e_correction *= 8.0 / 9.0 * np.pi * self.box.number_particles / volume * self.box.number_particles

        return e_correction
    
    def calculate_total_energy(self):
        tail_correction = self.calculate_tail_correction()
        pair_energy = self.calculate_total_pair_energy()
        return tail_correction + pair_energy
    
    def get_particle_energy(self, i_particle, particle_movement=None):
        """
        Calculate the interaction energy of a particle with its environment (all other particles in the system)

        Parameters
        ----------
        i_particle : int
            The particle index for which to calculate the energy
       
        Returns
        -------
        e_total : float 
            The pairwise interaction energy of he i_th particle with all other particles in the system.
        """

        # The C++ function always requires a displacement. So if it is None, just
        # displace by (0,0,0)

        if particle_movement is None:
            particle_movement = np.array([0,0,0])

        return mm_cpp.get_particle_energy(i_particle, self.box.box_length, self.box.particles, self.cutoff, particle_movement)

    
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

def _random(num_particles, box_length):
    if num_particles is None:
        raise ValueError('generate_initial_coordinates - "random" particle placement chosen, please input the number of particles using the num_particles argument.')
    if box_length is None:
        raise ValueError('generate_initial_coordinates - "random" particle placement chosen, please input the box length using the box_length argument.')
    # Randomly placing particles in a box
    coordinates = (0.5 - np.random.rand(num_particles, 3)) * box_length
    
    return coordinates, box_length

def _from_file(fname):
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

# Lennard Jones potential implementation

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

# Computation of the energy tail correction

def accept_or_reject(delta_e, beta):
    """
    Calculate if move is accepted or rejected based on the energy change, simulation temperature (beta), and (possibly) a random number.

    In Monte Carlo simuations, moves which result in a negative change in energy (delta_e) are accepted, while moves which result in a positive energy change are accepted or rejected depending on calculation of acceptance probability (p_acc = np.exp(-beta*delta_e)) compared to a random number generator.

    Parameters
    ----------
    delta_e : float
        Energy change caused by movement of a particle
    beta : float
        1 / (reduced temperature)
    
    Returns
    -------
    accept : bool
        True if move is accepted, False if move is rejected.
    """
    if delta_e <= 0.0:
        accept = True
    else:
        random_number = np.random.rand(1)
        p_acc = np.exp(-beta*delta_e)
        if random_number < p_acc:
            accept = True
        else:
            accept = False

    return accept


def adjust_displacement(n_trials, n_accept, max_displacement):
    """
    Adjust the max displacement to achieve ideal acceptance rate.

    Parameters
    ----------
    freq : int
        Step interval at which adjustment should be performed. 
    i_step : int
        The current step number
    n_trials : int
        Count of the number of trials since last adjustment
    n_accept : int
        Count of the number of trials which have resulted in a move acceptance.
    max_displacement : float
        The maximum distance a particle can be moved during a trial
    
    Return
    ------
    max_displacement : float

    """

    acc_rate = float(n_accept) / float(n_trials)
    if (acc_rate < 0.380):
        max_displacement *= 0.8
    elif (acc_rate > 0.42):
        max_displacement *= 1.2
    n_trials = 0
    n_accept = 0
    return max_displacement



if __name__ == "__main__":
        
    # ---------------
    # Parameter setup
    # ---------------

    reduced_temperature = 0.9
    max_displacement = 0.1
    n_steps = 50000
    freq = 1000
    tune_displacement = True
    simulation_cutoff = 3.0
    
    simulation_cutoff2 = np.power(simulation_cutoff, 2)
    beta = 1 / reduced_temperature

    # -------------------------
    # Simulation initialization
    # -------------------------

    # Coordinate initialization
    method = 'file'
    element = 'C'

    # Method = random
    # reduced_density = 0.9
    # num_particles = 100
    # box_length = np.cbrt(num_particles / reduced_density)
    # coordinates = generate_initial_coordinates(method=method, num_particles=num_particles, box_length=box_length)

    # Method = file

    file_name = os.path.join('..','nist_sample_config1.txt')
    coordinates, box_length = generate_initial_coordinates(method=method, fname=file_name)

    # Initialized objects
    box = Box(box_length=box_length, particles=coordinates)
    mc_system = MCState(box=box, cutoff=simulation_cutoff, max_displacement=max_displacement, reduced_temperature=reduced_temperature)
    num_particles = mc_system.box.number_particles

    n_trials = 0
    n_accept = 0
    energy_array = np.zeros(n_steps)

    total_pair_energy = mc_system.calculate_total_energy()
    print(mc_system.total_pair_energy)

    traj = open('traj.xyz', 'w') 

    # --------------------------------
    # Metropolis Monte Carlo algorithm
    # --------------------------------

    for i_step in range(n_steps):

        n_trials += 1

        # ---------------------------
        #  Propose a Monte Carlo Move
        # ---------------------------
        i_particle = np.random.randint(num_particles)
        random_displacement = (2.0 * np.random.rand(3) - 1.0)* max_displacement
    
        current_energy = mc_system.get_particle_energy(i_particle)
        
        # Make a copy before adding random displacement
        proposed_energy = mc_system.get_particle_energy(i_particle, random_displacement)

        delta_e = proposed_energy - current_energy
        
        accept = accept_or_reject(delta_e, beta)

        if accept:
            mc_system.total_pair_energy += delta_e
            n_accept += 1
            mc_system.box.particles[i_particle] += random_displacement

        total_energy = (mc_system.total_pair_energy + mc_system.tail_correction) / mc_system.box.number_particles

        energy_array[i_step] = total_energy

        if np.mod(i_step + 1, freq) == 0:
            # Update output file
            traj.write(str(num_particles) + '\n\n')
            for i_particle in range(mc_system.box.number_particles):
                particle = mc_system.box.particles[i_particle]
                traj.write("%s %10.5f %10.5f %10.5f \n" % (element, particle[0], particle[1], particle[2]))
            
            # Adjust displacement
            if tune_displacement:
                mc_system.max_displacement = adjust_displacement(n_trials, n_accept, mc_system.max_displacement)

            # Print info
            print (i_step + 1, energy_array[i_step])

    traj.close()
