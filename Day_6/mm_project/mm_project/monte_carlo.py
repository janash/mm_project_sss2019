"""
monte_carlo.py
A package for doing MC of a Lennard Jones fluid

Handles the primary functions
"""

import numpy as np
import os

from mm_project.coordinates import generate_initial_coordinates, Box
from mm_project.potential import lennard_jones_potential

from mm_project import mm_cpp

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

        # Call the fast c++ implementation
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

        e_total = mm_cpp.get_particle_energy(i_particle, self.box.box_length, self.box.particles, self.cutoff, particle_movement)
        return e_total

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


def _adjust_displacement(n_trials, n_accept, max_displacement):
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


def run(method, box_length=None, num_particles=None, file_name=None, reduced_temperature=0.9, simulation_cutoff=3.0, max_displacement=0.1, n_steps=5000, freq=1000, tune_displacement=True, output_trajectory='traj.xyz'):
    """
    Run Monte Carlo Simulation
    """
    # ---------------
    # Parameter setup
    # ---------------

    beta = 1 / reduced_temperature

    # -------------------------
    # Simulation initialization
    # -------------------------

    # Coordinate initialization
    element = 'C'

    coordinates, box_length = generate_initial_coordinates(method=method, box_length=box_length, num_particles=num_particles, fname=file_name)

    # Initialized objects
    box = Box(box_length=box_length, particles=coordinates)
    mc_system = MCState(box=box, cutoff=simulation_cutoff, max_displacement=max_displacement, reduced_temperature=reduced_temperature)
    num_particles = mc_system.box.number_particles

    n_trials = 0
    n_accept = 0
    energy_array = np.zeros(n_steps)

    total_pair_energy = mc_system.calculate_total_energy()
    print(mc_system.total_pair_energy)

    traj = open(output_trajectory, 'w') 

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
                mc_system.max_displacement = _adjust_displacement(n_trials, n_accept, mc_system.max_displacement)

            # Print info
            print (i_step + 1, energy_array[i_step])

    traj.close()

if __name__ == "__main__":
    # Do something if this file is invoked on its own
    current_directory = os.path.dirname(os.path.abspath(__file__))
    nist_file = os.path.join(current_directory, 'data', 'nist_sample_config1.txt')

    run(method="file", file_name=nist_file)
