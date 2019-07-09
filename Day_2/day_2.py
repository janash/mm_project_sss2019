import os
import numpy as np

# Generate initial state

def generate_initial_coordinates(method='random', fname=None, num_particles=None, box_length=None):
    """
    Generate initial coordinates for MC system.

    Parameters
    ----------
    method : str
        Method for creating initial configuration. Valid options are 'random' or 'file'. 'random' will place particles in the simulation box randomly, depending on the num_particles argument, while 'file' will read coordinates from an xyz file.
    fname : str
        File path of file to read coordinates from (for `method='file'`).
    num_particles : int
        Number of particles to place in box (for `method='random'`)
    box_length : float
        Length of the simulation box (for `method='file'`)

    Returns
    -------
    coordinates : numpy array
        Array of coordinates (x, y, z)
    """
    if method == 'random':
        if num_particles is None:
            raise ValueError('generate_initial_state - "random" particle placement chosen, please input the number of particles using the num_particles argument.')
        if box_length is None:
            raise ValueError('generate_initial_state - "random" particle placement chosen, please input the box length using the box_length argument.')
        # Randomly placing particles in a box
        coordinates = (0.5 - np.random.rand(num_particles, 3)) * box_length
    
    elif method == 'file':
        try:
            # Reading a reference configuration from NIST
            coordinates = np.loadtxt(fname, skiprows=2, usecols=(1,2,3))
        except ValueError:
            if fname is None:
                raise ValueError("generate_initial_state: Method set to 'file', but no filepath given. Please specify an input file")
            else:
                raise ValueError
        except OSError:
            raise OSError(F'File {fname} not found.')
        except Exception as e:
            print(e)
            raise
    
    return coordinates

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

# Minimum image distance implementation

def minimum_image_distance(r_i, r_j, box_length):
    """
    Calculate the distance between two particles using the minimum image convention.

    The minimum image convention is for calculating distances in boxes with periodic boundary conditions.

    Parameters
    ----------
    r_i : numpy array
        The coordinates of atom i np.array(x,y,z) 
    r_j : numpy array
        The coordinates of atom j np.array(x,y,z)
    box_length: float
        Length of simulation box
    
    Returns
    -------
    rij2 : float
        The square distance between two particles.
    """

    rij = r_i - r_j
    rij = rij - box_length * np.round(rij / box_length)
    rij2 = np.dot(rij, rij)
    return rij2

# Computation of the total system energy

def total_pair_energy(coordinates, box_length, cutoff):
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

    e_total = 0.0
    particle_count = len(coordinates)
    cutoff2 = cutoff ** 2

    for i_particle in range(particle_count):
        for j_particle in range(i_particle):
            r_i = coordinates[i_particle]
            r_j = coordinates[j_particle]
            rij2 = minimum_image_distance(r_i, r_j, box_length)
            if rij2 < cutoff2:
                e_pair = lennard_jones_potential(rij2)
                e_total += e_pair

    return e_total

# Computation of the energy tail correction

def tail_correction(box_length, cutoff, number_particles):
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

    volume = np.power(box_length, 3)
    sig_by_cutoff3 = np.power(1.0 / cutoff, 3)
    sig_by_cutoff9 = np.power(sig_by_cutoff3, 3)
    e_correction = sig_by_cutoff9 - 3.0 * sig_by_cutoff3
    e_correction *= 8.0 / 9.0 * np.pi * number_particles / volume * number_particles

    return e_correction

def get_particle_energy(coordinates, i_particle, cutoff):
    """
    Calculate the interaction energy of a particle with its environment (all other particles in the system)

    Parameters
    ----------
    coordinates : numpy array
        The coordinates for all particles in the system
    i_particle : int
        The particle number for which to calculate the energy
    cutoff : float
        The simulation cutoff. Beyond this distance, interactions are not calculated.
    
    Returns
    -------
    e_total : float 
        The pairwise interaction energy of he i_th particle with all other particles in the system.
    """

    e_total = 0.0
    i_position = coordinates[i_particle]

    cutoff2 = cutoff ** 2
    particle_count = len(coordinates)

    for j_particle in range(particle_count):
        if i_particle != j_particle:

            j_position = coordinates[j_particle]
            rij2 = minimum_image_distance(i_position, j_position, box_length)

            if rij2 < cutoff2:
                e_pair = lennard_jones_potential(rij2)
                e_total += e_pair

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

def adjust_displacement(freq, i_step, n_trials, n_accept, max_displacement):
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

    if np.mod(i_step + 1, freq) == 0:
        acc_rate = float(n_accept) / float(n_trials)
        if (acc_rate < 0.380):
            max_displacement *= 0.8
        elif (acc_rate > 0.42):
            max_displacement *= 1.2
        n_trials = 0
        n_accept = 0
    return max_displacement

def update_output_file(traj, i_step, freq, element, num_particles, coordinates):
    """
    Write coordinates to output file.

    Parameters
    ----------
    traj : file
        Open file object. Where data will be written.
    i_step : int
        The step number.
    freq : int
        The write frequency. Write to a file on when mod(i_step + 1, freq) == 0
    element: str
        Element of particles
    num_particles : int
        The number of particles in the simulation
    coordinates : numpy array
        The coordinates for all particles in the system. 
    """

    if np.mod(i_step + 1, freq) == 0:

        traj.write(str(num_particles) + '\n\n')
        for i_particle in range(num_particles):
            traj.write("%s %10.5f %10.5f %10.5f \n" % (element, coordinates[i_particle][0], coordinates[i_particle][1], coordinates[i_particle][2]))

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
    file_name = os.path.join('..', 'nist_sample_config1.txt')
    coordinates = generate_initial_coordinates(method=method, fname=file_name)
    num_particles = len(coordinates)
    with open(file_name) as f:
        f.readline()
        box_length = float(f.readline().split()[0])

    
    simulation_cutoff2 = np.power(simulation_cutoff, 2)
    n_trials = 0
    n_accept = 0
    energy_array = np.zeros(n_steps)

    total_pair_energy = total_pair_energy(coordinates, box_length, simulation_cutoff)
    tail_correction = tail_correction(box_length, simulation_cutoff, num_particles)
    print(total_pair_energy)

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
    
        current_energy = get_particle_energy(coordinates, i_particle, simulation_cutoff)
        
        # Make a copy before adding random displacement
        proposed_coordinates = coordinates.copy()
        proposed_coordinates[i_particle] += random_displacement
        proposed_energy = get_particle_energy(proposed_coordinates, i_particle, simulation_cutoff)

        delta_e = proposed_energy - current_energy
        
        accept = accept_or_reject(delta_e, beta)

        if accept:
            total_pair_energy += delta_e
            n_accept += 1
            coordinates[i_particle] += random_displacement

        if tune_displacement:
            max_displacement = adjust_displacement(freq, i_step, n_trials, n_accept, max_displacement)

        total_energy = (total_pair_energy + tail_correction) / num_particles

        energy_array[i_step] = total_energy

        if np.mod(i_step + 1, freq) == 0:
            # Update output file
            traj.write(str(num_particles) + '\n\n')
            for i_particle in range(num_particles):
                traj.write("%s %10.5f %10.5f %10.5f \n" % (element, coordinates[i_particle][0], coordinates[i_particle][1], coordinates[i_particle][2]))
            
            # Print info
            print (i_step + 1, energy_array[i_step])

    traj.close()