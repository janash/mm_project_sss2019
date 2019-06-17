import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#%matplotlib notebook

"""
This is the sample code for Day 2 of the software summer school MM project.

Day 2 covers
- PEP 8
- Numpy style docstrings
- Error and Exception handling
- testing using pytest

Student milestones:
Students will work with their teams on a common repository (MM_teamXX_2019) to fulfill the following milestones. Changes to the repo should be done using a Fork/PR model, where every change must be reviewed by one other person before merging.
1. Write numpy style docstrings for each function
    - generate_initial_state
    - lennard_jones_potential
    - minimum_image_distance
    - total_potential_energy
    - tail_correction
    - get_molecule_energy
    - move_particle
    - accept_or_reject
1. Use a linter to make sure your code adheres to PEP8 guidelines (yapf).
1. Add error handling to your functions.
    - For example, in the function `generate_initial_state`, your function should check that input parameters are compatible for each method.
        - if the method is "random", expected inputs are `num_particles` and `box_size`. Having additional or missing arguments should cause a TypeError.
        - if the method is "file", expected inputs are `fname`. Having additional or missing arguments should cause a TypeError.
        - You should check that the method is either "random" or "file". If not either of these, raise a ValueError.
        - use a `try except` clause to open the file for `method=file`
    *Instructor note* - you could walk them through writing this error checking - probably the most complicated of any function.
1. Make sure that each function has at least one unit test. Name the test file `test_mc.py`
    - create a fixture which returns a system - coordinates, atom names, box length.
    - Use @pytest.mark.parametrize to test the function `lennard_jones_potential` for a range of distances.
    - Identify one other function where you could use `parametrize` and write a test.
1. Test errors and exceptions 

"""

# Generate initial state

def generate_initial_state_2(method, **kwargs):

    if method != "random" and method != "file":
        raise ValueError(F'Unknown method {method} specified. Options are "random" or "file"')
    
    # Check inputs - this dictionary gives expected kwargs for each method.
    expected_kwargs = {
        "random": ['num_particles', 'box_length'],
        "file": ["fname"],
    }
    
    expected_keys = expected_kwargs[method]
    input_keys = list(kwargs.keys())
    
    ## Check that we only have expected arguments in kwargs
    list_diff = np.setdiff1d(input_keys, expected_keys)
    if list_diff:
        raise TypeError(F'Arguments {list_diff} not recognized for method="{method}"". Valid inputs are {expected_keys}')

    ## Check that we have all expected  arguments in kwargs
    list_diff = np.setdiff1d(expected_keys, input_keys)
    if list_diff:
        raise TypeError(F'Missing argument {list_diff} for input method {method}.')
    
    if method == "random":
        # Generate state based on inputs
        num_particles = kwargs['num_particles']
        box_length = kwargs['box_length']
        coordinates = (0.5 - np.random.rand(num_particles, 3)) * box_length
    else:
        # Else, method is file.
        fname = kwargs['fname']
        # Try reading the file
        try:
            # Reading a reference configuration from NIST
            coordinates = np.loadtxt(fname, skiprows=2, usecols=(1,2,3))
        except OSError:
            raise OSError(F'File {fname} not found.')
        except Exception as e:
            print(e)
            raise
    return coordinates
            

def generate_initial_state(method='random', fname=None, num_particles=None, box_length=None):

    if method is 'random':
        if num_particles is None:
            raise ValueError('generate_initial_state - "random" particle placement chosen, please input the number of particles using the num_particles argument.')
        if box_length is None:
            raise ValueError('generate_initial_state - "random" particle placement chosen, please input the box length using the box_length argument.')
        # Randomly placing particles in a box
        coordinates = (0.5 - np.random.rand(num_particles, 3)) * box_length
    
    elif method is 'file':
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

    sig_by_r6 = np.power(1 / rij2, 3)
    sig_by_r12 = np.power(sig_by_r6, 2)
    return 4.0 * (sig_by_r12 - sig_by_r6)

# Minimum image distance implementation

def minimum_image_distance(r_i, r_j, box_length):

    rij = r_i - r_j
    rij = rij - box_length * np.round(rij / box_length)
    rij2 = np.dot(rij, rij)
    return rij2

# Computation of the total system energy

def total_potential_energy(coordinates, box_length):

    e_total = 0.0

    for i_particle in range(num_particles):
        for j_particle in range(i_particle):
            r_i = coordinates[i_particle]
            r_j = coordinates[j_particle]
            rij2 = minimum_image_distance(r_i, r_j, box_length)
            if rij2 < cutoff2:
                e_pair = lennard_jones_potential(rij2)
                e_total += e_pair

    return e_total

# Computation of the energy tail correction

def tail_correction(box_length):

    volume = np.power(box_length, 3)
    sig_by_cutoff3 = np.power(1.0 / cutoff, 3)
    sig_by_cutoff9 = np.power(sig_by_cutoff3, 3)
    e_correction = sig_by_cutoff9 - 3.0 * sig_by_cutoff3
    e_correction *= 8.0 / 9.0 * np.pi * num_particles / volume * num_particles

    return e_correction

def get_molecule_energy(coordinates, i_particle):

    e_total = 0.0
    i_position = coordinates[i_particle]

    for j_particle in range(num_particles):
        if i_particle != j_particle:

            j_position = coordinates[j_particle]
            rij2 = minimum_image_distance(i_position, j_position, box_length)

            if rij2 < cutoff2:
                e_pair = lennard_jones_potential(rij2)
                e_total += e_pair

    return e_total


def move_particle(num_particles, max_displacement, coordinates):

    i_particle = np.random.randint(num_particles)
    random_displacement = (2.0 * np.random.rand(3) - 1.0)* max_displacement
    old_position = coordinates[i_particle].copy()
    old_energy = get_molecule_energy(coordinates, i_particle)
    coordinates[i_particle] += random_displacement
    new_energy = get_molecule_energy(coordinates, i_particle)

    return new_energy - old_energy, i_particle, old_position, coordinates

def accept_or_reject(delta_e, beta):

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

def restore_system(coordinates, i_particle, old_position):

    coordinates[i_particle] = old_position

    return coordinates

def adjust_displacement(freq, i_step, n_trials, n_accept, max_displacement):

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

    if np.mod(i_step + 1, freq) == 0:

        traj.write(str(num_particles) + '\n\n')
        for i_particle in range(num_particles):
            traj.write("%s %10.5f %10.5f %10.5f \n" % (element, coordinates[i_particle][0], coordinates[i_particle][1], coordinates[i_particle][2]))

if __name__ == "__main__":
        
    # ---------------
    # Parameter setup
    # ---------------

    reduced_density = 0.9
    reduced_temperature = 0.9
    max_displacement = 0.1
    n_steps = 50000
    freq = 1000
    tune_displacement = True

    beta = 1 / reduced_temperature

    # -------------------------
    # Simulation initialization
    # -------------------------

    # Coordinate initialization
    method = 'file'
    element = 'C'

    if method == 'random':
        num_particles = 100
        box_length = np.cbrt(num_particles / reduced_density)
        coordintes = generate_initial_state_2(method=method, num_particles=num_particles, box_length=box_length)
    else:
        file_name = os.path.join('..', 'nist_sample_config1.txt')
        coordinates = generate_initial_state_2(method=method, fname=file_name)
        num_particles = len(coordinates)
        with open(file_name) as f:
            f.readline()
            box_length = float(f.readline().split()[0])

    cutoff = 3.0
    cutoff2 = np.power(cutoff, 2)
    n_trials = 0
    n_accept = 0
    energy_array = np.zeros(n_steps)

    total_pair_energy = total_potential_energy(coordinates, box_length)
    tail_correction = tail_correction(box_length)
    print(total_pair_energy)

    traj = open('traj.xyz', 'w') 

    # --------------------------------
    # Metropolis Monte Carlo algorithm
    # --------------------------------

    for i_step in range(n_steps):

        n_trials += 1

        delta_e, i_particle, old_position, coordinates = move_particle(num_particles, max_displacement, coordinates)

        accept = accept_or_reject(delta_e, beta)

        if accept:
            total_pair_energy += delta_e
            n_accept += 1
        else:
            coordinates = restore_system(coordinates, i_particle, old_position)

        if tune_displacement:
            max_displacement = adjust_displacement(freq, i_step, n_trials, n_accept, max_displacement)

        total_energy = (total_pair_energy + tail_correction) / num_particles

        energy_array[i_step] = total_energy

        update_output_file(traj, i_step, freq, element, num_particles, coordinates)

        if np.mod(i_step + 1, 1000) == 0:
            print (i_step + 1, energy_array[i_step])

    traj.close()
