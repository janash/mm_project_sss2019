import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#%matplotlib notebook

# Generate initial state

def generate_initial_coordinates(method='random', file_name=None, num_particles=None, box_length=None):

    if method is 'random':
        coordinates = (0.5 - np.random.rand(num_particles, 3)) * box_length
    
    elif method is 'file':
        coordinates = np.loadtxt(file_name, skiprows=2, usecols=(1,2,3))
    
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

def calculate_total_pair_energy(coordinates, box_length, cutoff2):

    e_total = 0.0
    particle_count = len(coordinates)

    for i_particle in range(particle_count):
        for j_particle in range(i_particle):
            r_i = coordinates[i_particle]
            r_j = coordinates[j_particle]
            rij2 = minimum_image_distance(r_i, r_j, box_length)
            if rij2 < cutoff2:
                e_pair = lennard_jones_potential(rij2)
                e_total += e_pair

    return e_total

def calculate_tail_correction(box_length, cutoff, number_particles):


    volume = np.power(box_length, 3)
    sig_by_cutoff3 = np.power(1.0 / cutoff, 3)
    sig_by_cutoff9 = np.power(sig_by_cutoff3, 3)
    e_correction = sig_by_cutoff9 - 3.0 * sig_by_cutoff3
    e_correction *= 8.0 / 9.0 * np.pi * number_particles / volume * number_particles

    return e_correction

def get_particle_energy(coordinates, i_particle, cutoff2):

    e_total = 0.0
    i_position = coordinates[i_particle]

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
    n_steps = 500000
    freq = 1000
    tune_displacement = True
    simulation_cutoff = 3.0
    element = 'C'
    beta = 1 / reduced_temperature

    # -------------------------
    # Coordinate initialization
    # -------------------------

    # Method = random
    build_method = 'random'
    reduced_density = 0.9
    num_particles = 100
    box_length = np.cbrt(num_particles / reduced_density)
    coordinates = generate_initial_coordinates(method=build_method, num_particles=num_particles, box_length=box_length)

    # -------------------------
    # Simulation initialization
    # -------------------------

    # Method = file
#    build_method = 'file'
#    file_name = os.path.join('..', 'nist_sample_config1.txt')
#    coordinates = generate_initial_coordinates(method=build_method, file_name=file_name)
#    num_particles = len(coordinates)
#    with open(file_name) as f:
#        f.readline()
#        box_length = float(f.readline().split()[0])

    simulation_cutoff2 = np.power(simulation_cutoff, 2)
    n_trials = 0
    n_accept = 0
    energy_array = np.zeros(n_steps)

    total_pair_energy = calculate_total_pair_energy(coordinates, box_length, simulation_cutoff2)
    tail_correction = calculate_tail_correction(box_length, simulation_cutoff, num_particles)

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
    
        current_energy = get_particle_energy(coordinates, i_particle, simulation_cutoff2)
        
        # Make a copy before adding random displacement
        proposed_coordinates = coordinates.copy()
        proposed_coordinates[i_particle] += random_displacement
        proposed_energy = get_particle_energy(proposed_coordinates, i_particle, simulation_cutoff2)

        delta_e = proposed_energy - current_energy
        
        accept = accept_or_reject(delta_e, beta)

        if accept:
            total_pair_energy += delta_e
            n_accept += 1
            coordinates[i_particle] += random_displacement

        total_energy = (total_pair_energy + tail_correction) / num_particles

        energy_array[i_step] = total_energy

        if np.mod(i_step + 1, freq) == 0:
            # Update output file
            traj.write(str(num_particles) + '\n\n')
            for i_particle in range(num_particles):
                traj.write("%s %10.5f %10.5f %10.5f \n" % (element, coordinates[i_particle][0], coordinates[i_particle][1], coordinates[i_particle][2]))
            
            # Adjust displacement
            if tune_displacement:
                max_displacement = adjust_displacement(n_trials, n_accept, max_displacement)

            # Print info
            print (i_step + 1, energy_array[i_step])

    traj.close()

    #plt.plot(energy_array)
    #plt.xlabel('Monte Carlo steps')
    #plt.ylabel('Energy (reduced units)')
    #plt.grid(True)
    #plt.show()

    #plt.figure()
    #ax = plt.axes(projection='3d')
    #ax.plot3D(coordinates[:,0], coordinates[:,1], coordinates[:,2], 'o')
    #plt.show()
