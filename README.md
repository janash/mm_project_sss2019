# 2019 Software Summer School - MM Sample Repo
A sample repo for the software summer school project.

Each day will have milestones for the students to complete. 

## Day 1

### Summary
On the first day, students will implement a Monte Carlo script individually (following the MM lead, Eliseo). After a working script is implemented, students will work individually to complete the student milestones outlined here.

### Instructor Milestones
It has been decided that the Instructor Script will have the following qualities:

1. The instructor script will implement a Monte Carlo simulation of a Lennard Jones fluid.
1. The instructor script will be a flat script with no user defined functions.
1. The script will use the minimum image convention for periodic boundaries.
1. A cut-off will be used for LJ calculation - that is, if particles are separated by a distance greater than the cut-off distance, the LJ interaction will not be calculated.
1. The instructor script will implement a tail correction at the beginning of the simulation
1. The instructor script will perform a Monte Carlo simulation using the following user defined variables:
    - reduced temperature (reduced_temperature)
    - cut-off (cutoff) The cut-off for LJ energy calculation
    - number of moves (num_steps) - The number of displacements to try
1. The instructor script will read the initial configuration for simulation from an XYZ file.

### Student Milestones
The student will refactor the Instructor Script into functions.

1. The student will write the following functions and refactor the instructor script to use these.
    - read_xyz(filename):
        - Description: reads xyz file. Assume that if there is a periodic box, the box size will be on line two of the the xyz file (usually a comment for xyz formats.) Our program assumes a cubic box, meaning that the length of the x dimensions, y dimensions, and z dimensions are equal. If they are not equal on this line, return box length as the maximum of the three dimensions.
        - Parameters: 
            - filename - the file path to the xyz file
        - Returns: 
            - atom_id: the atom names from the xyz file (line)
            - coordinates: The coordinates of the atoms (numpy array)
            - box_length: The box length.  
    - calculate_distance(r_i, r_j, periodic=False, box_length=0):
        - Description: Calculate distance between two particles. If `periodic=True`, distance should be calculated based on the minimum image convention and box size. The function should raise an error if `periodic=True` and `box_length=0`.
        - Parameters:
            - r_i - the position vector of the ith particle (x_i, y_i, z_i)
            - r_j - the position vector of the jth particle (x_j, y_j, z_j)
        - Returns:
            - r_ij - the distance between particles (scalar value)
    - calculate_LJ(sigma, epsilon, r_ij): 
        - Description: This function calculates the Lennard Jones interaction energy of two particles based on a given distance between particles.
        - Parameters:
            - sigma: the sigma parameter of the LJ equation (ie the particle size)
            - epsilon: the epsilon parameter of the LJ equation (ie the depth of the potential well)
        - Returns:
            - E_LJ - The pairwise Lennard Jones energy for particles i and j
1. Adjust the maximum displacement of the simulation in order to achieve ~50% acceptance tranlation rates

