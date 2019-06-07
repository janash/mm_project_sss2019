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
1. The instructor script will use reduced units.
1. The instructor script will read an initial configuration from an xyz file from NIST.
1. The script will use the minimum image convention for periodic boundaries.
1. The instructor script will implement a tail correction at the beginning of the simulation
1. The instructor script will perform a Monte Carlo simulation using the following user defined variables:
    - reduced temperature (reduced_temperature)
    - cut-off (cutoff) The cut-off for LJ energy calculation
    - number of moves (num_steps) - The number of displacements to try
    - maximum displacement (max_displacement) - The maximum displacement of a particle in one dimension during a move.
1. A cut-off will be used for LJ calculation - that is, if particles are separated by a distance greater than the cut-off distance, the LJ interaction will not be calculated.
1. Accuracy Check: The instructor script will calculate the energy of the initial configuration of the system, and will match the energy given by NIST for the configuration. https://www.nist.gov/mml/csd/chemical-informatics-research-group/lennard-jones-fluid-reference-calculations


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
    - calculate_LJ(r_ij): 
        - Description: This function calculates the Lennard Jones interaction energy of two particles based on a given distance between particles.
        - Parameters:
            - r_ij: The scalar distance between particles i and j.
        - Returns:
            - E_LJ - The pairwise Lennard Jones energy for particles i and j.
    - calculate_total_energy(positions):
        - Description: This function calculates the total energy for a system configuration.
        - Parameters: 
            - positions - numpy array of system coordinates. Each particle has an x, y, and z position.
        - Returns:
            - system_energy - The energy of the system configuration.
1. Accuracy check - Verify that your MC simulation is accurate by comparint the energy to NIST benchmarks (https://mmlapps.nist.gov/srs/LJ_PURE/mc.htm). In other words, does your `U*` match watch NIST reports for your chosen reduced temperature (`T*`)
  - You will need a way to record the average system energy (after equilibration) for this. This can be computed as the simulation runs, or you may choose to write the system energy to a file after equilibration.
1. **Extension** Adjust the maximum displacement of the simulation in order to achieve ~50% acceptance tranlation rates.
1. **Extension**: Create a visualization of the system coordinates.
1. **Extension** Write a function to create an initial configuration on a lattice with a given box size and number of particles. 