# Day 3 - Object Oriented Programming and Design Patterns

## Student Milestones

1. Create a class called Particle which will hold the following attributes for each particle in the system.
    - name 
    - symbol
    - coordinate (x,y,z as numpy array)
    - total_pair_energy (interaction energy of particle with system)
    - **Question** - Do we need a particle class, or should we just keep track of particle coordinates and info with lists and numpy arrays in Box class?

1. Create a Box class.
    The Box class should have the the following attributes
        - box_length
        - particles - list of Particle objects in the box
        - number of particles
        - coordinates

    The Box class should have the following methods
        - minimum_image_distance
        - restore_system
            - set_particle_coordinates


1. Create a System class.
    The MCSystem class should have the following attributes
        - box - a Box object containing particles and box information
        - cutoff - the simulation cut off
        - max_displacement - the maximum displacement for a Monte Carlo move.
        - temperature
    The MCSystem class should have the following methods
        - lennard_jones_potential (a static method)
        - accept_or_reject (a static method)
        - total_potential_energy
        - tail_correction
        - get_molecule_energy
            - *Note* This should set the `total_pair_energy` attribute for a particle.
        - move_particle
        - accept_or_reject

1. Use the Factory Design patterns to create MCSystems
    