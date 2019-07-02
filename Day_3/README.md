# Day 3 - Object Oriented Programming and Design Patterns

OOP classes hold state.

## Student Milestones

1. Create a class called Particle which will hold the following attributes for each particle in the system.
    - name 
    - symbol
    - coordinate (x,y,z as numpy array)
    - total_pair_energy (interaction energy of particle with system)
    - **Question** - Do we need a particle class, or should we just keep track of particle coordinates and info with lists and numpy arrays in Box class?

1. Create a Box class.
    The Box class should have the the following attributes.
        - box dimensions
            - box_length (box dimensions)
            - could be extended to x,y,z and angles
                - Classes which inherit from box and have extensions.

    The Box class should have the following methods.
        - minimum_image_distance
        - volume
            - methods will return values rather than storing. Based on current state of class.     


1. Create a MCSystem class.
    The MCSystem class should have the following attributes
        - box - a Box object
        - particles - list of Particle objects
        - cutoff - the simulation cut off
        - max_displacement - the maximum displacement for a Monte Carlo move.
        - beta
    
    The MCSystem class should have the following methods
        - lennard_jones_potential (a static method (memory efficiency) - or function outside of class)
        - accept_or_reject (a static method (memory efficiency) - or function outside of class)
        - total_potential_energy
        - tail_correction
        - get_molecule_energy
            - *Note* This should set the `total_pair_energy` attribute for a particle.
        - move_particle

1. Use the Factory Design patterns to create MCSystems
    