# Day 3 - Object Oriented Programming and Design Patterns

OOP classes hold state.

## Student Milestones

1. Create a class called Particle which will hold the following attributes for each particle in the system.
    - symbol

1. Create a Box class.
    The Box class should be initialized with a box length, and should have the following.
    -  box dimensions
        *  box_length (box dimensions)
        * could be extended to x,y,z and angles
            * Classes which inherit from box and have extensions.

    The Box class should have the following methods.
    - minimum_image_distance
    - volume
        * methods will return values rather than storing. Based on current state of class.     


1. Create a MCSystem class.
    The MCSystem class should have the following attributes
    - box - a Box object
    - particles - list of Particle objects (the index of this list corresponds to index of the coordinate array.)
     - coordinates
        - a numpy array
    - cutoff - the simulation cut off
    - max_displacement - the maximum displacement for a Monte Carlo move.
    - beta
    
    The MCSystem class should have the following methods
        - lennard_jones_potential (a static method (memory efficiency) - or function outside of class)
        - accept_or_reject (a static method (memory efficiency) - or function outside of class)
        - total_potential_energy
        - tail_correction
        - get_particle_energy
            - *Note* This should set the `total_pair_energy` attribute for a particle.
        - move_particle -> rename `propose_move`
            - select random particle
            - apply random displacement based on self.max_displacement
            - calculate energy difference
            - return energy difference

1. Use the Factory Design patterns to create MCSystems
    - Specifiy method
        - associated arguments
    