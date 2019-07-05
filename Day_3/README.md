# Day 3 - Object Oriented Programming and Design Patterns

OOP classes hold state.

We do not have a particle class. Since the only attribute of a particle we care about is the coordinate, we do not have particle objects. Instead the MCSystem owns a list of coordinates which represent our particles.

## Student Milestones

1. Create a Box class.
    The Box class should be initialized with a box length, and should have the following attributes.
    -  box dimensions
        *  box_length (box dimensions)
        * could be extended to x,y,z and angles
            * Classes which inherit from box and have extensions.

    The Box class should have the following methods.
    - wrap
        Place points which are outside of box into box based on periodic boundaries.
    - minimum_image_distance
        - arguments are two points. Calculate the distance between these two points based on periodic boundaries.
    - volume
        * method will return values rather than storing. Based on current state of class.     


1. Create a MCState class.
    The MCState class should have the following attributes
    - box - a Box object
    - particles - numpy array
    - cutoff - the simulation cut off
    - max_displacement - the maximum displacement for a Monte Carlo move.
    - beta (1 / reduced_temperature)
    
    The MCState class should have the following methods
    - total_potential_energy
    - tail_correction
    - calculate_particle_energy
    

1. Functions
    - lennard_jones_potential
    - generate_initial_coordinates
    - accept_or_reject

1. Use the Factory Design patterns to create MCSystems
    - Specifiy method
        - associated arguments
        - return MCSystem class

1. **Extension** - How could you use inheritance with the box class to create other types of boxes? 
    