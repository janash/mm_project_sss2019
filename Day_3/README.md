# Day 3 - Object Oriented Programming and Design Patterns

The students should take about thirty minutes in their groups to discuss how to refactor their code using OOP principles. Where should OOP be used and why?

We do not have a particle class. Since the only attribute of a particle we care about is the coordinate, we do not have particle objects. Instead the MCSystem owns a list of coordinates which represent our particles.

## One Potential Solution

There are many ways you could sove this problem. The following is what the MolSSI team decided.

1. The Box class.
    The Box class is initialized with a box length and set of coordinates.
    -  box dimensions
        *  box_length
        *  coordinates - numpy array of particle coordinates

    The Box class has the following methods.
    - wrap
        Place points which are outside of box into box based on periodic boundaries.
    - minimum_image_distance
        - arguments are two points. Calculate the distance between these two points based on periodic boundaries.
    - volume
        * The box volume. Uses the `@property` decorator.


1. The MCState class.
Holds the state of the system (box with particles and thermodynamic quantities like reduced temperature)
    The MCState class has the following attributes.
    - box - a Box object
    - cutoff - the simulation cut off
    - max_displacement - the maximum displacement for a Monte Carlo move.
    - beta (1 / reduced_temperature)
    
    The MCState class has the following methods
    - total_potential_energy
    - tail_correction
    - calculate_particle_energy

1. The following functions are not associated with classes.
    - lennard_jones_potential
    - accept_or_reject
    - adjust_displacement
    - generate_initial_cooridinates

1. A Factory Design pattern is used with `generate_initial_coordinates` to specify method 

1. **Extension** - How could you use inheritance with the box class to create other types of boxes? 
    