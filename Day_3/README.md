# Day 3 - Object Oriented Programming and Design Patterns

The students should take about thirty minutes in their groups to discuss how to refactor their code using OOP principles. Where should OOP be used and why?

Considerations - Later in the week, you will be refactoring the energy calculation in this code to use C++. For this, you will still need all of the particle coordinates to be in a numpy array.

## Refactoring your code - the approach
For this day, you will have to refactor both the code and the associated tests. It may be more helpful to work on both at the same time. You should make a plan and assign tasks to team members before starting.

## One Potential Solution

There are many ways you could sove this problem. The following is what the MolSSI team decided.

We do not have a particle class. Since the only attribute of a particle we care about is the coordinate, we do not have particle objects. Instead the Box class (explained below) owns a list of coordinates which represent our particles.

1. The Box class.
    The Box class is initialized with a box length. A set of particle coordinates as a numpy array is an optional input. The default argument for `particles` is `None`. Meaning that you can have a box without particles. 

    The Box class has the following attributes.
    -  box dimensions (ie the `box_length`)
    - particles - `None` or a numpy array of particle coordinates

    The Box class has the following methods.
    - wrap  
        - Place particles which are outside of box into box based on periodic boundaries.
    - minimum_image_distance
        - arguments are two points. Calculate the distance between these two points based on periodic boundaries.
    - volume
        * The box volume. Uses the `@property` decorator.
    - number_particles
        * Uses the `@property` decorator. Returns 0 if there are no particles (ie `self.particles` is `None`) or the length of the numpy array if there are particles.


1. The MCState class.
    Holds the state of the system (box with particles and thermodynamic quantities like reduced temperature).

    The MCState class is initialized with a box object, cutoff distance, maximum displacement, and reduced temperature.

    The MCState class has the following attributes.
    - box - a Box object
    - cutoff - the simulation cut off
    - max_displacement - the maximum displacement for a Monte Carlo move.
    - beta (1 / reduced_temperature)
    - total_pair_energy
    - tail_correction
    - total_energy
    
    The MCState class has the following methods.
    - calculate_total_pair_energy
        - **returns** the total pair energy 
    - calculate_tail_correction
        - **returns** the tail correction
    - calculate_total_energy
        - **returns** the total energy based on calculation using `calculate_total_pair_energy` and `calculate_tail_correction`.
    - get_particle_energy
        - takes a particle index as an argument. Optional arugment of particle_movement (this is a vector describing out to move the specified particle).
        - **returns** the pairwise interaction energy of a with all other particles.

1. The following functions are not associated with classes.
    - lennard_jones_potential
    - accept_or_reject
    - adjust_displacement
    - generate_initial_cooridinates

1. A Factory Design pattern is used with `generate_initial_coordinates` to specify method. This function now returns a tuple, `coordinates, box_length`

1. You must also refactor your tests.
    1. You will need to change the tests to call the appropriate objects.
    1. You can now create additional fixtures to return `Box` objects or `MCState` objects.

    You should already have a `fixture` called `nist_file`. Here is our definition. 

    ~~~
    @pytest.fixture
    def nist_file():
        nist_file = os.path.join('..','nist_sample_config1.txt')
        coordinates = mc.generate_initial_coordinates(method='file', fname=nist_file)
        return coordinates, nist_file
    ~~~

    You can also write a fixture to create a box. You can build this from the nist file fixture by passing the `nist_file` fixture as an argument to the new function. For example, our code looks like this. 
    
    ~~~
    @pytest.fixture
    def mc_box(nist_file):
        coordinates = nist_file[0][0]
        box_length = nist_file[0][1]
        fname = nist_file[1]

        test_box = mc.Box(box_length, coordinates)

        return test_box
    ~~~