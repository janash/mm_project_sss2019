# Day 5 - Introduction to C++

The students should 

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
