# Introduction to C++

By the end of Day 6, the students will learn how to accelerate their python MM/QM code by re-writing the most (computationally) intensive functions in C++ and binding them to python via pybind11. The students will apply knowledge they acquired from the previous 4 days to:

- Install CMake and pybind11 with conda
- Compile a C++ program/library with CMake
- Profile python programs with cProfile and identify bottlenecks
- Bind C++11 functions with python using pybind11

The main objective is to show students how python can be slow for certain type of computations and how it can accelerated by extending the code base with a compiled, fast language such as C++.

## Code refactor
At the end of the day, the students will be asked to re-write certain functions in their MM/QM code in C++11 and profile the new code. 

## Extending python with C++
Students will 1st use cProfile to identify time-consuming functions in their python code. They will then write a new C++ file (`export.cpp`) that provides an interface for pybind11 to bind their C++ implementation (`day_6.QM/MM.cpp`) that uses eigen to provide equivalent numerical implementation for the following functions: 

### MM code
- `minimum_image_distance`
- `lennard_jones_potential`
- `get_particle_energy`
- `calculate_total_pair_energy`

### QM code
- `calculate_energy_mp2`
- `calculate_fock_matrix`: scales as O(N^4) 
- `calculate_fock_matrix_fast`: scales as O(N^2) but is much slower in python than `calculate_fock_matrix`

Students should observe that the C++ implementation of `calculate_fock_matrix_fast` is much more efficient than its python equivalent, and this leads to faster execution w.r.t `calculate_fock_matrix` in python when the number of particles is above 50 [?].
