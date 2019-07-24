# Day 5 - Introduction to C++11

By the end of Day 5, students will learn:

- The difference between compilation and interpretation
- The basics of C++: std input/output, functions, and classes
- How to compile a C++ program with g++/clang
- How to setup and use *conda* to install C++ packages
- How to install and use *Eigen* to perform numerical linear algebra operations
- Use cProfile to determine bottlenecks in their python MM/QM codes

Using Eigen, students will be asked at the end of the day to re-write the methods that calculate the potential energy (`get_particle_energy` for MM and `calculate_fock_matrix` for QM) in C++. They are also free to work on additional intensive functions such as: `calculate_total_pair_energy` for MM, and `calculate_energy_mp2` and `calculate_fock_matrix_fast` for QM. They should be able to compile the C++ code, which they will use on Day 6 to extend their python MM/QM code.
