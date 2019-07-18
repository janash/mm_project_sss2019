#include "day_6.hpp"
#include "pybind11/pybind11.h"
#include "pybind11/eigen.h"


PYBIND11_MODULE(mm_cpp, m)
{
    m.doc() = "Functions for the MM project implemented in C++";
    
    m.def("get_particle_energy", get_particle_energy,
          "Calculate the interaction energy of a particle with its environment");

    m.def("calculate_total_pair_energy", calculate_total_pair_energy,
          "Calculate the total potential energy of the system using a pairwise Lennard Jones potential");
}
