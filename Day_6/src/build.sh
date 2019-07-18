#!/bin/bash

######################################################################
# This build script should not be used and is only here as an example
# By this point, the students should be using CMake
######################################################################


# Exit on any command returning non-zero
# Exit when attempting to use undefined variables
set -eu

# CONDA_PREFIX is the location of the current conda environment

# We need to add paths to the header files for Eigen and Pybind11
# Eigen3 gets installed to ${CONDA_PREFIX}/include/eigen3
# Pybind11 gets installed to to ${CONDA_PREFIX}/include/python3.7m

EIGEN_PATH="${CONDA_PREFIX}/include/eigen3"
PYBIND11_PATH="${CONDA_PREFIX}/include/python3.7m"

# We also need the path to the location of the python development library
PYTHON_PATH="${CONDA_PREFIX}/lib"

g++ -fPIC -O2 -shared -I"${EIGEN_PATH}" -I"${PYBIND11_PATH}" -L"${PYTHON_PATH}" -o mm_cpp.so *.cpp -lpython3.7m
