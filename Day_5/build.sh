#!/bin/bash

# Exit on any command returning non-zero
# Exit when attempting to use undefined variables
set -eu

# CONDA_PREFIX is the location of the current conda environment

# We need to add paths to the header files for Eigen
# Eigen3 gets installed to ${CONDA_PREFIX}/include/eigen3
EIGEN_PATH="${CONDA_PREFIX}/include/eigen3"

g++ -Wall -O2 -std=c++11 -I"${EIGEN_PATH}" -o sss_cpp *.cpp
