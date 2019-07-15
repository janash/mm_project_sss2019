cmake_minimum_required(VERSION 3.6)
project(SSS_MM_CPP
        VERSION 0.1
        LANGUAGES CXX)


# Find the locations of pybind11 and Eigen
find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)


# The rest of the action happens in the mm_cpp subdirectory
add_subdirectory(mm_cpp)
