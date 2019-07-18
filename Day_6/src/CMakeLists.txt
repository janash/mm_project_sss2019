
# Create the library from the source files
add_library(mm_cpp MODULE day_6.cpp export.cpp)

# The output library should not have a prefix. By default,
# Modules/Libraries start with "lib". This removes that prefix.
set_target_properties(mm_cpp PROPERTIES PREFIX "")

# "Link" the library to pybind11 and eigen
# Link is not precisely the correct term. This also import
# paths to include files, etc
target_link_libraries(mm_cpp pybind11::module Eigen3::Eigen)
