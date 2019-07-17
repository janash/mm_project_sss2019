#!/bin/bash

################################
# A test script for Linux X86_64
################################

# Download conda installation script
wget https://repo.anaconda.com/archive/Anaconda3-2019.03-Linux-x86_64.sh

# Run installation script
bash Anaconda3-2019.03-Linux-x86_64.sh

# Initialize conda
source ~/.bashrc

# Install packages
conda install pybind11 cmake eigen

# Compile C++/python module (mm_cpp.so)
bash mm_cpp/build.sh

# Move module to current directory
mv mm_cpp/mm_cpp.so .
