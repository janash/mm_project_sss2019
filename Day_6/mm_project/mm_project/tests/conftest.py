"""
Fixtures for monte carlo tests
"""

# Import package, test suite, and other packages as needed
import mm_project as mc
import numpy as np
import os
import pytest
import sys

@pytest.fixture
def nist_file():
    current_directory = os.path.dirname(os.path.abspath(__file__))
    nist_file = os.path.join(current_directory,'..', 'data', 'nist_sample_config1.txt')
    coordinates = mc.generate_initial_coordinates(method='file', fname=nist_file)
    return coordinates, nist_file

@pytest.fixture
def mc_box(nist_file):
    coordinates = nist_file[0][0]
    box_length = nist_file[0][1]
    fname = nist_file[1]

    test_box = mc.Box(box_length, coordinates)

    return test_box