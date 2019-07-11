"""
Unit and regression test for the mm_project package.
"""

# Import package, test suite, and other packages as needed
import mm_project as mc
import numpy as np
import os
import pytest
import sys

def test_mm_project_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mm_project" in sys.modules

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

def test_generate_initial_coordinates_random():
    test_particles = 100
    box_length=10
    coords, box_length = mc.generate_initial_coordinates(method='random', num_particles=test_particles, box_length=box_length)
    assert test_particles == len(coords)
    assert box_length == 10
    
    # Check that coords are within bounds
    max_coord = np.max(coords)
    assert box_length/2 > max_coord, "Coordinates found outside box."
    
    min_coord = np.min(coords)
    assert -box_length/2 < min_coord

def test_generate_inital_coordinates_file(nist_file):
    
    coords = nist_file[0][0]
    box_length = nist_file[0][1]
    test_file = nist_file[1]

    with open(test_file) as f:
        data = f.readlines()
        num_particles = int(data[0])
        first_coord = [float(x) for x in data[2].split()][1:]
        last_coord = [float(x) for x in data[-1].split()][1:]

    assert len(coords) == num_particles
    assert box_length == 10
    assert np.all(np.equal(first_coord, coords[0]))
    assert np.all(np.equal(last_coord, coords[-1]))

@pytest.mark.parametrize("method, options, error", [
    ("random", {"num_particles": None, "box_length": 10}, ValueError),
    ("random", {"num_particles": 100, "box_length": None}, ValueError),
])
def test_generate_initial_coordinates_random_error(method, options, error):
    """
    Test failure for generate initial state - we are not giving specific directions for how to do this, so they may write several functions.
    """
    with pytest.raises(error):
        mc.generate_initial_coordinates(method=method, num_particles=options["num_particles"], box_length=options["box_length"])

@pytest.mark.parametrize("method, filename, error", [
    ("file", "this_file_doesnt_exist.txt", OSError),
    ("file", None, ValueError)
])
def test_generate_initial_state_file_error(method, filename, error):
    """
    Test failure for generate initial state - we are not giving specific directions for how to do this, so they may write several functions.
    """
    with pytest.raises(error):
        mc.generate_initial_coordinates(method=method, fname=filename)

@pytest.mark.parametrize("ri, rj, box_length, expected_distance", [
    (np.array([0,0,0]), np.array([0,1,0]), 10, 1),
    (np.array([0,0,4]), np.array([0,0,-4]), 10, 4),
]
)
def test_minimum_image_distance(ri, rj, box_length, expected_distance):
    test_box = mc.Box(box_length)
    assert expected_distance == test_box.minimum_image_distance(ri, rj)

