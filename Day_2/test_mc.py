"""
Unit tests for MC code
"""

import os

import pytest
import numpy as np

import day_2 as mc

@pytest.fixture
def nist_file():
    nist_file = os.path.join('..','nist_sample_config1.txt')
    coordinates = mc.generate_initial_coordinates(method='file', fname=nist_file)
    return coordinates, nist_file

def test_generate_initial_coordinates_random():
    test_particles = 100
    box_length=10
    coords = mc.generate_initial_coordinates(method='random', num_particles=test_particles, box_length=box_length)
    assert test_particles == len(coords)
    
    # Check that coords are within bounds
    max_coord = np.max(coords)
    assert box_length/2 > max_coord, "Coordinates found outside box."
    
    min_coord = np.min(coords)
    assert -box_length/2 < min_coord

def test_generate_inital_coordinates_file(nist_file):
    
    coords = nist_file[0]
    test_file = nist_file[1]

    with open(test_file) as f:
        data = f.readlines()
        num_particles = int(data[0])
        first_coord = [float(x) for x in data[2].split()][1:]
        last_coord = [float(x) for x in data[-1].split()][1:]

    assert len(coords) == num_particles
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

@pytest.mark.parametrize("distance2, expected_energy",[
    (0.5 , 4*((1/0.5)**6 - (1/0.5)**3) ),
    (1 , 0),
    (1.5 , 4*((1/1.5)**6 - (1/1.5)**3) ),
    (2.5 , 4*((1/2.5)**6 - (1/2.5)**3) ),
    (3 , 4*((1/3.0)**6 - (1/3.0)**3) ),
]
)
def test_lennard_jones(distance2, expected_energy):
    assert expected_energy == mc.lennard_jones_potential(distance2)

@pytest.mark.parametrize("ri, rj, box_length, expected_distance", [
    (np.array([0,0,0]), np.array([0,1,0]), 10, 1),
    (np.array([0,0,4]), np.array([0,0,-4]), 10, 4),
]
)
def test_minimum_image_distance(ri, rj, box_length, expected_distance):
    assert expected_distance == mc.minimum_image_distance(ri, rj, box_length)

@pytest.mark.parametrize("cutoff, nist_energy", [
    (3, -4.3515E+03),
    (4, -4.4675E+03)
])
def test_total_potential_energy(nist_file, cutoff, nist_energy):
    coords = nist_file[0]
    test_file = nist_file[1]

    calculated_energy = mc.total_potential_energy(coords, 10, cutoff)

    assert np.isclose(nist_energy, calculated_energy)

@pytest.mark.parametrize("cutoff, nist_reference", [
    (3, -1.9849E+02),
    (4, -8.3769E+01),
]
)
def test_tail_correction(cutoff, nist_reference):
    correction = mc.tail_correction(10, cutoff, 800)
    # Compare to reference from NIST
    assert np.isclose(correction, nist_reference)

def test_accept_or_reject_negative_E():
    delta_energy = -1
    beta = 1/0.9
    assert mc.accept_or_reject(delta_energy, beta) is True

def test_accept_or_reject_0():
    delta_energy = 0
    beta = 1/0.9
    assert mc.accept_or_reject(delta_energy, beta) is True

@pytest.mark.parametrize("seed, expected_response", [
    (0, True),
    (8, False)
    ])
def test_accept_or_reject_positive_energy(seed, expected_response):
    # Set the random seed
    np.random.seed(seed)
    delta_energy = 0.5
    beta = 1/0.9
    try:
        assert mc.accept_or_reject(delta_energy, beta) == expected_response
    finally:
        # Reset seed - a best practice in case other tests follow.
        np.random.seed()




