# Import package, test suite, and other packages as needed
import mm_project as mc
import numpy as np
import os
import pytest
import sys

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
