"""
Tests for Monte Carlo
"""

# Import package, test suite, and other packages as needed
import mm_project as mc
import numpy as np
import os
import pytest
import sys

@pytest.mark.parametrize("cutoff, nist_energy", [
    (3, -4.3515E+03),
    (4, -4.4675E+03)
])
def test_total_pair_energy(mc_box, cutoff, nist_energy):
    test_state = mc.MCState(mc_box, cutoff=cutoff, max_displacement=0.1, reduced_temperature=0.9)
    test_state.calculate_total_pair_energy()

    calculated_energy = test_state.total_pair_energy

    assert np.isclose(nist_energy, calculated_energy)

@pytest.mark.parametrize("cutoff, nist_energy", [
    (3, -4.3515E+03),
    (4, -4.4675E+03)
])
def test_total_pair_energy(mc_box, cutoff, nist_energy):
    test_state = mc.MCState(mc_box, cutoff=cutoff, max_displacement=0.1, reduced_temperature=0.9)
    test_state.calculate_total_pair_energy()

    calculated_energy = test_state.total_pair_energy

    assert np.isclose(nist_energy, calculated_energy)

@pytest.mark.parametrize("cutoff, nist_reference", [
    (3, -1.9849E+02),
    (4, -8.3769E+01),
]
)
def test_tail_correction(mc_box, cutoff, nist_reference):
    state = mc.MCState(mc_box, cutoff=cutoff, max_displacement=0.1, reduced_temperature=0.9)
    correction = state.calculate_tail_correction()
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




