"""Tests for force computations."""

import numpy as np
import pytest
from src.forces import compute_langevin_force


def test_langevin_force_shape():
    """Test that Langevin force has correct shape."""
    natoms = 100
    ndims = 3
    
    mass = np.ones(natoms)
    vels = np.random.randn(natoms, ndims)
    
    force = compute_langevin_force(mass, vels, 300, 1e-12, 1e-15)
    
    assert force.shape == (natoms, ndims)


def test_langevin_force_zero_velocity():
    """Test Langevin force with zero velocity."""
    mass = np.array([1.0])
    vels = np.zeros((1, 3))
    
    force = compute_langevin_force(mass, vels, 300, 1e-12, 1e-15)
    
    # With zero velocity, only random force remains
    # Check that force is non-zero (random term)
    assert not np.allclose(force, 0.0)