"""Tests for integration schemes."""

import numpy as np
import pytest
from src.integrators import euler_integrate


def test_euler_integrate_conservation():
    """Test that Euler integration conserves energy for zero force."""
    natoms = 10
    ndims = 3
    dt = 1e-15
    
    pos = np.random.rand(natoms, ndims)
    vels = np.random.rand(natoms, ndims)
    mass = np.ones(natoms)
    force = np.zeros((natoms, ndims))
    
    # Compute initial kinetic energy
    ke_initial = 0.5 * np.sum(mass[:, np.newaxis] * vels**2)
    
    # Integrate
    euler_integrate(pos, vels, force, mass, dt)
    
    # Compute final kinetic energy
    ke_final = 0.5 * np.sum(mass[:, np.newaxis] * vels**2)
    
    # For zero force, KE should be conserved
    assert np.isclose(ke_initial, ke_final)


def test_euler_integrate_motion():
    """Test that particles move correctly."""
    pos = np.array([[0.0, 0.0, 0.0]])
    vels = np.array([[1.0, 0.0, 0.0]])
    mass = np.array([1.0])
    force = np.zeros((1, 3))
    dt = 1.0
    
    euler_integrate(pos, vels, force, mass, dt)
    
    assert np.allclose(pos, [[1.0, 0.0, 0.0]])