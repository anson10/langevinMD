"""Tests for boundary conditions."""

import numpy as np
import pytest
from src.boundary import ReflectiveBoundary, PeriodicBoundary


def test_reflective_boundary():
    """Test reflective boundary conditions."""
    box = ((0, 10), (0, 10), (0, 10))
    boundary = ReflectiveBoundary(box)
    
    # Particle crossing lower boundary
    pos = np.array([[-1.0, 5.0, 5.0]])
    vels = np.array([[-1.0, 0.0, 0.0]])
    
    boundary.apply(pos, vels)
    
    assert pos[0, 0] == 0.0  # Position at boundary
    assert vels[0, 0] > 0  # Velocity reversed
    

def test_periodic_boundary():
    """Test periodic boundary conditions."""
    box = ((0, 10), (0, 10), (0, 10))
    boundary = PeriodicBoundary(box)
    
    # Particle crossing upper boundary
    pos = np.array([[11.0, 5.0, 5.0]])
    vels = np.array([[1.0, 0.0, 0.0]])
    
    boundary.apply(pos, vels)
    
    assert 0 <= pos[0, 0] < 10  # Wrapped back into box