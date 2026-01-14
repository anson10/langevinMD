"""
Force computation for Langevin dynamics.
"""

import numpy as np
from src.constants import BOLTZMANN


def compute_langevin_force(mass: np.ndarray, 
                          vels: np.ndarray, 
                          temperature: float,
                          relaxation_time: float, 
                          dt: float) -> np.ndarray:
    """
    Compute Langevin force with friction and random noise.
    
    The Langevin force consists of two terms:
    1. Friction force: F_friction = -γ * m * v
    2. Random force: F_random ~ N(0, σ²)
    
    where γ = 1/τ (friction coefficient) and σ² satisfies the
    fluctuation-dissipation theorem:
    σ² = 2 * m * k_B * T / (τ * dt)
    
    Parameters
    ----------
    mass : np.ndarray
        Particle masses, shape (natoms,)
    vels : np.ndarray
        Particle velocities, shape (natoms, ndims)
    temperature : float
        Target temperature in Kelvin
    relaxation_time : float
        Relaxation time τ in seconds
    dt : float
        Time step in seconds
    
    Returns
    -------
    np.ndarray
        Total force on each particle, shape (natoms, ndims)
    """
    natoms, ndims = vels.shape
    
    # Friction force: -γ * m * v = -(m/τ) * v
    friction_force = -(mass[:, np.newaxis] * vels) / relaxation_time
    
    # Random force standard deviation (fluctuation-dissipation theorem)
    sigma = np.sqrt(2.0 * mass * temperature * BOLTZMANN / (relaxation_time * dt))
    
    # Generate random forces
    random_force = np.random.randn(natoms, ndims) * sigma[:, np.newaxis]
    
    # Total Langevin force
    total_force = friction_force + random_force
    
    return total_force


def compute_lennard_jones_force(pos: np.ndarray, 
                               epsilon: float, 
                               sigma: float,
                               cutoff: float = None) -> np.ndarray:
    """
    Compute Lennard-Jones forces between particles (future feature).
    
    V(r) = 4ε[(σ/r)^12 - (σ/r)^6]
    
    Parameters
    ----------
    pos : np.ndarray
        Particle positions, shape (natoms, ndims)
    epsilon : float
        Energy parameter
    sigma : float
        Distance parameter
    cutoff : float, optional
        Cutoff distance for force calculation
    
    Returns
    -------
    np.ndarray
        Forces on each particle, shape (natoms, ndims)
    """
    # TODO: Implement Lennard-Jones interactions
    raise NotImplementedError("Lennard-Jones forces not yet implemented")