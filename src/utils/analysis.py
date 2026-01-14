"""
Analysis utilities for molecular dynamics simulations.
"""

import numpy as np
from src.constants import BOLTZMANN


def compute_temperature(mass: np.ndarray, vels: np.ndarray) -> float:
    """
    Compute instantaneous temperature from velocities.
    
    Uses the equipartition theorem:
    (1/2) * m * <v²> = (3/2) * k_B * T
    
    Parameters
    ----------
    mass : np.ndarray
        Particle masses, shape (natoms,)
    vels : np.ndarray
        Particle velocities, shape (natoms, ndims)
    
    Returns
    -------
    float
        Temperature in Kelvin
    """
    natoms, ndims = vels.shape
    
    # Center-of-mass velocity
    v_cm = np.mean(vels, axis=0)
    
    # Velocities relative to center of mass
    v_rel = vels - v_cm
    
    # Kinetic energy in center-of-mass frame
    ke = 0.5 * np.sum(mass[:, np.newaxis] * v_rel**2)
    
    # Temperature from equipartition: KE = (ndims * natoms / 2) * k_B * T
    temperature = 2.0 * ke / (ndims * natoms * BOLTZMANN)
    
    return temperature


def compute_kinetic_energy(mass: np.ndarray, vels: np.ndarray) -> float:
    """
    Compute total kinetic energy.
    
    Parameters
    ----------
    mass : np.ndarray
        Particle masses, shape (natoms,)
    vels : np.ndarray
        Particle velocities, shape (natoms, ndims)
    
    Returns
    -------
    float
        Total kinetic energy in Joules
    """
    return 0.5 * np.sum(mass[:, np.newaxis] * vels**2)


def compute_total_momentum(mass: np.ndarray, vels: np.ndarray) -> np.ndarray:
    """
    Compute total momentum vector.
    
    Parameters
    ----------
    mass : np.ndarray
        Particle masses, shape (natoms,)
    vels : np.ndarray
        Particle velocities, shape (natoms, ndims)
    
    Returns
    -------
    np.ndarray
        Total momentum vector, shape (ndims,)
    """
    return np.sum(mass[:, np.newaxis] * vels, axis=0)


def running_average(data: np.ndarray, window: int) -> np.ndarray:
    """
    Compute running average of data.
    
    Parameters
    ----------
    data : np.ndarray
        Input data
    window : int
        Window size for averaging
    
    Returns
    -------
    np.ndarray
        Running average
    """
    return np.convolve(data, np.ones(window)/window, mode='valid')