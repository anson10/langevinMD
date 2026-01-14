"""
Time integration schemes for molecular dynamics.
"""

import numpy as np


def euler_integrate(pos: np.ndarray, 
                   vels: np.ndarray, 
                   force: np.ndarray,
                   mass: np.ndarray, 
                   dt: float) -> None:
    """
    Forward Euler integration scheme.
    
    Updates positions and velocities in-place using:
    r(t+dt) = r(t) + v(t) * dt
    v(t+dt) = v(t) + F(t)/m * dt
    
    Parameters
    ----------
    pos : np.ndarray
        Particle positions, shape (natoms, ndims)
        Modified in-place
    vels : np.ndarray
        Particle velocities, shape (natoms, ndims)
        Modified in-place
    force : np.ndarray
        Forces on particles, shape (natoms, ndims)
    mass : np.ndarray
        Particle masses, shape (natoms,)
    dt : float
        Time step in seconds
    """
    # Update positions: r(t+dt) = r(t) + v(t) * dt
    pos += vels * dt
    
    # Update velocities: v(t+dt) = v(t) + a(t) * dt
    # where a(t) = F(t) / m
    vels += force * dt / mass[:, np.newaxis]


def velocity_verlet_integrate(pos: np.ndarray,
                              vels: np.ndarray,
                              force_func,
                              mass: np.ndarray,
                              dt: float,
                              *args) -> None:
    """
    Velocity Verlet integration scheme (future feature).
    
    More accurate than Euler for conservative forces:
    r(t+dt) = r(t) + v(t)*dt + 0.5*a(t)*dt²
    v(t+dt) = v(t) + 0.5*[a(t) + a(t+dt)]*dt
    
    Parameters
    ----------
    pos : np.ndarray
        Particle positions, shape (natoms, ndims)
    vels : np.ndarray
        Particle velocities, shape (natoms, ndims)
    force_func : callable
        Function to compute forces
    mass : np.ndarray
        Particle masses, shape (natoms,)
    dt : float
        Time step
    *args
        Additional arguments for force_func
    """
    # TODO: Implement velocity Verlet
    raise NotImplementedError("Velocity Verlet not yet implemented")