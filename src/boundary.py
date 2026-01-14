"""
Boundary condition implementations for molecular dynamics simulations.
"""

import numpy as np
from typing import Tuple


class BoundaryCondition:
    """Base class for boundary conditions."""
    
    def __init__(self, box: Tuple[Tuple[float, float], ...]):
        """
        Initialize boundary condition.
        
        Parameters
        ----------
        box : tuple of tuples
            Box dimensions as ((xmin, xmax), (ymin, ymax), (zmin, zmax))
        """
        self.box = box
        self.ndims = len(box)
    
    def apply(self, pos: np.ndarray, vels: np.ndarray) -> None:
        """
        Apply boundary conditions to positions and velocities.
        
        Parameters
        ----------
        pos : np.ndarray
            Particle positions, shape (natoms, ndims)
        vels : np.ndarray
            Particle velocities, shape (natoms, ndims)
        """
        raise NotImplementedError


class ReflectiveBoundary(BoundaryCondition):
    """
    Reflective boundary conditions.
    
    Particles that hit the boundary have their velocity component
    perpendicular to the wall reversed.
    """
    
    def apply(self, pos: np.ndarray, vels: np.ndarray) -> None:
        """
        Apply reflective boundary conditions.
        
        When a particle crosses a boundary, its velocity perpendicular
        to that boundary is reversed.
        
        Parameters
        ----------
        pos : np.ndarray
            Particle positions, shape (natoms, ndims)
        vels : np.ndarray
            Particle velocities, shape (natoms, ndims)
        """
        for i in range(self.ndims):
            # Check lower boundary
            mask_lower = pos[:, i] <= self.box[i][0]
            vels[mask_lower, i] = np.abs(vels[mask_lower, i])
            pos[mask_lower, i] = self.box[i][0]
            
            # Check upper boundary
            mask_upper = pos[:, i] >= self.box[i][1]
            vels[mask_upper, i] = -np.abs(vels[mask_upper, i])
            pos[mask_upper, i] = self.box[i][1]


class PeriodicBoundary(BoundaryCondition):
    """
    Periodic boundary conditions.
    
    Particles that cross one side of the box reappear on the opposite side.
    """
    
    def apply(self, pos: np.ndarray, vels: np.ndarray) -> None:
        """
        Apply periodic boundary conditions.
        
        When a particle crosses a boundary, it reappears on the opposite side.
        
        Parameters
        ----------
        pos : np.ndarray
            Particle positions, shape (natoms, ndims)
        vels : np.ndarray
            Particle velocities, shape (natoms, ndims)
        """
        for i in range(self.ndims):
            box_length = self.box[i][1] - self.box[i][0]
            
            # Wrap positions
            pos[:, i] = self.box[i][0] + (pos[:, i] - self.box[i][0]) % box_length


def create_boundary(boundary_type: str, box: Tuple[Tuple[float, float], ...]) -> BoundaryCondition:
    """
    Factory function to create boundary condition objects.
    
    Parameters
    ----------
    boundary_type : str
        Type of boundary condition ('reflective' or 'periodic')
    box : tuple of tuples
        Box dimensions
    
    Returns
    -------
    BoundaryCondition
        Boundary condition object
    
    Raises
    ------
    ValueError
        If boundary_type is not recognized
    """
    boundary_types = {
        'reflective': ReflectiveBoundary,
        'periodic': PeriodicBoundary,
    }
    
    if boundary_type.lower() not in boundary_types:
        raise ValueError(f"Unknown boundary type: {boundary_type}. "
                        f"Available types: {list(boundary_types.keys())}")
    
    return boundary_types[boundary_type.lower()](box)