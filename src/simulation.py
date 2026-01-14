"""
Main simulation class for Langevin dynamics.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import Tuple, Optional, Dict
import os

from src.constants import AVOGADRO
from src.boundary import create_boundary
from src.forces import compute_langevin_force
from src.integrators import euler_integrate
from src.io.lammps_writer import LAMMPSWriter
from src.utils.analysis import compute_temperature


class LangevinSimulation:
    """
    Langevin dynamics simulation for molecular systems.
    
    This class handles the setup and execution of molecular dynamics
    simulations with Langevin thermostat for temperature control.
    """
    
    def __init__(self,
                 natoms: int,
                 mass: float,
                 box: Tuple[Tuple[float, float], ...],
                 temperature: float,
                 dt: float,
                 relaxation_time: float,
                 boundary_type: str = 'reflective',
                 radius: float = 1e-10):
        """
        Initialize Langevin dynamics simulation.
        
        Parameters
        ----------
        natoms : int
            Number of atoms/particles
        mass : float
            Particle mass in kg/mol
        box : tuple of tuples
            Simulation box as ((xmin, xmax), (ymin, ymax), (zmin, zmax))
        temperature : float
            Target temperature in Kelvin
        dt : float
            Integration timestep in seconds
        relaxation_time : float
            Langevin thermostat relaxation time in seconds
        boundary_type : str, optional
            Type of boundary conditions ('reflective' or 'periodic')
        radius : float, optional
            Particle radius for visualization
        """
        self.natoms = natoms
        self.mass_per_particle = mass / AVOGADRO  # Convert to kg per particle
        self.mass = np.ones(natoms) * self.mass_per_particle
        self.box = box
        self.ndims = len(box)
        self.temperature = temperature
        self.dt = dt
        self.relaxation_time = relaxation_time
        self.radius = np.ones(natoms) * radius
        
        # Initialize boundary conditions
        self.boundary = create_boundary(boundary_type, box)
        
        # Initialize positions and velocities
        self._initialize_system()
        
        # Storage for trajectory data
        self.trajectory = {
            'time': [],
            'temperature': [],
            'kinetic_energy': [],
        }
    
    def _initialize_system(self) -> None:
        """Initialize particle positions and velocities."""
        # Random positions within box
        self.pos = np.random.rand(self.natoms, self.ndims)
        for i in range(self.ndims):
            self.pos[:, i] = (self.box[i][0] + 
                             (self.box[i][1] - self.box[i][0]) * self.pos[:, i])
        
        # Random velocities
        self.vels = np.random.randn(self.natoms, self.ndims)
        
        # Remove center-of-mass motion
        self.vels -= np.mean(self.vels, axis=0)
    
    def step(self) -> None:
        """Perform a single integration step."""
        # Compute Langevin forces
        force = compute_langevin_force(
            self.mass, 
            self.vels, 
            self.temperature,
            self.relaxation_time, 
            self.dt
        )
        
        # Integrate equations of motion
        euler_integrate(self.pos, self.vels, force, self.mass, self.dt)
        
        # Apply boundary conditions
        self.boundary.apply(self.pos, self.vels)
    
    def run(self,
            nsteps: int,
            output_file: Optional[str] = None,
            dump_frequency: int = 100,
            verbose: bool = True) -> Dict[str, np.ndarray]:
        """
        Run the simulation for a specified number of steps.
        
        Parameters
        ----------
        nsteps : int
            Number of timesteps to simulate
        output_file : str, optional
            Output trajectory file (LAMMPS format)
        dump_frequency : int, optional
            Frequency of trajectory dumps (every N steps)
        verbose : bool, optional
            Print progress information
        
        Returns
        -------
        dict
            Dictionary containing trajectory data (time, temperature, etc.)
        """
        # Initialize output writer
        writer = None
        if output_file:
            # Delete old file if exists
            if os.path.exists(output_file):
                os.remove(output_file)
            writer = LAMMPSWriter(output_file)
        
        # Main simulation loop
        for step in range(1, nsteps + 1):
            # Perform integration step
            self.step()
            
            # Compute observables
            current_temp = compute_temperature(self.mass, self.vels)
            current_time = step * self.dt
            
            # Store trajectory data
            self.trajectory['time'].append(current_time)
            self.trajectory['temperature'].append(current_temp)
            
            # Write trajectory frame
            if writer and step % dump_frequency == 0:
                writer.write_frame(
                    timestep=step,
                    natoms=self.natoms,
                    box=self.box,
                    radius=self.radius,
                    pos=self.pos,
                    v=self.vels
                )
            
            # Print progress
            if verbose and step % (nsteps // 10) == 0:
                print(f"Step {step}/{nsteps} | "
                      f"T = {current_temp:.2f} K | "
                      f"Time = {current_time*1e12:.2f} ps")
        
        # Convert trajectory lists to arrays
        for key in self.trajectory:
            self.trajectory[key] = np.array(self.trajectory[key])
        
        if verbose:
            avg_temp = np.mean(self.trajectory['temperature'])
            print(f"\nSimulation complete!")
            print(f"Average temperature: {avg_temp:.2f} K")
            if output_file:
                print(f"Trajectory saved to: {output_file}")
        
        return self.trajectory
    
    def plot_temperature(self, 
                    filename: Optional[str] = None,
                    show: bool = True) -> None:
        """
        Plot temperature evolution over time.
        
        Parameters
        ----------
        filename : str, optional
            Save plot to file
        show : bool, optional
            Display plot
        """
        if len(self.trajectory['time']) == 0:
            print("No trajectory data available. Run simulation first.")
            return
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        time_ps = np.array(self.trajectory['time']) * 1e12  # Convert to ps
        temp = np.array(self.trajectory['temperature'])
        
        ax.plot(time_ps, temp, 'b-', alpha=0.7, label='Instantaneous')
        ax.axhline(self.temperature, color='r', linestyle='--', 
                label=f'Target ({self.temperature} K)')
        
        # Running average
        window = min(100, len(temp) // 10)
        if window > 1:
            temp_avg = np.convolve(temp, np.ones(window)/window, mode='valid')
            time_avg = time_ps[window-1:]
            ax.plot(time_avg, temp_avg, 'g-', linewidth=2, label='Running average')
        
        ax.set_xlabel('Time (ps)', fontsize=12)
        ax.set_ylabel('Temperature (K)', fontsize=12)
        ax.set_title('Temperature Evolution', fontsize=14)
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if filename:
            plt.savefig(filename, dpi=300)
            print(f"Plot saved to: {filename}")
        
        if show:
            plt.show()
        else:
            plt.close()
    
    def plot_energy(self,
               filename: Optional[str] = None,
               show: bool = True) -> None:
        """
        Plot kinetic energy over time.
        
        Parameters
        ----------
        filename : str, optional
            Save plot to file
        show : bool, optional
            Display plot
        """
        if len(self.trajectory['time']) == 0:
            print("No trajectory data available. Run simulation first.")
            return
        
        # Compute kinetic energy from temperature
        from src.constants import BOLTZMANN
        ke = 1.5 * self.natoms * BOLTZMANN * np.array(self.trajectory['temperature'])
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        time_ps = np.array(self.trajectory['time']) * 1e12
        ax.plot(time_ps, ke * 1e18, 'b-', alpha=0.7)  # Convert to aJ
        
        ax.set_xlabel('Time (ps)', fontsize=12)
        ax.set_ylabel('Kinetic Energy (aJ)', fontsize=12)
        ax.set_title('Kinetic Energy Evolution', fontsize=14)
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        
        if filename:
            plt.savefig(filename, dpi=300)
            print(f"Plot saved to: {filename}")
        
        if show:
            plt.show()
        else:
            plt.close()