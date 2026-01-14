"""
Example: Hydrogen gas simulation at 300 K

This script demonstrates a basic Langevin dynamics simulation
of hydrogen atoms in a periodic box.
"""

import argparse
import os
from src.simulation import LangevinSimulation
from src.constants import HYDROGEN_MASS, HYDROGEN_RADIUS, FEMTOSECOND, PICOSECOND
from src.io.config_loader import load_config


def run_from_config(config_file: str):
    """Run simulation from configuration file."""
    config = load_config(config_file)
    
    sim = LangevinSimulation(
        natoms=config['natoms'],
        mass=config['mass'],
        box=tuple(tuple(b) for b in config['box']),
        temperature=config['temperature'],
        dt=config['dt'],
        relaxation_time=config['relaxation_time'],
        boundary_type=config.get('boundary_type', 'reflective'),
        radius=config.get('radius', HYDROGEN_RADIUS)
    )
    
    results = sim.run(
        nsteps=config['nsteps'],
        output_file=config.get('output_file'),
        dump_frequency=config.get('dump_frequency', 100),
        verbose=True
    )
    
    # Generate plots
    if config.get('plot_temperature', True):
        output_dir = os.path.dirname(config.get('output_file', 'outputs'))
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        sim.plot_temperature(
            filename=os.path.join(output_dir, 'temperature.png'),
            show=False
        )
    
    return sim, results


def run_default():
    """Run simulation with default parameters."""
    print("=" * 60)
    print("Hydrogen Gas Simulation")
    print("=" * 60)
    print(f"Number of atoms: 1000")
    print(f"Temperature: 300 K")
    print(f"Box size: 10 nm × 10 nm × 10 nm")
    print(f"Timestep: 1 fs")
    print(f"Total time: 10 ps")
    print("=" * 60)
    print()
    
    # Create simulation
    sim = LangevinSimulation(
        natoms=1000,
        mass=HYDROGEN_MASS,
        box=((0, 1e-8), (0, 1e-8), (0, 1e-8)),  # 10 nm box
        temperature=300,  # K
        dt=FEMTOSECOND,  # 1 fs
        relaxation_time=PICOSECOND,  # 1 ps
        boundary_type='reflective',
        radius=HYDROGEN_RADIUS
    )
    
    # Create output directory
    if not os.path.exists('outputs'):
        os.makedirs('outputs')
    
    # Run simulation
    results = sim.run(
        nsteps=10000,  # 10 ps
        output_file='outputs/hydrogen_300K.dump',
        dump_frequency=100,
        verbose=True
    )
    
    # Generate plots
    print("\nGenerating plots...")
    sim.plot_temperature(filename='outputs/temperature.png', show=False)
    sim.plot_energy(filename='outputs/energy.png', show=False)
    
    print("\nVisualize trajectory with OVITO:")
    print("  ovito outputs/hydrogen_300K.dump")
    
    return sim, results


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Run Langevin dynamics simulation of hydrogen gas'
    )
    parser.add_argument(
        '--config',
        type=str,
        help='Path to YAML configuration file'
    )
    
    args = parser.parse_args()
    
    if args.config:
        sim, results = run_from_config(args.config)
    else:
        sim, results = run_default()


if __name__ == '__main__':
    main()