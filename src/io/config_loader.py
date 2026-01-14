"""
Configuration file loader for simulation parameters.
"""

import yaml
from typing import Dict, Any


def load_config(filename: str) -> Dict[str, Any]:
    """
    Load simulation configuration from YAML file.
    
    Parameters
    ----------
    filename : str
        Path to YAML configuration file
    
    Returns
    -------
    dict
        Configuration dictionary
    """
    with open(filename, 'r') as f:
        config = yaml.safe_load(f)
    
    return config


def validate_config(config: Dict[str, Any]) -> None:
    """
    Validate configuration parameters.
    
    Parameters
    ----------
    config : dict
        Configuration dictionary
    
    Raises
    ------
    ValueError
        If required parameters are missing or invalid
    """
    required_keys = ['natoms', 'mass', 'temperature', 'dt', 'relaxation_time', 
                    'nsteps', 'box']
    
    for key in required_keys:
        if key not in config:
            raise ValueError(f"Required parameter '{key}' missing from config")
    
    if config['natoms'] <= 0:
        raise ValueError("natoms must be positive")
    
    if config['temperature'] < 0:
        raise ValueError("temperature cannot be negative")
    
    if config['dt'] <= 0:
        raise ValueError("dt must be positive")