"""
LAMMPS dump file writer for trajectory output.
"""

import numpy as np
from typing import Dict, Tuple, Any


class LAMMPSWriter:
    """
    Writer for LAMMPS custom dump format.
    
    Creates trajectory files compatible with OVITO and other
    visualization tools that read LAMMPS format.
    """
    
    def __init__(self, filename: str):
        """
        Initialize LAMMPS writer.
        
        Parameters
        ----------
        filename : str
            Output file path
        """
        self.filename = filename
        self.axis = ('x', 'y', 'z')
    
    def write_frame(self, 
                   timestep: int,
                   natoms: int, 
                   box: Tuple[Tuple[float, float], ...],
                   **data: Dict[str, np.ndarray]) -> None:
        """
        Write a single timestep frame to the dump file.
        
        Parameters
        ----------
        timestep : int
            Current timestep number
        natoms : int
            Number of atoms
        box : tuple of tuples
            Box dimensions as ((xmin, xmax), (ymin, ymax), (zmin, zmax))
        **data : dict
            Additional data arrays (pos, v, radius, etc.)
            Arrays can be 1D (natoms,) or 2D (natoms, ndims)
        """
        # Process data: split 2D arrays into separate columns
        processed_data = self._process_data(data)
        
        with open(self.filename, 'a') as f:
            # Write header
            self._write_header(f, timestep, natoms, box, processed_data.keys())
            
            # Write atom data
            self._write_atoms(f, natoms, processed_data)
    
    def _process_data(self, data: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """
        Process data arrays, splitting 2D arrays into column vectors.
        
        Parameters
        ----------
        data : dict
            Input data dictionary
        
        Returns
        -------
        dict
            Processed data with 1D arrays only
        """
        processed = {}
        
        for key, array in data.items():
            if len(array.shape) > 1:
                # 2D array: split into columns
                nrows, ncols = array.shape
                for i in range(ncols):
                    if key == 'pos':
                        # Special case: position becomes x, y, z
                        processed[self.axis[i]] = array[:, i]
                    else:
                        # Other arrays: key_x, key_y, key_z
                        processed[f'{key}_{self.axis[i]}'] = array[:, i]
            else:
                # 1D array: keep as is
                processed[key] = array
        
        return processed
    
    def _write_header(self, 
                     f, 
                     timestep: int,
                     natoms: int,
                     box: Tuple[Tuple[float, float], ...],
                     keys: Any) -> None:
        """Write the header section of a frame."""
        f.write('ITEM: TIMESTEP\n')
        f.write(f'{timestep}\n')
        
        f.write('ITEM: NUMBER OF ATOMS\n')
        f.write(f'{natoms}\n')
        
        f.write('ITEM: BOX BOUNDS pp pp pp\n')
        for box_bounds in box:
            f.write(f'{box_bounds[0]} {box_bounds[1]}\n')
        
        # Column headers
        header_str = 'id ' + ' '.join(keys)
        f.write(f'ITEM: ATOMS {header_str}\n')
    
    def _write_atoms(self, 
                    f,
                    natoms: int,
                    data: Dict[str, np.ndarray]) -> None:
        """Write atom data rows."""
        atom_ids = np.arange(1, natoms + 1)
        keys = list(data.keys())
        
        # Stack all columns
        output_data = np.column_stack([atom_ids] + [data[key] for key in keys])
        
        # Write rows
        for row in output_data:
            f.write(f'{int(row[0]):d}')  # Atom ID
            for val in row[1:]:
                f.write(f' {val:e}')  # Data in scientific notation
            f.write('\n')