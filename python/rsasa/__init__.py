"""
rsasa - Fast SASA (Solvent Accessible Surface Area) calculator

This module provides a fast implementation of the Shrake-Rupley algorithm
for calculating solvent accessible surface area of molecular structures.
"""

from typing import Optional
import numpy as np
from numpy.typing import NDArray

# Import the Rust extension module
from .rsasa import sasa as _sasa_rust

__version__ = "0.1.0"
__all__ = ["sasa"]


def sasa(
    n_frames: int,
    n_atoms: int,
    xyzlist: NDArray[np.float32],
    atom_radii: NDArray[np.float32],
    n_sphere_points: int,
    atom_mapping: NDArray[np.uintp],
    atom_selection_mask: NDArray[np.int32],
    n_groups: int,
) -> NDArray[np.float32]:
    """
    Calculate the solvent accessible surface area (SASA) for atoms in a trajectory.
    
    Parameters
    ----------
    n_frames : int
        Number of frames in the trajectory
    n_atoms : int
        Number of atoms in each frame
    xyzlist : ndarray of float32, shape (n_frames * n_atoms * 3,)
        Flattened array of atomic coordinates in Angstroms.
        Coordinates should be ordered as [x1, y1, z1, x2, y2, z2, ...]
        for all atoms in frame 1, then all atoms in frame 2, etc.
    atom_radii : ndarray of float32, shape (n_atoms,)
        Van der Waals radii plus probe radius for each atom in Angstroms.
        Typically, use VDW radius + 1.4 Å for water probe.
    n_sphere_points : int
        Number of points to generate on the unit sphere for sampling.
        Common values: 92 (fast), 960 (balanced), 9600 (high accuracy)
    atom_mapping : ndarray of uintp, shape (n_atoms,)
        Mapping from atoms to groups for accumulating SASA.
        Each value should be in range [0, n_groups).
    atom_selection_mask : ndarray of int32, shape (n_atoms,)
        Binary mask indicating which atoms to compute SASA for.
        1 = compute SASA, 0 = skip
    n_groups : int
        Number of groups in the output
    
    Returns
    -------
    ndarray of float32, shape (n_frames * n_groups,)
        SASA values in Ų (square Angstroms) for each group in each frame.
        Results are ordered as [group0_frame0, group1_frame0, ..., 
        group0_frame1, group1_frame1, ...]
    
    Raises
    ------
    ValueError
        If input arrays have incorrect shapes or invalid parameter values
    
    Examples
    --------
    Calculate SASA for a single atom:
    
    >>> import numpy as np
    >>> from rsasa import sasa
    >>> xyzlist = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    >>> atom_radii = np.array([1.5], dtype=np.float32)
    >>> atom_mapping = np.array([0], dtype=np.uintp)
    >>> atom_selection_mask = np.array([1], dtype=np.int32)
    >>> result = sasa(1, 1, xyzlist, atom_radii, 960, 
    ...               atom_mapping, atom_selection_mask, 1)
    >>> print(f"SASA: {result[0]:.2f} Ų")
    
    Calculate SASA for multiple atoms with grouping:
    
    >>> xyzlist = np.array([0.0, 0.0, 0.0, 1.8, 0.0, 0.0], dtype=np.float32)
    >>> atom_radii = np.array([1.5, 1.5], dtype=np.float32)
    >>> atom_mapping = np.array([0, 0], dtype=np.uintp)  # Both atoms in group 0
    >>> atom_selection_mask = np.array([1, 1], dtype=np.int32)
    >>> result = sasa(1, 2, xyzlist, atom_radii, 960,
    ...               atom_mapping, atom_selection_mask, 1)
    >>> print(f"Total SASA: {result[0]:.2f} Ų")
    """
    # Validate and convert inputs to correct types
    xyzlist = np.ascontiguousarray(xyzlist, dtype=np.float32)
    atom_radii = np.ascontiguousarray(atom_radii, dtype=np.float32)
    atom_mapping = np.ascontiguousarray(atom_mapping, dtype=np.uintp)
    atom_selection_mask = np.ascontiguousarray(atom_selection_mask, dtype=np.int32)
    
    # Validate shapes
    expected_xyz_len = n_frames * n_atoms * 3
    if xyzlist.shape[0] != expected_xyz_len:
        raise ValueError(
            f"xyzlist has length {xyzlist.shape[0]} but expected "
            f"{expected_xyz_len} (n_frames={n_frames} * n_atoms={n_atoms} * 3)"
        )
    
    if atom_radii.shape[0] != n_atoms:
        raise ValueError(
            f"atom_radii has length {atom_radii.shape[0]} but expected {n_atoms}"
        )
    
    if atom_mapping.shape[0] != n_atoms:
        raise ValueError(
            f"atom_mapping has length {atom_mapping.shape[0]} but expected {n_atoms}"
        )
    
    if atom_selection_mask.shape[0] != n_atoms:
        raise ValueError(
            f"atom_selection_mask has length {atom_selection_mask.shape[0]} "
            f"but expected {n_atoms}"
        )
    
    # Call the Rust implementation
    result = _sasa_rust(
        n_frames,
        n_atoms,
        xyzlist,
        atom_radii,
        n_sphere_points,
        atom_mapping,
        atom_selection_mask,
        n_groups,
    )
    
    return result