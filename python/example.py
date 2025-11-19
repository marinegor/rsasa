#!/usr/bin/env python3
"""
Example script demonstrating the use of rsasa Python module
"""

import numpy as np
from rsasa import sasa

def main():
    print("RSASA - Solvent Accessible Surface Area Calculator (Python)")
    print("=" * 60)
    print()

    # Example 1: Single isolated atom
    print("Example 1: Single isolated atom")
    print("-" * 40)
    n_frames = 1
    n_atoms = 1
    xyzlist = np.array([0.0, 0.0, 0.0], dtype=np.float32)
    atom_radii = np.array([1.5], dtype=np.float32)
    n_sphere_points = 960
    atom_mapping = np.array([0], dtype=np.uintp)
    atom_selection_mask = np.array([1], dtype=np.int32)
    n_groups = 1

    result = sasa(
        n_frames,
        n_atoms,
        xyzlist,
        atom_radii,
        n_sphere_points,
        atom_mapping,
        atom_selection_mask,
        n_groups,
    )

    expected_area = 4.0 * np.pi * atom_radii[0] ** 2
    print(f"  Calculated SASA: {result[0]:.2f} Ų")
    print(f"  Expected (full sphere): {expected_area:.2f} Ų")
    print(f"  Difference: {abs(result[0] - expected_area) / expected_area * 100:.2f}%")
    print()

    # Example 2: Two overlapping atoms
    print("Example 2: Two overlapping atoms")
    print("-" * 40)
    n_atoms = 2
    xyzlist = np.array([0.0, 0.0, 0.0, 1.8, 0.0, 0.0], dtype=np.float32)
    atom_radii = np.array([1.5, 1.5], dtype=np.float32)
    atom_mapping = np.array([0, 1], dtype=np.uintp)
    atom_selection_mask = np.array([1, 1], dtype=np.int32)
    n_groups = 2

    result = sasa(
        n_frames,
        n_atoms,
        xyzlist,
        atom_radii,
        n_sphere_points,
        atom_mapping,
        atom_selection_mask,
        n_groups,
    )

    full_sphere = 4.0 * np.pi * atom_radii[0] ** 2
    print(f"  Atom 1 SASA: {result[0]:.2f} Ų ({result[0] / full_sphere * 100:.1f}% of full sphere)")
    print(f"  Atom 2 SASA: {result[1]:.2f} Ų ({result[1] / full_sphere * 100:.1f}% of full sphere)")
    print(f"  Total SASA: {result[0] + result[1]:.2f} Ų")
    print()

    # Example 3: Three atoms with grouping
    print("Example 3: Three atoms with residue-level grouping")
    print("-" * 40)
    n_atoms = 3
    xyzlist = np.array(
        [0.0, 0.0, 0.0, 1.8, 0.0, 0.0, 0.9, 1.5, 0.0], dtype=np.float32
    )
    atom_radii = np.array([1.2, 1.4, 1.3], dtype=np.float32)
    atom_mapping = np.array([0, 0, 1], dtype=np.uintp)  # Atoms 1-2 in group 0, atom 3 in group 1
    atom_selection_mask = np.array([1, 1, 1], dtype=np.int32)
    n_groups = 2

    result = sasa(
        n_frames,
        n_atoms,
        xyzlist,
        atom_radii,
        n_sphere_points,
        atom_mapping,
        atom_selection_mask,
        n_groups,
    )

    print(f"  Group 0 (atoms 1-2) SASA: {result[0]:.2f} Ų")
    print(f"  Group 1 (atom 3) SASA: {result[1]:.2f} Ų")
    print(f"  Total SASA: {result[0] + result[1]:.2f} Ų")
    print()

    # Example 4: Selective computation
    print("Example 4: Selective SASA computation")
    print("-" * 40)
    n_atoms = 3
    xyzlist = np.array(
        [0.0, 0.0, 0.0, 1.8, 0.0, 0.0, 0.9, 1.5, 0.0], dtype=np.float32
    )
    atom_radii = np.array([1.2, 1.4, 1.3], dtype=np.float32)
    atom_mapping = np.array([0, 1, 2], dtype=np.uintp)
    atom_selection_mask = np.array([1, 0, 1], dtype=np.int32)  # Skip atom 2
    n_groups = 3

    result = sasa(
        n_frames,
        n_atoms,
        xyzlist,
        atom_radii,
        n_sphere_points,
        atom_mapping,
        atom_selection_mask,
        n_groups,
    )

    print(f"  Atom 1 SASA (computed): {result[0]:.2f} Ų")
    print(f"  Atom 2 SASA (skipped): {result[1]:.2f} Ų")
    print(f"  Atom 3 SASA (computed): {result[2]:.2f} Ų")
    print()
    print("Note: Atom 2 was skipped due to selection mask = 0")
    print()

    # Example 5: Multiple frames (trajectory)
    print("Example 5: Multiple frames (trajectory)")
    print("-" * 40)
    n_frames = 3
    n_atoms = 2
    # Three frames with atoms getting closer
    xyzlist = np.array(
        [
            # Frame 1: far apart
            0.0, 0.0, 0.0, 3.0, 0.0, 0.0,
            # Frame 2: closer
            0.0, 0.0, 0.0, 2.0, 0.0, 0.0,
            # Frame 3: overlapping
            0.0, 0.0, 0.0, 1.5, 0.0, 0.0,
        ],
        dtype=np.float32,
    )
    atom_radii = np.array([1.0, 1.0], dtype=np.float32)
    atom_mapping = np.array([0, 0], dtype=np.uintp)  # Both atoms in same group
    atom_selection_mask = np.array([1, 1], dtype=np.int32)
    n_groups = 1

    result = sasa(
        n_frames,
        n_atoms,
        xyzlist,
        atom_radii,
        n_sphere_points,
        atom_mapping,
        atom_selection_mask,
        n_groups,
    )

    for i in range(n_frames):
        print(f"  Frame {i + 1} total SASA: {result[i]:.2f} Ų")
    
    print()
    print("Note: SASA decreases as atoms get closer together")
    print()

    # Example 6: Using numpy arrays with different shapes
    print("Example 6: Working with 2D coordinate arrays")
    print("-" * 40)
    n_frames = 1
    n_atoms = 3
    
    # Start with a 2D array (n_atoms, 3) and flatten it
    coords_2d = np.array(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [1.0, 2.0, 0.0]], dtype=np.float32
    )
    xyzlist = coords_2d.flatten()
    
    atom_radii = np.array([1.0, 1.0, 1.0], dtype=np.float32)
    atom_mapping = np.array([0, 1, 2], dtype=np.uintp)
    atom_selection_mask = np.array([1, 1, 1], dtype=np.int32)
    n_groups = 3

    result = sasa(
        n_frames,
        n_atoms,
        xyzlist,
        atom_radii,
        n_sphere_points,
        atom_mapping,
        atom_selection_mask,
        n_groups,
    )

    print(f"  Input shape (before flatten): {coords_2d.shape}")
    print(f"  Flattened shape: {xyzlist.shape}")
    for i in range(n_atoms):
        print(f"  Atom {i + 1} SASA: {result[i]:.2f} Ų")
    print()


if __name__ == "__main__":
    main()