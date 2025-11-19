# rsasa

A Rust library for calculating Solvent Accessible Surface Area (SASA) of molecular structures using the Shrake-Rupley algorithm.

## Overview

This library implements the Shrake-Rupley algorithm for computing the solvent accessible surface area of atoms in molecular structures. SASA is an important metric in computational chemistry and structural biology for understanding protein folding, molecular interactions, and surface properties.

## Features

- Fast SASA calculation using the Shrake-Rupley algorithm
- Support for multiple frames/trajectories
- Configurable sphere point density for accuracy vs. performance trade-offs
- Atom selection masks for computing SASA on specific subsets
- Group-based accumulation for residue-level or custom grouping

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
rsasa = "0.1.0"
```

## Usage

```rust
use rsasa::sasa;

fn main() {
    // Define molecular structure
    let n_frames = 1;
    let n_atoms = 3;
    
    // Atom coordinates: [x1, y1, z1, x2, y2, z2, x3, y3, z3]
    let xyzlist = vec![
        0.0, 0.0, 0.0,  // Atom 1
        1.5, 0.0, 0.0,  // Atom 2
        0.0, 1.5, 0.0,  // Atom 3
    ];
    
    // Van der Waals radii + probe radius (typically 1.4 Å for water)
    let atom_radii = vec![1.0, 1.0, 1.0];
    
    // Number of sphere points (higher = more accurate but slower)
    let n_sphere_points = 960;
    
    // Map atoms to groups (e.g., for per-residue SASA)
    let atom_mapping = vec![0, 0, 1]; // Atoms 1-2 in group 0, atom 3 in group 1
    
    // Selection mask: 1 = compute SASA, 0 = skip
    let atom_selection_mask = vec![1, 1, 1];
    
    let n_groups = 2;
    
    // Output buffer
    let mut out = vec![0.0; n_frames * n_groups];
    
    // Calculate SASA
    sasa(
        n_frames,
        n_atoms,
        &xyzlist,
        &atom_radii,
        n_sphere_points,
        &atom_mapping,
        &atom_selection_mask,
        n_groups,
        &mut out,
    );
    
    println!("SASA for group 0: {:.2} Ų", out[0]);
    println!("SASA for group 1: {:.2} Ų", out[1]);
}
```

## API Documentation

### `sasa` Function

```rust
pub fn sasa(
    n_frames: usize,
    n_atoms: usize,
    xyzlist: &[f32],
    atom_radii: &[f32],
    n_sphere_points: usize,
    atom_mapping: &[usize],
    atom_selection_mask: &[i32],
    n_groups: usize,
    out: &mut [f32],
)
```

Calculate the solvent accessible surface area (SASA) for atoms in a trajectory.

#### Parameters

- `n_frames`: Number of frames in the trajectory
- `n_atoms`: Number of atoms in each frame
- `xyzlist`: Flattened 3D array of shape `[n_frames * n_atoms * 3]` containing atomic coordinates
- `atom_radii`: Array of length `[n_atoms]` containing van der Waals radii plus probe radius
- `n_sphere_points`: Number of points to generate sampling the unit sphere (recommended: 960)
- `atom_mapping`: Array of length `[n_atoms]` mapping atoms to groups for accumulation
- `atom_selection_mask`: Array of length `[n_atoms]` indicating which atoms to compute SASA for (1 = yes, 0 = no)
- `n_groups`: Number of groups in the output
- `out`: Output buffer of length `[n_frames * n_groups]` to store results

#### Notes

- Coordinates should be in Ångströms (Å)
- `atom_radii` should include the probe radius (typically 1.4 Å for water)
- Higher `n_sphere_points` values increase accuracy but decrease performance
- Common values: 92 (fast), 960 (balanced), 9600 (high accuracy)

## Testing

Run the test suite:

```bash
cargo test
```

Run with output:

```bash
cargo test -- --nocapture
```

## Algorithm

This implementation uses the Shrake-Rupley algorithm:

1. Generate uniformly distributed points on a unit sphere using the Golden Section Spiral
2. For each atom, scale the sphere points by the atom's radius
3. Check each point for occlusion by neighboring atoms
4. Count accessible points and compute surface area

The algorithm has O(n²) complexity for n atoms, with optimization for neighbor searching.

## Performance

The library is optimized for performance with:
- Efficient neighbor finding using radius cutoffs
- Minimal allocations with pre-allocated work buffers
- Smart iteration strategies for occlusion testing

## License

Licensed under either of:

- Apache License, Version 2.0
- MIT license

at your option.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## References

- Shrake, A.; Rupley, J. A. (1973). "Environment and exposure to solvent of protein atoms. Lysozyme and insulin". Journal of Molecular Biology. 79 (2): 351–371.