# rsasa

A high-performance library for calculating Solvent Accessible Surface Area (SASA) of molecular structures using the Shrake-Rupley algorithm.

Available as both a **Rust library** and a **Python package** with numpy integration.

## Overview

This library implements the Shrake-Rupley algorithm for computing the solvent accessible surface area of atoms in molecular structures. SASA is an important metric in computational chemistry and structural biology for understanding protein folding, molecular interactions, and surface properties.

The core implementation is written in Rust for maximum performance, with Python bindings via PyO3 for easy integration with scientific Python workflows.

## Features

- ⚡ **Fast**: Rust implementation with optimized algorithms
- 🐍 **Python Support**: Native numpy integration via PyO3
- 📊 **Accurate**: Standard Shrake-Rupley algorithm with configurable sphere point density
- 🔄 **Trajectory Support**: Process multiple frames efficiently
- 🎯 **Flexible**: Atom selection masks and group-based accumulation
- 🛡️ **Type Safe**: Strong typing and input validation with Result types
- 📦 **Zero Dependencies**: Rust library has no dependencies (Python package requires numpy)

---

## Python Usage

### Installation

#### From PyPI (when published)

```bash
pip install rsasa
```

#### From Source

Requirements:
- Python >= 3.8
- Rust toolchain (install from https://rustup.rs/)
- maturin (`pip install maturin`)

```bash
git clone <repository-url>
cd rsasa
maturin develop --features python
```

### Quick Start

```python
import numpy as np
from rsasa import sasa

# Define a simple system with two atoms
n_frames = 1
n_atoms = 2
xyzlist = np.array([0.0, 0.0, 0.0, 1.8, 0.0, 0.0], dtype=np.float32)
atom_radii = np.array([1.5, 1.5], dtype=np.float32)  # VDW + probe radius
n_sphere_points = 960
atom_mapping = np.array([0, 1], dtype=np.uintp)
atom_selection_mask = np.array([1, 1], dtype=np.int32)
n_groups = 2

result = sasa(
    n_frames, n_atoms, xyzlist, atom_radii,
    n_sphere_points, atom_mapping, atom_selection_mask, n_groups
)

print(f"Atom 1 SASA: {result[0]:.2f} Ų")
print(f"Atom 2 SASA: {result[1]:.2f} Ų")
```

### Python API

```python
def sasa(
    n_frames: int,
    n_atoms: int,
    xyzlist: NDArray[np.float32],      # shape: (n_frames * n_atoms * 3,)
    atom_radii: NDArray[np.float32],   # shape: (n_atoms,)
    n_sphere_points: int,
    atom_mapping: NDArray[np.uintp],   # shape: (n_atoms,)
    atom_selection_mask: NDArray[np.int32],  # shape: (n_atoms,)
    n_groups: int,
) -> NDArray[np.float32]:  # shape: (n_frames * n_groups,)
```

See `python/README.md` for detailed Python documentation and examples.

---

## Rust Usage

### Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
rsasa = "0.1.0"
```

### Basic Usage

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

### Rust API Documentation

The main function returns a `Result` for proper error handling:

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
) -> Result<Vec<f32>, SasaError>
```

#### Parameters

- `n_frames`: Number of frames in the trajectory
- `n_atoms`: Number of atoms in each frame
- `xyzlist`: Flattened array of coordinates `[n_frames * n_atoms * 3]` in Ångströms
- `atom_radii`: Van der Waals radii + probe radius `[n_atoms]`
- `n_sphere_points`: Sphere points for sampling (92/960/9600 recommended)
- `atom_mapping`: Atom-to-group mapping `[n_atoms]`
- `atom_selection_mask`: Binary mask `[n_atoms]` (1 = compute, 0 = skip)
- `n_groups`: Number of output groups

#### Returns

- `Ok(Vec<f32>)`: SASA values `[n_frames * n_groups]` in Ų
- `Err(SasaError)`: Shape mismatch or invalid parameters

#### Error Types

```rust
pub enum SasaError {
    ShapeMismatch(String),    // Array dimensions don't match
    InvalidParameter(String), // Invalid parameter values
}
```

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

## Building the Project

### Rust Library Only

```bash
cargo build --release
cargo test
```

### Python Package

```bash
# Install maturin
pip install maturin

# Development build (debug)
maturin develop --features python

# Release build
maturin build --release --features python

# Install the wheel
pip install target/wheels/rsasa-*.whl
```

### Run Examples

Rust example:
```bash
cargo run --example simple
```

Python example:
```bash
python python/example.py
```

## Performance

The Rust implementation provides significant performance improvements:
- **2-10x faster** than pure Python/NumPy implementations
- **O(n²)** complexity with neighbor search optimization
- **Minimal memory allocation** with pre-allocated work buffers
- **Parallel-friendly** design for trajectory processing

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

### Development Setup

1. Install Rust: https://rustup.rs/
2. Install Python development tools: `pip install maturin pytest numpy`
3. Run tests: `cargo test && pytest`

## References

- Shrake, A.; Rupley, J. A. (1973). "Environment and exposure to solvent of protein atoms. Lysozyme and insulin". Journal of Molecular Biology. 79 (2): 351–371.