# rsasa - Python Package

Fast SASA (Solvent Accessible Surface Area) calculator for molecular structures.

This is the Python binding for the rsasa library, which implements the Shrake-Rupley algorithm in Rust for high performance SASA calculations.

## Installation

### From PyPI (when published)

```bash
pip install rsasa
```

### From source

Requirements:
- Python >= 3.8
- Rust toolchain (install from https://rustup.rs/)
- maturin (`pip install maturin`)

```bash
# Clone the repository
git clone <repository-url>
cd rsasa

# Build and install the Python package
maturin develop --features python

# Or build a wheel
maturin build --release --features python
pip install target/wheels/rsasa-*.whl
```

## Quick Start

```python
import numpy as np
from rsasa import sasa

# Define a simple system with two atoms
n_frames = 1
n_atoms = 2
xyzlist = np.array([0.0, 0.0, 0.0, 1.8, 0.0, 0.0], dtype=np.float32)
atom_radii = np.array([1.5, 1.5], dtype=np.float32)  # VDW radius + probe radius
n_sphere_points = 960  # Number of sphere points for accuracy
atom_mapping = np.array([0, 1], dtype=np.uintp)  # Each atom in its own group
atom_selection_mask = np.array([1, 1], dtype=np.int32)  # Compute SASA for both atoms
n_groups = 2

# Calculate SASA
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

print(f"Atom 1 SASA: {result[0]:.2f} Ų")
print(f"Atom 2 SASA: {result[1]:.2f} Ų")
```

## API Reference

### `sasa` Function

```python
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
```

Calculate the solvent accessible surface area (SASA) for atoms in a trajectory.

#### Parameters

- **n_frames** (`int`): Number of frames in the trajectory
- **n_atoms** (`int`): Number of atoms in each frame
- **xyzlist** (`ndarray[float32]`): Flattened array of atomic coordinates in Ångströms
  - Shape: `(n_frames * n_atoms * 3,)`
  - Format: `[x1, y1, z1, x2, y2, z2, ...]` for all atoms in frame 1, then frame 2, etc.
- **atom_radii** (`ndarray[float32]`): Van der Waals radii plus probe radius in Ångströms
  - Shape: `(n_atoms,)`
  - Typically VDW radius + 1.4 Å for water probe
- **n_sphere_points** (`int`): Number of points on the unit sphere for sampling
  - Common values: 92 (fast), 960 (balanced), 9600 (high accuracy)
- **atom_mapping** (`ndarray[uintp]`): Mapping from atoms to groups
  - Shape: `(n_atoms,)`
  - Each value must be in range `[0, n_groups)`
- **atom_selection_mask** (`ndarray[int32]`): Binary mask for atom selection
  - Shape: `(n_atoms,)`
  - 1 = compute SASA, 0 = skip
- **n_groups** (`int`): Number of groups in the output

#### Returns

- **result** (`ndarray[float32]`): SASA values in Ų (square Ångströms)
  - Shape: `(n_frames * n_groups,)`
  - Format: `[group0_frame0, group1_frame0, ..., group0_frame1, ...]`

#### Raises

- **ValueError**: If input arrays have incorrect shapes or invalid parameter values

## Examples

### Example 1: Single Atom

Calculate SASA for a single isolated atom (should match theoretical sphere):

```python
import numpy as np
from rsasa import sasa

xyzlist = np.array([0.0, 0.0, 0.0], dtype=np.float32)
atom_radii = np.array([1.5], dtype=np.float32)
atom_mapping = np.array([0], dtype=np.uintp)
atom_selection_mask = np.array([1], dtype=np.int32)

result = sasa(1, 1, xyzlist, atom_radii, 960, 
              atom_mapping, atom_selection_mask, 1)

expected = 4 * np.pi * 1.5**2
print(f"Calculated: {result[0]:.2f} Ų")
print(f"Expected: {expected:.2f} Ų")
```

### Example 2: Residue-Level SASA

Group multiple atoms into residues:

```python
import numpy as np
from rsasa import sasa

# 6 atoms belonging to 2 residues
n_atoms = 6
xyzlist = np.array([
    0.0, 0.0, 0.0,  # Residue 1, atom 1
    1.0, 0.0, 0.0,  # Residue 1, atom 2
    2.0, 0.0, 0.0,  # Residue 1, atom 3
    5.0, 0.0, 0.0,  # Residue 2, atom 1
    6.0, 0.0, 0.0,  # Residue 2, atom 2
    7.0, 0.0, 0.0,  # Residue 2, atom 3
], dtype=np.float32)

atom_radii = np.array([1.2, 1.4, 1.3, 1.2, 1.4, 1.3], dtype=np.float32)
atom_mapping = np.array([0, 0, 0, 1, 1, 1], dtype=np.uintp)  # 3 atoms per residue
atom_selection_mask = np.ones(n_atoms, dtype=np.int32)

result = sasa(1, n_atoms, xyzlist, atom_radii, 960,
              atom_mapping, atom_selection_mask, 2)

print(f"Residue 1 SASA: {result[0]:.2f} Ų")
print(f"Residue 2 SASA: {result[1]:.2f} Ų")
```

### Example 3: Trajectory Analysis

Analyze SASA changes over a trajectory:

```python
import numpy as np
from rsasa import sasa

n_frames = 100
n_atoms = 10

# Generate random trajectory
np.random.seed(42)
xyzlist = np.random.randn(n_frames * n_atoms * 3).astype(np.float32) * 5.0
atom_radii = np.ones(n_atoms, dtype=np.float32) * 1.5
atom_mapping = np.arange(n_atoms, dtype=np.uintp)
atom_selection_mask = np.ones(n_atoms, dtype=np.int32)

result = sasa(n_frames, n_atoms, xyzlist, atom_radii, 960,
              atom_mapping, atom_selection_mask, n_atoms)

# Reshape to (n_frames, n_atoms)
result_2d = result.reshape(n_frames, n_atoms)

# Analyze
mean_sasa = result_2d.mean(axis=0)
std_sasa = result_2d.std(axis=0)

print("Per-atom SASA statistics:")
for i in range(n_atoms):
    print(f"  Atom {i}: {mean_sasa[i]:.2f} ± {std_sasa[i]:.2f} Ų")
```

### Example 4: Working with 2D Coordinate Arrays

Convert between 2D coordinate arrays and flat arrays:

```python
import numpy as np
from rsasa import sasa

# Start with structured coordinates
coords_2d = np.array([
    [0.0, 0.0, 0.0],
    [2.0, 0.0, 0.0],
    [1.0, 2.0, 0.0],
], dtype=np.float32)

# Flatten for rsasa
xyzlist = coords_2d.flatten()

atom_radii = np.array([1.0, 1.0, 1.0], dtype=np.float32)
atom_mapping = np.arange(3, dtype=np.uintp)
atom_selection_mask = np.ones(3, dtype=np.int32)

result = sasa(1, 3, xyzlist, atom_radii, 960,
              atom_mapping, atom_selection_mask, 3)

print(f"SASA values: {result}")
```

## Performance Tips

1. **Sphere Points**: Use fewer sphere points (92-240) for quick calculations, more (960-9600) for publication-quality results
2. **Batch Processing**: Process multiple frames together for better performance
3. **Memory**: Pre-allocate arrays with the correct dtype to avoid conversions
4. **Atom Selection**: Use `atom_selection_mask` to skip atoms you don't need

## Comparison with Other Tools

- **Faster**: Rust implementation is typically 2-10x faster than pure Python/NumPy implementations
- **Accurate**: Uses the standard Shrake-Rupley algorithm with customizable sphere point density
- **Memory Efficient**: Minimal memory overhead, streaming through frames
- **Type Safe**: Strong typing with numpy arrays prevents common errors

## Troubleshooting

### Import Error

If you get an import error, make sure the package is installed correctly:

```bash
pip install rsasa
# Or if building from source:
maturin develop --features python
```

### Shape Errors

Ensure your arrays have the correct shapes:
- `xyzlist`: must be flat with length `n_frames * n_atoms * 3`
- `atom_radii`: length `n_atoms`
- `atom_mapping`: length `n_atoms`, values in `[0, n_groups)`
- `atom_selection_mask`: length `n_atoms`, values 0 or 1

### Type Errors

Use the correct numpy dtypes:
- `xyzlist`: `np.float32`
- `atom_radii`: `np.float32`
- `atom_mapping`: `np.uintp` (platform-specific unsigned integer)
- `atom_selection_mask`: `np.int32`

## Contributing

Contributions are welcome! Please see the main repository for guidelines.

## License

Licensed under either of:
- Apache License, Version 2.0
- MIT license

at your option.

## Citation

If you use this software in your research, please cite:

- Shrake, A.; Rupley, J. A. (1973). "Environment and exposure to solvent of protein atoms. Lysozyme and insulin". Journal of Molecular Biology. 79 (2): 351–371.