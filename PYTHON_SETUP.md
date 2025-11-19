# Python Setup Guide for rsasa

This guide explains how to build, install, and use the rsasa Python extension.

## Overview

The rsasa library is written in Rust with Python bindings via PyO3. The Python package provides a native numpy interface to the high-performance Rust SASA calculator.

## Prerequisites

### Required

- **Python**: 3.8 or later
- **Rust**: Latest stable version (install from https://rustup.rs/)
- **pip**: Python package installer (usually comes with Python)

### Build Tools

- **maturin**: PyO3 build tool
  ```bash
  pip install maturin
  ```

## Installation Methods

### Method 1: Development Install (Recommended for Development)

This method builds the extension in debug mode and installs it in editable mode:

```bash
cd rsasa
maturin develop --features python
```

This is the fastest way to iterate during development. Changes to Python code are immediately reflected, but Rust changes require rebuilding.

### Method 2: Release Build

For production use or performance testing:

```bash
cd rsasa
maturin build --release --features python
```

This creates a wheel in `target/wheels/`. Install it with:

```bash
pip install target/wheels/rsasa-*.whl
```

### Method 3: Direct Install (All-in-One)

Build and install in one command:

```bash
cd rsasa
pip install . --features python
```

## Verifying Installation

Test that the package is correctly installed:

```python
import rsasa
print(rsasa.__version__)

# Run a simple calculation
import numpy as np
xyzlist = np.array([0.0, 0.0, 0.0], dtype=np.float32)
atom_radii = np.array([1.5], dtype=np.float32)
atom_mapping = np.array([0], dtype=np.uintp)
atom_selection_mask = np.array([1], dtype=np.int32)

result = rsasa.sasa(1, 1, xyzlist, atom_radii, 960, 
                    atom_mapping, atom_selection_mask, 1)
print(f"SASA: {result[0]:.2f} Ų")
```

Expected output:
```
0.1.0
SASA: 28.27 Ų
```

## Running Examples

### Python Example Script

```bash
cd rsasa
python python/example.py
```

This runs through several examples demonstrating different features.

### Interactive Python

```python
import numpy as np
from rsasa import sasa

# Simple two-atom example
xyzlist = np.array([0.0, 0.0, 0.0, 1.8, 0.0, 0.0], dtype=np.float32)
atom_radii = np.array([1.5, 1.5], dtype=np.float32)
atom_mapping = np.array([0, 1], dtype=np.uintp)
atom_selection_mask = np.array([1, 1], dtype=np.int32)

result = sasa(1, 2, xyzlist, atom_radii, 960,
              atom_mapping, atom_selection_mask, 2)

print(f"Atom 1: {result[0]:.2f} Ų")
print(f"Atom 2: {result[1]:.2f} Ų")
```

## Development Workflow

### Making Changes to Rust Code

1. Edit the Rust source in `src/lib.rs`
2. Rebuild the extension:
   ```bash
   maturin develop --features python
   ```
3. Test in Python (no need to restart Python interpreter for new builds)

### Making Changes to Python Code

Python files in `python/rsasa/` are used directly when installed with `maturin develop`, so changes are immediately visible.

## Troubleshooting

### Import Error: "No module named 'rsasa'"

**Solution**: Make sure you've run `maturin develop --features python` or installed the wheel.

### Import Error: "cannot import name 'sasa'"

**Solution**: Check that the module was built with the `--features python` flag:
```bash
maturin develop --features python
```

### Build Error: "pyo3 requires Python development headers"

**On Ubuntu/Debian**:
```bash
sudo apt-get install python3-dev
```

**On macOS**:
```bash
brew install python@3.11  # or your Python version
```

**On Windows**:
Download and install Python from python.org (includes dev headers).

### Build Error: "linking with `cc` failed"

This usually means Python libraries aren't being found. Make sure:
1. You're using a properly installed Python (not a system Python on macOS)
2. Python development headers are installed
3. Try setting `PYO3_PYTHON` environment variable:
   ```bash
   export PYO3_PYTHON=$(which python3)
   maturin develop --features python
   ```

### Runtime Error: "incompatible architecture"

If you get architecture errors (e.g., building on ARM but running on x86), rebuild:
```bash
maturin build --release --features python
```

### Type Errors: "array has wrong dtype"

Ensure you're using the correct numpy dtypes:
```python
xyzlist = np.array([...], dtype=np.float32)      # Not float64!
atom_radii = np.array([...], dtype=np.float32)   # Not float64!
atom_mapping = np.array([...], dtype=np.uintp)   # Platform-specific
atom_selection_mask = np.array([...], dtype=np.int32)  # Not int64!
```

### Shape Errors

Double-check array dimensions:
- `xyzlist`: length must be exactly `n_frames * n_atoms * 3`
- `atom_radii`: length must be exactly `n_atoms`
- `atom_mapping`: length must be exactly `n_atoms`
- `atom_selection_mask`: length must be exactly `n_atoms`

## Testing

### Run Rust Tests

```bash
cargo test --features python
```

### Run Python Tests (if pytest is set up)

```bash
pytest python/tests/
```

### Manual Testing

```bash
python python/example.py
```

## Performance Optimization

### Build with Optimizations

Always use `--release` for production:
```bash
maturin build --release --features python
```

Release builds are typically 10-20x faster than debug builds.

### Profile Python Code

```python
import time
import numpy as np
from rsasa import sasa

# Setup
n_frames = 100
n_atoms = 1000
xyzlist = np.random.randn(n_frames * n_atoms * 3).astype(np.float32)
atom_radii = np.ones(n_atoms, dtype=np.float32) * 1.5
atom_mapping = np.arange(n_atoms, dtype=np.uintp)
atom_selection_mask = np.ones(n_atoms, dtype=np.int32)

# Benchmark
start = time.time()
result = sasa(n_frames, n_atoms, xyzlist, atom_radii, 960,
              atom_mapping, atom_selection_mask, n_atoms)
elapsed = time.time() - start

print(f"Processed {n_frames} frames with {n_atoms} atoms in {elapsed:.2f}s")
print(f"Rate: {n_frames / elapsed:.1f} frames/sec")
```

## Distribution

### Build Wheels for Multiple Platforms

Using GitHub Actions or CI/CD:

```yaml
- uses: PyO3/maturin-action@v1
  with:
    command: build
    args: --release --features python
```

### Publish to PyPI (when ready)

```bash
maturin publish --features python
```

## API Quick Reference

```python
import numpy as np
from rsasa import sasa

result = sasa(
    n_frames=1,                                    # int
    n_atoms=10,                                    # int
    xyzlist=np.zeros(30, dtype=np.float32),       # (n_frames*n_atoms*3,)
    atom_radii=np.ones(10, dtype=np.float32),     # (n_atoms,)
    n_sphere_points=960,                           # int
    atom_mapping=np.arange(10, dtype=np.uintp),   # (n_atoms,)
    atom_selection_mask=np.ones(10, dtype=np.int32),  # (n_atoms,)
    n_groups=10,                                   # int
) -> np.ndarray  # shape: (n_frames*n_groups,), dtype: float32
```

## Additional Resources

- **Main README**: `README.md` - Overview of the entire project
- **Python README**: `python/README.md` - Detailed Python API documentation
- **Rust Documentation**: Run `cargo doc --open --features python`
- **Examples**: See `python/example.py` for comprehensive examples
- **PyO3 Guide**: https://pyo3.rs/ - Official PyO3 documentation

## Getting Help

If you encounter issues:

1. Check this guide first
2. Look at the examples in `python/example.py`
3. Check the main README.md
4. Review error messages carefully - they often indicate the exact problem
5. Open an issue on GitHub with:
   - Your OS and Python version
   - The exact command you ran
   - The complete error message
   - A minimal reproducible example

## Summary

**Quick start commands**:
```bash
# 1. Install maturin
pip install maturin

# 2. Build and install
cd rsasa
maturin develop --features python

# 3. Test
python python/example.py
```

That's it! You're ready to use rsasa in Python.