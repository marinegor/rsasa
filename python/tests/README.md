# rsasa Python Tests

This directory contains the test suite for the rsasa Python package.

## Test Structure

- `test_basic.py` - Core functionality tests (no external dependencies except numpy)
- `test_comparison_mdtraj.py` - Comparison tests against MDTraj's SASA implementation

## Installation

### Basic Tests Only

For basic tests, you only need pytest:

```bash
pip install pytest
```

### All Tests (Including MDTraj Comparison)

To run all tests including MDTraj comparisons:

```bash
pip install pytest mdtraj
# Or use the optional dependencies:
pip install -e ".[test]"
```

## Running Tests

### Run All Tests

```bash
# From the rsasa root directory
pytest

# Or with more verbose output
pytest -v

# Or with print statements visible
pytest -v -s
```

### Run Basic Tests Only (Skip MDTraj)

If you don't have MDTraj installed:

```bash
pytest python/tests/test_basic.py
```

### Run MDTraj Comparison Tests Only

```bash
pytest python/tests/test_comparison_mdtraj.py
```

### Run Specific Test Classes or Functions

```bash
# Run a specific test class
pytest python/tests/test_basic.py::TestBasicFunctionality

# Run a specific test function
pytest python/tests/test_basic.py::TestBasicFunctionality::test_single_atom
```

## Test Coverage

### test_basic.py

Tests core functionality without external dependencies:

- **TestBasicFunctionality**: Basic SASA calculations
  - Single atom SASA
  - Two atoms (overlapping and non-overlapping)
  - Multiple frames
  - Atom grouping
  - Atom selection masks

- **TestInputValidation**: Error handling
  - Shape mismatches
  - Invalid parameters
  - Invalid atom mappings

- **TestDataTypes**: Type handling
  - float64 to float32 conversion
  - int64 to uintp conversion

- **TestEdgeCases**: Corner cases
  - All atoms masked out
  - Different sphere point counts
  - Large number of atoms

### test_comparison_mdtraj.py

Compares rsasa results with MDTraj (requires MDTraj installation):

- Single atom comparison
- Simple trajectory comparison
- Multi-frame trajectory comparison
- Different sphere point counts
- Different probe radii
- Residue-level aggregation
- Performance comparison

## Expected Results

### Basic Tests

All basic tests should pass with rsasa correctly installed:

```bash
$ pytest python/tests/test_basic.py -v
============================== test session starts ===============================
collected 18 items

python/tests/test_basic.py::TestBasicFunctionality::test_single_atom PASSED
python/tests/test_basic.py::TestBasicFunctionality::test_two_atoms_non_overlapping PASSED
...
============================== 18 passed in 0.50s ================================
```

### MDTraj Comparison Tests

When MDTraj is installed, comparison tests should show close agreement:

```bash
$ pytest python/tests/test_comparison_mdtraj.py -v
============================== test session starts ===============================
collected 8 items

python/tests/test_comparison_mdtraj.py::test_single_atom_comparison PASSED
python/tests/test_comparison_mdtraj.py::test_simple_trajectory_comparison PASSED
...
============================== 8 passed in 2.35s =================================
```

The performance test will print timing information showing the speedup:

```
Performance comparison:
  MDTraj: 12.45 ms
  rsasa:  3.21 ms
  Speedup: 3.88x
```

## Tolerance Settings

The comparison tests use the following tolerances:

- **Standard**: 5% relative tolerance (`rtol=0.05`)
- **Low sphere points** (<500): 10% relative tolerance (`rtol=0.10`)

These tolerances account for:
- Numerical differences in sphere point generation algorithms
- Floating-point precision differences
- Minor implementation differences

## Troubleshooting

### Tests are skipped

If you see:
```
python/tests/test_basic.py SKIPPED (rsasa not installed)
```

**Solution**: Build and install rsasa:
```bash
maturin develop --features python
```

### MDTraj tests are skipped

If you see:
```
python/tests/test_comparison_mdtraj.py SKIPPED (MDTraj not installed)
```

**Solution**: Install MDTraj:
```bash
pip install mdtraj
# or
conda install -c conda-forge mdtraj
```

### Tests fail with shape errors

Make sure rsasa is built with the latest code:
```bash
maturin develop --features python --release
```

### Large differences from MDTraj

Small differences (within 5%) are expected due to:
- Different sphere point generation algorithms
- Different neighbor search optimizations
- Floating-point arithmetic order

If differences are larger:
1. Check that you're using the same probe radius
2. Check that you're using the same number of sphere points
3. Verify coordinate units (rsasa: Å, MDTraj: nm)

## Contributing Tests

When adding new tests:

1. Add to `test_basic.py` if testing core functionality
2. Add to `test_comparison_mdtraj.py` if comparing with MDTraj
3. Use descriptive test names
4. Include docstrings explaining what's being tested
5. Use appropriate assertions with clear error messages
6. Mark slow tests with `@pytest.mark.slow`

Example:

```python
def test_my_feature(self):
    """Test description of what this tests."""
    # Setup
    xyzlist = np.array([...], dtype=np.float32)
    
    # Execute
    result = sasa(...)
    
    # Assert
    assert result.shape == expected_shape, "Result has wrong shape"
    assert_allclose(result, expected, rtol=0.05)
```

## Running Tests During Development

For rapid iteration during development:

```bash
# Build in debug mode for faster compilation
maturin develop --features python

# Run tests
pytest python/tests/test_basic.py -v

# Run a specific test while debugging
pytest python/tests/test_basic.py::TestBasicFunctionality::test_single_atom -v -s
```

## Continuous Integration

These tests are designed to run in CI environments:

```yaml
# Example GitHub Actions workflow
- name: Install dependencies
  run: |
    pip install maturin pytest numpy
    pip install mdtraj  # Optional

- name: Build rsasa
  run: maturin develop --features python

- name: Run tests
  run: pytest -v
```

## Performance Benchmarking

To get detailed performance information:

```bash
pytest python/tests/test_comparison_mdtraj.py::test_performance_comparison -v -s
```

This will show:
- Execution time for MDTraj
- Execution time for rsasa
- Speedup factor

## Additional Resources

- [pytest documentation](https://docs.pytest.org/)
- [MDTraj SASA documentation](https://mdtraj.org/latest/api/generated/mdtraj.shrake_rupley.html)
- Main project README: `../README.md`
- Python package README: `../python/README.md`
