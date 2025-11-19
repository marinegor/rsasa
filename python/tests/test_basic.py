"""
Basic unit tests for rsasa Python interface.

These tests don't require MDTraj and test the core functionality.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose, assert_array_equal

try:
    from rsasa import sasa
    HAS_RSASA = True
except ImportError:
    HAS_RSASA = False


pytestmark = pytest.mark.skipif(not HAS_RSASA, reason="rsasa not installed")


class TestBasicFunctionality:
    """Test basic SASA calculations."""
    
    def test_single_atom(self):
        """Test SASA for a single isolated atom."""
        n_frames = 1
        n_atoms = 1
        xyzlist = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        atom_radii = np.array([1.5], dtype=np.float32)
        n_sphere_points = 960
        atom_mapping = np.array([0], dtype=np.uintp)
        atom_selection_mask = np.array([1], dtype=np.int32)
        n_groups = 1
        
        result = sasa(
            n_frames, n_atoms, xyzlist, atom_radii,
            n_sphere_points, atom_mapping, atom_selection_mask, n_groups
        )
        
        # Expected: surface area of sphere with radius 1.5 Angstroms
        expected = 4.0 * np.pi * 1.5**2
        
        assert result.shape == (1,)
        assert result.dtype == np.float32
        assert_allclose(result[0], expected, rtol=0.05)
    
    def test_two_atoms_non_overlapping(self):
        """Test SASA for two atoms that don't overlap."""
        n_frames = 1
        n_atoms = 2
        # Atoms 5 Angstroms apart (radii 1.0 each, so no overlap)
        xyzlist = np.array([0.0, 0.0, 0.0, 5.0, 0.0, 0.0], dtype=np.float32)
        atom_radii = np.array([1.0, 1.0], dtype=np.float32)
        n_sphere_points = 960
        atom_mapping = np.array([0, 1], dtype=np.uintp)
        atom_selection_mask = np.array([1, 1], dtype=np.int32)
        n_groups = 2
        
        result = sasa(
            n_frames, n_atoms, xyzlist, atom_radii,
            n_sphere_points, atom_mapping, atom_selection_mask, n_groups
        )
        
        # Each atom should have nearly full sphere area
        expected = 4.0 * np.pi * 1.0**2
        
        assert result.shape == (2,)
        assert_allclose(result[0], expected, rtol=0.05)
        assert_allclose(result[1], expected, rtol=0.05)
    
    def test_two_atoms_overlapping(self):
        """Test SASA for two overlapping atoms."""
        n_frames = 1
        n_atoms = 2
        # Atoms 1.5 Angstroms apart with radii 1.0 each (overlapping)
        xyzlist = np.array([0.0, 0.0, 0.0, 1.5, 0.0, 0.0], dtype=np.float32)
        atom_radii = np.array([1.0, 1.0], dtype=np.float32)
        n_sphere_points = 960
        atom_mapping = np.array([0, 1], dtype=np.uintp)
        atom_selection_mask = np.array([1, 1], dtype=np.int32)
        n_groups = 2
        
        result = sasa(
            n_frames, n_atoms, xyzlist, atom_radii,
            n_sphere_points, atom_mapping, atom_selection_mask, n_groups
        )
        
        # Each atom should have less than full sphere area due to overlap
        expected_full = 4.0 * np.pi * 1.0**2
        
        assert result.shape == (2,)
        assert result[0] < expected_full
        assert result[1] < expected_full
        assert result[0] > 0
        assert result[1] > 0
    
    def test_multiple_frames(self):
        """Test SASA calculation over multiple frames."""
        n_frames = 3
        n_atoms = 2
        # Three frames with atoms getting closer
        xyzlist = np.array([
            # Frame 1: far apart
            0.0, 0.0, 0.0, 5.0, 0.0, 0.0,
            # Frame 2: closer
            0.0, 0.0, 0.0, 2.0, 0.0, 0.0,
            # Frame 3: very close
            0.0, 0.0, 0.0, 1.5, 0.0, 0.0,
        ], dtype=np.float32)
        atom_radii = np.array([1.0, 1.0], dtype=np.float32)
        n_sphere_points = 960
        atom_mapping = np.array([0, 1], dtype=np.uintp)
        atom_selection_mask = np.array([1, 1], dtype=np.int32)
        n_groups = 2
        
        result = sasa(
            n_frames, n_atoms, xyzlist, atom_radii,
            n_sphere_points, atom_mapping, atom_selection_mask, n_groups
        )
        
        # Result should have shape (n_frames * n_groups,)
        assert result.shape == (6,)
        
        # Reshape to (n_frames, n_groups)
        result_2d = result.reshape(n_frames, n_groups)
        
        # Total SASA should decrease as atoms get closer
        total_frame1 = result_2d[0].sum()
        total_frame2 = result_2d[1].sum()
        total_frame3 = result_2d[2].sum()
        
        assert total_frame1 > total_frame2 > total_frame3
    
    def test_atom_grouping(self):
        """Test grouping multiple atoms together."""
        n_frames = 1
        n_atoms = 4
        xyzlist = np.array([
            0.0, 0.0, 0.0,
            1.5, 0.0, 0.0,
            5.0, 0.0, 0.0,
            6.5, 0.0, 0.0,
        ], dtype=np.float32)
        atom_radii = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float32)
        n_sphere_points = 960
        # Group atoms 0,1 into group 0 and atoms 2,3 into group 1
        atom_mapping = np.array([0, 0, 1, 1], dtype=np.uintp)
        atom_selection_mask = np.array([1, 1, 1, 1], dtype=np.int32)
        n_groups = 2
        
        result = sasa(
            n_frames, n_atoms, xyzlist, atom_radii,
            n_sphere_points, atom_mapping, atom_selection_mask, n_groups
        )
        
        assert result.shape == (2,)
        # Both groups should have positive SASA
        assert result[0] > 0
        assert result[1] > 0
    
    def test_atom_selection_mask(self):
        """Test selective computation with atom_selection_mask."""
        n_frames = 1
        n_atoms = 3
        xyzlist = np.array([
            0.0, 0.0, 0.0,
            2.0, 0.0, 0.0,
            4.0, 0.0, 0.0,
        ], dtype=np.float32)
        atom_radii = np.array([1.0, 1.0, 1.0], dtype=np.float32)
        n_sphere_points = 960
        atom_mapping = np.array([0, 1, 2], dtype=np.uintp)
        # Only compute for atoms 0 and 2 (skip atom 1)
        atom_selection_mask = np.array([1, 0, 1], dtype=np.int32)
        n_groups = 3
        
        result = sasa(
            n_frames, n_atoms, xyzlist, atom_radii,
            n_sphere_points, atom_mapping, atom_selection_mask, n_groups
        )
        
        assert result.shape == (3,)
        # Atoms 0 and 2 should have positive SASA
        assert result[0] > 0
        assert result[2] > 0
        # Atom 1 should have zero SASA (skipped)
        assert result[1] == 0.0


class TestInputValidation:
    """Test error handling and input validation."""
    
    def test_shape_mismatch_xyzlist(self):
        """Test error when xyzlist has wrong shape."""
        n_frames = 1
        n_atoms = 2
        # Wrong size: should be 6 but is 3
        xyzlist = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        atom_radii = np.array([1.0, 1.0], dtype=np.float32)
        n_sphere_points = 960
        atom_mapping = np.array([0, 1], dtype=np.uintp)
        atom_selection_mask = np.array([1, 1], dtype=np.int32)
        n_groups = 2
        
        with pytest.raises(ValueError, match="xyzlist"):
            sasa(n_frames, n_atoms, xyzlist, atom_radii,
                 n_sphere_points, atom_mapping, atom_selection_mask, n_groups)
    
    def test_shape_mismatch_atom_radii(self):
        """Test error when atom_radii has wrong shape."""
        n_frames = 1
        n_atoms = 2
        xyzlist = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0], dtype=np.float32)
        # Wrong size: should be 2 but is 1
        atom_radii = np.array([1.0], dtype=np.float32)
        n_sphere_points = 960
        atom_mapping = np.array([0, 1], dtype=np.uintp)
        atom_selection_mask = np.array([1, 1], dtype=np.int32)
        n_groups = 2
        
        with pytest.raises(ValueError, match="atom_radii"):
            sasa(n_frames, n_atoms, xyzlist, atom_radii,
                 n_sphere_points, atom_mapping, atom_selection_mask, n_groups)
    
    def test_shape_mismatch_atom_mapping(self):
        """Test error when atom_mapping has wrong shape."""
        n_frames = 1
        n_atoms = 2
        xyzlist = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0], dtype=np.float32)
        atom_radii = np.array([1.0, 1.0], dtype=np.float32)
        n_sphere_points = 960
        # Wrong size
        atom_mapping = np.array([0], dtype=np.uintp)
        atom_selection_mask = np.array([1, 1], dtype=np.int32)
        n_groups = 2
        
        with pytest.raises(ValueError, match="atom_mapping"):
            sasa(n_frames, n_atoms, xyzlist, atom_radii,
                 n_sphere_points, atom_mapping, atom_selection_mask, n_groups)
    
    def test_invalid_atom_mapping_value(self):
        """Test error when atom_mapping has invalid group index."""
        n_frames = 1
        n_atoms = 2
        xyzlist = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0], dtype=np.float32)
        atom_radii = np.array([1.0, 1.0], dtype=np.float32)
        n_sphere_points = 960
        # Invalid: group 5 doesn't exist (n_groups=2)
        atom_mapping = np.array([0, 5], dtype=np.uintp)
        atom_selection_mask = np.array([1, 1], dtype=np.int32)
        n_groups = 2
        
        with pytest.raises(ValueError, match="atom_mapping"):
            sasa(n_frames, n_atoms, xyzlist, atom_radii,
                 n_sphere_points, atom_mapping, atom_selection_mask, n_groups)
    
    def test_invalid_n_sphere_points(self):
        """Test error when n_sphere_points is invalid."""
        n_frames = 1
        n_atoms = 1
        xyzlist = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        atom_radii = np.array([1.0], dtype=np.float32)
        n_sphere_points = 0  # Invalid
        atom_mapping = np.array([0], dtype=np.uintp)
        atom_selection_mask = np.array([1], dtype=np.int32)
        n_groups = 1
        
        with pytest.raises(ValueError, match="n_sphere_points"):
            sasa(n_frames, n_atoms, xyzlist, atom_radii,
                 n_sphere_points, atom_mapping, atom_selection_mask, n_groups)


class TestDataTypes:
    """Test handling of different data types."""
    
    def test_float64_conversion(self):
        """Test that float64 arrays are converted to float32."""
        n_frames = 1
        n_atoms = 1
        # Input as float64
        xyzlist = np.array([0.0, 0.0, 0.0], dtype=np.float64)
        atom_radii = np.array([1.5], dtype=np.float64)
        n_sphere_points = 960
        atom_mapping = np.array([0], dtype=np.uintp)
        atom_selection_mask = np.array([1], dtype=np.int32)
        n_groups = 1
        
        # Should work (converted internally)
        result = sasa(
            n_frames, n_atoms, xyzlist, atom_radii,
            n_sphere_points, atom_mapping, atom_selection_mask, n_groups
        )
        
        # Result should be float32
        assert result.dtype == np.float32
        expected = 4.0 * np.pi * 1.5**2
        assert_allclose(result[0], expected, rtol=0.05)
    
    def test_int64_mapping_conversion(self):
        """Test that int64 atom_mapping is converted."""
        n_frames = 1
        n_atoms = 1
        xyzlist = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        atom_radii = np.array([1.5], dtype=np.float32)
        n_sphere_points = 960
        # Input as int64
        atom_mapping = np.array([0], dtype=np.int64)
        atom_selection_mask = np.array([1], dtype=np.int32)
        n_groups = 1
        
        # Should work (converted internally)
        result = sasa(
            n_frames, n_atoms, xyzlist, atom_radii,
            n_sphere_points, atom_mapping, atom_selection_mask, n_groups
        )
        
        assert result.shape == (1,)
        assert result[0] > 0


class TestEdgeCases:
    """Test edge cases and corner scenarios."""
    
    def test_all_atoms_masked_out(self):
        """Test when all atoms are masked out."""
        n_frames = 1
        n_atoms = 2
        xyzlist = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0], dtype=np.float32)
        atom_radii = np.array([1.0, 1.0], dtype=np.float32)
        n_sphere_points = 960
        atom_mapping = np.array([0, 1], dtype=np.uintp)
        # All atoms masked out
        atom_selection_mask = np.array([0, 0], dtype=np.int32)
        n_groups = 2
        
        result = sasa(
            n_frames, n_atoms, xyzlist, atom_radii,
            n_sphere_points, atom_mapping, atom_selection_mask, n_groups
        )
        
        # All SASA values should be zero
        assert_array_equal(result, np.zeros(2, dtype=np.float32))
    
    def test_different_sphere_point_counts(self):
        """Test with different sphere point counts."""
        n_frames = 1
        n_atoms = 1
        xyzlist = np.array([0.0, 0.0, 0.0], dtype=np.float32)
        atom_radii = np.array([1.5], dtype=np.float32)
        atom_mapping = np.array([0], dtype=np.uintp)
        atom_selection_mask = np.array([1], dtype=np.int32)
        n_groups = 1
        
        expected = 4.0 * np.pi * 1.5**2
        
        # Test different sphere point counts
        for n_points in [92, 240, 960, 9600]:
            result = sasa(
                n_frames, n_atoms, xyzlist, atom_radii,
                n_points, atom_mapping, atom_selection_mask, n_groups
            )
            
            # More points should give more accurate results
            tolerance = 0.15 if n_points < 500 else 0.05
            assert_allclose(result[0], expected, rtol=tolerance)
    
    def test_large_number_of_atoms(self):
        """Test with a larger number of atoms."""
        n_frames = 1
        n_atoms = 50
        
        # Create a grid of atoms
        np.random.seed(42)
        xyzlist = np.random.randn(n_atoms * 3).astype(np.float32) * 10.0
        atom_radii = np.ones(n_atoms, dtype=np.float32) * 1.5
        n_sphere_points = 960
        atom_mapping = np.arange(n_atoms, dtype=np.uintp)
        atom_selection_mask = np.ones(n_atoms, dtype=np.int32)
        n_groups = n_atoms
        
        result = sasa(
            n_frames, n_atoms, xyzlist, atom_radii,
            n_sphere_points, atom_mapping, atom_selection_mask, n_groups
        )
        
        assert result.shape == (n_atoms,)
        # All atoms should have some SASA
        assert np.all(result > 0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])