"""
Tests comparing rsasa with MDTraj's SASA implementation.

These tests require MDTraj to be installed. If MDTraj is not available,
the tests will be skipped.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

# Try to import rsasa
try:
    from rsasa import sasa
except ImportError:
    pytest.skip("rsasa not installed", allow_module_level=True)

# Try to import MDTraj
try:
    import mdtraj as md
    HAS_MDTRAJ = True
except ImportError:
    HAS_MDTRAJ = False


pytestmark = pytest.mark.skipif(not HAS_MDTRAJ, reason="MDTraj not installed")


@pytest.fixture
def simple_trajectory():
    """Create a simple test trajectory with a few atoms."""
    # Create a simple topology with 3 atoms
    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue("ALA", chain)
    
    # Add three atoms (simplified - just using carbon atoms)
    atom1 = topology.add_atom("CA", md.element.carbon, residue)
    atom2 = topology.add_atom("CB", md.element.carbon, residue)
    atom3 = topology.add_atom("CG", md.element.carbon, residue)
    
    # Create coordinates (in nanometers for MDTraj)
    # Frame 1: atoms spread out
    coords = np.array([
        [[0.0, 0.0, 0.0],
         [0.20, 0.0, 0.0],
         [0.10, 0.20, 0.0]]
    ], dtype=np.float32)
    
    # Create trajectory
    traj = md.Trajectory(coords, topology)
    return traj


@pytest.fixture
def multi_frame_trajectory():
    """Create a trajectory with multiple frames."""
    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue("ALA", chain)
    
    atom1 = topology.add_atom("CA", md.element.carbon, residue)
    atom2 = topology.add_atom("CB", md.element.carbon, residue)
    
    # Create 5 frames with varying distances
    coords = np.array([
        [[0.0, 0.0, 0.0], [0.30, 0.0, 0.0]],  # Far apart
        [[0.0, 0.0, 0.0], [0.25, 0.0, 0.0]],  # Closer
        [[0.0, 0.0, 0.0], [0.20, 0.0, 0.0]],  # Even closer
        [[0.0, 0.0, 0.0], [0.18, 0.0, 0.0]],  # Overlapping
        [[0.0, 0.0, 0.0], [0.15, 0.0, 0.0]],  # More overlap
    ], dtype=np.float32)
    
    traj = md.Trajectory(coords, topology)
    return traj


def get_vdw_radii_for_mdtraj(topology):
    """
    Get van der Waals radii for atoms in the topology.
    MDTraj uses radii from the element.
    """
    radii = []
    for atom in topology.atoms:
        # Get VDW radius in nanometers
        # Carbon: ~0.17 nm, use a standard value
        if atom.element.symbol == 'C':
            radii.append(0.17)  # nm
        elif atom.element.symbol == 'N':
            radii.append(0.155)  # nm
        elif atom.element.symbol == 'O':
            radii.append(0.152)  # nm
        elif atom.element.symbol == 'S':
            radii.append(0.180)  # nm
        else:
            radii.append(0.17)  # Default
    return np.array(radii, dtype=np.float32)


def mdtraj_to_rsasa_format(traj, probe_radius=0.14):
    """
    Convert MDTraj trajectory to rsasa input format.
    
    Parameters
    ----------
    traj : mdtraj.Trajectory
        Input trajectory
    probe_radius : float
        Probe radius in nanometers (default: 0.14 nm = 1.4 Å for water)
    
    Returns
    -------
    dict with rsasa parameters
    """
    n_frames = traj.n_frames
    n_atoms = traj.n_atoms
    
    # Convert coordinates from nm to Angstroms and flatten
    # MDTraj uses nm, rsasa uses Angstroms
    xyzlist = (traj.xyz * 10.0).flatten().astype(np.float32)
    
    # Get VDW radii and add probe radius, convert to Angstroms
    vdw_radii = get_vdw_radii_for_mdtraj(traj.topology)
    atom_radii = ((vdw_radii + probe_radius) * 10.0).astype(np.float32)
    
    # Create per-atom mapping and selection
    atom_mapping = np.arange(n_atoms, dtype=np.uintp)
    atom_selection_mask = np.ones(n_atoms, dtype=np.int32)
    
    return {
        'n_frames': n_frames,
        'n_atoms': n_atoms,
        'xyzlist': xyzlist,
        'atom_radii': atom_radii,
        'atom_mapping': atom_mapping,
        'atom_selection_mask': atom_selection_mask,
        'n_groups': n_atoms,
    }


def test_single_atom_comparison():
    """Test SASA calculation for a single isolated atom."""
    # Create a single atom trajectory
    topology = md.Topology()
    chain = topology.add_chain()
    residue = topology.add_residue("ALA", chain)
    atom = topology.add_atom("CA", md.element.carbon, residue)
    
    coords = np.array([[[0.0, 0.0, 0.0]]], dtype=np.float32)
    traj = md.Trajectory(coords, topology)
    
    # Calculate with MDTraj
    probe_radius = 0.14  # nm (1.4 Angstroms)
    n_sphere_points = 960
    mdtraj_sasa = md.shrake_rupley(traj, probe_radius=probe_radius, 
                                    n_sphere_points=n_sphere_points)
    
    # Convert to rsasa format and calculate
    rsasa_params = mdtraj_to_rsasa_format(traj, probe_radius=probe_radius)
    rsasa_result = sasa(
        rsasa_params['n_frames'],
        rsasa_params['n_atoms'],
        rsasa_params['xyzlist'],
        rsasa_params['atom_radii'],
        n_sphere_points,
        rsasa_params['atom_mapping'],
        rsasa_params['atom_selection_mask'],
        rsasa_params['n_groups'],
    )
    
    # MDTraj returns nm^2, rsasa returns Angstrom^2
    # 1 nm^2 = 100 Angstrom^2
    mdtraj_sasa_angstrom = mdtraj_sasa * 100.0
    
    # For a single isolated atom, both should give the full sphere surface area
    # Allow some tolerance due to sphere point sampling differences
    assert_allclose(rsasa_result, mdtraj_sasa_angstrom.flatten(), rtol=0.05)


def test_simple_trajectory_comparison(simple_trajectory):
    """Test SASA calculation for a simple trajectory."""
    probe_radius = 0.14
    n_sphere_points = 960
    
    # Calculate with MDTraj
    mdtraj_sasa = md.shrake_rupley(simple_trajectory, 
                                    probe_radius=probe_radius,
                                    n_sphere_points=n_sphere_points)
    
    # Calculate with rsasa
    rsasa_params = mdtraj_to_rsasa_format(simple_trajectory, 
                                          probe_radius=probe_radius)
    rsasa_result = sasa(
        rsasa_params['n_frames'],
        rsasa_params['n_atoms'],
        rsasa_params['xyzlist'],
        rsasa_params['atom_radii'],
        n_sphere_points,
        rsasa_params['atom_mapping'],
        rsasa_params['atom_selection_mask'],
        rsasa_params['n_groups'],
    )
    
    # Convert units
    mdtraj_sasa_angstrom = mdtraj_sasa * 100.0
    
    # Compare results (allow 5% tolerance)
    assert_allclose(rsasa_result, mdtraj_sasa_angstrom.flatten(), rtol=0.05)


def test_multi_frame_comparison(multi_frame_trajectory):
    """Test SASA calculation for multiple frames."""
    probe_radius = 0.14
    n_sphere_points = 960
    
    # Calculate with MDTraj
    mdtraj_sasa = md.shrake_rupley(multi_frame_trajectory,
                                    probe_radius=probe_radius,
                                    n_sphere_points=n_sphere_points)
    
    # Calculate with rsasa
    rsasa_params = mdtraj_to_rsasa_format(multi_frame_trajectory,
                                          probe_radius=probe_radius)
    rsasa_result = sasa(
        rsasa_params['n_frames'],
        rsasa_params['n_atoms'],
        rsasa_params['xyzlist'],
        rsasa_params['atom_radii'],
        n_sphere_points,
        rsasa_params['atom_mapping'],
        rsasa_params['atom_selection_mask'],
        rsasa_params['n_groups'],
    )
    
    # Convert units and reshape
    mdtraj_sasa_angstrom = mdtraj_sasa * 100.0
    n_frames = multi_frame_trajectory.n_frames
    n_atoms = multi_frame_trajectory.n_atoms
    rsasa_result_2d = rsasa_result.reshape(n_frames, n_atoms)
    
    # Compare results frame by frame
    assert_allclose(rsasa_result_2d, mdtraj_sasa_angstrom, rtol=0.05)
    
    # Check that SASA decreases as atoms get closer (qualitative check)
    total_sasa_per_frame = rsasa_result_2d.sum(axis=1)
    assert total_sasa_per_frame[0] > total_sasa_per_frame[-1], \
        "SASA should decrease as atoms get closer"


def test_different_sphere_points_comparison(simple_trajectory):
    """Test with different numbers of sphere points."""
    probe_radius = 0.14
    
    for n_sphere_points in [92, 240, 960]:
        # Calculate with MDTraj
        mdtraj_sasa = md.shrake_rupley(simple_trajectory,
                                        probe_radius=probe_radius,
                                        n_sphere_points=n_sphere_points)
        
        # Calculate with rsasa
        rsasa_params = mdtraj_to_rsasa_format(simple_trajectory,
                                              probe_radius=probe_radius)
        rsasa_result = sasa(
            rsasa_params['n_frames'],
            rsasa_params['n_atoms'],
            rsasa_params['xyzlist'],
            rsasa_params['atom_radii'],
            n_sphere_points,
            rsasa_params['atom_mapping'],
            rsasa_params['atom_selection_mask'],
            rsasa_params['n_groups'],
        )
        
        # Convert units
        mdtraj_sasa_angstrom = mdtraj_sasa * 100.0
        
        # Allow more tolerance for fewer sphere points
        tolerance = 0.10 if n_sphere_points < 500 else 0.05
        assert_allclose(rsasa_result, mdtraj_sasa_angstrom.flatten(), 
                       rtol=tolerance)


def test_different_probe_radii_comparison(simple_trajectory):
    """Test with different probe radii."""
    n_sphere_points = 960
    
    # Test different probe radii (in nm)
    probe_radii = [0.10, 0.14, 0.20]  # 1.0, 1.4, 2.0 Angstroms
    
    for probe_radius in probe_radii:
        # Calculate with MDTraj
        mdtraj_sasa = md.shrake_rupley(simple_trajectory,
                                        probe_radius=probe_radius,
                                        n_sphere_points=n_sphere_points)
        
        # Calculate with rsasa
        rsasa_params = mdtraj_to_rsasa_format(simple_trajectory,
                                              probe_radius=probe_radius)
        rsasa_result = sasa(
            rsasa_params['n_frames'],
            rsasa_params['n_atoms'],
            rsasa_params['xyzlist'],
            rsasa_params['atom_radii'],
            n_sphere_points,
            rsasa_params['atom_mapping'],
            rsasa_params['atom_selection_mask'],
            rsasa_params['n_groups'],
        )
        
        # Convert units
        mdtraj_sasa_angstrom = mdtraj_sasa * 100.0
        
        # Compare
        assert_allclose(rsasa_result, mdtraj_sasa_angstrom.flatten(), rtol=0.05)
        
        # Larger probe radius should give larger SASA
        if probe_radius > 0.10:
            assert rsasa_result.sum() > 0, "SASA should be positive"


def test_residue_level_comparison(simple_trajectory):
    """Test residue-level SASA aggregation."""
    probe_radius = 0.14
    n_sphere_points = 960
    
    # Calculate per-atom SASA with MDTraj
    mdtraj_sasa_atom = md.shrake_rupley(simple_trajectory,
                                         probe_radius=probe_radius,
                                         n_sphere_points=n_sphere_points)
    
    # Calculate residue-level SASA by summing
    # (simple_trajectory has 1 residue with 3 atoms)
    mdtraj_sasa_residue = mdtraj_sasa_atom.sum(axis=1, keepdims=True)
    
    # Calculate with rsasa at residue level
    rsasa_params = mdtraj_to_rsasa_format(simple_trajectory,
                                          probe_radius=probe_radius)
    # Map all atoms to group 0 (single residue)
    rsasa_params['atom_mapping'] = np.zeros(rsasa_params['n_atoms'], 
                                            dtype=np.uintp)
    rsasa_params['n_groups'] = 1
    
    rsasa_result = sasa(
        rsasa_params['n_frames'],
        rsasa_params['n_atoms'],
        rsasa_params['xyzlist'],
        rsasa_params['atom_radii'],
        n_sphere_points,
        rsasa_params['atom_mapping'],
        rsasa_params['atom_selection_mask'],
        rsasa_params['n_groups'],
    )
    
    # Convert units
    mdtraj_sasa_angstrom = mdtraj_sasa_residue * 100.0
    
    # Compare
    assert_allclose(rsasa_result, mdtraj_sasa_angstrom.flatten(), rtol=0.05)


@pytest.mark.skipif(not HAS_MDTRAJ, reason="MDTraj not installed")
def test_performance_comparison(multi_frame_trajectory):
    """
    Rough performance comparison between rsasa and MDTraj.
    This is not a precise benchmark, just a sanity check.
    """
    import time
    
    probe_radius = 0.14
    n_sphere_points = 960
    
    # Time MDTraj
    start = time.time()
    mdtraj_sasa = md.shrake_rupley(multi_frame_trajectory,
                                    probe_radius=probe_radius,
                                    n_sphere_points=n_sphere_points)
    mdtraj_time = time.time() - start
    
    # Time rsasa
    rsasa_params = mdtraj_to_rsasa_format(multi_frame_trajectory,
                                          probe_radius=probe_radius)
    start = time.time()
    rsasa_result = sasa(
        rsasa_params['n_frames'],
        rsasa_params['n_atoms'],
        rsasa_params['xyzlist'],
        rsasa_params['atom_radii'],
        n_sphere_points,
        rsasa_params['atom_mapping'],
        rsasa_params['atom_selection_mask'],
        rsasa_params['n_groups'],
    )
    rsasa_time = time.time() - start
    
    print(f"\nPerformance comparison:")
    print(f"  MDTraj: {mdtraj_time*1000:.2f} ms")
    print(f"  rsasa:  {rsasa_time*1000:.2f} ms")
    print(f"  Speedup: {mdtraj_time/rsasa_time:.2f}x")
    
    # Just verify results are close
    mdtraj_sasa_angstrom = mdtraj_sasa * 100.0
    n_frames = multi_frame_trajectory.n_frames
    n_atoms = multi_frame_trajectory.n_atoms
    rsasa_result_2d = rsasa_result.reshape(n_frames, n_atoms)
    assert_allclose(rsasa_result_2d, mdtraj_sasa_angstrom, rtol=0.05)


if __name__ == "__main__":
    # Allow running tests directly
    pytest.main([__file__, "-v", "-s"])