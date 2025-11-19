use std::f32::consts::PI;
use std::fmt;

#[cfg(feature = "python")]
use numpy::{PyArray1, PyReadonlyArray1};
#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Error types for SASA calculation
#[derive(Debug, Clone)]
pub enum SasaError {
    /// Input arrays have incompatible shapes
    ShapeMismatch(String),
    /// Invalid parameter values
    InvalidParameter(String),
}

impl fmt::Display for SasaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SasaError::ShapeMismatch(msg) => write!(f, "Shape mismatch: {}", msg),
            SasaError::InvalidParameter(msg) => write!(f, "Invalid parameter: {}", msg),
        }
    }
}

impl std::error::Error for SasaError {}

#[cfg(feature = "python")]
impl From<SasaError> for PyErr {
    fn from(err: SasaError) -> PyErr {
        pyo3::exceptions::PyValueError::new_err(err.to_string())
    }
}

/// Calculate the accessible surface area of each atom in a single snapshot
///
/// # Parameters
/// - `frame`: 2D array of shape `[n_atoms, 3]` containing the coordinates of the nuclei
/// - `n_atoms`: Number of atoms in the frame
/// - `atom_radii`: 1D array of shape `[n_atoms]` containing the van der Waals radii of the atoms plus the probe radius
/// - `sphere_points`: 2D array of shape `[n_sphere_points, 3]` containing uniformly distributed points on a sphere
/// - `n_sphere_points`: Number of sphere points
/// - `neighbor_indices`: Work buffer 2D array of shape `[n_atoms]` for intermediate calculations
/// - `centered_sphere_points`: Work buffer 2D array of shape `[n_sphere_points, 3]` for intermediate calculations
/// - `atom_selection_mask`: 1D array of shape `[n_atoms]` indicating whether the SASA should be computed for each atom
/// - `areas`: Output buffer 1D array of shape `[n_atoms]` to place the results in
fn asa_frame(
    frame: &[f32],
    n_atoms: usize,
    atom_radii: &[f32],
    sphere_points: &[f32],
    n_sphere_points: usize,
    neighbor_indices: &mut [usize],
    centered_sphere_points: &mut [f32],
    atom_selection_mask: &[i32],
    areas: &mut [f32],
) {
    let constant = 4.0 * PI / n_sphere_points as f32;

    for i in 0..n_atoms {
        // Skip atom if not in selection
        let in_selection = atom_selection_mask[i];
        if in_selection == 0 {
            continue;
        }

        let atom_radius_i = atom_radii[i];
        let r_i = [frame[i * 3], frame[i * 3 + 1], frame[i * 3 + 2]];

        // Get all the atoms close to atom `i`
        let mut n_neighbor_indices = 0;
        for j in 0..n_atoms {
            if i == j {
                continue;
            }

            let r_j = [frame[j * 3], frame[j * 3 + 1], frame[j * 3 + 2]];
            let r_ij = [r_i[0] - r_j[0], r_i[1] - r_j[1], r_i[2] - r_j[2]];
            let atom_radius_j = atom_radii[j];

            // Look for atoms `j` that are nearby atom `i`
            let radius_cutoff = atom_radius_i + atom_radius_j;
            let radius_cutoff2 = radius_cutoff * radius_cutoff;
            let r2 = r_ij[0] * r_ij[0] + r_ij[1] * r_ij[1] + r_ij[2] * r_ij[2];
            if r2 < radius_cutoff2 {
                neighbor_indices[n_neighbor_indices] = j;
                n_neighbor_indices += 1;
            }
            if r2 < 1e-10f32 {
                panic!(
                    "ERROR: THIS CODE IS KNOWN TO FAIL WHEN ATOMS ARE VIRTUALLY ON TOP OF ONE ANOTHER. YOU SUPPLIED TWO ATOMS {} APART. QUITTING NOW",
                    r2.sqrt()
                );
            }
        }

        // Center the sphere points on atom i
        for j in 0..n_sphere_points {
            centered_sphere_points[3 * j] = frame[3 * i] + atom_radius_i * sphere_points[3 * j];
            centered_sphere_points[3 * j + 1] =
                frame[3 * i + 1] + atom_radius_i * sphere_points[3 * j + 1];
            centered_sphere_points[3 * j + 2] =
                frame[3 * i + 2] + atom_radius_i * sphere_points[3 * j + 2];
        }

        // Check if each of these points is accessible
        let mut k_closest_neighbor = 0;
        for j in 0..n_sphere_points {
            let mut is_accessible = true;
            let r_j = [
                centered_sphere_points[3 * j],
                centered_sphere_points[3 * j + 1],
                centered_sphere_points[3 * j + 2],
            ];

            // Iterate through the sphere points by cycling through them
            // in a circle, starting with k_closest_neighbor and then wrapping
            // around
            for k in k_closest_neighbor..n_neighbor_indices + k_closest_neighbor {
                let k_prime = k % n_neighbor_indices;
                let r = atom_radii[neighbor_indices[k_prime]];

                let index = neighbor_indices[k_prime];
                let r_jk = [
                    r_j[0] - frame[3 * index],
                    r_j[1] - frame[3 * index + 1],
                    r_j[2] - frame[3 * index + 2],
                ];
                if r_jk[0] * r_jk[0] + r_jk[1] * r_jk[1] + r_jk[2] * r_jk[2] < r * r {
                    k_closest_neighbor = k;
                    is_accessible = false;
                    break;
                }
            }

            if is_accessible {
                areas[i] += 1.0;
            }
        }

        areas[i] *= constant * atom_radii[i] * atom_radii[i];
    }
}

/// Generate points on a sphere using the Golden Section Spiral algorithm
///
/// # Parameters
/// - `sphere_points`: Empty array of length `n_points * 3` to be filled with the points
/// - `n_points`: Number of points to generate on the sphere
fn generate_sphere_points(sphere_points: &mut [f32], n_points: usize) {
    let inc = PI * (3.0 - 5.0f32.sqrt());
    let offset = 2.0 / n_points as f32;

    for i in 0..n_points {
        let y = i as f32 * offset - 1.0 + (offset / 2.0);
        let r = (1.0 - y * y).sqrt();
        let phi = i as f32 * inc;

        sphere_points[3 * i] = phi.cos() * r;
        sphere_points[3 * i + 1] = y;
        sphere_points[3 * i + 2] = phi.sin() * r;
    }
}

/// Validate input parameters for SASA calculation
fn validate_inputs(
    n_frames: usize,
    n_atoms: usize,
    xyzlist: &[f32],
    atom_radii: &[f32],
    n_sphere_points: usize,
    atom_mapping: &[usize],
    atom_selection_mask: &[i32],
    n_groups: usize,
) -> Result<(), SasaError> {
    // Validate n_sphere_points
    if n_sphere_points == 0 {
        return Err(SasaError::InvalidParameter(
            "n_sphere_points must be greater than 0".to_string(),
        ));
    }

    // Validate n_atoms
    if n_atoms == 0 {
        return Err(SasaError::InvalidParameter(
            "n_atoms must be greater than 0".to_string(),
        ));
    }

    // Validate n_frames
    if n_frames == 0 {
        return Err(SasaError::InvalidParameter(
            "n_frames must be greater than 0".to_string(),
        ));
    }

    // Validate n_groups
    if n_groups == 0 {
        return Err(SasaError::InvalidParameter(
            "n_groups must be greater than 0".to_string(),
        ));
    }

    // Validate xyzlist shape
    let expected_xyz_len = n_frames * n_atoms * 3;
    if xyzlist.len() != expected_xyz_len {
        return Err(SasaError::ShapeMismatch(format!(
            "xyzlist has length {} but expected {} (n_frames={} * n_atoms={} * 3)",
            xyzlist.len(),
            expected_xyz_len,
            n_frames,
            n_atoms
        )));
    }

    // Validate atom_radii shape
    if atom_radii.len() != n_atoms {
        return Err(SasaError::ShapeMismatch(format!(
            "atom_radii has length {} but expected {} (n_atoms)",
            atom_radii.len(),
            n_atoms
        )));
    }

    // Validate atom_mapping shape
    if atom_mapping.len() != n_atoms {
        return Err(SasaError::ShapeMismatch(format!(
            "atom_mapping has length {} but expected {} (n_atoms)",
            atom_mapping.len(),
            n_atoms
        )));
    }

    // Validate atom_selection_mask shape
    if atom_selection_mask.len() != n_atoms {
        return Err(SasaError::ShapeMismatch(format!(
            "atom_selection_mask has length {} but expected {} (n_atoms)",
            atom_selection_mask.len(),
            n_atoms
        )));
    }

    // Validate atom_mapping values
    for (i, &group_idx) in atom_mapping.iter().enumerate() {
        if group_idx >= n_groups {
            return Err(SasaError::InvalidParameter(format!(
                "atom_mapping[{}] = {} is >= n_groups = {}",
                i, group_idx, n_groups
            )));
        }
    }

    // Validate atom_radii values (should be positive)
    for (i, &radius) in atom_radii.iter().enumerate() {
        if radius <= 0.0 {
            return Err(SasaError::InvalidParameter(format!(
                "atom_radii[{}] = {} must be positive",
                i, radius
            )));
        }
    }

    Ok(())
}

/// Calculate the solvent accessible surface area (SASA) for atoms in a trajectory
///
/// # Parameters
/// - `n_frames`: Number of frames in the trajectory
/// - `n_atoms`: Number of atoms in each frame
/// - `xyzlist`: Flattened array of coordinates with shape `[n_frames * n_atoms * 3]`
/// - `atom_radii`: Array of van der Waals radii plus probe radius with shape `[n_atoms]`
/// - `n_sphere_points`: Number of points to generate sampling the unit sphere
/// - `atom_mapping`: Mapping from atoms onto groups with shape `[n_atoms]`
/// - `atom_selection_mask`: Array indicating which atoms to compute SASA for with shape `[n_atoms]`
/// - `n_groups`: Number of groups
///
/// # Returns
/// - `Result<Vec<f32>, SasaError>`: Vector of SASA values with shape `[n_frames * n_groups]`
///
/// # Errors
/// - Returns `SasaError::ShapeMismatch` if input arrays have incorrect shapes
/// - Returns `SasaError::InvalidParameter` if parameters are invalid
pub fn sasa(
    n_frames: usize,
    n_atoms: usize,
    xyzlist: &[f32],
    atom_radii: &[f32],
    n_sphere_points: usize,
    atom_mapping: &[usize],
    atom_selection_mask: &[i32],
    n_groups: usize,
) -> Result<Vec<f32>, SasaError> {
    // Validate all inputs
    validate_inputs(
        n_frames,
        n_atoms,
        xyzlist,
        atom_radii,
        n_sphere_points,
        atom_mapping,
        atom_selection_mask,
        n_groups,
    )?;

    // Allocate output buffer
    let mut out = vec![0.0f32; n_frames * n_groups];

    // Generate the sphere points
    let mut sphere_points = vec![0.0f32; n_sphere_points * 3];
    generate_sphere_points(&mut sphere_points, n_sphere_points);

    // Work buffers
    let mut wb1 = vec![0usize; n_atoms];
    let mut wb2 = vec![0.0f32; 3 * n_sphere_points];
    let mut outframebuffer = vec![0.0f32; n_atoms];

    for i in 0..n_frames {
        // Reset output frame buffer
        outframebuffer.fill(0.0);

        asa_frame(
            &xyzlist[i * n_atoms * 3..(i + 1) * n_atoms * 3],
            n_atoms,
            atom_radii,
            &sphere_points,
            n_sphere_points,
            &mut wb1,
            &mut wb2,
            atom_selection_mask,
            &mut outframebuffer,
        );

        let outframe = &mut out[i * n_groups..(i + 1) * n_groups];
        for j in 0..n_atoms {
            outframe[atom_mapping[j]] += outframebuffer[j];
        }
    }

    Ok(out)
}

#[cfg(feature = "python")]
#[pyfunction]
#[pyo3(name = "sasa")]
fn py_sasa<'py>(
    py: Python<'py>,
    n_frames: usize,
    n_atoms: usize,
    xyzlist: PyReadonlyArray1<f32>,
    atom_radii: PyReadonlyArray1<f32>,
    n_sphere_points: usize,
    atom_mapping: PyReadonlyArray1<usize>,
    atom_selection_mask: PyReadonlyArray1<i32>,
    n_groups: usize,
) -> PyResult<pyo3::Bound<'py, PyArray1<f32>>> {
    // Convert numpy arrays to slices
    let xyzlist_slice = xyzlist.as_slice()?;
    let atom_radii_slice = atom_radii.as_slice()?;
    let atom_mapping_slice = atom_mapping.as_slice()?;
    let atom_selection_mask_slice = atom_selection_mask.as_slice()?;

    // Call the Rust function
    let result = sasa(
        n_frames,
        n_atoms,
        xyzlist_slice,
        atom_radii_slice,
        n_sphere_points,
        atom_mapping_slice,
        atom_selection_mask_slice,
        n_groups,
    )?;

    // Convert result to numpy array
    Ok(PyArray1::from_vec_bound(py, result))
}

#[cfg(feature = "python")]
#[pymodule]
fn rsasa(m: &pyo3::Bound<'_, pyo3::types::PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(py_sasa, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_atom() {
        // Test with a single atom
        let frame = vec![0.0, 0.0, 0.0];
        let n_atoms = 1;
        let atom_radii = vec![1.0];
        let n_sphere_points = 100;
        let atom_mapping = vec![0];
        let atom_selection_mask = vec![1];
        let n_groups = 1;

        let result = sasa(
            1, // n_frames
            n_atoms,
            &frame,
            &atom_radii,
            n_sphere_points,
            &atom_mapping,
            &atom_selection_mask,
            n_groups,
        )
        .unwrap();

        // The area should be close to the surface area of a sphere with radius 1.0
        let expected_area = 4.0 * PI * atom_radii[0] * atom_radii[0];
        assert!((result[0] - expected_area).abs() < 0.1 * expected_area);
    }

    #[test]
    fn test_two_atoms() {
        // Test with two atoms
        let frame = vec![0.0, 0.0, 0.0, 1.5, 0.0, 0.0];
        let n_atoms = 2;
        let atom_radii = vec![1.0, 1.0];
        let n_sphere_points = 100;
        let atom_mapping = vec![0, 1];
        let atom_selection_mask = vec![1, 1];
        let n_groups = 2;

        let result = sasa(
            1, // n_frames
            n_atoms,
            &frame,
            &atom_radii,
            n_sphere_points,
            &atom_mapping,
            &atom_selection_mask,
            n_groups,
        )
        .unwrap();

        // The areas should be less than the surface area of a sphere with radius 1.0
        // because the atoms are close to each other
        let expected_area = 4.0 * PI * atom_radii[0] * atom_radii[0];
        assert!(result[0] < expected_area);
        assert!(result[1] < expected_area);
    }

    #[test]
    fn test_three_atoms() {
        // Test with three atoms
        let frame = vec![0.0, 0.0, 0.0, 1.5, 0.0, 0.0, 0.0, 1.5, 0.0];
        let n_atoms = 3;
        let atom_radii = vec![1.0, 1.0, 1.0];
        let n_sphere_points = 100;
        let atom_mapping = vec![0, 1, 2]; // Each atom in its own group
        let atom_selection_mask = vec![1, 1, 1];
        let n_groups = 3;

        let result = sasa(
            1, // n_frames
            n_atoms,
            &frame,
            &atom_radii,
            n_sphere_points,
            &atom_mapping,
            &atom_selection_mask,
            n_groups,
        )
        .unwrap();

        // The areas should be less than the surface area of a sphere with radius 1.0
        // because the atoms are close to each other
        let expected_area = 4.0 * PI * atom_radii[0] * atom_radii[0];
        assert!(result[0] < expected_area);
        assert!(result[1] < expected_area);
        assert!(result[2] < expected_area);
    }

    #[test]
    fn test_shape_mismatch() {
        let n_frames = 1;
        let n_atoms = 2;
        let xyzlist = vec![0.0, 0.0, 0.0]; // Wrong size - should be 6
        let atom_radii = vec![1.0, 1.0];
        let n_sphere_points = 100;
        let atom_mapping = vec![0, 0];
        let atom_selection_mask = vec![1, 1];
        let n_groups = 1;

        let result = sasa(
            n_frames,
            n_atoms,
            &xyzlist,
            &atom_radii,
            n_sphere_points,
            &atom_mapping,
            &atom_selection_mask,
            n_groups,
        );

        assert!(result.is_err());
        if let Err(SasaError::ShapeMismatch(msg)) = result {
            assert!(msg.contains("xyzlist"));
        } else {
            panic!("Expected ShapeMismatch error");
        }
    }

    #[test]
    fn test_invalid_parameters() {
        let n_frames = 1;
        let n_atoms = 2;
        let xyzlist = vec![0.0; 6];
        let atom_radii = vec![1.0, 1.0];
        let n_sphere_points = 0; // Invalid
        let atom_mapping = vec![0, 0];
        let atom_selection_mask = vec![1, 1];
        let n_groups = 1;

        let result = sasa(
            n_frames,
            n_atoms,
            &xyzlist,
            &atom_radii,
            n_sphere_points,
            &atom_mapping,
            &atom_selection_mask,
            n_groups,
        );

        assert!(result.is_err());
        if let Err(SasaError::InvalidParameter(msg)) = result {
            assert!(msg.contains("n_sphere_points"));
        } else {
            panic!("Expected InvalidParameter error");
        }
    }

    #[test]
    fn test_invalid_atom_mapping() {
        let n_frames = 1;
        let n_atoms = 2;
        let xyzlist = vec![0.0; 6];
        let atom_radii = vec![1.0, 1.0];
        let n_sphere_points = 100;
        let atom_mapping = vec![0, 5]; // 5 is >= n_groups
        let atom_selection_mask = vec![1, 1];
        let n_groups = 2;

        let result = sasa(
            n_frames,
            n_atoms,
            &xyzlist,
            &atom_radii,
            n_sphere_points,
            &atom_mapping,
            &atom_selection_mask,
            n_groups,
        );

        assert!(result.is_err());
        if let Err(SasaError::InvalidParameter(msg)) = result {
            assert!(msg.contains("atom_mapping"));
        } else {
            panic!("Expected InvalidParameter error");
        }
    }

    #[test]
    fn test_negative_radius() {
        let n_frames = 1;
        let n_atoms = 1;
        let xyzlist = vec![0.0; 3];
        let atom_radii = vec![-1.0]; // Invalid
        let n_sphere_points = 100;
        let atom_mapping = vec![0];
        let atom_selection_mask = vec![1];
        let n_groups = 1;

        let result = sasa(
            n_frames,
            n_atoms,
            &xyzlist,
            &atom_radii,
            n_sphere_points,
            &atom_mapping,
            &atom_selection_mask,
            n_groups,
        );

        assert!(result.is_err());
        if let Err(SasaError::InvalidParameter(msg)) = result {
            assert!(msg.contains("atom_radii"));
        } else {
            panic!("Expected InvalidParameter error");
        }
    }

    #[test]
    fn test_multiple_frames() {
        let n_frames = 2;
        let n_atoms = 2;
        let xyzlist = vec![
            // Frame 1
            0.0, 0.0, 0.0, 2.0, 0.0, 0.0, // Frame 2
            0.0, 0.0, 0.0, 1.5, 0.0, 0.0,
        ];
        let atom_radii = vec![1.0, 1.0];
        let n_sphere_points = 100;
        let atom_mapping = vec![0, 1];
        let atom_selection_mask = vec![1, 1];
        let n_groups = 2;

        let result = sasa(
            n_frames,
            n_atoms,
            &xyzlist,
            &atom_radii,
            n_sphere_points,
            &atom_mapping,
            &atom_selection_mask,
            n_groups,
        )
        .unwrap();

        assert_eq!(result.len(), n_frames * n_groups);

        // Frame 1: atoms are touching (distance = 2.0, sum of radii = 2.0)
        // Frame 2: atoms are overlapping (distance = 1.5, sum of radii = 2.0)
        // So frame 2 should have smaller SASA than frame 1
        let frame1_total = result[0] + result[1];
        let frame2_total = result[2] + result[3];
        assert!(frame2_total < frame1_total);
    }
}
