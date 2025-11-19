use std::f32::consts::PI;

/// Calculate the accessible surface area of each atom in a single snapshot
///
/// # Parameters
/// - `frame`: 2D array of shape `[n_atoms, 3]` containing the coordinates of the nuclei
/// - `n_atoms`: Number of atoms in the frame
/// - `atom_radii`: 1D array of shape `[n_atoms]` containing the van der Waals radii of the atoms plus the probe radius
/// - `sphere_points`: 2D array of shape `[n_sphere_points, 3]` containing uniformly distributed points on a sphere
/// - `n_sphere_points`: Number of sphere points
/// - `atom_selection_mask`: 1D array of shape `[n_atoms]` indicating whether the SASA should be computed for each atom
/// - `centered_sphere_points`: Work buffer 2D array of shape `[n_sphere_points, 3]` for intermediate calculations
/// - `neighbor_indices`: Work buffer 2D array of shape `[n_atoms]` for intermediate calculations
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

/// Calculate the accessible surface area of each atom in each frame of a trajectory
///
/// # Parameters
/// - `n_frames`: Number of frames in the trajectory
/// - `n_atoms`: Number of atoms in each frame
/// - `xyzlist`: 3D array of shape `[n_frames, n_atoms, 3]` containing the coordinates of the nuclei
/// - `atom_radii`: 1D array of shape `[n_atoms]` containing the van der Waals radii of the atoms plus the probe radius
/// - `n_sphere_points`: Number of points to generate sampling the unit sphere
/// - `atom_mapping`: Mapping from atoms onto groups, over which to accumulate the SASA
/// - `atom_selection_mask`: 1D array of shape `[n_atoms]` indicating whether the SASA should be computed for each atom
/// - `n_groups`: Number of groups
/// - `out`: Output buffer 2D array of shape `[n_frames, n_groups]` to place the results in
/// # Parameters
/// - `n_frames`: Number of frames in the trajectory
/// - `n_atoms`: Number of atoms in each frame
/// - `xyzlist`: 3D array of shape `[n_frames, n_atoms, 3]` containing the coordinates of the nuclei
/// - `atom_radii`: 1D array of shape `[n_atoms]` containing the van der Waals radii of the atoms plus the probe radius
/// - `n_sphere_points`: Number of points to generate sampling the unit sphere
/// - `atom_mapping`: Mapping from atoms onto groups, over which to accumulate the SASA
/// - `atom_selection_mask`: 1D array of shape `[n_atoms]` indicating whether the SASA should be computed for each atom
/// - `n_groups`: Number of groups
/// - `out`: Output buffer 2D array of shape `[n_frames, n_groups]` to place the results in
/// Calculate the accessible surface area of each atom in a single snapshot
///
/// # Parameters
/// - `frame`: 2D array of shape `[n_atoms, 3]` containing the coordinates of the nuclei
/// - `n_atoms`: Number of atoms in the frame
/// - `atom_radii`: 1D array of shape `[n_atoms]` containing the van der Waals radii of the atoms plus the probe radius
/// - `sphere_points`: 2D array of shape `[n_sphere_points, 3]` containing uniformly distributed points on a sphere
/// - `n_sphere_points`: Number of sphere points
/// - `atom_selection_mask`: 1D array of shape `[n_atoms]` indicating whether the SASA should be computed for each atom
/// - `centered_sphere_points`: Work buffer 2D array of shape `[n_sphere_points, 3]` for intermediate calculations
/// - `neighbor_indices`: Work buffer 2D array of shape `[n_atoms]` for intermediate calculations
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
) {
    // Generate the sphere points
    let mut sphere_points = vec![0.0f32; n_sphere_points * 3];
    generate_sphere_points(&mut sphere_points, n_sphere_points);

    // Work buffers
    let mut wb1 = vec![0usize; n_atoms];
    let mut wb2 = vec![0.0f32; 3 * n_sphere_points];
    let mut outframebuffer = vec![0.0f32; n_atoms];

    for i in 0..n_frames {
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
}
    // Generate the sphere points
    let mut sphere_points = vec![0.0f32; n_sphere_points * 3];
    generate_sphere_points(&mut sphere_points, n_sphere_points);

    // Work buffers
    let mut wb1 = vec![0usize; n_atoms];
    let mut wb2 = vec![0.0f32; 3 * n_sphere_points];
    let mut outframebuffer = vec![0.0f32; n_atoms];

    for i in 0..n_frames {
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
        let mut sphere_points = vec![0.0; n_sphere_points * 3];
        generate_sphere_points(&mut sphere_points, n_sphere_points);
        let mut neighbor_indices = vec![0; n_atoms];
        let mut centered_sphere_points = vec![0.0; n_sphere_points * 3];
        let atom_selection_mask = vec![1];
        let mut areas = vec![0.0; n_atoms];

        asa_frame(
            &frame,
            n_atoms,
            &atom_radii,
            &sphere_points,
            n_sphere_points,
            &mut neighbor_indices,
            &mut centered_sphere_points,
            &atom_selection_mask,
            &mut areas,
        );

        // The area should be close to the surface area of a sphere with radius 1.0
        let expected_area = 4.0 * PI * atom_radii[0] * atom_radii[0];
        assert!((areas[0] - expected_area).abs() < 0.1 * expected_area);
    }

    #[test]
    fn test_two_atoms() {
        // Test with two atoms
        let frame = vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0];
        let n_atoms = 2;
        let atom_radii = vec![1.0, 1.0];
        let n_sphere_points = 100;
        let mut sphere_points = vec![0.0; n_sphere_points * 3];
        generate_sphere_points(&mut sphere_points, n_sphere_points);
        let mut neighbor_indices = vec![0; n_atoms];
        let mut centered_sphere_points = vec![0.0; n_sphere_points * 3];
        let atom_selection_mask = vec![1, 1];
        let mut areas = vec![0.0; n_atoms];

        asa_frame(
            &frame,
            n_atoms,
            &atom_radii,
            &sphere_points,
            n_sphere_points,
            &mut neighbor_indices,
            &mut centered_sphere_points,
            &atom_selection_mask,
            &mut areas,
        );

        // The areas should be less than the surface area of a sphere with radius 1.0
        // because the atoms are close to each other
        let expected_area = 4.0 * PI * atom_radii[0] * atom_radii[0];
        assert!(areas[0] < expected_area);
        assert!(areas[1] < expected_area);
    }

    #[test]
    fn test_three_atoms() {
        // Test with three atoms
        let frame = vec![0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0, 0.0];
        let n_atoms = 3;
        let atom_radii = vec![1.0, 1.0, 1.0];
        let n_sphere_points = 100;
        let mut sphere_points = vec![0.0; n_sphere_points * 3];
        generate_sphere_points(&mut sphere_points, n_sphere_points);
        let mut neighbor_indices = vec![0; n_atoms];
        let mut centered_sphere_points = vec![0.0; n_sphere_points * 3];
        let atom_selection_mask = vec![1, 1, 1];
        let mut areas = vec![0.0; n_atoms];

        asa_frame(
            &frame,
            n_atoms,
            &atom_radii,
            &sphere_points,
            n_sphere_points,
            &mut neighbor_indices,
            &mut centered_sphere_points,
            &atom_selection_mask,
            &mut areas,
        );

        // The areas should be less than the surface area of a sphere with radius 1.0
        // because the atoms are close to each other
        let expected_area = 4.0 * PI * atom_radii[0] * atom_radii[0];
        assert!(areas[0] < expected_area);
        assert!(areas[1] < expected_area);
        assert!(areas[2] < expected_area);
    }
}
