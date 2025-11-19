//! Simple example demonstrating SASA calculation for a small molecular system

use rsasa::sasa;

fn main() {
    println!("RSASA - Solvent Accessible Surface Area Calculator");
    println!("==================================================\n");

    // Example 1: Single isolated atom
    println!("Example 1: Single isolated atom");
    println!("--------------------------------");
    let n_frames = 1;
    let n_atoms = 1;
    let xyzlist = vec![0.0, 0.0, 0.0]; // Single atom at origin
    let atom_radii = vec![1.5]; // Radius of 1.5 Å
    let n_sphere_points = 960; // Good balance of accuracy and speed
    let atom_mapping = vec![0]; // Single group
    let atom_selection_mask = vec![1]; // Compute SASA for this atom
    let n_groups = 1;

    let out = sasa(
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

    let expected_area = 4.0 * std::f32::consts::PI * atom_radii[0] * atom_radii[0];
    println!("  Calculated SASA: {:.2} Ų", out[0]);
    println!("  Expected (full sphere): {:.2} Ų", expected_area);
    println!(
        "  Difference: {:.2}%\n",
        (out[0] - expected_area).abs() / expected_area * 100.0
    );

    // Example 2: Two overlapping atoms
    println!("Example 2: Two overlapping atoms");
    println!("----------------------------------");
    let n_atoms = 2;
    let xyzlist = vec![
        0.0, 0.0, 0.0, // First atom at origin
        1.8, 0.0, 0.0, // Second atom 1.8 Å away (overlapping)
    ];
    let atom_radii = vec![1.5, 1.5];
    let atom_mapping = vec![0, 1]; // Each atom in its own group
    let atom_selection_mask = vec![1, 1];
    let n_groups = 2;

    let out = sasa(
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

    let full_sphere = 4.0 * std::f32::consts::PI * atom_radii[0] * atom_radii[0];
    println!(
        "  Atom 1 SASA: {:.2} Ų ({:.1}% of full sphere)",
        out[0],
        out[0] / full_sphere * 100.0
    );
    println!(
        "  Atom 2 SASA: {:.2} Ų ({:.1}% of full sphere)",
        out[1],
        out[1] / full_sphere * 100.0
    );
    println!("  Total SASA: {:.2} Ų\n", out[0] + out[1]);

    // Example 3: Three atoms with grouping
    println!("Example 3: Three atoms with residue-level grouping");
    println!("---------------------------------------------------");
    let n_atoms = 3;
    let xyzlist = vec![
        0.0, 0.0, 0.0, // Atom 1
        1.8, 0.0, 0.0, // Atom 2
        0.9, 1.5, 0.0, // Atom 3
    ];
    let atom_radii = vec![1.2, 1.4, 1.3];
    let atom_mapping = vec![0, 0, 1]; // Atoms 1-2 in group 0, atom 3 in group 1
    let atom_selection_mask = vec![1, 1, 1];
    let n_groups = 2;

    let out = sasa(
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

    println!("  Group 0 (atoms 1-2) SASA: {:.2} Ų", out[0]);
    println!("  Group 1 (atom 3) SASA: {:.2} Ų", out[1]);
    println!("  Total SASA: {:.2} Ų\n", out[0] + out[1]);

    // Example 4: Selective computation
    println!("Example 4: Selective SASA computation");
    println!("--------------------------------------");
    let n_atoms = 3;
    let xyzlist = vec![0.0, 0.0, 0.0, 1.8, 0.0, 0.0, 0.9, 1.5, 0.0];
    let atom_radii = vec![1.2, 1.4, 1.3];
    let atom_mapping = vec![0, 1, 2];
    let atom_selection_mask = vec![1, 0, 1]; // Only compute for atoms 1 and 3
    let n_groups = 3;

    let out = sasa(
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

    println!("  Atom 1 SASA (computed): {:.2} Ų", out[0]);
    println!("  Atom 2 SASA (skipped): {:.2} Ų", out[1]);
    println!("  Atom 3 SASA (computed): {:.2} Ų", out[2]);
    println!("\nNote: Atom 2 was skipped due to selection mask = 0\n");
}
