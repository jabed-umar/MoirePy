use std::collections::BTreeMap;
use wasm_bindgen::prelude::*;

// ---------------------------------------------------------------------------
// Internal helpers (mirrors of the JS helper functions)
// ---------------------------------------------------------------------------

fn angle_from_x(p: [f64; 2]) -> f64 {
    p[1].atan2(p[0]).to_degrees()
}

/// Generate all lattice points within `radius` for the lattice spanned by lv1, lv2.
fn generate_lattice_points(lv1: [f64; 2], lv2: [f64; 2], radius: f64) -> Vec<[f64; 2]> {
    let mut points = Vec::new();
    // lv2[1] is the y-component; safe guard against zero
    let max_grid = if lv2[1].abs() > 1e-12 {
        (radius / lv2[1].abs()).floor() as i64 + 5
    } else {
        (radius + 5.0) as i64
    };

    for i in -max_grid..=max_grid {
        for j in -max_grid..=max_grid {
            let point = [
                i as f64 * lv1[0] + j as f64 * lv2[0],
                i as f64 * lv1[1] + j as f64 * lv2[1],
            ];
            let dist_sq = point[0] * point[0] + point[1] * point[1];
            if dist_sq <= radius * radius {
                points.push(point);
            }
        }
    }
    points
}

/// Group lattice points by their distance from the origin, rounded to `tol` decimal places.
/// Key = round(distance * 10^tol) as i64, Value = list of points at that distance.
fn process_lattice(points: &[[f64; 2]], tol: u32) -> BTreeMap<i64, Vec<[f64; 2]>> {
    let scale = 10f64.powi(tol as i32);
    let mut dist_map: BTreeMap<i64, Vec<[f64; 2]>> = BTreeMap::new();
    for &p in points {
        let dist = (p[0] * p[0] + p[1] * p[1]).sqrt();
        let key = (dist * scale).round() as i64;
        dist_map.entry(key).or_default().push(p);
    }
    dist_map
}

/// Solve p = n * lv1 + m * lv2 for integers [n, m].
fn calc_indices(p: [f64; 2], lv1: [f64; 2], lv2: [f64; 2]) -> [i64; 2] {
    let [a, b] = lv1;
    let [c, d] = lv2;
    let [x, y] = p;
    let det = a * d - b * c;
    let nx = (d * x - c * y) / det;
    let ny = (a * y - b * x) / det;
    [nx.round() as i64, ny.round() as i64]
}

// ---------------------------------------------------------------------------
// Public WASM export
// ---------------------------------------------------------------------------

/// Find commensurate twist angles for a bilayer system.
///
/// # Arguments
/// * `radius`           – truncation radius for both lattices
/// * `layer1_vectors`   – flat [lv1x, lv1y, lv2x, lv2y] for layer 1
/// * `layer2_vectors`   – flat [lv1x, lv1y, lv2x, lv2y] for layer 2 (currently unused –
///                        both lattices are generated from layer1_vectors, matching the JS)
/// * `tol`              – decimal precision for distance grouping
/// * `angle_start`      – lower bound of the angle search window (degrees, exclusive)
/// * `angle_end`        – upper bound of the angle search window (degrees, exclusive)
///
/// # Returns
/// A flat `Vec<f64>` with **7 values per result row**, in this order:
///   `angle_deg | angle_rad | ll1 | ll2 | ul1 | ul2 | cells`
/// Total length = 7 × N results.
#[wasm_bindgen]
pub fn find_values(
    radius: f64,
    layer1_vectors: &[f64],
    _layer2_vectors: &[f64], // reserved – both lattices use layer1 vectors (matches JS behaviour)
    tol: u32,
    angle_start: f64,
    angle_end: f64,
) -> Vec<f64> {
    let lv1 = [layer1_vectors[0], layer1_vectors[1]];
    let lv2 = [layer1_vectors[2], layer1_vectors[3]];

    // Beta = angle between the two basis vectors (used for supercell orientation)
    let beta = lv2[1].atan2(lv2[0]) - lv1[1].atan2(lv1[0]);
    let cos_b = beta.cos();
    let sin_b = beta.sin();

    let lattice1 = generate_lattice_points(lv1, lv2, radius);
    let lattice2 = generate_lattice_points(lv1, lv2, radius); // same as lattice1

    let dict1 = process_lattice(&lattice1, tol);
    let dict2 = process_lattice(&lattice2, tol);

    let threshold = 10f64.powi(-(tol as i32));

    // Common distance keys, sorted ascending, skipping the zero-distance shell (first key)
    let mut common_keys: Vec<i64> = dict1
        .keys()
        .filter(|k| dict2.contains_key(k))
        .copied()
        .collect();
    common_keys.sort_unstable();
    if !common_keys.is_empty() {
        common_keys.remove(0); // remove key for distance ≈ 0
    }

    // ------------------------------------------------------------------
    // 1. Collect all valid candidates
    // ------------------------------------------------------------------
    struct Candidate {
        angle: f64,
        p1: [f64; 2],
        p2: [f64; 2],
        shell_radius: f64, // |p1|, used for clash resolution
        phi_res: f64,      // canonical orientation, used for clash resolution
    }

    let mut candidates: Vec<Candidate> = Vec::new();

    for key in &common_keys {
        let pts1 = &dict1[key];
        let pts2 = &dict2[key];

        for &p1 in pts1 {
            if p1[1] < 0.0 { continue; } // below x-axis → mirror duplicate, skip
            let theta1 = angle_from_x(p1);
            let shell_r = (p1[0] * p1[0] + p1[1] * p1[1]).sqrt();

            // V_sum = mlv1 + mlv2, where mlv2 = p1 rotated by beta
            let v_sum = [
                p1[0] + (p1[0] * cos_b - p1[1] * sin_b),
                p1[1] + (p1[0] * sin_b + p1[1] * cos_b),
            ];
            let phi_res = (v_sum[1].atan2(v_sum[0]).to_degrees() + 360.0).rem_euclid(360.0);

            for &p2 in pts2 {
                if p2[1] < 0.0 { continue; } // below x-axis → mirror duplicate, skip
                let theta2 = angle_from_x(p2);
                let angle = (theta2 - theta1 + 360.0).rem_euclid(360.0);

                // Step A: range filter (exclusive bounds, matching JS `<= angle_start || >= angle_end`)
                if angle <= angle_start || angle >= angle_end {
                    continue;
                }

                candidates.push(Candidate {
                    angle,
                    p1,
                    p2,
                    shell_radius: shell_r,
                    phi_res,
                });
            }
        }
    }

    // ------------------------------------------------------------------
    // 2. Sort candidates by twist angle
    // ------------------------------------------------------------------
    candidates.sort_by(|a, b| a.angle.partial_cmp(&b.angle).unwrap_or(std::cmp::Ordering::Equal));

    // ------------------------------------------------------------------
    // 3. Linear scan: keep one winner per unique angle (clash resolution)
    // ------------------------------------------------------------------
    let mut winning_indices: Vec<usize> = Vec::new();

    if !candidates.is_empty() {
        let mut best = 0usize;

        for i in 1..candidates.len() {
            if (candidates[i].angle - candidates[best].angle).abs() < threshold {
                // Clash: smaller shell radius wins; tie-break on phi_res
                if candidates[i].shell_radius < candidates[best].shell_radius {
                    best = i;
                } else if (candidates[i].shell_radius - candidates[best].shell_radius).abs()
                    < threshold
                    && candidates[i].phi_res < candidates[best].phi_res
                {
                    best = i;
                }
            } else {
                winning_indices.push(best);
                best = i;
            }
        }
        winning_indices.push(best); // final winner
    }

    // ------------------------------------------------------------------
    // 4. Build flat output: 7 f64 values per row
    //    [ angle_deg, angle_rad, ll1, ll2, ul1, ul2, cells ]
    // ------------------------------------------------------------------
    let mut output = Vec::with_capacity(winning_indices.len() * 7);

    for idx in winning_indices {
        let c = &candidates[idx];
        let angle_rad = c.angle.to_radians();

        // calc_indices returns [n, m] s.t. p = n*lv1 + m*lv2
        let [i1, j1] = calc_indices(c.p1, lv1, lv2);
        let [i2, j2] = calc_indices(c.p2, lv1, lv2);

        // Estimated cells in the moiré supercell (both layers)
        let cells = (2.0 * (c.p1[0] * c.p1[0] + c.p1[1] * c.p1[1])).round();

        output.push(c.angle);  // angle_deg
        output.push(angle_rad); // angle_rad
        output.push(i2 as f64); // ll1  (lower layer index 1)
        output.push(j2 as f64); // ll2  (lower layer index 2)
        output.push(i1 as f64); // ul1  (upper layer index 1)
        output.push(j1 as f64); // ul2  (upper layer index 2)
        output.push(cells);     // cells
    }

    output
}
