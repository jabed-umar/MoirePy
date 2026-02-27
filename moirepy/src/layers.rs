use kiddo::immutable::float::kdtree::ImmutableKdTree;
use ndarray::{array, Array1, Array2, ArrayView2};
use numpy::{PyArray1, PyArray2, PyReadonlyArray2, ToPyArray};
use pyo3::prelude::*;
use std::collections::HashMap;

// A: f64 (float type), T: u64 (ID type), K: 2 (Dimensions), B: 32 (Bucket size)
type LatticeTree = ImmutableKdTree<f64, u64, 2, 32>;

#[pyclass]
pub struct Layer {
    // --- 1. Fundamental Unit Cell Geometry ---
    // #[getter]
    pub lv1: Array1<f64>, // Primary lattice vector 1 (rotated if applicable)
    // #[getter]
    pub lv2: Array1<f64>, // Primary lattice vector 2 (rotated if applicable)

    // #[getter]
    pub basis_points: Array2<f64>, // Positions of atoms within the unit cell

    pub basis_type_maps: HashMap<String, usize>, // Type name → type ID (Internal use only)

    #[pyo3(get)]
    pub basis_types: Vec<String>, // Reverse mapping by index: ["A", "B", ...]

    // #[getter]
    pub neighbours: Vec<Array2<f64>>, // Each index corresponds to a type ID (e.g. 0 for "A")

    // --- 2. Global Moiré Lattice Parameters ---
    // #[getter]
    pub mlv1: Array1<f64>, // Moiré lattice vector 1 (defines the supercell)
    // #[getter]
    pub mlv2: Array1<f64>, // Moiré lattice vector 2
    #[pyo3(get)]
    pub mln1: usize, // Number of Moiré cells repeated along mlv1
    #[pyo3(get)]
    pub mln2: usize, // Number of Moiré cells repeated along mlv2
    #[pyo3(get)]
    pub pbc: bool, // Periodic vs Open Boundary Conditions
    #[pyo3(get)]
    pub study_proximity: i32, // Buffer size for neighbor searching beyond the first shell

    pub toll_scale: f64, // Dynamic tolerance for boundary & coordinate checks (Internal)

    // --- 3. Generated Real-Space Points ---
    // #[getter]
    pub points: Option<Array2<f64>>, // Final (N, 2) array of coordinates for primary lattice
    #[pyo3(get)]
    pub point_types: Option<Vec<usize>>, // Atom type ID for every point in 'points'
    // #[getter]
    pub bigger_points: Option<Array2<f64>>, // Padded set of points for PBC search logic (Internal)
    pub mappings: Option<Vec<usize>>, // Map from 'bigger_points' indices back to primary 'points' (Internal)

    // --- 4. Search & Transformation State ---
    pub kdtree: Option<LatticeTree>, // Spatial index (Kiddo) for fast neighbor queries (Internal)

    // #[getter]
    pub rot_m: Array2<f64>, // Current 2x2 rotation matrix applied to the layer

    // #[getter]
    pub translation: Array1<f64>, // Current [dx, dy] translation applied to the layer
}

#[pymethods]
impl Layer {
    #[new]
    #[pyo3(signature = (lv1, lv2, basis_data, neighbours, pbc, study_proximity))]
    pub fn new(
        lv1: [f64; 2],
        lv2: [f64; 2],
        basis_data: Vec<(f64, f64, String)>,
        neighbours: HashMap<String, Vec<[f64; 2]>>,
        pbc: bool,
        study_proximity: i32,
    ) -> Self {
        let v1 = array![lv1[0], lv1[1]];
        let v2 = array![lv2[0], lv2[1]];

        let norm1 = (v1[0].powi(2) + v1[1].powi(2)).sqrt();
        let norm2 = (v2[0].powi(2) + v2[1].powi(2)).sqrt();
        let toll_scale = f64::max(norm1, norm2);

        let mut basis_coords = Vec::new();
        let mut basis_type_maps: HashMap<String, usize> = HashMap::new();
        let mut basis_types: Vec<String> = Vec::new();

        for (x, y, name) in basis_data {
            // Fail early on duplicate key
            if basis_type_maps.contains_key(&name) {
                panic!("Duplicate key in the line 'lattice_points'");
                // If using PyO3 properly, you'd return a PyErr instead.
            }

            // Safe to insert now
            let id = basis_types.len();
            basis_types.push(name.clone());
            basis_type_maps.insert(name, id);

            basis_coords.push(x);
            basis_coords.push(y);
        }

        let n_basis = basis_coords.len() / 2;
        let neighbours_by_id: Vec<Array2<f64>> = basis_types
            .iter()
            .map(|name| {
                let deltas = neighbours.get(name).cloned().unwrap_or_default();
                let n_neighs = deltas.len();

                // Flatten Vec<[f64; 2]> into a flat Vec<f64> for ndarray consumption
                let flat_deltas: Vec<f64> = deltas.into_iter().flatten().collect();

                Array2::from_shape_vec((n_neighs, 2), flat_deltas)
                    .unwrap_or_else(|_| Array2::zeros((0, 2)))
            })
            .collect();

        Layer {
            lv1: v1,
            lv2: v2,
            basis_points: Array2::from_shape_vec((n_basis, 2), basis_coords).unwrap(),
            basis_type_maps,
            basis_types,
            neighbours: neighbours_by_id,

            mlv1: array![0.0, 0.0],
            mlv2: array![0.0, 0.0],
            mln1: 1,
            mln2: 1,
            pbc,
            study_proximity,
            toll_scale,

            points: None,
            point_types: None,
            bigger_points: None,
            mappings: None,

            kdtree: None,
            rot_m: Array2::eye(2),
            translation: array![0.0, 0.0],
        }
    }

    pub fn perform_rotation_translation(
        &mut self,
        rot: f64,
        translation: [f64; 2],
    ) -> PyResult<()> {
        // Use the Option pattern to check if points already exist
        if self.points.is_some() {
            return Err(pyo3::exceptions::PyAssertionError::new_err(
                "Cannot perform rotation and translation after points have been generated.",
            ));
        }

        let rot_m = crate::utils::get_rotation_matrix(rot);
        let trans_v = array![translation[0], translation[1]];

        // 1. Update internal transformation state
        self.rot_m = rot_m.clone();
        self.translation = trans_v.clone();

        // 2. Rotate Lattice Vectors: lv = R @ lv
        self.lv1 = rot_m.dot(&self.lv1);
        self.lv2 = rot_m.dot(&self.lv2);

        // 3. Rotate and Translate Basis Points
        // Transformation: P_new = (P_old + V) @ R^T
        // We use .dot() with the transpose because basis_points are in rows (N, 2)
        let shifted = &self.basis_points + &trans_v;
        self.basis_points = shifted.dot(&rot_m.t());

        // 4. Rotate Neighbours (Deltas)
        // Neighbours are displacements (deltas), so translation does not apply.
        // Each entry in the Vec is an Array2<f64> representing neighbors for a specific type.
        for type_matrix in self.neighbours.iter_mut() {
            let rotated = type_matrix.dot(&rot_m.t());
            type_matrix.assign(&rotated);
        }

        Ok(())
    }

    #[pyo3(signature = (mlv1, mlv2, mln1=1, mln2=1, test=false))]
    pub fn generate_points(
        &mut self,
        mlv1: [f64; 2],
        mlv2: [f64; 2],
        mln1: usize,
        mln2: usize,
        test: bool,
    ) -> PyResult<()> {
        // Validation logic
        if mln1 == 0 || mln2 == 0 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "mln1 and mln2 must be positive integers.",
            ));
        }
        if test && (mln1 != 1 || mln2 != 1) {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "If test is True, mln1 and mln2 must be 1.",
            ));
        }

        self.mlv1 = array![mlv1[0], mlv1[1]];
        self.mlv2 = array![mlv2[0], mlv2[1]];
        self.mln1 = mln1;
        self.mln2 = mln2;

        // Step 1: Determine grid size 'n' based on maximum Moiré distance
        let p_m1 = &self.mlv1;
        let p_m2 = &self.mlv2;
        let p_m3 = p_m1 + p_m2;

        let max_dist = [
            (p_m1[0].powi(2) + p_m1[1].powi(2)).sqrt(),
            (p_m2[0].powi(2) + p_m2[1].powi(2)).sqrt(),
            (p_m3[0].powi(2) + p_m3[1].powi(2)).sqrt(),
        ]
        .into_iter()
        .fold(0.0, f64::max);

        let norm_lv1 = (self.lv1[0].powi(2) + self.lv1[1].powi(2)).sqrt();
        let norm_lv2 = (self.lv2[0].powi(2) + self.lv2[1].powi(2)).sqrt();
        let min_lv_norm = f64::min(norm_lv1, norm_lv2);

        let n = (max_dist / min_lv_norm).ceil() as i32 * 2;

        // Step 2: Generate points for a single Moiré unit cell
        let mut step1_coords = Vec::new();
        let mut step1_types = Vec::new();

        for i in -n..=n {
            for j in -n..=n {
                let point_o = &self.lv1 * i as f64 + &self.lv2 * j as f64;
                for (k, basis_pt) in self.basis_points.rows().into_iter().enumerate() {
                    let pt = &point_o + &basis_pt;
                    step1_coords.push([pt[0], pt[1]]);
                    step1_types.push(self.point_types.as_ref().map_or(k, |_| 0));
                    // Placeholder or map to type ID
                }
            }
        }

        // Filter step1 points using the boundary check
        let step1_ndarray = Array2::from_shape_vec(
            (step1_coords.len(), 2),
            step1_coords.iter().flatten().cloned().collect(),
        )
        .unwrap();

        let mask = self.inside_boundaries(&step1_ndarray, Some(1), Some(1));

        let mut filtered_step1_coords = Vec::new();
        let mut filtered_step1_types = Vec::new();
        for (idx, &is_inside) in mask.iter().enumerate() {
            if is_inside {
                filtered_step1_coords.push(step1_coords[idx]);
                // Map back to the original type index based on basis_points row order
                filtered_step1_types.push(idx % self.basis_points.nrows());
            }
        }

        // Step 3: Tile the unit cell to create the full lattice
        let mut final_coords = Vec::with_capacity(filtered_step1_coords.len() * mln1 * mln2);
        let mut final_types = Vec::with_capacity(filtered_step1_types.len() * mln1 * mln2);

        for i in 0..mln1 {
            for j in 0..mln2 {
                let translation = &self.mlv1 * i as f64 + &self.mlv2 * j as f64;
                for (coord, &ptype) in filtered_step1_coords.iter().zip(&filtered_step1_types) {
                    final_coords.push([coord[0] + translation[0], coord[1] + translation[1]]);
                    final_types.push(ptype);
                }
            }
        }

        // Update state
        self.points = Some(
            Array2::from_shape_vec(
                (final_coords.len(), 2),
                final_coords.iter().flatten().cloned().collect(),
            )
            .unwrap(),
        );
        self.point_types = Some(final_types);

        if !test {
            // Calling internal Result-based kdtree generation
            self.generate_kdtree()
                .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e))?;
        }

        Ok(())
    }

    pub fn first_nearest_neighbours<'py>(
        &self,
        py: Python<'py>,                    // Explicitly add this to manage lifetimes
        points: PyReadonlyArray2<'py, f64>, // Link points to the same lifetime
        types: Vec<String>,
    ) -> PyResult<(
        Bound<'py, PyArray2<i64>>,
        Bound<'py, PyArray2<i64>>,
        Bound<'py, PyArray2<f64>>,
    )> {
        let pts = points.as_array();

        let type_ids: Vec<usize> = types
            .iter()
            .map(|name| {
                self.basis_type_maps.get(name).cloned().ok_or_else(|| {
                    pyo3::exceptions::PyValueError::new_err(format!("Unknown basis type: {}", name))
                })
            })
            .collect::<Result<Vec<_>, _>>()?;

        let (bigger, smaller, dists) = self
            .first_nearest_neighbours_internal(&pts, &type_ids)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e))?;

        // Use the 'py' token passed in the arguments
        Ok((
            bigger.to_pyarray(py),
            smaller.to_pyarray(py),
            dists.to_pyarray(py),
        ))
    }

    pub fn get_neighbors_within_radius<'py>(
        &self,
        py: Python<'py>,
        query_points: PyReadonlyArray2<'py, f64>,
        radius: f64,
    ) -> PyResult<(
        Bound<'py, PyArray1<usize>>,
        Bound<'py, PyArray1<usize>>,
        Bound<'py, PyArray2<f64>>,
    )> {
        let q_pts = query_points.as_array();

        let (q_idx, l_idx, coords) = self
            .get_neighbors_within_radius_internal(&q_pts, radius)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e))?;

        Ok((
            q_idx.to_pyarray(py),
            l_idx.to_pyarray(py),
            coords.to_pyarray(py),
        ))
    }

    // --- 1D Array Getters (Size 2) ---
    #[getter]
    pub fn lv1<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.lv1.to_pyarray(py) //
    }

    #[getter]
    pub fn lv2<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.lv2.to_pyarray(py) //
    }

    #[getter]
    pub fn mlv1<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mlv1.to_pyarray(py) //
    }

    #[getter]
    pub fn mlv2<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mlv2.to_pyarray(py) //
    }

    #[getter]
    pub fn translation<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.translation.to_pyarray(py) //
    }

    // --- 2D Array Getters ---
    #[getter]
    pub fn basis_points<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.basis_points.to_pyarray(py) //
    }

    #[getter]
    pub fn rot_m<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.rot_m.to_pyarray(py) //
    }

    // --- Optional 2D Array ---
    #[getter]
    pub fn points<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        self.points.as_ref().map(|p| p.to_pyarray(py)) //
    }

    #[getter]
    pub fn bigger_points<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        self.bigger_points.as_ref().map(|p| p.to_pyarray(py)) //
    }

    #[getter]
    pub fn point_types(&self) -> Option<Vec<String>> {
        self.point_types.as_ref().map(|indices| {
            indices
                .iter()
                .map(|&i| self.basis_types[i].clone())
                .collect()
        })
    }

    // --- Complex / Nested Types ---
    #[getter]
    pub fn neighbours<'py>(&self, py: Python<'py>) -> Vec<Bound<'py, PyArray2<f64>>> {
        self.neighbours.iter().map(|n| n.to_pyarray(py)).collect()
    }
}

// Internal Rust-only methods (not exported to Python)
impl Layer {
    pub fn inside_polygon(&self, points: &Array2<f64>, polygon: &Array2<f64>) -> Vec<bool> {
        let n_points = points.nrows();
        let n_poly = polygon.nrows();
        let mut results = Vec::with_capacity(n_points);

        // Pre-extract polygon coordinates for speed
        let px = polygon.column(0);
        let py = polygon.column(1);

        for i in 0..n_points {
            let x = points[[i, 0]];
            let y = points[[i, 1]];
            let mut inside = false;
            let mut j = n_poly - 1;

            for k in 0..n_poly {
                // Ray-casting condition:
                // 1. Is the point's Y between the edge's Y-coordinates?
                // 2. Is the point to the left of the intersection of the ray with the edge?
                if ((py[k] > y) != (py[j] > y))
                    && (x < (px[j] - px[k]) * (y - py[k]) / (py[j] - py[k]) + px[k])
                {
                    inside = !inside;
                }
                j = k;
            }
            results.push(inside);
        }
        results
    }

    pub fn generate_kdtree(&mut self) -> Result<(), String> {
        // 1. Ensure primary points exist
        let points = self
            .points
            .as_ref()
            .ok_or("Cannot generate KDTree: primary points are missing.")?;

        // Handle Open Boundary Conditions (OBC)
        if !self.pbc {
            let points_slice: Vec<[f64; 2]> =
                points.rows().into_iter().map(|r| [r[0], r[1]]).collect();
            self.kdtree = Some(ImmutableKdTree::new_from_slice(&points_slice));
            return Ok(());
        }

        // Handle Periodic Boundary Conditions (PBC)
        // 2. Tiling logic: generate 9 copies of the lattice (-1, 0, 1 shifts)
        let mut all_points = Vec::with_capacity(points.nrows() * 9);
        let mut all_types = Vec::with_capacity(points.nrows() * 9);

        let v1_total = &self.mlv1 * self.mln1 as f64;
        let v2_total = &self.mlv2 * self.mln2 as f64;

        for i in -1..=1 {
            for j in -1..=1 {
                let shift = &v1_total * i as f64 + &v2_total * j as f64;
                for (row, &ptype) in points
                    .rows()
                    .into_iter()
                    .zip(self.point_types.as_ref().unwrap())
                {
                    all_points.push([row[0] + shift[0], row[1] + shift[1]]);
                    all_types.push(ptype);
                }
            }
        }

        // 3. Padding calculation for filtering
        let norm_lv1 = (self.lv1[0].powi(2) + self.lv1[1].powi(2)).sqrt();
        let norm_lv2 = (self.lv2[0].powi(2) + self.lv2[1].powi(2)).sqrt();
        let norm_v1 = (v1_total[0].powi(2) + v1_total[1].powi(2)).sqrt();
        let norm_v2 = (v2_total[0].powi(2) + v2_total[1].powi(2)).sqrt();

        let neigh_pad_1 = (1.0 + self.study_proximity as f64) * norm_lv1 / norm_v1;
        let neigh_pad_2 = (1.0 + self.study_proximity as f64) * norm_lv2 / norm_v2;

        // 4. Create the padding polygon
        let p1 = &v1_total * -neigh_pad_1 + &v2_total * -neigh_pad_2;
        let p2 = &v1_total * (1.0 + neigh_pad_1) + &v2_total * -neigh_pad_2;
        let p3 = &v1_total * (1.0 + neigh_pad_1) + &v2_total * (1.0 + neigh_pad_2);
        let p4 = &v1_total * -neigh_pad_1 + &v2_total * (1.0 + neigh_pad_2);

        let pad_polygon = array![
            [p1[0], p1[1]],
            [p2[0], p2[1]],
            [p3[0], p3[1]],
            [p4[0], p4[1]]
        ];

        // 5. Filter points inside the padded boundary
        let all_points_ndarray = Array2::from_shape_vec(
            (all_points.len(), 2),
            all_points.iter().flatten().cloned().collect(),
        )
        .unwrap();

        let mask = self.inside_polygon(&all_points_ndarray, &pad_polygon);

        // Collect filtered results
        let mut filtered_points = Vec::new();
        let mut filtered_types = Vec::new();

        for (idx, is_inside) in mask.into_iter().enumerate() {
            if is_inside {
                filtered_points.push(all_points[idx]);
                filtered_types.push(all_types[idx]);
            }
        }

        // 6. Update internal state
        // We need to convert filtered_points back to an Array2 for storage
        let n_filtered = filtered_points.len();
        self.bigger_points = Some(
            Array2::from_shape_vec(
                (n_filtered, 2),
                filtered_points.iter().flatten().cloned().collect(),
            )
            .unwrap(),
        );

        // Build the final KDTree from the "bigger" slice
        self.kdtree = Some(ImmutableKdTree::new_from_slice(&filtered_points));

        // 7. Finally, generate the mapping between Bigger -> Original
        self.generate_mapping()
    }

    /// Calculates Moiré cell boundaries and filters points
    pub fn inside_boundaries(
        &self,
        points: &Array2<f64>,
        mln1: Option<usize>,
        mln2: Option<usize>,
    ) -> Vec<bool> {
        let n1 = mln1.unwrap_or(self.mln1) as f64;
        let n2 = mln2.unwrap_or(self.mln2) as f64;

        let v1 = &self.mlv1 * n1;
        let v2 = &self.mlv2 * n2;

        // Define the 4 corners of the parallelogram
        let p1 = array![0.0, 0.0];
        let p2 = &v1;
        let p3 = &v1 + &v2;
        let p4 = &v2;

        // Create polygon array (N, 2)
        let polygon = array![
            [p1[0], p1[1]],
            [p2[0], p2[1]],
            [p3[0], p3[1]],
            [p4[0], p4[1]]
        ];

        // Shift the polygon slightly inward to avoid edge-case "ghost" points
        // that belong to the neighboring cell (matching your Python logic)
        let shift_dir = -(&v1 + &v2);
        let shift_mag = (shift_dir[0].powi(2) + shift_dir[1].powi(2)).sqrt();
        let shift = if shift_mag > 0.0 {
            (shift_dir / shift_mag) * self.toll_scale * 1e-4
        } else {
            array![0.0, 0.0]
        };

        let mut shifted_poly = polygon;
        for mut row in shifted_poly.rows_mut() {
            row += &shift;
        }

        self.inside_polygon(points, &shifted_poly)
    }

    fn generate_mapping(&mut self) -> Result<(), String> {
        let points = self.points.as_ref().ok_or("Points not generated")?;
        let bigger_points = self
            .bigger_points
            .as_ref()
            .ok_or("Bigger points not generated")?;

        let points_slice: Vec<[f64; 2]> = points.rows().into_iter().map(|r| [r[0], r[1]]).collect();

        let tree: LatticeTree = ImmutableKdTree::new_from_slice(&points_slice);

        let v1_total = &self.mlv1 * self.mln1 as f64;
        let v2_total = &self.mlv2 * self.mln2 as f64;
        let translations = self.get_point_positions(bigger_points, &v1_total, &v2_total);

        let mut mapped_coords = bigger_points.clone();
        for (mut coord, trans) in mapped_coords
            .rows_mut()
            .into_iter()
            .zip(translations.rows().into_iter())
        {
            coord[0] -= trans[0] * v1_total[0] + trans[1] * v2_total[0];
            coord[1] -= trans[0] * v1_total[1] + trans[1] * v2_total[1];
        }

        let mut mapping_vec = Vec::with_capacity(bigger_points.nrows());
        let tol_sq = (self.toll_scale * 1e-3).powi(2);

        for (i, coord) in mapped_coords.rows().into_iter().enumerate() {
            let query_point = [coord[0], coord[1]];

            // Explicitly specify the metric to resolve the ambiguity
            let nearest = tree.nearest_one::<kiddo::SquaredEuclidean>(&query_point);

            if nearest.distance > tol_sq {
                return Err(format!(
                    "Distance {} exceeds tolerance for point index {}",
                    nearest.distance.sqrt(),
                    i
                ));
            }

            mapping_vec.push(nearest.item as usize);
        }

        self.mappings = Some(mapping_vec);
        Ok(())
    }

    fn get_point_positions(
        &self,
        points: &Array2<f64>,
        a: &Array1<f64>,
        b: &Array1<f64>,
    ) -> Array2<f64> {
        let n = points.nrows();
        let mut results = Array2::<f64>::zeros((n, 2));
        let tol = self.toll_scale * 1e-2;

        // Extract columns for easier reading (no-copy views)
        let x = points.column(0);
        let y = points.column(1);

        // Vectorized determinants for Y-position (relative to vectors A and BC)
        // det = x * Ay - y * Ax
        let det_oa = &x * a[1] - &y * a[0];
        let det_bc = (&x - b[0]) * a[1] - (&y - b[1]) * a[0];

        // Vectorized determinants for X-position (relative to vectors B and AC)
        // det = x * By - y * Bx
        let det_ob = &x * b[1] - &y * b[0];
        let det_ac = (&x - a[0]) * b[1] - (&y - a[1]) * b[0];

        // Zip through the results and assign based on the tolerance checks
        // We use azip! here for the best balance of readability and speed.
        ndarray::azip!((
            mut res in results.rows_mut(),
            &d_oa in &det_oa,
            &d_bc in &det_bc,
            &d_ob in &det_ob,
            &d_ac in &det_ac
        ) {
            let pos_y = if d_oa <= tol { 1.0 } else { 0.0 } + if d_bc <= tol { 1.0 } else { 0.0 };
            let pos_x = if d_ob > -tol { 1.0 } else { 0.0 } + if d_ac > -tol { 1.0 } else { 0.0 };

            res[0] = pos_x - 1.0;
            res[1] = pos_y - 1.0;
        });

        results
    }

    pub fn first_nearest_neighbours_internal(
        &self,
        points: &ArrayView2<f64>,
        types: &[usize],
    ) -> Result<(Array2<i64>, Array2<i64>, Array2<f64>), String> {
        let tree = self.kdtree.as_ref().ok_or("KDTree not initialized.")?;
        let n_points = points.nrows();

        // Find unique types in the input to determine the width of the output arrays
        let mut unique_types: Vec<usize> = types.to_vec();
        unique_types.sort_unstable();
        unique_types.dedup();

        let max_neighs = unique_types
            .iter()
            .map(|&t| self.neighbours[t].nrows())
            .max()
            .unwrap_or(0);

        // Initialize arrays with default values (-1 for indices, NaN for distances)
        let mut bigger_indices = Array2::<i64>::from_elem((n_points, max_neighs), -1);
        let mut smaller_indices = Array2::<i64>::from_elem((n_points, max_neighs), -1);
        let mut distances_arr = Array2::<f64>::from_elem((n_points, max_neighs), f64::NAN);

        let tol_sq = (self.toll_scale * 1e-2).powi(2);

        // Iterate through unique types to perform bulk-like neighbor searches
        for &t_id in &unique_types {
            let rel_neighs = &self.neighbours[t_id];

            for (i, &current_type_id) in types.iter().enumerate() {
                if current_type_id != t_id {
                    continue;
                }

                let point = points.row(i);

                for (n_idx, delta) in rel_neighs.rows().into_iter().enumerate() {
                    let target = [point[0] + delta[0], point[1] + delta[1]];

                    // KDTree query using the squared distance metric
                    let nearest = tree.nearest_one::<kiddo::SquaredEuclidean>(&target);

                    if nearest.distance > tol_sq {
                        return Err(format!(
                            "Distance {} exceeds tolerance for type {} at point {}",
                            nearest.distance.sqrt(),
                            self.basis_types[t_id],
                            i
                        ));
                    }

                    let b_idx = nearest.item as i64;
                    bigger_indices[[i, n_idx]] = b_idx;
                    distances_arr[[i, n_idx]] = nearest.distance.sqrt();

                    if let Some(ref map) = self.mappings {
                        smaller_indices[[i, n_idx]] = map[b_idx as usize] as i64;
                    } else {
                        smaller_indices[[i, n_idx]] = b_idx;
                    }
                }
            }
        }

        Ok((bigger_indices, smaller_indices, distances_arr))
    }

    pub fn get_neighbors_within_radius_internal(
        &self,
        query_points: &ArrayView2<f64>,
        radius: f64,
    ) -> Result<(Array1<usize>, Array1<usize>, Array2<f64>), String> {
        let tree = self.kdtree.as_ref().ok_or("KDTree not initialized.")?;
        let radius_sq = radius.powi(2); // Kiddo uses squared distance

        let mut q_indices = Vec::new();
        let mut l_indices = Vec::new();
        let mut coords_vec = Vec::new();

        // Use bigger_points for PBC, otherwise use primary points
        let source_points = if self.pbc {
            self.bigger_points
                .as_ref()
                .ok_or("Bigger points missing for PBC.")?
        } else {
            self.points.as_ref().ok_or("Points not generated.")?
        };

        for (i, q_pt) in query_points.rows().into_iter().enumerate() {
            let matches = tree.within::<kiddo::SquaredEuclidean>(&[q_pt[0], q_pt[1]], radius_sq);

            for m in matches {
                let b_idx = m.item as usize;
                q_indices.push(i);

                // Map to primary unit cell if PBC is active
                if let Some(ref maps) = self.mappings {
                    l_indices.push(maps[b_idx]);
                } else {
                    l_indices.push(b_idx);
                }

                let row = source_points.row(b_idx);
                coords_vec.push(row[0]);
                coords_vec.push(row[1]);
            }
        }

        let n_found = q_indices.len();
        let coords_arr = Array2::from_shape_vec((n_found, 2), coords_vec)
            .map_err(|_| "Failed to build coordinate array.")?;

        Ok((Array1::from(q_indices), Array1::from(l_indices), coords_arr))
    }
}
