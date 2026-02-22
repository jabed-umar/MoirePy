use ndarray::{array, Array1, Array2, Axis};
use kiddo::immutable::float::kdtree::ImmutableKdTree;
use pyo3::prelude::*;
use numpy::{PyArray1, PyArray2, PyReadonlyArray1, ToPyArray};
use std::collections::HashMap;

// A: f64 (float type), T: u64 (ID type), K: 2 (Dimensions), B: 32 (Bucket size)
type LatticeTree = ImmutableKdTree<f64, u64, 2, 32>;

#[pyclass]
pub struct Layer {
    // --- 1. Fundamental Unit Cell Geometry ---
    pub lv1: Array1<f64>,             // Primary lattice vector 1 (rotated if applicable)
    pub lv2: Array1<f64>,             // Primary lattice vector 2 (rotated if applicable)
    pub basis_points: Array2<f64>,    // Positions of atoms within the unit cell
    #[pyo3(get)]
    pub basis_types: Vec<usize>,      // Integer mapping of atom names (e.g., 0 for 'A', 1 for 'B')
    #[pyo3(get)]
    pub point_type_names: Vec<String>, // The original string names (e.g., ["A", "B"])
    #[pyo3(get)]
    pub neighbours: HashMap<String, Vec<[f64; 2]>>, // Relative displacements for hopping logic

    // --- 2. Global Moiré Lattice Parameters ---
    pub mlv1: Array1<f64>,            // Moiré lattice vector 1 (defines the supercell)
    pub mlv2: Array1<f64>,            // Moiré lattice vector 2
    #[pyo3(get)]
    pub mln1: usize,                  // Number of Moiré cells repeated along mlv1
    #[pyo3(get)]
    pub mln2: usize,                  // Number of Moiré cells repeated along mlv2
    #[pyo3(get)]
    pub pbc: bool,                    // Periodic vs Open Boundary Conditions
    #[pyo3(get)]
    pub study_proximity: i32,         // Buffer size for neighbor searching beyond the first shell
    #[pyo3(get)]
    pub toll_scale: f64,              // Dynamic tolerance for boundary & coordinate checks

    // --- 3. Generated Real-Space Points ---
    pub points: Array2<f64>,          // The final (N, 2) array of coordinates for the primary lattice
    #[pyo3(get)]
    pub point_types: Vec<usize>,      // Atom type ID for every point in 'points'
    pub bigger_points: Option<Array2<f64>>, // Padded set of points for PBC search logic
    #[pyo3(get)]
    pub mappings: Vec<usize>,         // Map from 'bigger_points' indices back to primary 'points'

    // --- 4. Search & Transformation State ---
    pub kdtree: Option<LatticeTree>,  // Spatial index (Kiddo) for fast neighbor queries
    pub rot_m: Array2<f64>,           // Current 2x2 rotation matrix applied to the layer
    pub translation: Array1<f64>,     // Current [dx, dy] translation applied to the layer
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
        
        // Calculate toll_scale based on the norm of lattice vectors
        let norm1 = (v1[0].powi(2) + v1[1].powi(2)).sqrt();
        let norm2 = (v2[0].powi(2) + v2[1].powi(2)).sqrt();
        let toll_scale = f64::max(norm1, norm2);

        let mut basis_coords = Vec::new();
        let mut basis_types = Vec::new();
        let mut point_type_names = Vec::new();

        for (x, y, name) in basis_data {
            basis_coords.push(x);
            basis_coords.push(y);
            
            // Map string types (A, B...) to integer indices for performance
            let type_id = point_type_names.iter().position(|r| r == &name)
                .unwrap_or_else(|| {
                    point_type_names.push(name);
                    point_type_names.len() - 1
                });
            basis_types.push(type_id);
        }

        Layer {
            // Unit Cell Geometry
            lv1: v1,
            lv2: v2,
            basis_points: Array2::from_shape_vec((basis_types.len(), 2), basis_coords).unwrap(),
            basis_types,
            point_type_names,
            neighbours,

            // Moiré Parameters
            mlv1: array![0.0, 0.0],
            mlv2: array![0.0, 0.0],
            mln1: 1,
            mln2: 1,
            pbc,
            study_proximity,
            toll_scale,

            // Generated Data
            points: Array2::default((0, 2)),
            point_types: Vec::new(),
            bigger_points: None,
            mappings: Vec::new(),

            // Search & Transformations
            kdtree: None,
            rot_m: Array2::eye(2),
            translation: array![0.0, 0.0],
        }
    }


    pub fn perform_rotation_translation(&mut self, rot: f64, translation: PyReadonlyArray1<f64>) {
        // 1. Calculate and store the rotation matrix
        let rot_m = crate::utils::get_rotation_matrix(rot);
        self.rot_m = rot_m.clone();

        // 2. Convert the input NumPy array view to an ndarray for math
        let trans_view = translation.as_array();
        let trans_v = array![trans_view[0], trans_view[1]];
        self.translation = trans_v.clone();

        // 3. Rotate Lattice Vectors
        self.lv1 = self.rot_m.dot(&self.lv1);
        self.lv2 = self.rot_m.dot(&self.lv2);

        // 4. Rotate and Translate Basis Points
        // New position = rot_m @ (old_pos + translation)
        for mut row in self.basis_points.axis_iter_mut(Axis(0)) {
            let old_pos = array![row[0], row[1]];
            let new_pos = self.rot_m.dot(&(old_pos + &trans_v));
            row[0] = new_pos[0];
            row[1] = new_pos[1];
        }

        // 5. Rotate Neighbours Dictionary
        for neighbour_list in self.neighbours.values_mut() {
            for n in neighbour_list.iter_mut() {
                let n_vec = array![n[0], n[1]];
                let rotated_n = self.rot_m.dot(&n_vec);
                n[0] = rotated_n[0];
                n[1] = rotated_n[1];
            }
        }
    }

    pub fn generate_points(
        &mut self,
        mlv1: PyReadonlyArray1<f64>,
        mlv2: PyReadonlyArray1<f64>,
        mln1: usize,
        mln2: usize,
    ) {
        let v1_view = mlv1.as_array();
        let v2_view = mlv2.as_array();
        self.mlv1 = array![v1_view[0], v1_view[1]];
        self.mlv2 = array![v2_view[0], v2_view[1]];

        // Step 1: Find the maximum distance to determine grid range
        let moire_corners = [
            array![0.0, 0.0],
            self.mlv1.clone(),
            self.mlv2.clone(),
            &self.mlv1 + &self.mlv2,
        ];

        let max_distance = moire_corners
            .iter()
            .map(|p| (p[0].powi(2) + p[1].powi(2)).sqrt())
            .fold(0.0, f64::max);

        let min_lv_norm = f64::min(
            (self.lv1[0].powi(2) + self.lv1[1].powi(2)).sqrt(),
            (self.lv2[0].powi(2) + self.lv2[1].powi(2)).sqrt(),
        );

        let n = (max_distance / min_lv_norm).ceil() as i32 * 2;

        // Step 2: Generate points inside one moire unit cell
        let mut step1_coords = Vec::new();
        let mut step1_types = Vec::new();

        for i in -n..=n {
            for j in -n..=n {
                // Calculate the origin point for this (i, j) cell
                let point_o = &self.lv1 * i as f64 + &self.lv2 * j as f64;
                for idx in 0..self.basis_types.len() {
                    let basis_pt = self.basis_points.row(idx);
                    let point = array![point_o[0] + basis_pt[0], point_o[1] + basis_pt[1]];
                    
                    // Filter points using the boundary check
                    if self.check_boundaries(&[point[0], point[1]], 1, 1) {
                        step1_coords.push(point[0]);
                        step1_coords.push(point[1]);
                        step1_types.push(self.basis_types[idx]);
                    }
                }
            }
        }

        // Step 3: Copy and translate the unit cell to create the full lattice
        let mut final_coords = Vec::new();
        let mut final_types = Vec::new();

        for i in 0..mln1 {
            for j in 0..mln2 {
                let translation = &self.mlv1 * i as f64 + &self.mlv2 * j as f64;
                for p_idx in 0..step1_types.len() {
                    final_coords.push(step1_coords[p_idx * 2] + translation[0]);
                    final_coords.push(step1_coords[p_idx * 2 + 1] + translation[1]);
                    final_types.push(step1_types[p_idx]);
                }
            }
        }

        self.points = Array2::from_shape_vec((final_types.len(), 2), final_coords).unwrap();
        self.point_types = final_types;
        
        self.generate_kdtree();
    }

    // // Port of first_nearest_neighbours
    // // PLAN: Take query points, run bulk kiddo queries, resolve via self.mappings if PBC.
    // pub fn first_nearest_neighbours(&self, query_points: Array2<f64>) -> (Vec<usize>, Vec<f64>) { 
    //     todo!() 
    // }

    // Explicit Getters to expose ndarray as Read-Only NumPy views in Python
    #[getter]
    fn lv1<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.lv1.to_pyarray(py)
    }

    #[getter]
    fn lv2<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.lv2.to_pyarray(py)
    }

    #[getter]
    fn mlv1<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mlv1.to_pyarray(py)
    }

    #[getter]
    fn mlv2<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mlv2.to_pyarray(py)
    }

    #[getter]
    fn points<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.points.to_pyarray(py)
    }

    #[getter]
    fn basis_points<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.basis_points.to_pyarray(py)
    }

    #[getter]
    fn rot_m<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.rot_m.to_pyarray(py)
    }

    #[getter]
    fn translation<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.translation.to_pyarray(py)
    }

    #[getter]
    fn point_types(&self) -> Vec<usize> {
        self.point_types.clone()
    }

    #[getter]
    fn mappings(&self) -> Vec<usize> {
        self.mappings.clone()
    }

    #[getter]
    fn bigger_points<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        self.bigger_points.as_ref().map(|p| p.to_pyarray(py))
    }

}

// Internal Rust-only methods (not exported to Python)
impl Layer {
    fn generate_kdtree(&mut self) {
        if !self.pbc {
            // OBC: Just build the tree from the primary points
            let tree_points: Vec<[f64; 2]> = self.points
                .axis_iter(ndarray::Axis(0))
                .map(|row| [row[0], row[1]])
                .collect();
            
            self.kdtree = Some(LatticeTree::new_from_slice(&tree_points));
            return;
        }

        // PBC: Generate a larger set of points (3x3 grid of moire cells)
        let mut all_points_coords = Vec::new();
        let mut all_point_types = Vec::new();

        let total_mlv1 = &self.mlv1 * self.mln1 as f64;
        let total_mlv2 = &self.mlv2 * self.mln2 as f64;

        for i in -1..=1 {
            for j in -1..=1 {
                let shift = &total_mlv1 * i as f64 + &total_mlv2 * j as f64;
                for (idx, row) in self.points.axis_iter(ndarray::Axis(0)).enumerate() {
                    all_points_coords.push(row[0] + shift[0]);
                    all_points_coords.push(row[1] + shift[1]);
                    all_point_types.push(self.point_types[idx]);
                }
            }
        }

        // Calculate padding for the boundary mask
        let norm_v1 = (total_mlv1[0].powi(2) + total_mlv1[1].powi(2)).sqrt();
        let norm_v2 = (total_mlv2[0].powi(2) + total_mlv2[1].powi(2)).sqrt();
        let norm_lv1 = (self.lv1[0].powi(2) + self.lv1[1].powi(2)).sqrt();
        let norm_lv2 = (self.lv2[0].powi(2) + self.lv2[1].powi(2)).sqrt();

        let neigh_pad_1 = (1.0 + self.study_proximity as f64) * norm_lv1 / norm_v1;
        let neigh_pad_2 = (1.0 + self.study_proximity as f64) * norm_lv2 / norm_v2;

        // Construct the "bigger" polygon for filtering
        let p_corners = [
            &total_mlv1 * (-neigh_pad_1) + &total_mlv2 * (-neigh_pad_2),
            &total_mlv1 * (1.0 + neigh_pad_1) + &total_mlv2 * (-neigh_pad_2),
            &total_mlv1 * (1.0 + neigh_pad_1) + &total_mlv2 * (1.0 + neigh_pad_2),
            &total_mlv1 * (-neigh_pad_1) + &total_mlv2 * (1.0 + neigh_pad_2),
        ];

        let polygon = [
            [p_corners[0][0], p_corners[0][1]],
            [p_corners[1][0], p_corners[1][1]],
            [p_corners[2][0], p_corners[2][1]],
            [p_corners[3][0], p_corners[3][1]],
        ];

        // Filter points using the polygon
        let mut bigger_coords = Vec::new();
        let mut bigger_types = Vec::new();

        for (idx, chunk) in all_points_coords.chunks(2).enumerate() {
            let pt = [chunk[0], chunk[1]];
            if self.is_inside_polygon(&pt, &polygon) {
                bigger_coords.push(pt[0]);
                bigger_coords.push(pt[1]);
                bigger_types.push(all_point_types[idx]);
            }
        }

        // Store the bigger points and build the final tree
        let n_bigger = bigger_types.len();
        let bigger_points_arr = Array2::from_shape_vec((n_bigger, 2), bigger_coords).unwrap();
        
        let tree_points: Vec<[f64; 2]> = bigger_points_arr
            .axis_iter(ndarray::Axis(0))
            .map(|row| [row[0], row[1]])
            .collect();

        self.bigger_points = Some(bigger_points_arr);
        // Note: You'll need to add 'bigger_point_types' to your struct if you want to keep them
        self.kdtree = Some(LatticeTree::new_from_slice(&tree_points));

        // Generate the index mapping
        self.generate_mapping();
    }


    // Port of _inside_polygon
    // PLAN: Implementation of ray-casting algorithm for a single point vs polygon.
    fn is_inside_polygon(&self, point: &[f64; 2], polygon: &[[f64; 2]]) -> bool {
        let (x, y) = (point[0], point[1]);
        let mut inside = false;
        let n = polygon.len();

        for i in 0..n {
            let j = (i + 1) % n;
            let (pi_x, pi_y) = (polygon[i][0], polygon[i][1]);
            let (pj_x, pj_y) = (polygon[j][0], polygon[j][1]);

            // Ray-casting condition
            let intersect = ((pi_y > y) != (pj_y > y)) &&
                (x < (pj_x - pi_x) * (y - pi_y) / (pj_y - pi_y) + pi_x);
            
            if intersect {
                inside = !inside;
            }
        }
        inside
    }

    // Port of _inside_boundaries
    // PLAN: Construct the polygon vertices based on mlv and mln, then call is_inside_polygon.
    fn check_boundaries(&self, point: &[f64; 2], mln1: usize, mln2: usize) -> bool {
        let v1 = &self.mlv1 * mln1 as f64;
        let v2 = &self.mlv2 * mln2 as f64;

        let p1 = [0.0, 0.0];
        let p2 = [v1[0], v1[1]];
        let p3 = [v2[0], v2[1]];
        let p4 = [v1[0] + v2[0], v1[1] + v2[1]];

        // Calculate the shift (same as Python logic)
        let shift_dir = -(&v1 + &v2);
        let norm = (shift_dir[0].powi(2) + shift_dir[1].powi(2)).sqrt();
        let shift = if norm > 0.0 {
            let s = self.toll_scale * 1e-4 / norm;
            [shift_dir[0] * s, shift_dir[1] * s]
        } else {
            [0.0, 0.0]
        };

        // Define polygon vertices (p1, p2, p4, p3 order for the parallelogram)
        let polygon = [
            [p1[0] + shift[0], p1[1] + shift[1]],
            [p2[0] + shift[0], p2[1] + shift[1]],
            [p4[0] + shift[0], p4[1] + shift[1]],
            [p3[0] + shift[0], p3[1] + shift[1]],
        ];

        self.is_inside_polygon(point, &polygon)
    }

    // Port of _generate_mapping
    // PLAN: For each 'bigger_point', find its primary counterpart using point_positions math.
    fn generate_mapping(&mut self) {
        // 1. Create a temporary KD-Tree of the small (primary) point set
        // We use this to find which primary point corresponds to a translated "bigger" point.
        let tree_points: Vec<[f64; 2]> = self.points
            .axis_iter(ndarray::Axis(0))
            .map(|row| [row[0], row[1]])
            .collect();
        let tree = LatticeTree::new_from_slice(&tree_points);

        let bigger_pts = self.bigger_points.as_ref().unwrap();
        let mut mappings = Vec::with_capacity(bigger_pts.nrows());

        let total_mlv1 = &self.mlv1 * self.mln1 as f64;
        let total_mlv2 = &self.mlv2 * self.mln2 as f64;

        // point positions... for each point in self.point, point position is a array of length 2 (x, y)
        // where the elements are -1, 0 and 1... this is what their value mean about their position
        // (-1, 1) | (0, 1) | (1, 1)
        // -----------------------------
        // (-1, 0) | (0, 0) | (1, 0)
        // -----------------------------
        // (-1,-1) | (0,-1) | (1,-1)
        // -----------------------------
        // (0, 0) is our actual lattice part...
        for row in bigger_pts.axis_iter(ndarray::Axis(0)) {
            let pt = [row[0], row[1]];
            
            // Get (dx, dy) which tells us which neighboring moire cell the point is in
            let (dx, dy) = self.get_point_positions(&pt, &total_mlv1, &total_mlv2);

            // all point with point_positions = (x, y) need to be translated by
            // (-x*self.mlv1*self.mln1 - y*self.mlv2*self.mln2) to get the corresponding point inside the lattice
            let mapped_pt = [
                pt[0] - (dx as f64 * total_mlv1[0] + dy as f64 * total_mlv2[0]),
                pt[1] - (dx as f64 * total_mlv1[1] + dy as f64 * total_mlv2[1]),
            ];

            // query on the kdtree of the smaller points to get the index 
            // of the corresponding point inside the lattice
            let neighbor = tree.nearest_one::<kiddo::SquaredEuclidean>(&mapped_pt);
            
            if neighbor.distance >= (self.toll_scale * 1e-3).powi(2) {
                panic!(
                    "FATAL ERROR: Distance {} exceeds tolerance for point {:?} mapped at {:?}",
                    neighbor.distance.sqrt(), pt, mapped_pt
                );
            }

            // store that in `self.mappings`: keys are indices in bigger_points, values are indices in points
            mappings.push(neighbor.item as usize);
        }

        self.mappings = mappings;
    }

    // Port of _point_positions
    // PLAN: Use Cramer's rule/determinants to find which "cell" a point belongs to.
    fn get_point_positions(&self, point: &[f64; 2], a: &Array1<f64>, b: &Array1<f64>) -> (isize, isize) {
        let toll = self.toll_scale * 1e-2;

        // Positions relative to OA and BC (Y direction)
        let det_oa = (point[0] * a[1] - point[1] * a[0]) <= toll;
        let det_bc = ((point[0] - b[0]) * a[1] - (point[1] - b[1]) * a[0]) <= toll;
        let pos_y = (det_oa as isize) + (det_bc as isize) - 1;

        // Positions relative to OB and AC (X direction)
        let det_ob = (point[0] * b[1] - point[1] * b[0]) > -toll;
        let det_ac = ((point[0] - a[0]) * b[1] - (point[1] - a[1]) * b[0]) > -toll;
        let pos_x = (det_ob as isize) + (det_ac as isize) - 1;

        (pos_x, pos_y)
    }
}