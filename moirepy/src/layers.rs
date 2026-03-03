use kiddo::immutable::float::kdtree::ImmutableKdTree;
use ndarray::{array, Array1, Array2, ArrayView2};
use numpy::{PyArray1, PyArray2, PyReadonlyArray2, ToPyArray};
use pyo3::prelude::*;
use std::collections::HashMap;

/// Concrete type alias for the immutable KD-tree used throughout this module.
///
/// Parameters: `f64` coordinate type, `u64` item index (cast to `usize` on use),
/// `2` spatial dimensions, `32` bucket size (leaf capacity before splitting).
///
/// Using a type alias keeps the rest of the code readable and makes it trivial
/// to change these parameters in one place if needed.
type LatticeTree = ImmutableKdTree<f64, u64, 2, 32>;

/// Represents a simple 2D crystal layer, the fundamental building block of a Moiré system.
///
/// A `Layer` encapsulates everything needed to describe one crystal plane:
/// its geometry (lattice vectors + basis), the generated real-space point
/// cloud, nearest-neighbour connectivity, and the spatial index (KD-tree)
/// used for fast radius and nearest neighbour (NN) searches.
///
/// # Lifecycle
/// 1. Construct with [`Layer::new`]: sets geometry, does NOT generate points.
/// 2. (Optional) Rotate/translate with [`Layer::perform_rotation_translation`].
/// 3. Generate points with [`Layer::generate_points`]: populates `points`,
///    `point_types`, and (if `pbc = true`) `bigger_points`, `mappings`, `kdtree`.
/// 4. Query with [`Layer::first_nearest_neighbours`] or
///    [`Layer::get_neighbors_within_radius`].
///
/// Fields marked `#[pyo3(get)]` are directly readable from Python. Other fields
/// either have hand-written `#[getter]` methods (ndarray -> NumPy conversion)
/// or are private to Rust.
#[pyclass]
pub struct Layer {
    // -------------------------------------------------------------------------
    // 1. Fundamental Unit Cell Geometry
    // These define the crystal's intrinsic structure before any replication.
    // -------------------------------------------------------------------------

    /// Primary lattice vector **a1** in Cartesian coordinates `[x, y]`.
    ///
    /// By convention `lv1` must lie along the x-axis (`lv1[1] == 0`).
    /// Things may (or may not) break if this is not followed.
    pub lv1: Array1<f64>,  // #[getter]

    /// Primary lattice vector **a2** in Cartesian coordinates `[x, y]`.
    ///
    /// Must have a positive y-component (`lv2[1] > 0`).
    pub lv2: Array1<f64>,  // #[getter]

    /// Positions of all atoms inside one unit cell, shape `(B, 2)`.
    ///
    /// ```text
    /// basis_points = [
    ///     [x_A, y_A],   // basis_types[0] == "A"
    ///     [x_B, y_B],   // basis_types[1] == "B"
    ///     ...
    /// ]
    /// ```
    ///
    /// Updated in place when the layer is rotated/translated.
    pub basis_points: Array2<f64>,  // #[getter]

    /// Name-to-ID mapping for basis types. Inverse of `basis_types`. Private.
    ///
    /// ```text
    /// basis_type_maps = { "A" => 0, "B" => 1, ... }
    /// ```
    pub basis_type_maps: HashMap<String, usize>,    // private

    /// Ordered list of basis-type names. Used as reverse search of `basis_type_maps`.
    ///
    /// ```text
    /// basis_types = ["A", "B", ...]
    /// or u can think of it like:
    /// basis_types = [0 => "A", 1 => "B", ...]
    /// ```
    #[pyo3(get)]
    pub basis_types: Vec<String>,

    /// First-neighbour displacement vectors per basis type. `Vec` index == type-ID.
    ///
    /// ```text
    /// neighbours[k] =            // shape (N_k, 2), N_k neighbours of type-k atom
    ///     [[dx_0, dy_0],         // displacement to neighbour 0
    ///      [dx_1, dy_1],         // displacement to neighbour 1
    ///      ...]                  // N_k may differ per type (e.g. honeycomb: 3, Kagome: 4)
    /// ```
    ///
    /// Rotated in place by [`Layer::perform_rotation_translation`].
    pub neighbours: Vec<Array2<f64>>,  // #[getter]

    // -------------------------------------------------------------------------
    // 2. Global Moiré / Supercell Parameters
    // Set during generate_points(); describe the full tiled supercell.
    // -------------------------------------------------------------------------

    /// First Moiré supercell lattice vector **M1**, shape `[x, y]`.
    ///
    /// Defines the supercell boundary in one direction. Supplied by the caller
    /// in [`Layer::generate_points`] and stored for PBC tiling and boundary polygon.
    pub mlv1: Array1<f64>,  // #[getter]

    /// Second Moiré supercell lattice vector **M2**, shape `[x, y]`.
    pub mlv2: Array1<f64>,  // #[getter]

    /// Number of Moiré unit cells tiled along M1.
    /// The final supercell spans `mln1 * mlv1` in that direction.
    #[pyo3(get)]
    pub mln1: usize,

    /// Number of Moiré unit cells tiled along M2.
    #[pyo3(get)]
    pub mln2: usize,

    /// Whether this layer uses Periodic (`true`) or Open (`false`) Boundary Conditions.
    ///
    /// PBC: the KD-tree is built from `bigger_points` (a 3x3-tiled supercell padded
    /// to `study_proximity` shells); every matched index is remapped via `mappings`.
    /// OBC: the KD-tree is built directly from `points`; no remapping.
    #[pyo3(get)]
    pub pbc: bool,

    /// Number of additional neighbour shells to include in the PBC padding.
    ///
    /// A value of `1` means atleast the first shell of neighbours outside the
    /// primary cell is kept in `bigger_points`. Increase for systems where
    /// you need to search further than one lattice constant beyond the
    /// supercell boundary.
    #[pyo3(get)]
    pub study_proximity: i32,

    /// Dynamic tolerance scale: `toll_scale = max(|lv1|, |lv2|)`.
    ///
    /// Used as the base for all length-scale-dependent epsilons:
    /// boundary polygon inward-shift (`* 1e-4`), KD-tree match check (`* 1e-2`)^2,
    /// and mapping distance check (`* 1e-3`)^2.
    /// Keeping it proportional to the lattice constant makes the code
    /// lattice-agnostic (same thresholds work for angstrom-scale and nm-scale systems).
    pub toll_scale: f64,  // private

    // -------------------------------------------------------------------------
    // 3. Generated Real-Space Points
    // Populated by generate_points(); None before that call.
    // -------------------------------------------------------------------------

    /// All atom positions in the primary supercell, shape `(N, 2)` where N is the total atom count.
    /// Elements have a one to one correspondence with elements of point_types.
    ///
    /// ```text
    /// points = [[x_0, y_0],
    ///           [x_1, y_1],
    ///           ...           // N rows total
    ///           [x_N, y_N]]
    /// ```
    ///
    /// `None` until [`Layer::generate_points`] succeeds.
    pub points: Option<Array2<f64>>,  // #[getter]

    /// Basis-type ID for every point in `points`, length N.
    ///
    /// ```text
    /// point_types = [0, 1, 0, 1, 1, 1, 0 ...]   // point_types[k] indexes into basis_types
    /// ```
    #[pyo3(get)]
    pub point_types: Option<Vec<usize>>,

    /// Padded, tiled point set used only when `pbc = true`, shape `(M, 2)`.
    ///
    /// Built inside [`Layer::generate_kdtree`] by replicating `points` into a
    /// 3x3 super-tiling and keeping only copies inside the padded boundary polygon
    /// (controlled by `study_proximity`). The KD-tree is built from this array
    /// so that radius searches near the supercell edge pick up periodic images.
    /// Internal use only.
    pub bigger_points: Option<Array2<f64>>,  // #[getter]

    /// Maps each `bigger_points` index back to its primary-cell `points` index. Length M.
    /// As M > N, this is a many to one mapping.
    ///
    /// ```text
    /// mappings = [0, 1, 0, 3, ...]   // mappings[b] == p  =>  bigger_points[b] is a periodic image of points[p]
    /// ```
    ///
    /// `None` for OBC layers (no periodic images needed).
    pub mappings: Option<Vec<usize>>,  // private

    // -------------------------------------------------------------------------
    // 4. Search & Transformation State
    // -------------------------------------------------------------------------

    /// The spatial index powering all neighbour queries.
    ///
    /// Built from `bigger_points` (PBC) or `points` (OBC) at the end of
    /// [`Layer::generate_kdtree`]. `None` before `generate_points` is called,
    /// or in `test` mode where the tree is intentionally skipped.
    pub kdtree: Option<LatticeTree>,                // private

    /// 2x2 rotation matrix applied to the layer, shape `(2, 2)`.
    ///
    /// ```text
    /// rot_m = [[ cos, -sin],
    ///          [ sin,  cos]]   // standard 2D CCW rotation
    /// ```
    ///
    /// Starts as the identity (theta = 0). Updated by [`Layer::perform_rotation_translation`].
    pub rot_m: Array2<f64>,  // #[getter]

    /// Translation vector `[dx, dy]` applied to the basis before replication.
    ///
    /// Starts as `[0, 0]`.  Updated by [`Layer::perform_rotation_translation`].
    pub translation: Array1<f64>,  // #[getter]
}

/// PyO3-exported methods, callable from Python.
#[pymethods]
impl Layer {
    /// Constructs a new `Layer` without generating any lattice points.
    ///
    /// Records geometry and connectivity data but does NOT build the point cloud
    /// or KD-tree. Call [`Layer::generate_points`] for that.
    ///
    /// # Arguments
    /// * `lv1`: First lattice vector `[x, y]`. Must lie along the x-axis (`lv1[1] == 0`).
    /// * `lv2`: Second lattice vector `[x, y]`. Must have `lv2[1] > 0`.
    /// * `basis_data`: List of `(x, y, type_name)` tuples for each atom in the
    ///   unit cell. Duplicate `type_name` values will panic.
    /// * `neighbours`: Map from each `type_name` to a list of `[dx, dy]` displacement
    ///   vectors pointing to its first neighbours. Must be in the same Cartesian
    ///   frame as `lv1`/`lv2`.
    /// * `pbc`: `true` for Periodic Boundary Conditions, `false` for Open.
    /// * `study_proximity`: Extra neighbour shells retained in the PBC padding region.
    ///   `1` is sufficient for first-NN searches.
    ///
    /// # Panics
    /// Panics if `basis_data` contains a duplicate type name.
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

    /// Applies a rigid-body rotation and translation to the layer before points are generated.
    ///
    /// Must be called before [`Layer::generate_points`]; returns `PyAssertionError` otherwise.
    ///
    /// Steps: (1) builds 2x2 rotation matrix R for the given angle,
    /// (2) rotates `lv1`/`lv2` in-place (`lv = R @ lv`),
    /// (3) shifts then rotates basis points (`p_new = (p_old + translation) @ R^T`),
    /// (4) rotates neighbour displacement vectors (translation has no effect on deltas).
    ///
    /// # Arguments
    /// * `rot`: Rotation angle in radians (counter-clockwise positive).
    /// * `translation`: `[dx, dy]` shift applied to basis points before rotation.
    ///
    /// # Errors
    /// Returns `PyAssertionError` if `self.points` is already populated.
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
    /// Generates all real-space atom positions for the full Moiré supercell.
    ///
    /// Proceeds in three stages:
    /// 1. Single-cell generation: sweeps a grid of integer multiples of `lv1`/`lv2`
    ///    large enough to cover one Moiré cell `mlv1 x mlv2`, places every basis atom
    ///    at each grid point, keeps only those inside the boundary
    ///    (via [`Layer::inside_boundaries`]).
    /// 2. Tiling: replicates the single-cell result `mln1 x mln2` times.
    /// 3. KD-tree construction: calls [`Layer::generate_kdtree`], which also builds
    ///    `bigger_points` and `mappings` when `pbc = true`.
    ///
    /// After this call, `self.points`, `self.point_types`, and `self.kdtree` are
    /// all populated.
    ///
    /// # Arguments
    /// * `mlv1`: First Moiré super-cell vector `[x, y]`.
    /// * `mlv2`: Second Moiré super-cell vector `[x, y]`.
    /// * `mln1`: Number of cells to tile along `mlv1` (must be >= 1).
    /// * `mln2`: Number of cells to tile along `mlv2` (must be >= 1).
    /// * `test`: If `true`, skip KD-tree construction (for tests that only need coordinates).
    ///
    /// # Errors
    /// Returns `PyValueError` if `mln1`/`mln2` is zero, or if `test = true` with `mln > 1`.
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

    /// Finds the first nearest neighbours of every point in `points`.
    ///
    /// For each input point, looks up its basis type, retrieves the corresponding
    /// neighbour displacement vectors from `self.neighbours`, and queries the
    /// KD-tree at `point + delta` for each delta. Thin wrapper around
    /// [`Layer::first_nearest_neighbours_internal`].
    ///
    /// # Arguments
    /// * `py`: Python interpreter token (provided automatically by PyO3).
    /// * `points`: `(N, 2)` NumPy array of Cartesian coordinates to search from.
    /// * `types`: List of N basis-type strings (e.g. `["A", "B", "A", ...]`).
    ///   Determines which neighbour deltas are used per point.
    ///
    /// # Returns
    /// A 3-tuple of 2-D NumPy arrays, each shape `(N, max_neighbours)`:
    /// 1. `bigger_indices` (`i64`): raw KD-tree item index of each neighbour.
    ///    Unfilled entries are `-1`.
    /// 2. `smaller_indices` (`i64`): neighbour index mapped to the primary cell
    ///    via `self.mappings`. Unfilled entries are `-1`.
    /// 3. `distances` (`f64`): Euclidean distance to each neighbour. Unfilled: `NaN`.
    ///
    /// # Errors
    /// Returns `PyValueError` if an unknown type string is passed, or if any
    /// KD-tree hit exceeds the geometry tolerance.
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

    /// Finds all neighbours within a specified Euclidean distance of every query point.
    ///
    /// Queries the internal KD-tree to find every point (including periodic images
    /// if PBC is enabled) that falls within `radius` of each point in `query_points`.
    ///
    /// Because different query points may have different neighbour counts, results
    /// are returned as flat 1-D arrays. The `q_indices` array encodes which query
    /// point each hit belongs to (like a CSR row-pointer array). Thin wrapper
    /// around [`Layer::get_neighbors_within_radius_internal`].
    ///
    /// # Arguments
    /// * `py`: Python interpreter token (provided automatically by PyO3).
    /// * `query_points`: `(M, 2)` NumPy array of coordinates to search from.
    /// * `radius`: Maximum Euclidean distance for a neighbour to be included.
    ///
    /// # Returns
    /// A 4-tuple of 1-D NumPy arrays, all length K (total hits across all query points):
    /// 1. `q_indices` (`usize`): index of the query point for this hit, e.g. `[0,0,0,1,1,2,...]`.
    /// 2. `bigger_indices` (`usize`): raw KD-tree item index of the matched neighbour.
    /// 3. `smaller_indices` (`usize`): matched index mapped back to the primary unit cell.
    /// 4. `distances` (`f64`): actual Euclidean distance to the neighbour.
    ///
    /// # Errors
    /// Returns `PyRuntimeError` if the KD-tree has not been initialized.
    pub fn get_neighbors_within_radius<'py>(
        &self,
        py: Python<'py>,
        query_points: PyReadonlyArray2<'py, f64>,
        radius: f64,
    ) -> PyResult<(
        Bound<'py, PyArray1<usize>>,
        Bound<'py, PyArray1<usize>>,
        Bound<'py, PyArray1<usize>>,
        Bound<'py, PyArray1<f64>>,
    )> {
        let q_pts = query_points.as_array();

        let (q_indices, bigger_indices, smaller_indices, distances) = self
            .get_neighbors_within_radius_internal(&q_pts, radius)
            .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(e))?;

        Ok((
            q_indices.to_pyarray(py),
            bigger_indices.to_pyarray(py),
            smaller_indices.to_pyarray(py),
            distances.to_pyarray(py),
        ))
    }

    // Read-only NumPy getters. Convert internal ndarray fields to NumPy arrays on the fly.

    /// Returns `lv1` as a 1-D NumPy array of shape `(2,)`.
    #[getter]
    pub fn lv1<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.lv1.to_pyarray(py) //
    }

    /// Returns `lv2` as a 1-D NumPy array of shape `(2,)`.
    #[getter]
    pub fn lv2<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.lv2.to_pyarray(py) //
    }

    /// Returns the first Moiré supercell vector `mlv1` as a NumPy array `(2,)`.
    #[getter]
    pub fn mlv1<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mlv1.to_pyarray(py) //
    }

    /// Returns the second Moiré supercell vector `mlv2` as a NumPy array `(2,)`.
    #[getter]
    pub fn mlv2<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.mlv2.to_pyarray(py) //
    }

    /// Returns the `[dx, dy]` translation vector as a NumPy array `(2,)`.
    #[getter]
    pub fn translation<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        self.translation.to_pyarray(py) //
    }

    /// Returns `basis_points` (atom positions in one unit cell) as a NumPy
    /// array of shape `(B, 2)` where B is the number of basis atoms.
    #[getter]
    pub fn basis_points<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.basis_points.to_pyarray(py) //
    }

    /// Returns the current 2x2 rotation matrix as a NumPy array of shape `(2, 2)`.
    #[getter]
    pub fn rot_m<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray2<f64>> {
        self.rot_m.to_pyarray(py) //
    }

    /// Returns the primary-cell atom positions as an `(N, 2)` NumPy array,
    /// or `None` if `generate_points` has not been called yet.
    #[getter]
    pub fn points<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        self.points.as_ref().map(|p| p.to_pyarray(py)) //
    }

    /// Returns the PBC-padded point set (`bigger_points`) as an `(M, 2)` NumPy
    /// array, or `None` for OBC layers or before `generate_points` is called.
    #[getter]
    pub fn bigger_points<'py>(&self, py: Python<'py>) -> Option<Bound<'py, PyArray2<f64>>> {
        self.bigger_points.as_ref().map(|p| p.to_pyarray(py)) //
    }

    /// Returns the basis-type *name* for every point in `points` as a Python
    /// list of strings, or `None` before `generate_points` is called.
    ///
    /// Note: internally types are stored as integer IDs; this getter converts
    /// them back to the original string names via `basis_types`.
    #[getter]
    pub fn point_types(&self) -> Option<Vec<String>> {
        self.point_types.as_ref().map(|indices| {
            indices
                .iter()
                .map(|&i| self.basis_types[i].clone())
                .collect()
        })
    }

    /// Returns the neighbour displacement arrays as a Python list of NumPy
    /// arrays.  Index `k` is an `(N_k, 2)` array of `[dx, dy]` vectors for
    /// basis type `k`.
    #[getter]
    pub fn neighbours<'py>(&self, py: Python<'py>) -> Vec<Bound<'py, PyArray2<f64>>> {
        self.neighbours.iter().map(|n| n.to_pyarray(py)).collect()
    }
}

/// Internal Rust-only methods, not exported to Python.
///
/// All methods here are `pub` (so `moire.rs` can call them) but are deliberately
/// NOT inside the `#[pymethods]` block, so PyO3 never exposes them to Python.
impl Layer {
    /// Tests which of the given `points` lie strictly inside a convex `polygon`
    /// using the ray-casting algorithm.
    ///
    /// A point is inside if a horizontal ray cast rightward from it crosses an
    /// odd number of polygon edges. Points exactly on an edge are treated as outside.
    ///
    /// # Arguments
    /// * `points`: `(N, 2)` array of `[x, y]` coordinates to test.
    /// * `polygon`: `(V, 2)` array of polygon vertices in order. The edge list
    ///   wraps from vertex `V-1` back to vertex `0`.
    ///
    /// # Returns
    /// `Vec<bool>` of length N; `true` if the corresponding point is inside.
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

    /// Builds the KD-tree (and PBC padding structures if needed).
    ///
    /// Called automatically at the end of [`Layer::generate_points`].
    ///
    /// OBC path (`pbc = false`): builds `self.kdtree` directly from `self.points`.
    ///
    /// PBC path (`pbc = true`):
    /// 1. Creates a 3x3 tiling (shifts -1, 0, +1 along both supercell vectors).
    /// 2. Computes a padded boundary polygon extending `study_proximity` shells beyond the cell.
    /// 3. Filters tiled points inside that polygon -> stored as `self.bigger_points`.
    /// 4. Builds `self.kdtree` from `bigger_points`.
    /// 5. Calls [`Layer::generate_mapping`] to build the bigger-to-primary index map.
    ///
    /// # Errors
    /// Returns `Err(String)` if `self.points` is `None`, or if [`Layer::generate_mapping`] fails.
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
        // 2. Tiling logic: generate 9-1 copies of the lattice (-1, 0, 1 shifts on two axes)
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

    /// Tests which `points` lie inside the Moiré supercell boundary parallelogram.
    ///
    /// The boundary is defined by `mlv1 * mln1` and `mlv2 * mln2`. The polygon
    /// is shifted slightly inward (by `toll_scale * 1e-4` along the inward diagonal)
    /// so that atoms on a shared edge are assigned to only one cell.
    ///
    /// # Arguments
    /// * `points`: `(N, 2)` array of Cartesian coordinates to test.
    /// * `mln1`: Override for the `mln1` field (`None` uses `self.mln1`).
    /// * `mln2`: Override for the `mln2` field (`None` uses `self.mln2`).
    ///
    /// # Returns
    /// `Vec<bool>` of length N; `true` means the point is inside.
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

    /// Builds the `bigger_points -> primary points` index map (`self.mappings`).
    ///
    /// For each point in `bigger_points`, finds the closest point in `points`
    /// (after folding the bigger point back into the primary cell via
    /// [`Layer::get_point_positions`]). Stores the result as `self.mappings`:
    /// a `Vec<usize>` of length `bigger_points.nrows()` where `mappings[b]`
    /// is the primary-cell index for the periodic image at `bigger_points[b]`.
    ///
    /// # Errors
    /// Returns `Err(String)` if `self.points` or `self.bigger_points` is `None`,
    /// or if any entry cannot be matched within `toll_scale * 1e-3` (geometry mismatch).
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

    /// Computes the supercell fractional coordinates of each point w.r.t. vectors `a` and `b`.
    ///
    /// Returns an `(N, 2)` integer-valued array where each row `[n_a, n_b]` is the
    /// number of full supercell lengths the point is offset from the origin.
    /// Used by [`Layer::generate_mapping`] to fold `bigger_points` back into the
    /// primary cell before nearest-neighbour matching.
    ///
    /// Implemented via 2x2 determinant tests (parallelogram intersection), which
    /// is equivalent to solving `point = n_a * a + n_b * b` with integer rounding.
    ///
    /// # Arguments
    /// * `points`: `(N, 2)` array of Cartesian coordinates.
    /// * `a`: First supercell vector `(2,)`.
    /// * `b`: Second supercell vector `(2,)`.
    ///
    /// # Returns
    /// `(N, 2)` array; column 0 is the integer offset along `a`, column 1 along `b`.
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

    /// Core nearest-neighbour search, the computational heart of the NN pipeline.
    ///
    /// For each of the N input points (with associated basis types), looks up the
    /// pre-stored neighbour displacement vectors for that type, queries the KD-tree
    /// at `point + delta` for every delta, and records the best match.
    /// Results are stored in fixed-width 2-D arrays (`N x max_neighbours`) so that
    /// the output is a single contiguous block of memory regardless of how many
    /// neighbours each type has. Unfilled slots are filled with sentinel values:
    /// `-1` for indices, `NaN` for distances.
    ///
    /// This is kept separate from the `#[pymethods]` wrapper so that `moire.rs`
    /// can call it directly without Python/NumPy conversion overhead.
    ///
    /// # Arguments
    /// * `points`: `(N, 2)` view of Cartesian coordinates.
    /// * `types`: Slice of N integer type-IDs (indices into `self.neighbours`).
    ///
    /// # Returns
    /// A 3-tuple of owned arrays:
    /// 1. `bigger_indices`: `(N, max_neighbours)` `i64` array of raw KD-tree item indices.
    ///    Unfilled slots are `-1`.
    /// 2. `smaller_indices`: `(N, max_neighbours)` `i64` array of primary-cell indices
    ///    (after applying `self.mappings`). Unfilled slots are `-1`.
    /// 3. `distances_arr`: `(N, max_neighbours)` `f64` Euclidean distances.
    ///    Unfilled slots are `NaN`.
    ///
    /// # Errors
    /// Returns `Err(String)` if the KD-tree is not initialized, or if any KD-tree hit
    /// is further from the expected position than `toll_scale * 1e-2`.
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

    /// Core radius search: finds all points within `radius` of every query point.
    ///
    /// Unlike [`Layer::first_nearest_neighbours_internal`] (which uses pre-defined
    /// displacement vectors), this search is purely geometric: queries the KD-tree
    /// with `.within()` (squared-distance metric) and collects every hit.
    ///
    /// Because different query points may have different hit counts, the output is
    /// four flat 1-D arrays. The `q_indices` array acts as a row-pointer:
    /// entry k belongs to query point `q_indices[k]`.
    ///
    /// # Arguments
    /// * `query_points`: `(M, 2)` view of Cartesian query coordinates.
    /// * `radius`: Maximum Euclidean search radius. KD-tree is queried with `radius^2`
    ///   (Kiddo uses squared distances internally).
    ///
    /// # Returns
    /// A 4-tuple of `Array1`s, all length K (total number of hits):
    /// 1. `q_indices` (`usize`): which query point this hit belongs to.
    ///    Monotonically non-decreasing, e.g. `[0,0,1,1,1,2,...]`.
    /// 2. `bigger_indices` (`usize`): raw KD-tree item index of the matched point.
    ///    For OBC equals the primary-cell index; for PBC may point into `bigger_points`.
    /// 3. `smaller_indices` (`usize`): the matched point remapped to the primary cell
    ///    via `self.mappings` (same as `bigger_indices` for OBC).
    /// 4. `distances_arr` (`f64`): `sqrt(squared_distance)` for each hit.
    ///
    /// # Errors
    /// Returns `Err(String)` if `self.kdtree` is `None`.
    pub fn get_neighbors_within_radius_internal(
        &self,
        query_points: &ArrayView2<f64>,
        radius: f64,
    ) -> Result<(Array1<usize>, Array1<usize>, Array1<usize>, Array1<f64>), String> {
        let tree = self.kdtree.as_ref().ok_or("KDTree not initialized.")?;
        let radius_sq = radius.powi(2);

        let mut q_indices = Vec::new();
        let mut bigger_indices = Vec::new();
        let mut smaller_indices = Vec::new();
        let mut distances_vec = Vec::new();

        for (i, q_pt) in query_points.rows().into_iter().enumerate() {
            let matches = tree.within::<kiddo::SquaredEuclidean>(&[q_pt[0], q_pt[1]], radius_sq);

            for m in matches {
                let b_idx = m.item as usize;
                q_indices.push(i);
                bigger_indices.push(b_idx);

                if let Some(ref maps) = self.mappings {
                    smaller_indices.push(maps[b_idx]);
                } else {
                    smaller_indices.push(b_idx);
                }

                distances_vec.push(m.distance.sqrt());
            }
        }

        Ok((
            Array1::from(q_indices),
            Array1::from(bigger_indices),
            Array1::from(smaller_indices),
            Array1::from(distances_vec),
        ))
    }
}
