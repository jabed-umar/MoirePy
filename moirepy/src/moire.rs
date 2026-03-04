use num_complex::Complex;
use ndarray::{Array1};
use pyo3::prelude::*;
use numpy::{PyReadonlyArray1,  IntoPyArray, Element, PyArray1};
use crate::layers::Layer;
use crate::utils::{COOBuilder, HoppingInstruction, SelfInstruction, ShiftInstruction};

/// A trait that abstracts over scalar and vector hopping amplitudes.
///
/// When assembling the Hamiltonian it is common to have either:
/// - A single number (e.g. `-2.7 eV`) that is the same for every bond, or
/// - A per-bond list of numbers (one entry per hopping pair).
///
/// Rather than branching on this distinction everywhere, this trait is
/// implemented for both `T` (scalar) and `Vec<T>` (vector) and Rust dispatches
/// automatically. The `Input<T>` enum delegates to whichever branch
/// matches the Python argument.
///
/// # Type parameters
/// * `T`: the numeric type of the amplitude (e.g. `f64` or `Complex<f64>`).
trait Pushable<T> {
    /// Fills the COO builder with on-site (self-energy) entries.
    ///
    /// # Arguments
    /// * `builder`: the COO triplet accumulator to write into.
    /// * `inst`: the `SelfInstruction` describing which sites to fill and their types.
    /// * `off`: row/column offset for this layer block (0 for lower, `n_lower * orbitals` for upper).
    /// * `k`: number of orbitals per site (block size).
    fn push_self(&self, builder: &mut COOBuilder<T>, inst: &SelfInstruction, off: i32, k: i32);

    /// Fills the COO builder with hopping (off-diagonal) entries.
    ///
    /// # Arguments
    /// * `builder`: the COO triplet accumulator to write into.
    /// * `inst`: the `HoppingInstruction` listing `(i, j)` site pairs and their basis types.
    /// * `r_off`: row-block offset for the source layer.
    /// * `c_off`: column-block offset for the target layer.
    /// * `k`: number of orbitals per site (block size).
    fn push_hop(&self, builder: &mut COOBuilder<T>, inst: &HoppingInstruction, r_off: i32, c_off: i32, k: i32);
}

/// Pushable impl for a scalar amplitude `T`.
///
/// Each hopping/self-energy entry gets the same value (`*self`) regardless of
/// which site pair it belongs to. Fills an entire `k x k` orbital block uniformly.
impl<T: Copy> Pushable<T> for T {
    fn push_self(&self, builder: &mut COOBuilder<T>, inst: &SelfInstruction, off: i32, k: i32) {
        for &i in &inst.site_i {
            // Fill the entire k x k block for this site
            for row in 0..k {
                for col in 0..k {
                    builder.add(i * k + row + off, i * k + col + off, *self);
                }
            }
        }
    }

    fn push_hop(&self, builder: &mut COOBuilder<T>, inst: &HoppingInstruction, r_off: i32, c_off: i32, k: i32) {
        for (&i, &j) in inst.site_i.iter().zip(inst.site_j.iter()) {
            // Fill the entire k x k block for this hopping pair
            for row in 0..k {
                for col in 0..k {
                    builder.add(i * k + row + r_off, j * k + col + c_off, *self);
                }
            }
        }
    }
}

/// Pushable impl for a vector of amplitudes `Vec<T>`.
///
/// The vector must be a flat concatenation of `k x k` blocks, one block per
/// hopping/self-energy pair: `self[pair_idx * k^2 .. (pair_idx+1) * k^2]` is
/// the orbital block for pair `pair_idx`.
/// This is the path taken when the caller passes a Python list/array of per-bond amplitudes.
impl<T: Copy> Pushable<T> for Vec<T> {
    fn push_self(&self, builder: &mut COOBuilder<T>, inst: &SelfInstruction, off: i32, k: i32) {
        let k_sq = (k * k) as usize;
        // Zip the site indices with k*k chunks of the flattened data
        for (&i, block) in inst.site_i.iter().zip(self.chunks_exact(k_sq)) {
            for row in 0..k {
                for col in 0..k {
                    let val = block[(row * k + col) as usize];
                    // Self-energy can have off-diagonal orbital terms too!
                    builder.add(i * k + row + off, i * k + col + off, val);
                }
            }
        }
    }

    fn push_hop(&self, builder: &mut COOBuilder<T>, inst: &HoppingInstruction, r_off: i32, c_off: i32, k: i32) {
        let k_sq = (k * k) as usize;
        // Zip (i, j) pairs with k*k chunks of the flattened data
        for ((&i, &j), block) in inst.site_i.iter().zip(inst.site_j.iter()).zip(self.chunks_exact(k_sq)) {
            for row in 0..k {
                for col in 0..k {
                    let val = block[(row * k + col) as usize];
                    builder.add(i * k + row + r_off, j * k + col + c_off, val);
                }
            }
        }
    }
}


/// Represents a hopping amplitude that is either a single value or per-bond.
///
/// When calling [`BilayerMoire::build_ham`] from Python you can pass:
/// - A scalar `float` or `complex`: the same amplitude for all bonds.
/// - A list/NumPy array: one amplitude per bond, in the order produced
///   by [`BilayerMoire::generate_connections`].
///
/// PyO3's `FromPyObject` derive handles the conversion automatically.
#[derive(FromPyObject)]
pub enum Input<T> {
    /// A single amplitude applied uniformly to every entry.
    Scalar(T),
    /// Per-bond amplitudes; must have as many elements as hopping pairs times `orbitals^2`.
    Vector(Vec<T>),
}

// Now we implement the Pushable trait FOR the Input enum itself.
// This is where the routing happens.
/// Pushable impl for the `Input<T>` enum.
///
/// Matches on whether the Python caller passed a scalar or a list and forwards
/// to the correct concrete impl. All methods on `BilayerMoire` accept `Input<T>`.
impl<T: Copy> Pushable<T> for Input<T> {
    fn push_self(&self, builder: &mut COOBuilder<T>, inst: &SelfInstruction, off: i32, k: i32) {
        match self {
            // If it's a Scalar, call Definition 1 (The Scalar Worker)
            Input::Scalar(val) => val.push_self(builder, inst, off, k),
            // If it's a Vector, call Definition 2 (The Vector Worker)
            Input::Vector(vec) => vec.push_self(builder, inst, off, k),
        }
    }

    fn push_hop(&self, builder: &mut COOBuilder<T>, inst: &HoppingInstruction, r_off: i32, c_off: i32, k: i32) {
        match self {
            Input::Scalar(val) => val.push_hop(builder, inst, r_off, c_off, k),
            Input::Vector(vec) => vec.push_hop(builder, inst, r_off, c_off, k),
        }
    }
}



/// Represents a twisted bilayer Moiré system built from two [`Layer`]s.
///
/// Owns references to two pre-initialised `Layer` objects (lower + upper), stores
/// the geometric parameters of the Moiré supercell, and caches the hopping
/// instructions computed in [`BilayerMoire::generate_connections`]. The
/// Hamiltonian is assembled on demand in [`BilayerMoire::build_ham`] or
/// [`BilayerMoire::build_ham_complex`].
#[pyclass]
pub struct BilayerMoire {
    // -------------------------------------------------------------------------
    // 1. Layer Ownership
    // -------------------------------------------------------------------------

    /// Reference to the lower crystal layer.
    ///
    /// Stored as `Py<Layer>` (a Python-owned smart pointer) so the same `Layer`
    /// object is shared between Rust and Python without copying.
    /// Use `.bind(py).borrow()` to get an immutable Rust reference inside a method.
    pub lower_lattice: Py<Layer>,                   // private

    /// Reference to the upper crystal layer.
    pub upper_lattice: Py<Layer>,                   // private

    // -------------------------------------------------------------------------
    // 2. Geometric (AVC) Coefficients
    // These integers parameterise the commensurate twist angle and come from
    // the AVC (Angle-Vector-Commensurate) construction.
    // -------------------------------------------------------------------------

    /// Integer coefficient `l1` for the lower lattice supercell vector.
    #[pyo3(get)]
    pub ll1: i32,

    /// Integer coefficient `l2` for the lower lattice supercell vector.
    #[pyo3(get)]
    pub ll2: i32,

    /// Integer coefficient `u1` for the upper lattice supercell vector.
    #[pyo3(get)]
    pub ul1: i32,

    /// Integer coefficient `u2` for the upper lattice supercell vector.
    #[pyo3(get)]
    pub ul2: i32,

    // -------------------------------------------------------------------------
    // 3. Supercell Parameters
    // -------------------------------------------------------------------------

    /// Number of Moiré super-cells tiled along the first Moiré vector.
    #[pyo3(get)]
    pub n1: usize,

    /// Number of Moiré super-cells tiled along the second Moiré vector.
    #[pyo3(get)]
    pub n2: usize,

    /// Twist angle in **radians** (the physical parameter of interest).
    #[pyo3(get)]
    pub theta: f64,

    /// First Moiré supercell vector M1, shape `(2,)`.
    ///
    /// Kept as an owned `Array1<f64>`; exposed to Python via the hand-written `#[getter]`.
    pub mlv1: Array1<f64>,  // #[getter]

    /// Second Moiré supercell vector M2, shape `(2,)`.
    pub mlv2: Array1<f64>,  // #[getter]

    /// Translation vector `[dx, dy]` applied to the upper layer basis.
    ///
    /// Needed when the two layers are not on the same origin (AA vs AB stacking).
    /// Exposed to Python via `#[getter]`.
    pub translate_upper: Array1<f64>,  // #[getter]

    // -------------------------------------------------------------------------
    // 4. Global Simulation State
    // -------------------------------------------------------------------------

    /// Whether Periodic Boundary Conditions are active for both layers.
    #[pyo3(get)]
    pub pbc: bool,

    /// Number of orbitals per lattice site (block size of the Hamiltonian).
    ///
    /// `orbitals = 1` produces a scalar tight-binding model; `orbitals = 2`
    /// produces a 2x2 block Hamiltonian (e.g. two sublattice orbitals per site).
    #[pyo3(get)]
    pub orbitals: usize,

    // -------------------------------------------------------------------------
    // 5. Hopping Instructions (cached after generate_connections)
    // Each Option is None before generate_connections is called.
    // -------------------------------------------------------------------------

    /// Intra-layer hopping pairs for the lower layer.
    ///
    /// Each entry is a `(site_i, site_j, type_i, type_j)` tuple describing one directed
    /// hopping. Built from the nearest-neighbour search in [`BilayerMoire::generate_connections`].
    pub hop_ll: Option<HoppingInstruction>,  // #[getter]

    /// Intra-layer hopping pairs for the upper layer.
    pub hop_uu: Option<HoppingInstruction>,  // #[getter]

    /// Inter-layer hopping from lower to upper (H_lu block).
    pub hop_lu: Option<HoppingInstruction>,  // #[getter]

    /// Inter-layer hopping from upper to lower (H_ul block).
    pub hop_ul: Option<HoppingInstruction>,  // #[getter]

    /// On-site (self-energy) entries for the lower layer sites.
    pub hop_l: Option<SelfInstruction>,  // #[getter]

    /// On-site (self-energy) entries for the upper layer sites.
    pub hop_u: Option<SelfInstruction>,  // #[getter]

    /// Cached count of non-zero entries for pre-allocating the COO builder.
    ///
    /// Computed at the end of `generate_connections` via `calculate_total_nnz`.
    /// `None` before that call.
    pub nnz: Option<usize>,                         // private

    /// Cached COO shift entries used for k-space phase factors.
    ///
    /// Each COO entry stores a row, column, and 2D shift vector `(dx, dy)`.
    /// `None` until `generate_connections` is called.
    shifts: Option<ShiftInstruction>,               // private
}



/// PyO3-exported methods, callable from Python.
#[pymethods]
impl BilayerMoire {
    /// Constructs a new `BilayerMoire` system from two pre-initialised layers.
    ///
    /// The caller must generate points on both layers (`layer.generate_points(...)`) before
    /// creating the `BilayerMoire`; this constructor only stores references and geometry.
    ///
    /// # Arguments
    /// * `lower`: the lower crystal layer (Python `Layer` object).
    /// * `upper`: the upper crystal layer (rotated by `theta`).
    /// * `ll1`, `ll2`: AVC integer coefficients for the lower lattice.
    /// * `ul1`, `ul2`: AVC integer coefficients for the upper lattice.
    /// * `n1`, `n2`: tiling counts for the Moiré supercell.
    /// * `theta`: twist angle in radians.
    /// * `mlv1`, `mlv2`: Moiré supercell vectors as 1-D NumPy arrays.
    /// * `translate_upper`: `[dx, dy]` stacking offset for the upper layer.
    /// * `pbc`: whether to use periodic boundary conditions.
    /// * `orbitals`: number of orbitals per site (Hamiltonian block size).
    #[new]
    #[pyo3(signature = (lower, upper, ll1, ll2, ul1, ul2, n1, n2, theta, mlv1, mlv2, translate_upper, pbc, orbitals))]
    pub fn new(
        lower: Bound<'_, Layer>,
        upper: Bound<'_, Layer>,
        ll1: i32,
        ll2: i32,
        ul1: i32,
        ul2: i32,
        n1: usize,
        n2: usize,
        theta: f64,
        mlv1: PyReadonlyArray1<f64>,
        mlv2: PyReadonlyArray1<f64>,
        translate_upper: PyReadonlyArray1<f64>,
        pbc: bool,
        orbitals: usize,
    ) -> Self {
        BilayerMoire {
            // unbind() captures the Python-owned Layer into a Py<Layer> smart pointer
            lower_lattice: lower.unbind(),
            upper_lattice: upper.unbind(),
            ll1,
            ll2,
            ul1,
            ul2,
            n1,
            n2,
            theta,
            // Converting NumPy views to owned ndarray Array1s
            mlv1: mlv1.as_array().to_owned(),
            mlv2: mlv2.as_array().to_owned(),
            translate_upper: translate_upper.as_array().to_owned(),
            pbc,
            orbitals,
            hop_ll: None,
            hop_uu: None,
            hop_lu: None,
            hop_ul: None,
            hop_l: None,
            hop_u: None,
            nnz: None,
            shifts: None,
        }
    }

    /// Computes and caches all hopping and self-energy instructions.
    ///
    /// Must be called once before [`BilayerMoire::build_ham`]. Calling again is
    /// safe: it clears and rebuilds all instructions.
    ///
    /// Steps:
    /// 1. Self-energies (`hop_l`, `hop_u`): records every site index and its type.
    /// 2. Intra-layer LL (`hop_ll`): calls [`Layer::first_nearest_neighbours_internal`]
    ///    on the lower layer; stores all `(i, j)` pairs where `j >= 0`.
    /// 3. Intra-layer UU (`hop_uu`): same for the upper layer.
    /// 4. Inter-layer LU/UL (`hop_lu`, `hop_ul`): radius search on the lower KD-tree
    ///    using each upper-layer site as the query point. For PBC, raw indices are
    ///    remapped via `lower.mappings`.
    ///
    /// Also updates `self.nnz` for memory pre-allocation and refreshes the cached
    /// shift COO entries used by `get_phase`.
    ///
    /// # Arguments
    /// * `py`: Python interpreter token.
    /// * `inter_layer_radius`: maximum Euclidean distance for an inter-layer bond
    ///   to be included (e.g. ~3.35 A for graphene interlayer distance).
    pub fn generate_connections(&mut self, py: Python<'_>, inter_layer_radius: f64) {
        let lower = self.lower_lattice.bind(py).borrow();
        let upper = self.upper_lattice.bind(py).borrow();

        let points_l = lower.points.as_ref().expect("Lower points missing");
        let points_u = upper.points.as_ref().expect("Upper points missing");

        // Use internal usize type IDs directly
        let internal_types_l = lower.point_types.as_ref().expect("Lower types missing");
        let internal_types_u = upper.point_types.as_ref().expect("Upper types missing");

        let n_l = points_l.nrows();
        let n_u = points_u.nrows();

        // Efficiently reuse or initialize instruction buffers
        if let Some(ref mut h) = self.hop_ll { h.clear(); } else { self.hop_ll = Some(HoppingInstruction::new(n_l)); }
        if let Some(ref mut h) = self.hop_uu { h.clear(); } else { self.hop_uu = Some(HoppingInstruction::new(n_u)); }
        if let Some(ref mut h) = self.hop_ul { h.clear(); } else { self.hop_ul = Some(HoppingInstruction::new(n_u)); }
        if let Some(ref mut h) = self.hop_lu { h.clear(); } else { self.hop_lu = Some(HoppingInstruction::new(n_l)); }
        if let Some(ref mut h) = self.hop_l { h.clear(); } else { self.hop_l = Some(SelfInstruction::new(n_l)); }
        if let Some(ref mut h) = self.hop_u { h.clear(); } else { self.hop_u = Some(SelfInstruction::new(n_u)); }
        if let Some(ref mut s) = self.shifts { s.clear(); } else { self.shifts = Some(ShiftInstruction::new((n_l + n_u) * 8)); }

        let h_ll = self.hop_ll.as_mut().unwrap();
        let h_uu = self.hop_uu.as_mut().unwrap();
        let h_ul = self.hop_ul.as_mut().unwrap();
        let h_lu = self.hop_lu.as_mut().unwrap();
        let h_l = self.hop_l.as_mut().unwrap();
        let h_u = self.hop_u.as_mut().unwrap();
        let shifts = self.shifts.as_mut().unwrap();
        let offset = n_l as i32;

        let lower_big = if lower.pbc {
            lower.bigger_points.as_ref().expect("Lower bigger_points missing under PBC")
        } else {
            points_l
        };
        let upper_big = if upper.pbc {
            upper.bigger_points.as_ref().expect("Upper bigger_points missing under PBC")
        } else {
            points_u
        };

        // 1. Self-Energies
        for i in 0..n_l {
            h_l.add(i as i32, internal_types_l[i] as i32);
            shifts.add(i as i32, i as i32, 0.0, 0.0);
        }
        for i in 0..n_u {
            h_u.add(i as i32, internal_types_u[i] as i32);
            shifts.add(i as i32 + offset, i as i32 + offset, 0.0, 0.0);
        }

        // 2. Intra-layer (Lower) - Direct Internal Call
        let (big_idx_l, small_idx_l, _) = lower.first_nearest_neighbours_internal(
            &points_l.view(),
            internal_types_l
        ).expect("Lower neighbors failed");

        for i in 0..n_l {
            for (&j, &b) in small_idx_l[i].iter().zip(big_idx_l[i].iter()) {
                let j_us = j as usize;
                let b_us = b as usize;
                h_ll.add(i as i32, j as i32, internal_types_l[i] as i32, internal_types_l[j_us] as i32);
                let dx = lower_big[[b_us, 0]] - points_l[[j_us, 0]];
                let dy = lower_big[[b_us, 1]] - points_l[[j_us, 1]];
                shifts.add(i as i32, j as i32, dx, dy);
            }
        }

        // 3. Intra-layer (Upper) - Direct Internal Call
        let (big_idx_u, small_idx_u, _) = upper.first_nearest_neighbours_internal(
            &points_u.view(),
            internal_types_u
        ).expect("Upper neighbors failed");

        for i in 0..n_u {
            for (&j, &b) in small_idx_u[i].iter().zip(big_idx_u[i].iter()) {
                let j_us = j as usize;
                let b_us = b as usize;
                h_uu.add(i as i32, j as i32, internal_types_u[i] as i32, internal_types_u[j_us] as i32);
                let dx = upper_big[[b_us, 0]] - points_u[[j_us, 0]];
                let dy = upper_big[[b_us, 1]] - points_u[[j_us, 1]];
                shifts.add(i as i32 + offset, j as i32 + offset, dx, dy);
            }
        }

        // 4. Inter-layer
        if let Some(ref tree) = lower.kdtree {
            let r_sq = inter_layer_radius.powi(2);
            for i in 0..n_u {
                let pos = [points_u[[i, 0]], points_u[[i, 1]]];
                let results = tree.within::<kiddo::SquaredEuclidean>(&pos, r_sq);
                for res in results {
                    let j_big = res.item as usize;
                    let j_small = if lower.pbc { lower.mappings.as_ref().unwrap()[j_big] } else { j_big };
            
                    h_ul.add(i as i32, j_small as i32, internal_types_u[i] as i32, internal_types_l[j_small] as i32);
                    h_lu.add(j_small as i32, i as i32, internal_types_l[j_small] as i32, internal_types_u[i] as i32);
                    let dx = lower_big[[j_big, 0]] - points_l[[j_small, 0]];
                    let dy = lower_big[[j_big, 1]] - points_l[[j_small, 1]];
                    shifts.add(i as i32 + offset, j_small as i32, dx, dy);
                    shifts.add(j_small as i32, i as i32 + offset, -dx, -dy);
                }
            }
        }
        self.calculate_total_nnz();
    }

    /// Assembles the real-valued tight-binding Hamiltonian in COO format.
    ///
    /// Thin wrapper around [`BilayerMoire::build_ham_internal`] accepting `f64` amplitudes.
    ///
    /// # Arguments
    /// * `tll`: intra-layer hopping amplitude(s) for the lower layer.
    /// * `tuu`: intra-layer hopping amplitude(s) for the upper layer.
    /// * `tlu`: inter-layer hopping amplitude(s) from lower to upper.
    /// * `tul`: inter-layer hopping amplitude(s) from upper to lower.
    /// * `tuself`: on-site energy for upper-layer sites.
    /// * `tlself`: on-site energy for lower-layer sites.
    ///
    /// Each argument can be a scalar `float` or a list/array of per-bond values.
    ///
    /// # Returns
    /// `(data, rows, cols)` as 1-D NumPy arrays, ready for `scipy.sparse.coo_matrix`.
    ///
    /// # Errors
    /// Returns `PyRuntimeError` if `generate_connections` was not called first.
    pub fn build_ham<'py>(
        &self,
        py: Python<'py>,
        tll: Input<f64>,
        tuu: Input<f64>,
        tlu: Input<f64>,
        tul: Input<f64>,
        tuself: Input<f64>,
        tlself: Input<f64>,
    ) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<i32>>, Bound<'py, PyArray1<i32>>)> {
        self.build_ham_internal(py, tll, tuu, tlu, tul, tuself, tlself)
    }

    /// Assembles the complex-valued tight-binding Hamiltonian in COO format.
    ///
    /// Same as [`BilayerMoire::build_ham`] but accepts `Complex<f64>` amplitudes,
    /// enabling Hamiltonians with complex hopping phases (e.g. Haldane model,
    /// magnetic flux, spin-orbit coupling).
    ///
    /// # Returns
    /// `(data, rows, cols)` where `data` is a `Complex<f64>` array.
    ///
    /// # Errors
    /// Same as [`BilayerMoire::build_ham`].
    pub fn build_ham_complex<'py>(
        &self,
        py: Python<'py>,
        tll: Input<Complex<f64>>,
        tuu: Input<Complex<f64>>,
        tlu: Input<Complex<f64>>,
        tul: Input<Complex<f64>>,
        tuself: Input<Complex<f64>>,
        tlself: Input<Complex<f64>>,
    ) -> PyResult<(Bound<'py, PyArray1<Complex<f64>>>, Bound<'py, PyArray1<i32>>, Bound<'py, PyArray1<i32>>)> {
        self.build_ham_internal(py, tll, tuu, tlu, tul, tuself, tlself)
    }

    /// Builds COO phase data from cached shifts for a given k-vector.
    ///
    /// For each cached entry `(row, col, dx, dy)`, writes
    /// `data = exp(-i * (dx * kx + dy * ky))` into the COO output.
    ///
    /// # Errors
    /// Returns `PyRuntimeError` if `generate_connections` has not been called.
    pub fn get_phase<'py>(
        &self,
        py: Python<'py>,
        kx: f64,
        ky: f64,
    ) -> PyResult<(Bound<'py, PyArray1<Complex<f64>>>, Bound<'py, PyArray1<i32>>, Bound<'py, PyArray1<i32>>)> {
        let shifts = self
            .shifts
            .as_ref()
            .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>("shifts missing. Run generate_connections first"))?;

        let k_orb = self.orbitals as i32;
        let mut phase = COOBuilder::<Complex<f64>>::new(shifts.rows.len() * (k_orb * k_orb) as usize);
        for (((&row, &col), &dx), &dy) in shifts
            .rows
            .iter()
            .zip(shifts.cols.iter())
            .zip(shifts.dx.iter())
            .zip(shifts.dy.iter())
        {
            let theta = dx * kx + dy * ky;
            let (sin_t, cos_t) = theta.sin_cos();
            let val = Complex::new(cos_t, -sin_t);
            for r_orb in 0..k_orb {
                for c_orb in 0..k_orb {
                    phase.add(row * k_orb + r_orb, col * k_orb + c_orb, val);
                }
            }
        }

        Ok((
            phase.data.into_pyarray(py),
            phase.rows.into_pyarray(py),
            phase.cols.into_pyarray(py),
        ))
    }



    // -------------------------------------------------------------------------
    // Read-only getters for the ndarray fields that can't use #[pyo3(get)].
    // -------------------------------------------------------------------------

    /// Returns the first Moiré supercell vector `mlv1` as a 1-D NumPy array `(2,)`.
    #[getter]
    fn mlv1<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<f64>> {
        use numpy::ToPyArray;
        self.mlv1.to_pyarray(py)
    }

    /// Returns the second Moiré supercell vector `mlv2` as a 1-D NumPy array `(2,)`.
    #[getter]
    fn mlv2<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<f64>> {
        use numpy::ToPyArray;
        self.mlv2.to_pyarray(py)
    }

    /// Returns the stacking translation vector `[dx, dy]` as a 1-D NumPy array `(2,)`.
    #[getter]
    fn translate_upper<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<f64>> {
        use numpy::ToPyArray;
        self.translate_upper.to_pyarray(py)
    }


    /// Returns the cached intra-layer lower hopping instruction, or `None`
    /// if `generate_connections` has not been called yet.
    #[getter]
    fn hop_ll(&self) -> Option<HoppingInstruction> {
        self.hop_ll.as_ref().map(|h| HoppingInstruction {
            site_i: h.site_i.clone(),
            site_j: h.site_j.clone(),
            type1: h.type1.clone(),
            type2: h.type2.clone(),
        }) //
    }

    /// Returns the cached intra-layer upper hopping instruction, or `None`.
    #[getter]
    fn hop_uu(&self) -> Option<HoppingInstruction> {
        self.hop_uu.as_ref().map(|h| HoppingInstruction {
            site_i: h.site_i.clone(),
            site_j: h.site_j.clone(),
            type1: h.type1.clone(),
            type2: h.type2.clone(),
        }) //
    }

    /// Returns the cached inter-layer upper -> lower hopping instruction, or `None`.
    #[getter]
    fn hop_ul(&self) -> Option<HoppingInstruction> {
        self.hop_ul.as_ref().map(|h| HoppingInstruction {
            site_i: h.site_i.clone(),
            site_j: h.site_j.clone(),
            type1: h.type1.clone(),
            type2: h.type2.clone(),
        }) //
    }

    /// Returns the cached inter-layer lower -> upper hopping instruction, or `None`.
    #[getter]
    fn hop_lu(&self) -> Option<HoppingInstruction> {
        self.hop_lu.as_ref().map(|h| HoppingInstruction {
            site_i: h.site_i.clone(),
            site_j: h.site_j.clone(),
            type1: h.type1.clone(),
            type2: h.type2.clone(),
        }) //
    }

    /// Returns the cached self-energy instruction for lower-layer sites, or `None`.
    #[getter]
    fn hop_l(&self) -> Option<SelfInstruction> {
        self.hop_l.as_ref().map(|h| SelfInstruction {
            site_i: h.site_i.clone(),
            ptype: h.ptype.clone(),
        }) //
    }

    /// Returns the cached self-energy instruction for upper-layer sites, or `None`.
    #[getter]
    fn hop_u(&self) -> Option<SelfInstruction> {
        self.hop_u.as_ref().map(|h| SelfInstruction {
            site_i: h.site_i.clone(),
            ptype: h.ptype.clone(),
        }) //
    }

}


/// Private Rust-only helpers, not exported to Python.
impl BilayerMoire {
    /// Computes and caches the **total number of non-zero** Hamiltonian entries.
    ///
    /// Called automatically at the end of [`BilayerMoire::generate_connections`].
    /// The count is a conservative upper bound: every instruction entry is
    /// assumed to contribute a full `k^2 = orbitals^2` block regardless of
    /// whether the amplitude is zero.  This avoids reallocation during the
    /// COO assembly loop.
    ///
    /// The result is stored in `self.nnz` and used to pre-allocate the
    /// `COOBuilder` inside [`BilayerMoire::build_ham_internal`].
    fn calculate_total_nnz(&mut self) {
        let k_sq = (self.orbitals * self.orbitals) as usize; // k^2 entries per instruction
        let mut count = 0;

        // We assume worst-case (all blocks are full) for safety and speed
        if let Some(ref h) = self.hop_l { count += h.site_i.len() * k_sq; }
        if let Some(ref h) = self.hop_u { count += h.site_i.len() * k_sq; }
        if let Some(ref h) = self.hop_ll { count += h.site_i.len() * k_sq; }
        if let Some(ref h) = self.hop_uu { count += h.site_i.len() * k_sq; }
        if let Some(ref h) = self.hop_lu { count += h.site_i.len() * k_sq; }
        if let Some(ref h) = self.hop_ul { count += h.site_i.len() * k_sq; }

        self.nnz = Some(count);
    }


    /// Generic Hamiltonian assembly in COO (coordinate) sparse format.
    ///
    /// This is the single implementation behind both [`BilayerMoire::build_ham`]
    /// (real `f64`) and [`BilayerMoire::build_ham_complex`] (`Complex<f64>`).
    /// Generics let us avoid code duplication while keeping Rust's type safety.
    ///
    /// ## Assembly order
    /// The COO triplets are added in this order (important for understanding
    /// the row/column layout of the sparse matrix):
    /// 1. Lower self-energies (rows/cols in `[0, n_lower * k)`).
    /// 2. Upper self-energies (rows/cols in `[n_lower * k, (n_lower + n_upper) * k)`).
    /// 3. Lower intra-layer hoppings (LL block, top-left quadrant).
    /// 4. Upper intra-layer hoppings (UU block, bottom-right quadrant).
    /// 5. Lower -> Upper inter-layer hoppings (LU block, top-right quadrant).
    /// 6. Upper -> Lower inter-layer hoppings (UL block, bottom-left quadrant).
    ///
    /// ## Orbital blocking
    /// If `orbitals = k > 1`, each site expands into a `k x k` sub-block.
    /// Scalar amplitudes fill the full block uniformly; vector amplitudes
    /// must provide `k^2` values per pair.
    ///
    /// # Type parameters
    /// * `T`: numeric type of the amplitude (`f64` or `Complex<f64>`).
    /// * `S1..S6`: concrete `Pushable<T>` types inferred from arguments.
    ///
    /// # Returns
    /// `(data, rows, cols)` as 1-D PyO3-bound NumPy arrays.
    ///
    /// # Errors
    /// Returns `PyRuntimeError` if any of the six hopping instructions is `None`
    /// (meaning `generate_connections` was not called first).
    fn build_ham_internal<'py, T, S1, S2, S3, S4, S5, S6>(
        &self,
        py: Python<'py>,
        tll: S1,
        tuu: S2,
        tlu: S3,
        tul: S4,
        tuself: S5,
        tlself: S6,
    ) -> PyResult<(Bound<'py, PyArray1<T>>, Bound<'py, PyArray1<i32>>, Bound<'py, PyArray1<i32>>)>
    where
        T: Copy + Element,
        S1: Pushable<T>, S2: Pushable<T>, S3: Pushable<T>,
        S4: Pushable<T>, S5: Pushable<T>, S6: Pushable<T>,
    {
        let h_l = self.hop_l.as_ref().ok_or_else(|| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>("hop_l missing. Run generate_connections first"))?;
        let h_u = self.hop_u.as_ref().ok_or_else(|| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>("hop_u missing. Run generate_connections first"))?;
        let h_ll = self.hop_ll.as_ref().ok_or_else(|| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>("hop_ll missing. Run generate_connections first"))?;
        let h_uu = self.hop_uu.as_ref().ok_or_else(|| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>("hop_uu missing. Run generate_connections first"))?;
        let h_lu = self.hop_lu.as_ref().ok_or_else(|| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>("hop_lu missing. Run generate_connections first"))?;
        let h_ul = self.hop_ul.as_ref().ok_or_else(|| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>("hop_ul missing. Run generate_connections first"))?;

        let offset = {
            let lower = self.lower_lattice.bind(py).borrow();
            (lower.points.as_ref().map_or(0, |p| p.nrows()) as i32) * (self.orbitals as i32)
        };
        let k_orb = self.orbitals as i32;

        let mut builder = COOBuilder::<T>::new(self.nnz.unwrap_or(64));

        // The logic is now "attached" to the data.
        // This calls either the Scalar loop or the Vector loop automatically.
        tlself.push_self(&mut builder, h_l, 0, k_orb);
        tuself.push_self(&mut builder, h_u, offset, k_orb);

        tll.push_hop(&mut builder, h_ll, 0, 0, k_orb);
        tuu.push_hop(&mut builder, h_uu, offset, offset, k_orb);

        tlu.push_hop(&mut builder, h_lu, 0, offset, k_orb);
        tul.push_hop(&mut builder, h_ul, offset, 0, k_orb);

        Ok((
            builder.data.into_pyarray(py),
            builder.rows.into_pyarray(py),
            builder.cols.into_pyarray(py),
        ))
    }


}
