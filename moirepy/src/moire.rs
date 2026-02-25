use num_complex::Complex;
use ndarray::{Array1};
use pyo3::prelude::*;
use numpy::{PyReadonlyArray1,  IntoPyArray, Element, PyArray1};
use crate::layers::Layer;
use crate::utils::{COOBuilder, HoppingInstruction, SelfInstruction};

trait Pushable<T> {
    fn push_self(&self, builder: &mut COOBuilder<T>, inst: &SelfInstruction, off: i32, k: i32);
    fn push_hop(&self, builder: &mut COOBuilder<T>, inst: &HoppingInstruction, r_off: i32, c_off: i32, k: i32);
}

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


#[derive(FromPyObject)]
pub enum Input<T> {
    Scalar(T),      // This variant catches a single value (f64 or Complex)
    Vector(Vec<T>), // This variant catches a Python list or NumPy array
}

// Now we implement the Pushable trait FOR the Input enum itself.
// This is where the routing happens.
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



#[pyclass]
pub struct BilayerMoire {
    // --- 1. Layer Ownership (Shared with Python) ---
    // Using Py<Layer> allows Rust to point to the same memory Python uses.
    pub lower_lattice: Py<Layer>,
    pub upper_lattice: Py<Layer>,

    // --- 2. Geometric Coefficients (AVC Tool) ---
    #[pyo3(get)]
    pub ll1: i32,
    #[pyo3(get)]
    pub ll2: i32,
    #[pyo3(get)]
    pub ul1: i32,
    #[pyo3(get)]
    pub ul2: i32,

    // --- 3. Supercell Parameters ---
    #[pyo3(get)]
    pub n1: usize,
    #[pyo3(get)]
    pub n2: usize,
    #[pyo3(get)]
    pub theta: f64,
    // #[getter]
    pub mlv1: Array1<f64>,
    // #[getter]
    pub mlv2: Array1<f64>,
    // #[getter]
    pub translate_upper: Array1<f64>,

    // --- 4. Global Simulation State ---
    #[pyo3(get)]
    pub pbc: bool,
    #[pyo3(get)]
    pub orbitals: usize,

    // --- 5. Hopping Instructions ---
    // #[getter] for all the 6 below
    pub hop_ll: Option<HoppingInstruction>,
    pub hop_uu: Option<HoppingInstruction>,
    pub hop_lu: Option<HoppingInstruction>,
    pub hop_ul: Option<HoppingInstruction>,
    pub hop_l: Option<SelfInstruction>,
    pub hop_u: Option<SelfInstruction>,

    pub nnz: Option<usize>,
}



#[pymethods]
impl BilayerMoire {
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
        }
    }

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

        let h_ll = self.hop_ll.as_mut().unwrap();
        let h_uu = self.hop_uu.as_mut().unwrap();
        let h_ul = self.hop_ul.as_mut().unwrap();
        let h_lu = self.hop_lu.as_mut().unwrap();
        let h_l = self.hop_l.as_mut().unwrap();
        let h_u = self.hop_u.as_mut().unwrap();

        // 1. Self-Energies
        for i in 0..n_l { h_l.add(i as i32, internal_types_l[i] as i32); }
        for i in 0..n_u { h_u.add(i as i32, internal_types_u[i] as i32); }

        // 2. Intra-layer (Lower) - Direct Internal Call
        let (_, small_idx_l, _) = lower.first_nearest_neighbours_internal(
            &points_l.view(), 
            internal_types_l
        ).expect("Lower neighbors failed");

        for i in 0..n_l {
            for &j in small_idx_l.row(i) {
                if j >= 0 {
                    h_ll.add(i as i32, j as i32, internal_types_l[i] as i32, internal_types_l[j as usize] as i32);
                }
            }
        }

        // 3. Intra-layer (Upper) - Direct Internal Call
        let (_, small_idx_u, _) = upper.first_nearest_neighbours_internal(
            &points_u.view(), 
            internal_types_u
        ).expect("Upper neighbors failed");

        for i in 0..n_u {
            for &j in small_idx_u.row(i) {
                if j >= 0 {
                    h_uu.add(i as i32, j as i32, internal_types_u[i] as i32, internal_types_u[j as usize] as i32);
                }
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
                }
            }
        }
        self.calculate_total_nnz();
    }

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



    // --- Read-Only Getters for Python ---
    
    #[getter]
    fn mlv1<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<f64>> {
        use numpy::ToPyArray;
        self.mlv1.to_pyarray(py)
    }

    #[getter]
    fn mlv2<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<f64>> {
        use numpy::ToPyArray;
        self.mlv2.to_pyarray(py)
    }

    #[getter]
    fn translate_upper<'py>(&self, py: Python<'py>) -> Bound<'py, numpy::PyArray1<f64>> {
        use numpy::ToPyArray;
        self.translate_upper.to_pyarray(py)
    }


    #[getter]
    fn hop_ll(&self) -> Option<HoppingInstruction> {
        self.hop_ll.as_ref().map(|h| HoppingInstruction {
            site_i: h.site_i.clone(),
            site_j: h.site_j.clone(),
            type1: h.type1.clone(),
            type2: h.type2.clone(),
        }) //
    }

    #[getter]
    fn hop_uu(&self) -> Option<HoppingInstruction> {
        self.hop_uu.as_ref().map(|h| HoppingInstruction {
            site_i: h.site_i.clone(),
            site_j: h.site_j.clone(),
            type1: h.type1.clone(),
            type2: h.type2.clone(),
        }) //
    }

    #[getter]
    fn hop_ul(&self) -> Option<HoppingInstruction> {
        self.hop_ul.as_ref().map(|h| HoppingInstruction {
            site_i: h.site_i.clone(),
            site_j: h.site_j.clone(),
            type1: h.type1.clone(),
            type2: h.type2.clone(),
        }) //
    }

    #[getter]
    fn hop_lu(&self) -> Option<HoppingInstruction> {
        self.hop_lu.as_ref().map(|h| HoppingInstruction {
            site_i: h.site_i.clone(),
            site_j: h.site_j.clone(),
            type1: h.type1.clone(),
            type2: h.type2.clone(),
        }) //
    }

    #[getter]
    fn hop_l(&self) -> Option<SelfInstruction> {
        self.hop_l.as_ref().map(|h| SelfInstruction {
            site_i: h.site_i.clone(),
            ptype: h.ptype.clone(),
        }) //
    }

    #[getter]
    fn hop_u(&self) -> Option<SelfInstruction> {
        self.hop_u.as_ref().map(|h| SelfInstruction {
            site_i: h.site_i.clone(),
            ptype: h.ptype.clone(),
        }) //
    }
}


impl BilayerMoire {
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
