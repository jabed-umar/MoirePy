use std::ops::{Mul, AddAssign};
use num_complex::Complex;
use ndarray::{Array1, array};
use pyo3::prelude::*;
use numpy::{PyReadonlyArray1, ToPyArray, PyArrayMethods, IntoPyArray, Element, PyArray1};
use crate::layers::Layer;
use crate::utils::{COOBuilder, HoppingInstruction, SelfInstruction};


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

    pub fn build_ham_from_scalars<'py>(
        &self,
        py: Python<'py>,
        tll: f64,
        tuu: f64,
        tlu: f64,
        tul: f64,
        tuself: f64,
        tlself: f64,
    ) -> PyResult<(Bound<'py, PyArray1<f64>>, Bound<'py, PyArray1<i32>>, Bound<'py, PyArray1<i32>>)> {
        self.build_ham_from_scalars_internal(py, tll, tuu, tlu, tul, tuself, tlself)
    }

    // Exposes the complex128 version to Python
    pub fn build_ham_from_complex_scalars<'py>(
        &self,
        py: Python<'py>,
        tll: Complex<f64>,
        tuu: Complex<f64>,
        tlu: Complex<f64>,
        tul: Complex<f64>,
        tuself: Complex<f64>,
        tlself: Complex<f64>,
    ) -> PyResult<(Bound<'py, PyArray1<Complex<f64>>>, Bound<'py, PyArray1<i32>>, Bound<'py, PyArray1<i32>>)> {
        self.build_ham_from_scalars_internal(py, tll, tuu, tlu, tul, tuself, tlself)
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
        let k = self.orbitals;
        let mut count = 0;

        if let Some(ref h) = self.hop_l { count += h.site_i.len() * k; }
        if let Some(ref h) = self.hop_u { count += h.site_i.len() * k; }
        if let Some(ref h) = self.hop_ll { count += h.site_i.len() * k; }
        if let Some(ref h) = self.hop_uu { count += h.site_i.len() * k; }
        if let Some(ref h) = self.hop_lu { count += h.site_i.len() * k; }
        if let Some(ref h) = self.hop_ul { count += h.site_i.len() * k; }

        self.nnz = Some(count);
    }

    fn build_ham_from_scalars_internal<'py, T>(
        &self,
        py: Python<'py>,
        tll: T,
        tuu: T,
        tlu: T,
        tul: T,
        tuself: T,
        tlself: T,
    ) -> PyResult<(Bound<'py, PyArray1<T>>, Bound<'py, PyArray1<i32>>, Bound<'py, PyArray1<i32>>)>
   // -> PyResult<(Bound<'py, PyArray1<T>>, Bound<'py, Bound<'py, PyArray1<i32>>, Bound<'py, PyArray1<i32>>)>
    where 
        T: Copy + Element 
    {
        // 1. Pre-flight check: Ensure all instructions exist before starting the work
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

        // 2. High-speed, un-checked pushes
        Self::push_self_data(&mut builder, h_l, tlself, 0, k_orb);
        Self::push_self_data(&mut builder, h_u, tuself, offset, k_orb);
        Self::push_hop_data(&mut builder, h_ll, tll, 0, 0, k_orb);
        Self::push_hop_data(&mut builder, h_uu, tuu, offset, offset, k_orb);
        Self::push_hop_data(&mut builder, h_lu, tlu, 0, offset, k_orb);
        Self::push_hop_data(&mut builder, h_ul, tul, offset, 0, k_orb);
        // 4. Return zero-copy NumPy arrays
        Ok((
            builder.data.into_pyarray(py),
            builder.rows.into_pyarray(py),
            builder.cols.into_pyarray(py),
        ))
    }


    


    fn push_self_data<T>(
        builder: &mut COOBuilder<T>, 
        s: &SelfInstruction, // Changed from &Option<SelfInstruction>
        val: T, 
        off: i32, 
        k: i32
    ) where T: Copy {
        // No 'if let Some' needed anymore. We know it exists.
        for &i in &s.site_i {
            for orb in 0..k {
                let idx = i * k + orb + off;
                builder.add(idx, idx, val);
            }
        }
    }

    fn push_hop_data<T>(
        builder: &mut COOBuilder<T>, 
        h: &HoppingInstruction, // Changed from &Option<HoppingInstruction>
        val: T, 
        row_off: i32, 
        col_off: i32, 
        k: i32
    ) where T: Copy {
        for (&i, &j) in h.site_i.iter().zip(h.site_j.iter()) {
            for orb in 0..k {
                builder.add(i * k + orb + row_off, j * k + orb + col_off, val);
            }
        }
    }
}
