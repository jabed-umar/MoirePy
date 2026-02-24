use ndarray::{Array1, array};
use pyo3::prelude::*;
use numpy::{PyReadonlyArray1, ToPyArray, PyArrayMethods};
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

    // --- 6. Hamiltonian Buffers ---
    #[pyo3(get)]
	pub ham_builder: Option<COOBuilder>,
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
            ham_builder: None,
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
    }

    pub fn build_coo_from_scalars(
        &mut self,
        py: Python<'_>,
        tll: f64,
        tuu: f64,
        tlu: f64,
        tul: f64,
        tuself: f64,
        tlself: f64,
    ) {
        // 1. Initialize or clear the persistent COOBuilder
        if self.ham_builder.is_none() {
            self.ham_builder = Some(COOBuilder::new(1024));
        } else {
            self.ham_builder.as_mut().unwrap().clear();
        }

        // 2. Extract metadata in a separate block to release 'lower' borrow
        let (n_l, k) = {
            let lower = self.lower_lattice.bind(py).borrow();
            let n = lower.points.as_ref().map_or(0, |p| p.nrows());
            (n, self.orbitals as i32)
        };
        let offset = (n_l as i32) * k;

        // 3. Obtain mutable reference to the builder
        let builder = self.ham_builder.as_mut().unwrap();

        // 4. Populate with data
        Self::push_self_data(builder, &self.hop_l, tlself, 0, k);
        Self::push_self_data(builder, &self.hop_u, tuself, offset, k);
        
        Self::push_hop_data(builder, &self.hop_ll, tll, 0, 0, k);
        Self::push_hop_data(builder, &self.hop_uu, tuu, offset, offset, k);
        
        Self::push_hop_data(builder, &self.hop_lu, tlu, 0, offset, k);
        Self::push_hop_data(builder, &self.hop_ul, tul, offset, 0, k);
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
    fn push_hop_data(
        builder: &mut COOBuilder,
        h_inst: &Option<HoppingInstruction>,
        val: f64,
        row_off: i32,
        col_off: i32,
        k: i32,
    ) {
        if let Some(ref h) = h_inst {
            // Zip the site indices and iterate
            for (&i, &j) in h.site_i.iter().zip(h.site_j.iter()) {
                for orb in 0..k {
                    builder.add(i * k + orb + row_off, j * k + orb + col_off, val);
                }
            }
        }
    }

    fn push_self_data(
        builder: &mut COOBuilder,
        s_inst: &Option<SelfInstruction>,
        val: f64,
        off: i32,
        k: i32,
    ) {
        if let Some(ref s) = s_inst {
            for &i in &s.site_i {
                for orb in 0..k {
                    builder.add(i * k + orb + off, i * k + orb + off, val);
                }
            }
        }
    }
}
