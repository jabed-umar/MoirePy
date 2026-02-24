use ndarray::{array, Array2, Array1};
use pyo3::prelude::*;

// Returns a 2D Array (Matrix)
pub fn get_rotation_matrix(theta: f64) -> Array2<f64> {
    let (sin_t, cos_t) = theta.sin_cos();
    array![
        [cos_t, -sin_t],
        [sin_t,  cos_t]
    ]
}

// Accepts 1D Arrays (Vectors)
pub fn are_coeffs_integers(v1: Array1<f64>, v2: Array1<f64>, v3: Array1<f64>, tol: f64) -> bool {
    let det = v1[0] * v2[1] - v1[1] * v2[0];
    
    if det.abs() < tol {
        return false;
    }

    let a = (v3[0] * v2[1] - v3[1] * v2[0]) / det;
    let b = (v1[0] * v3[1] - v1[1] * v3[0]) / det;

    (a - a.round()).abs() < tol && (b - b.round()).abs() < tol
}

#[pyclass]
#[derive(Clone)]
pub struct COOBuilder {
    #[pyo3(get)]
    pub data: Vec<f64>,
    #[pyo3(get)]
    pub rows: Vec<i32>,
    #[pyo3(get)]
    pub cols: Vec<i32>,
}

impl COOBuilder {
    pub fn new(capacity: usize) -> Self {
        COOBuilder {
            data: Vec::with_capacity(capacity),
            rows: Vec::with_capacity(capacity),
            cols: Vec::with_capacity(capacity),
        }
    }

    pub fn add(&mut self, row: i32, col: i32, val: f64) {
        self.data.push(val);
        self.rows.push(row);
        self.cols.push(col);
    }

    pub fn clear(&mut self) {
        self.data.clear();
        self.rows.clear();
        self.cols.clear();
    }
}




// In utils.rs
#[pyclass]
pub struct HoppingInstruction {
    #[pyo3(get)]
    pub site_i: Vec<i32>,
    #[pyo3(get)]
    pub site_j: Vec<i32>,
    #[pyo3(get)]
    pub type1: Vec<i32>,
    #[pyo3(get)]
    pub type2: Vec<i32>,
}



impl HoppingInstruction {
    pub fn new(n: usize) -> Self {
        let capacity = 4 * n;
        HoppingInstruction {
            site_i: Vec::with_capacity(capacity),
            site_j: Vec::with_capacity(capacity),
            type1: Vec::with_capacity(capacity),
            type2: Vec::with_capacity(capacity),
        }
    }

    pub fn add(&mut self, i: i32, j: i32, t1: i32, t2: i32) {
        self.site_i.push(i);
        self.site_j.push(j);
        self.type1.push(t1);
        self.type2.push(t2);
    }
    
    pub fn clear(&mut self) {
        self.site_i.clear();
        self.site_j.clear();
        self.type1.clear();
        self.type2.clear();
    }
}

#[pyclass]
pub struct SelfInstruction {
    #[pyo3(get)]
    pub site_i: Vec<i32>,
    #[pyo3(get)]
    pub ptype: Vec<i32>,
}

impl SelfInstruction {
    pub fn new(n: usize) -> Self {
        let capacity = 4 * n;
        SelfInstruction {
            site_i: Vec::with_capacity(capacity),
            ptype: Vec::with_capacity(capacity),
        }
    }

    pub fn add(&mut self, i: i32, ptype: i32) {
        self.site_i.push(i);
        self.ptype.push(ptype);
    }
    
    pub fn clear(&mut self) {
        self.site_i.clear();
        self.ptype.clear();
    }
}
