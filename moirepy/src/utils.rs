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

pub struct COOBuilder<T> {
    pub data: Vec<T>,
    pub rows: Vec<i32>,
    pub cols: Vec<i32>,
}

impl<T> COOBuilder<T> {
    pub fn new(capacity: usize) -> Self {
        COOBuilder {
            data: Vec::with_capacity(capacity),
            rows: Vec::with_capacity(capacity),
            cols: Vec::with_capacity(capacity),
        }
    }

    pub fn add(&mut self, row: i32, col: i32, val: T) {
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



// pub fn coo_to_csr<T>(
//     rows: Vec<i32>,
//     cols: Vec<i32>,
//     data: Vec<T>,
//     num_rows: usize,
// ) -> (Vec<T>, Vec<i32>, Vec<i32>) 
// where 
//     T: Copy + AddAssign, 
// {
//     let m = rows.len();
//     if m == 0 {
//         return (Vec::new(), Vec::new(), vec![0; num_rows + 1]);
//     }

//     let mut indices: Vec<usize> = (0..m).collect();
//     indices.sort_unstable_by(|&a, &b| {
//         rows[a].cmp(&rows[b]).then(cols[a].cmp(&cols[b]))
//     });

//     let mut csr_data = Vec::with_capacity(m);
//     let mut csr_indices = Vec::with_capacity(m);
//     let mut csr_indptr = vec![0i32; num_rows + 1];

//     let first_idx = indices[0];
//     let mut last_row = rows[first_idx] as usize;
//     let mut last_col = cols[first_idx] as usize;
//     let mut current_sum = data[first_idx];

//     // CRITICAL: Initialize indptr for any empty rows at the very beginning
//     for j in 0..=last_row {
//         csr_indptr[j] = 0;
//     }

//     for &idx in indices.iter().skip(1) {
//         let r = rows[idx] as usize;
//         let c = cols[idx] as usize;

//         if r == last_row && c == last_col {
//             current_sum += data[idx];
//         } else {
//             csr_data.push(current_sum);
//             csr_indices.push(last_col as i32);
            
//             if r != last_row {
//                 for j in (last_row + 1)..=r {
//                     csr_indptr[j] = csr_indices.len() as i32;
//                 }
//             }
//             current_sum = data[idx];
//             last_row = r;
//             last_col = c;
//         }
//     }
    
//     csr_data.push(current_sum);
//     csr_indices.push(last_col as i32);

//     for j in (last_row + 1)..=num_rows {
//         csr_indptr[j] = csr_indices.len() as i32;
//     }

//     (csr_data, csr_indices, csr_indptr)
// }