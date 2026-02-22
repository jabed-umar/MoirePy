use ndarray::{array, Array2, Array1};

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
