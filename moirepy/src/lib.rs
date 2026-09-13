use pyo3::prelude::*;

mod kdtree;
mod hamiltonian;

#[pymodule]
fn moirepy_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(kdtree::find_first_nearest_neighbors, m)?)?;
    m.add_function(wrap_pyfunction!(kdtree::find_neighbors_within_radius, m)?)?;
    m.add_function(wrap_pyfunction!(hamiltonian::assemble_hamiltonian_sparse, m)?)?;
    Ok(())
}
