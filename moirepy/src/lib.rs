use pyo3::prelude::*;

#[pyfunction]
fn hello_world() -> PyResult<()> {
    println!("Hello World from Rust!");
    Ok(())
}

#[pymodule]
fn moirepy_rust(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(hello_world, m)?)?;
    Ok(())
}
