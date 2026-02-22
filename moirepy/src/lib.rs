pub mod utils;
pub mod layers;
pub mod moire;

use pyo3::prelude::*;

#[pymodule]
fn moirepy_rust(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<layers::Layer>()?;
    Ok(())
}
