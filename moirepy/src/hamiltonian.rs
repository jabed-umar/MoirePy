use pyo3::prelude::*;
use pyo3::types::PyTuple;
use numpy::{PyReadonlyArray1, PyReadonlyArray2, IntoPyArray};
use num_complex::Complex;

enum BlockValues {
    Real(Vec<f64>),
    Complex(Vec<Complex<f64>>),
}

struct EvaluatedBlock<'py> {
    rows: PyReadonlyArray1<'py, i64>,
    cols: PyReadonlyArray1<'py, i64>,
    r_shifts: PyReadonlyArray2<'py, f64>,
    values: BlockValues,
    row_offset: i64,
    col_offset: i64,
    is_conj: bool,
    n_entries: usize,
}

fn extract_block_values<'py>(
    py: Python<'py>,
    val_obj: Bound<'py, PyAny>,
    expected_len: usize,
    is_conj: bool,
) -> PyResult<BlockValues> {
    // 1. Try PyReadonlyArray1<f64>
    if let Ok(arr) = val_obj.extract::<PyReadonlyArray1<f64>>() {
        let s = arr.as_slice()?;
        let v = s.to_vec();
        return Ok(BlockValues::Real(v));
    }
    // 2. Try PyReadonlyArray1<Complex<f64>>
    if let Ok(arr) = val_obj.extract::<PyReadonlyArray1<Complex<f64>>>() {
        let s = arr.as_slice()?;
        let mut v = s.to_vec();
        if is_conj {
            for c in v.iter_mut() {
                *c = c.conj();
            }
        }
        return Ok(BlockValues::Complex(v));
    }
    // 3. Try scalar f64
    if let Ok(s) = val_obj.extract::<f64>() {
        return Ok(BlockValues::Real(vec![s; expected_len]));
    }
    // 4. Try scalar Complex<f64>
    if let Ok(mut c) = val_obj.extract::<Complex<f64>>() {
        if is_conj {
            c = c.conj();
        }
        return Ok(BlockValues::Complex(vec![c; expected_len]));
    }
    // 5. Fallback via np.asarray
    let np = py.import("numpy")?;
    let np_arr = np.getattr("asarray")?.call1((&val_obj,))?;
    if let Ok(arr) = np_arr.extract::<PyReadonlyArray1<f64>>() {
        let s = arr.as_slice()?;
        let v = s.to_vec();
        return Ok(BlockValues::Real(v));
    }
    if let Ok(arr) = np_arr.extract::<PyReadonlyArray1<Complex<f64>>>() {
        let s = arr.as_slice()?;
        let mut v = s.to_vec();
        if is_conj {
            for c in v.iter_mut() {
                *c = c.conj();
            }
        }
        return Ok(BlockValues::Complex(v));
    }

    Err(pyo3::exceptions::PyTypeError::new_err(
        "Hopping function returned an unsupported type; expected 1D float or complex array",
    ))
}

/// Assembles Sparse COO triplet arrays across blocks with a single FFI crossing per block.
///
/// blocks: list of tuples:
///   (call_args, py_func, c_i, c_j, r_shifts, row_offset, col_offset, is_conj)
///
/// Returns:
///   (H_data, H_row, H_col, Rx_data, Ry_data, Rz_data)
#[pyfunction]
pub fn assemble_hamiltonian_sparse<'py>(
    py: Python<'py>,
    blocks: Vec<Bound<'py, PyTuple>>,
) -> PyResult<(
    PyObject,
    Bound<'py, numpy::PyArray1<i64>>,
    Bound<'py, numpy::PyArray1<i64>>,
    Bound<'py, numpy::PyArray1<f64>>,
    Bound<'py, numpy::PyArray1<f64>>,
    Bound<'py, numpy::PyArray1<f64>>,
)> {
    let mut evaluated_blocks = Vec::with_capacity(blocks.len());

    for block in blocks {
        let call_args = block.get_item(0)?.downcast_into::<PyTuple>()?;
        let py_func = block.get_item(1)?;
        let row_indices = block.get_item(2)?.extract::<PyReadonlyArray1<i64>>()?;
        let col_indices = block.get_item(3)?.extract::<PyReadonlyArray1<i64>>()?;
        let r_shifts = block.get_item(4)?.extract::<PyReadonlyArray2<f64>>()?;
        let row_offset = block.get_item(5)?.extract::<i64>()?;
        let col_offset = block.get_item(6)?.extract::<i64>()?;
        let is_conj = block.get_item(7)?.extract::<bool>()?;

        let n_entries = row_indices.as_slice()?.len();
        if n_entries == 0 {
            continue;
        }

        // CRITICAL SPEED STRATEGY: Cross FFI boundary ONCE per block!
        let val_obj = py_func.call1(&call_args)?;
        let values = extract_block_values(py, val_obj, n_entries, is_conj)?;

        let v_len = match &values {
            BlockValues::Real(v) => v.len(),
            BlockValues::Complex(v) => v.len(),
        };

        if v_len != n_entries {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "Hopping function returned {} elements, expected {}",
                v_len, n_entries
            )));
        }

        evaluated_blocks.push(EvaluatedBlock {
            rows: row_indices,
            cols: col_indices,
            r_shifts,
            values,
            row_offset,
            col_offset,
            is_conj,
            n_entries,
        });
    }

    let has_complex = evaluated_blocks
        .iter()
        .any(|b| matches!(b.values, BlockValues::Complex(_)));
    let total_nnz: usize = evaluated_blocks.iter().map(|b| b.n_entries).sum();

    let mut h_row: Vec<i64> = Vec::with_capacity(total_nnz);
    let mut h_col: Vec<i64> = Vec::with_capacity(total_nnz);
    let mut rx: Vec<f64> = Vec::with_capacity(total_nnz);
    let mut ry: Vec<f64> = Vec::with_capacity(total_nnz);
    let mut rz: Vec<f64> = Vec::with_capacity(total_nnz);

    for b in &evaluated_blocks {
        let r_slice = b.rows.as_slice()?;
        let c_slice = b.cols.as_slice()?;
        let s_shift = b.r_shifts.as_slice()?;
        let r_off = b.row_offset;
        let c_off = b.col_offset;
        let is_conj = b.is_conj;

        if r_off == 0 && c_off == 0 && !is_conj {
            h_row.extend_from_slice(r_slice);
            h_col.extend_from_slice(c_slice);
            for k in 0..b.n_entries {
                rx.push(s_shift[3 * k]);
                ry.push(s_shift[3 * k + 1]);
                rz.push(s_shift[3 * k + 2]);
            }
        } else {
            for k in 0..b.n_entries {
                h_row.push(r_slice[k] + r_off);
                h_col.push(c_slice[k] + c_off);
                if is_conj {
                    rx.push(-s_shift[3 * k]);
                    ry.push(-s_shift[3 * k + 1]);
                    rz.push(-s_shift[3 * k + 2]);
                } else {
                    rx.push(s_shift[3 * k]);
                    ry.push(s_shift[3 * k + 1]);
                    rz.push(s_shift[3 * k + 2]);
                }
            }
        }
    }

    let h_data_obj = if has_complex {
        let mut h_data: Vec<Complex<f64>> = Vec::with_capacity(total_nnz);
        for b in evaluated_blocks {
            match b.values {
                BlockValues::Complex(v) => h_data.extend(v),
                BlockValues::Real(v) => {
                    for x in v {
                        h_data.push(Complex::new(x, 0.0));
                    }
                }
            }
        }
        h_data.into_pyarray(py).into_any().unbind()
    } else {
        let mut h_data: Vec<f64> = Vec::with_capacity(total_nnz);
        for b in evaluated_blocks {
            match b.values {
                BlockValues::Real(v) => h_data.extend(v),
                BlockValues::Complex(_) => unreachable!(),
            }
        }
        h_data.into_pyarray(py).into_any().unbind()
    };

    Ok((
        h_data_obj,
        h_row.into_pyarray(py),
        h_col.into_pyarray(py),
        rx.into_pyarray(py),
        ry.into_pyarray(py),
        rz.into_pyarray(py),
    ))
}
