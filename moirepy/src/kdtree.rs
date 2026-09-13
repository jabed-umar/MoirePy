use pyo3::prelude::*;
use numpy::{PyReadonlyArray1, PyReadonlyArray2, PyArray1, PyArray2, IntoPyArray};
use kiddo::{ImmutableKdTree, SquaredEuclidean};
use rayon::prelude::*;

/// Finds first nearest neighbours within an intra-layer point set.
///
/// Returns (primary_indices, bigger_indices) as 1D int64 NumPy arrays.
#[pyfunction]
#[pyo3(signature = (points, bigger_points, point_type_indices, deltas_per_type, dist_cutoff=1e-5))]
pub fn find_first_nearest_neighbors<'py>(
    py: Python<'py>,
    points: PyReadonlyArray2<'py, f64>,
    bigger_points: PyReadonlyArray2<'py, f64>,
    point_type_indices: PyReadonlyArray1<'py, i64>,
    deltas_per_type: Vec<PyReadonlyArray2<'py, f64>>,
    dist_cutoff: f64,
) -> PyResult<(Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<i64>>)> {
    let pts_slice = points.as_slice()?;
    let big_slice = bigger_points.as_slice()?;
    let type_ids = point_type_indices.as_slice()?;

    let n_pts = pts_slice.len() / 3;
    let n_big = big_slice.len() / 3;

    if n_pts == 0 || n_big == 0 {
        return Ok((
            Vec::<i64>::new().into_pyarray(py),
            Vec::<i64>::new().into_pyarray(py),
        ));
    }

    let pts_3d: &[[f64; 3]] = unsafe {
        std::slice::from_raw_parts(pts_slice.as_ptr() as *const [f64; 3], n_pts)
    };
    let big_3d: &[[f64; 3]] = unsafe {
        std::slice::from_raw_parts(big_slice.as_ptr() as *const [f64; 3], n_big)
    };

    let deltas_by_type: Vec<Vec<[f64; 3]>> = deltas_per_type
        .iter()
        .map(|d| {
            let s = d.as_slice().unwrap_or(&[]);
            let count = s.len() / 3;
            if count == 0 {
                Vec::new()
            } else {
                unsafe {
                    std::slice::from_raw_parts(s.as_ptr() as *const [f64; 3], count).to_vec()
                }
            }
        })
        .collect();

    let tree: ImmutableKdTree<f64, 3> = ImmutableKdTree::new_from_slice(big_3d);
    let cutoff_sq = dist_cutoff * dist_cutoff;

    let chunk_size = 512;
    let chunks: Vec<(usize, &[[f64; 3]])> = pts_3d
        .chunks(chunk_size)
        .enumerate()
        .map(|(c_idx, c)| (c_idx * chunk_size, c))
        .collect();

    let thread_results: Vec<(Vec<i64>, Vec<i64>)> = chunks
        .into_par_iter()
        .map(|(start_idx, chunk)| {
            let mut local_i = Vec::with_capacity(chunk.len() * 4);
            let mut local_j = Vec::with_capacity(chunk.len() * 4);

            for (offset, pt) in chunk.iter().enumerate() {
                let p_idx = start_idx + offset;
                if p_idx >= type_ids.len() {
                    continue;
                }
                let type_id = type_ids[p_idx];
                if type_id < 0 || (type_id as usize) >= deltas_by_type.len() {
                    continue;
                }
                let deltas = &deltas_by_type[type_id as usize];
                for delta in deltas {
                    let target = [pt[0] + delta[0], pt[1] + delta[1], pt[2] + delta[2]];
                    let nearest = tree.nearest_one::<SquaredEuclidean>(&target);
                    if nearest.distance <= cutoff_sq {
                        local_i.push(p_idx as i64);
                        local_j.push(nearest.item as i64);
                    }
                }
            }
            (local_i, local_j)
        })
        .collect();

    let total_matches: usize = thread_results.iter().map(|(i, _)| i.len()).sum();
    let mut all_i = Vec::with_capacity(total_matches);
    let mut all_j = Vec::with_capacity(total_matches);

    for (i_v, j_v) in thread_results {
        all_i.extend(i_v);
        all_j.extend(j_v);
    }

    Ok((all_i.into_pyarray(py), all_j.into_pyarray(py)))
}

/// Finds all neighbours within radius between two layers.
///
/// If primary_pts_j and bigger_to_primary are provided:
///   Returns (c_i, c_j, R_shifts) where c_j is the primary cell index in layer j,
///   and R_shifts is the supercell offset vector (b_pt_j - pt_j).
/// Else:
///   Returns (c_i, b_j, R_shifts) where b_j is the index into bigger_pts_j,
///   and R_shifts is (b_pt_j - pt_i).
#[pyfunction]
#[pyo3(signature = (pts_i, bigger_pts_j, primary_pts_j=None, bigger_to_primary=None, radius=0.0))]
pub fn find_neighbors_within_radius<'py>(
    py: Python<'py>,
    pts_i: PyReadonlyArray2<'py, f64>,
    bigger_pts_j: PyReadonlyArray2<'py, f64>,
    primary_pts_j: Option<PyReadonlyArray2<'py, f64>>,
    bigger_to_primary: Option<PyReadonlyArray1<'py, i64>>,
    radius: f64,
) -> PyResult<(Bound<'py, PyArray1<i64>>, Bound<'py, PyArray1<i64>>, Bound<'py, PyArray2<f64>>)> {
    let s_i = pts_i.as_slice()?;
    let s_bj = bigger_pts_j.as_slice()?;

    let n_i = s_i.len() / 3;
    let n_bj = s_bj.len() / 3;

    if n_i == 0 || n_bj == 0 || radius <= 0.0 {
        let empty_i = Vec::<i64>::new().into_pyarray(py);
        let empty_j = Vec::<i64>::new().into_pyarray(py);
        let empty_r = ndarray::Array2::<f64>::zeros((0, 3)).into_pyarray(py);
        return Ok((empty_i, empty_j, empty_r));
    }

    let pts_i_3d: &[[f64; 3]] = unsafe {
        std::slice::from_raw_parts(s_i.as_ptr() as *const [f64; 3], n_i)
    };
    let bigger_j_3d: &[[f64; 3]] = unsafe {
        std::slice::from_raw_parts(s_bj.as_ptr() as *const [f64; 3], n_bj)
    };

    let tree: ImmutableKdTree<f64, 3> = ImmutableKdTree::new_from_slice(bigger_j_3d);
    let r_sq = radius * radius;

    let chunk_size = 512;
    let chunks: Vec<(usize, &[[f64; 3]])> = pts_i_3d
        .chunks(chunk_size)
        .enumerate()
        .map(|(c_idx, c)| (c_idx * chunk_size, c))
        .collect();

    let thread_results: Vec<(Vec<i64>, Vec<i64>, Vec<[f64; 3]>)> = if let (Some(prim_pts), Some(b2p_arr)) = (primary_pts_j, bigger_to_primary) {
        let s_pj = prim_pts.as_slice()?;
        let n_pj = s_pj.len() / 3;
        let prim_j_3d: &[[f64; 3]] = unsafe {
            std::slice::from_raw_parts(s_pj.as_ptr() as *const [f64; 3], n_pj)
        };
        let b2p = b2p_arr.as_slice()?;

        chunks
            .into_par_iter()
            .map(|(start_idx, chunk)| {
                let mut local_i = Vec::new();
                let mut local_j = Vec::new();
                let mut local_r = Vec::new();

                for (offset, pt) in chunk.iter().enumerate() {
                    let idx_i = (start_idx + offset) as i64;
                    let matches = tree.within_unsorted::<SquaredEuclidean>(pt, r_sq);
                    for m in matches {
                        let b_idx = m.item as usize;
                        if b_idx < b2p.len() {
                            let prim_j = b2p[b_idx];
                            if (prim_j as usize) < prim_j_3d.len() {
                                let b_pt = bigger_j_3d[b_idx];
                                let p_pt = prim_j_3d[prim_j as usize];
                                let shift = [b_pt[0] - p_pt[0], b_pt[1] - p_pt[1], b_pt[2] - p_pt[2]];

                                local_i.push(idx_i);
                                local_j.push(prim_j);
                                local_r.push(shift);
                            }
                        }
                    }
                }
                (local_i, local_j, local_r)
            })
            .collect()
    } else {
        chunks
            .into_par_iter()
            .map(|(start_idx, chunk)| {
                let mut local_i = Vec::new();
                let mut local_j = Vec::new();
                let mut local_r = Vec::new();

                for (offset, pt) in chunk.iter().enumerate() {
                    let idx_i = (start_idx + offset) as i64;
                    let matches = tree.within_unsorted::<SquaredEuclidean>(pt, r_sq);
                    for m in matches {
                        let b_idx = m.item as usize;
                        let b_pt = bigger_j_3d[b_idx];
                        let shift = [b_pt[0] - pt[0], b_pt[1] - pt[1], b_pt[2] - pt[2]];

                        local_i.push(idx_i);
                        local_j.push(b_idx as i64);
                        local_r.push(shift);
                    }
                }
                (local_i, local_j, local_r)
            })
            .collect()
    };

    let total: usize = thread_results.iter().map(|(i, _, _)| i.len()).sum();
    let mut out_i = Vec::with_capacity(total);
    let mut out_j = Vec::with_capacity(total);
    let mut out_r = Vec::with_capacity(total * 3);

    for (i_v, j_v, r_v) in thread_results {
        out_i.extend(i_v);
        out_j.extend(j_v);
        for shift in r_v {
            out_r.push(shift[0]);
            out_r.push(shift[1]);
            out_r.push(shift[2]);
        }
    }

    let r_arr = ndarray::Array2::from_shape_vec((total, 3), out_r)
        .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

    Ok((
        out_i.into_pyarray(py),
        out_j.into_pyarray(py),
        r_arr.into_pyarray(py),
    ))
}
