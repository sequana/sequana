//! Rust accelerators for sequana_coverage, exposed to Python via PyO3.
mod bed;
mod csvfmt;
mod em;
mod roi;

use numpy::{IntoPyArray, PyReadonlyArray1};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

/// merge_rois(pos, cov, rm, zscore, low, high, step=1) -> dict of numpy arrays
///
/// Columns: start, end, size, mean_cov, mean_rm, mean_zscore, log2_ratio,
/// max_zscore, max_cov. The caller adds the `chr` column.
#[pyfunction]
#[pyo3(signature = (pos, cov, rm, zscore, low, high, step=1))]
#[allow(clippy::too_many_arguments)]
fn merge_rois<'py>(
    py: Python<'py>,
    pos: PyReadonlyArray1<'py, i64>,
    cov: PyReadonlyArray1<'py, f64>,
    rm: PyReadonlyArray1<'py, f64>,
    zscore: PyReadonlyArray1<'py, f64>,
    low: f64,
    high: f64,
    step: i64,
) -> PyResult<Bound<'py, PyDict>> {
    let rois = roi::merge_rois(
        pos.as_slice()?,
        cov.as_slice()?,
        rm.as_slice()?,
        zscore.as_slice()?,
        low,
        high,
        step,
    );
    let d = PyDict::new(py);
    d.set_item("start", rois.iter().map(|r| r.start).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("end", rois.iter().map(|r| r.end).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("size", rois.iter().map(|r| r.size).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("mean_cov", rois.iter().map(|r| r.mean_cov).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("mean_rm", rois.iter().map(|r| r.mean_rm).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("mean_zscore", rois.iter().map(|r| r.mean_zscore).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("log2_ratio", rois.iter().map(|r| r.log2_ratio).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("max_zscore", rois.iter().map(|r| r.max_zscore).collect::<Vec<_>>().into_pyarray(py))?;
    d.set_item("max_cov", rois.iter().map(|r| r.max_cov).collect::<Vec<_>>().into_pyarray(py))?;
    Ok(d)
}

/// total_roi_size(pos, zscore, low, high, step=1) -> int
#[pyfunction]
#[pyo3(signature = (pos, zscore, low, high, step=1))]
fn total_roi_size(
    pos: PyReadonlyArray1<'_, i64>,
    zscore: PyReadonlyArray1<'_, f64>,
    low: f64,
    high: f64,
    step: i64,
) -> PyResult<i64> {
    Ok(roi::total_roi_size(pos.as_slice()?, zscore.as_slice()?, low, high, step))
}

/// to_csv(names, columns, precision=3) -> str
///
/// Equivalent of ``DataFrame.to_csv(index=False, float_format="%.<precision>g")``
/// for int64 / float64 numpy columns.
#[pyfunction]
#[pyo3(signature = (names, columns, precision=3))]
fn to_csv(names: Vec<String>, columns: &Bound<'_, PyList>, precision: usize) -> PyResult<String> {
    if names.len() != columns.len() {
        return Err(pyo3::exceptions::PyValueError::new_err("names and columns differ in length"));
    }
    let mut ints: Vec<PyReadonlyArray1<i64>> = Vec::new();
    let mut floats: Vec<PyReadonlyArray1<f64>> = Vec::new();
    // kind[i] = (is_float, index in ints/floats)
    let mut kind = Vec::with_capacity(names.len());
    for item in columns.iter() {
        if let Ok(a) = item.extract::<PyReadonlyArray1<i64>>() {
            kind.push((false, ints.len()));
            ints.push(a);
        } else if let Ok(a) = item.extract::<PyReadonlyArray1<f64>>() {
            kind.push((true, floats.len()));
            floats.push(a);
        } else {
            return Err(pyo3::exceptions::PyTypeError::new_err("columns must be int64 or float64 arrays"));
        }
    }
    let cols: Vec<csvfmt::Column> = kind
        .iter()
        .map(|&(is_float, k)| {
            if is_float {
                csvfmt::Column::Float(floats[k].as_slice().unwrap())
            } else {
                csvfmt::Column::Int(ints[k].as_slice().unwrap())
            }
        })
        .collect();
    Ok(csvfmt::to_csv(&names, &cols, precision))
}

/// em_fit(data, mu, sigma, pi, max_iter=100) -> dict
///
/// Gaussian-mixture EM loop of ``sequana.mixture.EM``. Keys: mu, sigma, pi
/// (final parameters), x (interleaved parameters to report), n_k, counter, ok,
/// gamma (flat k*n array of responsibilities).
#[pyfunction]
#[pyo3(signature = (data, mu, sigma, pi, max_iter=100))]
fn em_fit<'py>(
    py: Python<'py>,
    data: PyReadonlyArray1<'py, f64>,
    mu: Vec<f64>,
    sigma: Vec<f64>,
    pi: Vec<f64>,
    max_iter: usize,
) -> PyResult<Bound<'py, PyDict>> {
    if mu.len() != sigma.len() || mu.len() != pi.len() {
        return Err(pyo3::exceptions::PyValueError::new_err("mu, sigma and pi differ in length"));
    }
    let r = em::em_fit(data.as_slice()?, &mu, &sigma, &pi, max_iter);
    let d = PyDict::new(py);
    d.set_item("mu", r.mu)?;
    d.set_item("sigma", r.sigma)?;
    d.set_item("pi", r.pi)?;
    d.set_item("x", r.x)?;
    d.set_item("n_k", r.n_k)?;
    d.set_item("counter", r.counter)?;
    d.set_item("ok", r.ok)?;
    d.set_item("gamma", r.gamma.into_pyarray(py))?;
    Ok(d)
}

/// neg_log_likelihood(data, x) -> float
///
/// ``-sum(log(pdf(data)))`` of the Gaussian mixture ``x`` (mu, sigma, pi
/// interleaved, pis normalised), as ``GaussianMixtureModel.log_likelihood``.
#[pyfunction]
fn neg_log_likelihood(data: PyReadonlyArray1<'_, f64>, x: Vec<f64>) -> PyResult<f64> {
    if x.len() % 3 != 0 {
        return Err(pyo3::exceptions::PyValueError::new_err("x must be made of (mu, sigma, pi) triplets"));
    }
    Ok(em::neg_log_likelihood(data.as_slice()?, &x))
}

/// scan_bed(path) -> (names, first_rows, last_rows, byte_offsets, total_rows)
#[pyfunction]
fn scan_bed(path: &str) -> PyResult<(Vec<String>, Vec<u64>, Vec<u64>, Vec<u64>, u64)> {
    let s = bed::scan_bed(std::path::Path::new(path))
        .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))?;
    Ok((s.names, s.first, s.last, s.offsets, s.total))
}

/// read_bed(path, offset, nrows) -> dict(pos, cov, cov_is_int, extra, extra_is_int, next_offset)
///
/// Raises ValueError for anything the fast reader does not handle (the caller
/// then falls back to pandas).
#[pyfunction]
fn read_bed<'py>(py: Python<'py>, path: &str, offset: u64, nrows: usize) -> PyResult<Bound<'py, PyDict>> {
    let chunk = py
        .detach(|| bed::read_bed(std::path::Path::new(path), offset, nrows))
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    let d = PyDict::new(py);
    d.set_item("pos", chunk.pos.into_pyarray(py))?;
    d.set_item("cov", chunk.cov.values.into_pyarray(py))?;
    d.set_item("cov_is_int", chunk.cov.is_int)?;
    match chunk.extra {
        Some(e) => {
            d.set_item("extra", e.values.into_pyarray(py))?;
            d.set_item("extra_is_int", e.is_int)?;
        }
        None => {
            d.set_item("extra", py.None())?;
            d.set_item("extra_is_int", false)?;
        }
    }
    d.set_item("next_offset", chunk.next_offset)?;
    Ok(d)
}

#[pymodule]
fn sequana_cov_rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(merge_rois, m)?)?;
    m.add_function(wrap_pyfunction!(total_roi_size, m)?)?;
    m.add_function(wrap_pyfunction!(to_csv, m)?)?;
    m.add_function(wrap_pyfunction!(em_fit, m)?)?;
    m.add_function(wrap_pyfunction!(neg_log_likelihood, m)?)?;
    m.add_function(wrap_pyfunction!(scan_bed, m)?)?;
    m.add_function(wrap_pyfunction!(read_bed, m)?)?;
    Ok(())
}
