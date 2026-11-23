//! Port of `FilteredGenomeCov._merge_region` / `_merge_row` (sequana/bedtools.py).
//!
//! Input arrays are the *filtered* rows (zscore outside the double thresholds),
//! sorted by position. In Python the dataframe index is `pos`, so the label slice
//! `.loc[start:stop]` selects the contiguous rows between the region's first and
//! last row; here we simply keep row indices.

/// One merged region of interest (struct-of-arrays is built by the caller).
#[derive(Debug, Clone, PartialEq)]
pub struct Roi {
    pub start: i64,
    /// exclusive end (stop + 1), as in the Python implementation
    pub end: i64,
    pub size: i64,
    pub mean_cov: f64,
    pub mean_rm: f64,
    pub mean_zscore: f64,
    pub log2_ratio: f64,
    pub max_zscore: f64,
    pub max_cov: f64,
}

fn mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}

fn summarize(pos: &[i64], cov: &[f64], rm: &[f64], z: &[f64], i0: usize, i1: usize) -> Roi {
    let (c, r, zs) = (&cov[i0..=i1], &rm[i0..=i1], &z[i0..=i1]);
    let (start, stop) = (pos[i0], pos[i1]);
    let mean_cov = mean(c);
    let mean_rm = mean(r);
    let mean_zscore = mean(zs);
    let max_cov = c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    // extreme value in the direction of the mean
    let max_zscore = if mean_zscore >= 0.0 {
        zs.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
    } else {
        zs.iter().cloned().fold(f64::INFINITY, f64::min)
    };
    // Python returned None here (missing `return`); NaN is the intended value.
    let log2_ratio = if mean_rm != 0.0 && mean_cov != 0.0 {
        (mean_cov / mean_rm).log2()
    } else {
        f64::NAN
    };
    Roi {
        start,
        end: stop + 1,
        size: stop - start + 1,
        mean_cov,
        mean_rm,
        mean_zscore,
        log2_ratio,
        max_zscore,
        max_cov,
    }
}

/// Scan sorted positions and call `emit(first_row, last_row)` for every region.
///
/// A region is closed when positions are not `step` apart, or when the zscore
/// changes sign (a -5 followed by a +5 are two different regions). Rows are
/// part of a region when `z > high` or `z < low`.
fn scan_regions<F: FnMut(usize, usize)>(pos: &[i64], zscore: &[f64], low: f64, high: f64, step: i64, mut emit: F) {
    let n = pos.len();
    assert!(zscore.len() == n, "length mismatch");

    let mut region: Option<(usize, usize)> = None; // (first row, last row)
    let (mut start, mut stop, mut prev) = (1_i64, 1_i64, 1_i64);
    let mut region_z = 0.0_f64;

    for i in 0..n {
        let z = zscore[i];
        stop = pos[i];
        if stop - step == prev && z * region_z >= 0.0 {
            prev = stop;
        } else {
            if let Some((a, b)) = region.take() {
                emit(a, b);
            }
            start = stop;
            prev = stop;
            region_z = z;
        }

        if (z > 0.0 && z > high) || (z < 0.0 && z < low) {
            region = Some(match region {
                None => (i, i),
                Some((a, _)) => (a, i),
            });
        }
    }
    if start < stop {
        if let Some((a, b)) = region {
            emit(a, b);
        }
    }
}

/// Merge consecutive positions whose zscore exceeds `high` or falls below `low`.
pub fn merge_rois(
    pos: &[i64],
    cov: &[f64],
    rm: &[f64],
    zscore: &[f64],
    low: f64,
    high: f64,
    step: i64,
) -> Vec<Roi> {
    assert!(cov.len() == pos.len() && rm.len() == pos.len(), "length mismatch");
    let mut out = Vec::new();
    scan_regions(pos, zscore, low, high, step, |a, b| out.push(summarize(pos, cov, rm, zscore, a, b)));
    out
}

/// Sum of the sizes (`stop - start + 1`) of all regions; what
/// `get_centralness` needs, without building any statistics.
pub fn total_roi_size(pos: &[i64], zscore: &[f64], low: f64, high: f64, step: i64) -> i64 {
    let mut total = 0;
    scan_regions(pos, zscore, low, high, step, |a, b| total += pos[b] - pos[a] + 1);
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merges_consecutive_and_splits_on_gap() {
        let pos = [10, 11, 12, 20, 21];
        let cov = [1.0, 2.0, 3.0, 10.0, 12.0];
        let rm = [2.0; 5];
        let z = [-5.0, -6.0, -5.0, 5.0, 6.0];
        let r = merge_rois(&pos, &cov, &rm, &z, -4.0, 4.0, 1);
        assert_eq!(r.len(), 2);
        assert_eq!((r[0].start, r[0].end, r[0].size), (10, 13, 3));
        assert_eq!(r[0].mean_cov, 2.0);
        assert_eq!(r[0].max_zscore, -6.0);
        assert_eq!((r[1].start, r[1].end, r[1].size), (20, 22, 2));
        assert_eq!(r[1].max_zscore, 6.0);
    }

    #[test]
    fn sign_flip_splits_adjacent_regions() {
        let pos = [5, 6, 7, 8];
        let cov = [1.0; 4];
        let rm = [1.0; 4];
        let z = [-5.0, -5.0, 5.0, 5.0];
        let r = merge_rois(&pos, &cov, &rm, &z, -4.0, 4.0, 1);
        assert_eq!(r.len(), 2);
        assert_eq!((r[0].start, r[0].size), (5, 2));
        assert_eq!((r[1].start, r[1].size), (7, 2));
    }

    #[test]
    fn zero_rm_gives_nan_log2() {
        let r = merge_rois(&[3, 4], &[1.0, 1.0], &[0.0, 0.0], &[-5.0, -5.0], -4.0, 4.0, 1);
        assert!(r[0].log2_ratio.is_nan());
    }

    #[test]
    fn total_size_matches_regions() {
        let pos = [10, 11, 12, 20, 21, 30];
        let z = [-5.0, -6.0, -5.0, 5.0, 6.0, 0.1];
        assert_eq!(total_roi_size(&pos, &z, -4.0, 4.0, 1), 5);
    }

    #[test]
    fn nothing_below_threshold() {
        let r = merge_rois(&[1, 2, 3], &[1.0; 3], &[1.0; 3], &[0.1, -0.2, 0.3], -4.0, 4.0, 1);
        assert!(r.is_empty());
    }
}
