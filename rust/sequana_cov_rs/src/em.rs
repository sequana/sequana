//! Port of the EM loop of `sequana.mixture.EM.estimate` (Gaussian mixture).
//!
//! Same algorithm and same quirks as the Python code: fixed number of
//! iterations (no convergence test), the denominator of the responsibilities
//! uses the parameters of the previous iteration, and the loop stops early when
//! the sanity checks on `sum(N_k)` and `sum(pi_k)` fail (NaN included).

const SQRT_2PI: f64 = 2.506_628_274_631_000_7;

#[derive(Debug)]
pub struct EmResult {
    pub mu: Vec<f64>,
    pub sigma: Vec<f64>,
    pub pi: Vec<f64>,
    /// parameters (mu, sigma, pi interleaved) at the start of the last iteration
    /// when the loop failed, or at the end when it succeeded
    pub x: Vec<f64>,
    pub n_k: Vec<f64>,
    pub counter: usize,
    pub ok: bool,
    /// responsibilities, row-major (k, n)
    pub gamma: Vec<f64>,
}

/// scipy.stats.norm.pdf(x, loc, scale)
#[inline]
fn norm_pdf(x: f64, loc: f64, scale: f64) -> f64 {
    if scale.is_nan() || scale <= 0.0 {
        return f64::NAN;
    }
    let z = (x - loc) / scale;
    (-z * z / 2.0).exp() / SQRT_2PI / scale
}

fn interleave(mu: &[f64], sigma: &[f64], pi: &[f64]) -> Vec<f64> {
    (0..mu.len()).flat_map(|k| [mu[k], sigma[k], pi[k]]).collect()
}

pub fn em_fit(data: &[f64], mu0: &[f64], sigma0: &[f64], pi0: &[f64], max_iter: usize) -> EmResult {
    let k = mu0.len();
    let n = data.len();
    assert!(sigma0.len() == k && pi0.len() == k, "length mismatch");
    let size = n as f64;

    let (mut mu, mut sigma, mut pi) = (mu0.to_vec(), sigma0.to_vec(), pi0.to_vec());
    let mut gamma = vec![0.0; k * n];
    let mut n_k = vec![0.0; k];
    let mut denom = vec![0.0; n];
    let mut x = interleave(&mu, &sigma, &pi);
    let mut counter = 0;
    let mut ok = true;

    loop {
        // numerators and denominators use the parameters of the previous iteration
        let any_negative_pi = pi.iter().any(|&p| p < 0.0);
        denom.iter_mut().for_each(|d| *d = 0.0);
        for j in 0..k {
            let row = &mut gamma[j * n..(j + 1) * n];
            for (i, g) in row.iter_mut().enumerate() {
                *g = pi[j] * norm_pdf(data[i], mu[j], sigma[j]);
            }
            // GaussianMixtureModel.pdf skips null sigma and returns 0 if a pi < 0
            if sigma[j] != 0.0 {
                for (d, g) in denom.iter_mut().zip(row.iter()) {
                    *d += *g;
                }
            }
        }
        if any_negative_pi {
            denom.iter_mut().for_each(|d| *d = 0.0);
        }

        for j in 0..k {
            let row = &mut gamma[j * n..(j + 1) * n];
            for (g, d) in row.iter_mut().zip(denom.iter()) {
                *g /= *d;
            }
            let nj: f64 = row.iter().sum();
            n_k[j] = nj;
            mu[j] = row.iter().zip(data).map(|(g, d)| g * d).sum::<f64>() / nj;
            let var: f64 = row.iter().zip(data).map(|(g, d)| g * (d - mu[j]).powi(2)).sum::<f64>() / nj;
            sigma[j] = var.sqrt();
            pi[j] = nj / size;
        }

        let sane = (n_k.iter().sum::<f64>() - size).abs() / size < 1e-6 && (pi.iter().sum::<f64>() - 1.0).abs() < 1e-6;
        if !sane {
            ok = false;
            break; // `x` still holds the parameters before this iteration
        }
        x = interleave(&mu, &sigma, &pi);
        counter += 1;
        if counter >= max_iter {
            break;
        }
    }

    EmResult { mu, sigma, pi, x, n_k, counter, ok, gamma }
}

/// `GaussianMixtureModel.log_likelihood(params, sample)`: -sum(log(pdf(sample))),
/// with the pis normalised to sum to 1. `x` is (mu, sigma, pi) interleaved.
pub fn neg_log_likelihood(data: &[f64], x: &[f64]) -> f64 {
    assert!(x.len() % 3 == 0, "x must be made of (mu, sigma, pi) triplets");
    let k = x.len() / 3;
    let pis: Vec<f64> = (0..k).map(|j| x[3 * j + 2]).collect();
    if pis.iter().any(|&p| p < 0.0) {
        // pdf() returns 0: log(0) = -inf for every sample
        return if data.is_empty() { 0.0 } else { f64::INFINITY };
    }
    let total: f64 = pis.iter().sum();
    let mut sum = 0.0;
    for &v in data {
        let mut pdf = 0.0;
        for j in 0..k {
            let (mu, sigma) = (x[3 * j], x[3 * j + 1]);
            if sigma != 0.0 {
                pdf += pis[j] / total * norm_pdf(v, mu, sigma);
            }
        }
        sum += pdf.ln();
    }
    -sum
}

#[cfg(test)]
mod tests {
    use super::*;

    // deterministic pseudo data: two well separated groups
    fn data() -> Vec<f64> {
        let mut v = Vec::new();
        for i in 0..700 {
            v.push(1.0 + ((i * 37 % 101) as f64 / 101.0 - 0.5) * 0.4);
        }
        for i in 0..300 {
            v.push(3.0 + ((i * 53 % 97) as f64 / 97.0 - 0.5) * 0.4);
        }
        v
    }

    #[test]
    fn recovers_two_groups() {
        let r = em_fit(&data(), &[1.5, 2.5], &[0.5, 0.5], &[0.5, 0.5], 100);
        assert!(r.ok);
        assert_eq!(r.counter, 100);
        assert!((r.mu[0] - 1.0).abs() < 0.05 && (r.mu[1] - 3.0).abs() < 0.05);
        assert!((r.pi[0] - 0.7).abs() < 0.01 && (r.pi[1] - 0.3).abs() < 0.01);
    }

    #[test]
    fn log_likelihood_of_single_standard_normal() {
        // -sum(log N(0,1)(0)) = n * 0.5 * ln(2 pi)
        let ll = neg_log_likelihood(&[0.0, 0.0], &[0.0, 1.0, 1.0]);
        assert!((ll - 2.0 * 0.5 * (2.0 * std::f64::consts::PI).ln()).abs() < 1e-12);
    }

    #[test]
    fn stops_on_nan() {
        let r = em_fit(&data(), &[1.5, 2.5], &[0.0, 0.5], &[0.5, 0.5], 100);
        assert!(!r.ok);
    }
}
