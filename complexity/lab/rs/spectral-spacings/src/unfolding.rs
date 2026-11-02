use crate::numerics::mean;
use faer::prelude::SolveLstsq;
use faer::Mat;

pub struct Spacings {
    pub values: Vec<f64>,
    pub negatives: usize,
}

fn spacings(unfolded: &[f64]) -> Spacings {
    let raw: Vec<f64> = unfolded.windows(2).map(|pair| pair[1] - pair[0]).collect();
    let negatives = raw.iter().filter(|s| **s < 0.0).count();
    let clamped: Vec<f64> = raw.iter().map(|s| s.max(0.0)).collect();
    let scale = mean(&clamped);
    Spacings {
        values: clamped.iter().map(|s| s / scale).collect(),
        negatives,
    }
}

fn chebyshev(x: f64, degree: usize) -> Vec<f64> {
    let mut out = vec![1.0, x];
    for k in 2..=degree {
        out.push(2.0 * x * out[k - 1] - out[k - 2]);
    }
    out.truncate(degree + 1);
    out
}

pub fn polynomial(values: &[f64], degree: usize) -> Spacings {
    let n = values.len();
    let (low, high) = (values[0], values[n - 1]);
    let scaled: Vec<f64> = values
        .iter()
        .map(|v| 2.0 * (v - low) / (high - low) - 1.0)
        .collect();
    let mut design = Mat::<f64>::zeros(n, degree + 1);
    let mut target = Mat::<f64>::zeros(n, 1);
    for (i, x) in scaled.iter().enumerate() {
        for (k, basis) in chebyshev(*x, degree).iter().enumerate() {
            design[(i, k)] = *basis;
        }
        target[(i, 0)] = i as f64 + 0.5;
    }
    let coefficients = design.qr().solve_lstsq(&target);
    let unfolded: Vec<f64> = scaled
        .iter()
        .map(|x| {
            chebyshev(*x, degree)
                .iter()
                .enumerate()
                .map(|(k, basis)| coefficients[(k, 0)] * basis)
                .sum()
        })
        .collect();
    spacings(&unfolded)
}

pub fn window(values: &[f64], half: usize) -> Spacings {
    let n = values.len();
    let mut unfolded = vec![0.0; n];
    for i in 1..n {
        let lo = i.saturating_sub(half);
        let hi = (i + half).min(n - 1);
        let width = values[hi] - values[lo];
        let gap = values[i] - values[i - 1];
        let step = if gap > 0.0 {
            gap * (hi - lo) as f64 / width
        } else {
            0.0
        };
        unfolded[i] = unfolded[i - 1] + step;
    }
    spacings(&unfolded)
}
