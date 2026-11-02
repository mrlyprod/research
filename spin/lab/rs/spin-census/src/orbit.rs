use crate::design::{floats, plane, BASE};
use mrly_math::spin::harmonics;

pub fn spectrum(code: u128, level: usize, rings: usize, orders: usize) -> Vec<f64> {
    let grid = plane(code, BASE, level);
    let side = grid.shape[0];
    harmonics(&floats(&grid), side, rings, orders).unwrap_or_default()
}

pub fn agree(left: &[f64], right: &[f64], tolerance: f64) -> bool {
    let scale = left[0].max(right[0]).max(1e-30);
    left.iter()
        .zip(right)
        .all(|(a, b)| (a - b).abs() / scale < tolerance)
}

pub fn gap(left: &[f64], right: &[f64]) -> f64 {
    let scale = left[0].max(right[0]).max(1e-30);
    left.iter()
        .zip(right)
        .map(|(a, b)| (a - b).abs() / scale)
        .fold(0.0f64, f64::max)
}

pub fn census(level: usize, rings: usize, orders: usize) -> Vec<Vec<f64>> {
    let mut out = vec![vec![0.0; orders + 1]];
    for code in 1..512u128 {
        out.push(spectrum(code, level, rings, orders));
    }
    out
}
