use census::save;
use mrly_core::error::Result;
use mrly_core::{Json, Map};

const NAME: &str = "research-weights";
const LEVEL: u32 = 8;
const SAMPLES: usize = 720;
const CELLS: [(i64, i64); 3] = [(0, 0), (2, 0), (0, 2)];
const CURVES: [&str; 4] = ["lattice", "smooth", "thin", "fat"];

// FOLDS

fn root(costs: [f64; 3]) -> f64 {
    let (mut lo, mut hi) = (0.0f64, 8.0f64);
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        let mass: f64 = costs.iter().map(|c| (-mid * c).exp()).sum();
        if mass > 1.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

fn pascal(n: usize) -> Vec<Vec<u128>> {
    let mut rows: Vec<Vec<u128>> = vec![vec![1]];
    for i in 1..=n {
        let mut row = vec![1u128; i + 1];
        for k in 1..i {
            row[k] = rows[i - 1][k - 1] + rows[i - 1][k];
        }
        rows.push(row);
    }
    rows
}

fn stopping(budget: f64, cheap: f64, dear: f64, binom: &[Vec<u128>]) -> u128 {
    let mut total = 0u128;
    let mut p = 0usize;
    while p as f64 * cheap <= budget {
        let mut q = 0usize;
        while p as f64 * cheap + q as f64 * dear <= budget {
            total += binom[p + q][q] << p;
            q += 1;
        }
        p += 1;
    }
    total
}

fn mass_side(dear: f64, span: (f64, f64)) -> Vec<f64> {
    let cheap = 2.0f64.ln();
    let delta = root([cheap, cheap, dear]);
    let binom = pascal(2 + (span.1 / cheap) as usize);
    (0..SAMPLES)
        .map(|i| {
            let a = span.0 + (span.1 - span.0) * i as f64 / (SAMPLES - 1) as f64;
            (stopping(a, cheap, dear, &binom) as f64).ln() - delta * a
        })
        .collect()
}

fn gasket(nums: [u64; 3]) -> Vec<(i64, i64, u128)> {
    let mut pts = vec![(0i64, 0i64, 1u128)];
    for _ in 0..LEVEL {
        let mut next = Vec::with_capacity(pts.len() * 3);
        for &(x, y, m) in &pts {
            for (k, &(dx, dy)) in CELLS.iter().enumerate() {
                next.push((x * 3 + dx, y * 3 + dy, m * nums[k] as u128));
            }
        }
        pts = next;
    }
    pts
}

fn length_side(nums: [u64; 3], den: u64) -> (f64, Vec<f64>) {
    let pts = gasket(nums);
    assert_eq!(pts.len(), 3usize.pow(LEVEL));
    let whole: u128 = pts.iter().map(|p| p.2).sum();
    assert_eq!(whole, (den as u128).pow(LEVEL));
    let mut rows: Vec<(u128, u128)> = pts
        .iter()
        .map(|&(x, y, m)| ((x * x + y * y) as u128, m))
        .collect();
    rows.sort_by_key(|row| row.0);
    let mut running = 0u128;
    let shells: Vec<(f64, f64)> = rows
        .iter()
        .map(|&(square, m)| {
            running += m;
            (square as f64, running as f64)
        })
        .collect();
    let base = 3.0f64.ln();
    let alpha = -(nums[0] as f64 / den as f64).ln() / base;
    let ys = (0..SAMPLES)
        .map(|i| {
            let u = 4.0 + 4.0 * i as f64 / (SAMPLES - 1) as f64;
            let radius = (u * base).exp();
            let reach = radius * radius;
            let cut = shells.partition_point(|shell| shell.0 <= reach);
            let held = if cut == 0 { 0.0 } else { shells[cut - 1].1 };
            held.ln() - LEVEL as f64 * (den as f64).ln() - alpha * u * base
        })
        .collect();
    (alpha, ys)
}

// WRITE

fn main() -> Result<()> {
    let window = (48.0, 48.0 + 6.0 * 2.0f64.ln());
    let lattice = mass_side(4.0f64.ln(), window);
    let smooth = mass_side(3.0f64.ln(), window);
    assert!((root([2.0f64.ln(), 2.0f64.ln(), 4.0f64.ln()]) - 1.2715533032).abs() < 1e-9);
    assert!((root([2.0f64.ln(), 2.0f64.ln(), 3.0f64.ln()]) - 1.3646005647).abs() < 1e-9);
    let (thin, ring_thin) = length_side([2, 2, 1], 5);
    let (fat, ring_fat) = length_side([3, 3, 2], 8);
    assert!((thin - 0.834043767).abs() < 1e-9);
    assert!((fat - 0.892789261).abs() < 1e-9);
    let curves = [lattice, smooth, ring_thin, ring_fat];
    assert!(curves.iter().flatten().all(|y| y.is_finite()));
    assert!(curves.iter().all(|ys| ys.len() == SAMPLES));
    let mut body = Map::new();
    for (key, ys) in CURVES.iter().zip(curves) {
        body.insert(key.to_string(), Json::from(ys));
    }
    save(NAME, &Json::Object(body))
}
