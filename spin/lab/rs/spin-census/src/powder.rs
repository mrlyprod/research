use mrly_core::Tensor;
use mrly_num::fft::fft2;

pub struct Powder {
    pub slope: f64,
    pub low: f64,
    pub high: f64,
    pub swing: f64,
}

fn rings(grid: &Tensor, pad: usize, bins: usize) -> Vec<(f64, f64, f64)> {
    let side = grid.shape[0];
    let mut re = vec![0.0f64; pad * pad];
    let mut im = vec![0.0f64; pad * pad];
    for row in 0..side {
        for col in 0..side {
            re[row * pad + col] = grid.at(row * side + col) as f64;
        }
    }
    fft2(&mut re, &mut im, pad, false).expect("the pad square is pad by pad");
    let half = (pad / 2) as i64;
    let top = (half as f64 * 2f64.sqrt()).ln();
    let mut sums = vec![0.0f64; bins];
    let mut logs = vec![0.0f64; bins];
    let mut hits = vec![0.0f64; bins];
    for row in 0..pad {
        for col in 0..pad {
            let u = if (row as i64) <= half {
                row as i64
            } else {
                row as i64 - pad as i64
            };
            let v = if (col as i64) <= half {
                col as i64
            } else {
                col as i64 - pad as i64
            };
            let norm = ((u * u + v * v) as f64).sqrt();
            if norm < 1.0 {
                continue;
            }
            let bin = ((norm.ln() / top) * bins as f64) as usize;
            if bin >= bins {
                continue;
            }
            let flat = row * pad + col;
            let power = re[flat] * re[flat] + im[flat] * im[flat];
            sums[bin] += power;
            logs[bin] += power.max(1e-300).ln();
            hits[bin] += 1.0;
        }
    }
    let mut out = Vec::new();
    for bin in 0..bins {
        if hits[bin] == 0.0 {
            continue;
        }
        let centre = ((bin as f64 + 0.5) / bins as f64 * top).exp();
        out.push((centre, sums[bin] / hits[bin], (logs[bin] / hits[bin]).exp()));
    }
    out
}

fn fit(points: &[(f64, f64)]) -> (f64, f64) {
    let n = points.len() as f64;
    let sx: f64 = points.iter().map(|p| p.0.ln()).sum();
    let sy: f64 = points.iter().map(|p| p.1.ln()).sum();
    let sxx: f64 = points.iter().map(|p| p.0.ln() * p.0.ln()).sum();
    let sxy: f64 = points.iter().map(|p| p.0.ln() * p.1.ln()).sum();
    let slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    (slope, (sy - slope * sx) / n)
}

fn slide(points: &[(f64, f64)], width: f64, step: f64) -> (f64, f64) {
    let first = points.first().map(|p| p.0.ln()).unwrap_or(0.0);
    let last = points.last().map(|p| p.0.ln()).unwrap_or(0.0);
    let mut low = f64::INFINITY;
    let mut high = f64::NEG_INFINITY;
    let mut start = first;
    while start + width <= last + 1e-12 {
        let window: Vec<(f64, f64)> = points
            .iter()
            .copied()
            .filter(|(k, _)| k.ln() >= start && k.ln() <= start + width)
            .collect();
        if window.len() >= 8 {
            let (value, _) = fit(&window);
            low = low.min(value);
            high = high.max(value);
        }
        start += step;
    }
    (low, high)
}

pub fn powder(
    grid: &Tensor,
    pad: usize,
    low: f64,
    high: f64,
    bins: usize,
    phases: usize,
) -> Powder {
    let all = rings(grid, pad, bins);
    let band: Vec<(f64, f64, f64)> = all
        .iter()
        .copied()
        .filter(|(k, p, g)| *k >= low && *k <= high && *p > 0.0 && *g > 0.0)
        .collect();
    let arithmetic: Vec<(f64, f64)> = band.iter().map(|(k, p, _)| (*k, *p)).collect();
    let (slope, intercept) = fit(&arithmetic);
    let period = (crate::design::BASE as f64).ln();
    let (least, most) = slide(&arithmetic, 3.0 * period, period / 4.0);
    let mut sums = vec![0.0f64; phases];
    let mut hits = vec![0.0f64; phases];
    for (k, power) in &arithmetic {
        let residual = power.ln() - slope * k.ln() - intercept;
        let phase = (k.ln() / period).rem_euclid(1.0);
        let bin = ((phase * phases as f64) as usize).min(phases - 1);
        sums[bin] += residual;
        hits[bin] += 1.0;
    }
    let curve: Vec<f64> = sums
        .iter()
        .zip(&hits)
        .map(|(sum, hit)| if *hit == 0.0 { f64::NAN } else { sum / hit })
        .collect();
    let swing = crate::mass::spread(&curve);
    Powder {
        slope,
        low: least,
        high: most,
        swing,
    }
}
