use crate::design::BASE;
use mrly_core::Tensor;

pub struct Mass {
    pub radii: Vec<f64>,
    pub cumulative: Vec<u64>,
}

impl Mass {
    pub fn at(&self, radius: f64) -> u64 {
        match self
            .radii
            .binary_search_by(|probe| probe.partial_cmp(&radius).expect("finite radii"))
        {
            Ok(index) => self.cumulative[index],
            Err(0) => 0,
            Err(index) => self.cumulative[index - 1],
        }
    }

    pub fn total(&self) -> u64 {
        *self.cumulative.last().unwrap_or(&0)
    }
}

pub fn shells(grid: &Tensor, digit: (usize, usize)) -> Mass {
    let side = grid.shape[0];
    let anchor = (side as i64 * digit.0 as i64, side as i64 * digit.1 as i64);
    let mut keys: Vec<u32> = Vec::new();
    for row in 0..side {
        for col in 0..side {
            if grid.at(row * side + col) == 0 {
                continue;
            }
            let dr = 2 * row as i64 + 1 - anchor.0;
            let dc = 2 * col as i64 + 1 - anchor.1;
            keys.push((dr * dr + dc * dc) as u32);
        }
    }
    keys.sort_unstable();
    let mut radii = Vec::new();
    let mut cumulative = Vec::new();
    let mut running = 0u64;
    for (index, key) in keys.iter().enumerate() {
        running += 1;
        if index + 1 == keys.len() || keys[index + 1] != *key {
            radii.push((*key as f64).sqrt() / 2.0);
            cumulative.push(running);
        }
    }
    Mass { radii, cumulative }
}

/// The nearest filled cell to the fixed point of the digit, as four times its squared distance, an exact integer.
pub fn nearest_cell(grid: &Tensor, digit: (usize, usize)) -> u64 {
    let side = grid.shape[0];
    let anchor = (side as i64 * digit.0 as i64, side as i64 * digit.1 as i64);
    let reach = |at: i64, low: i64| -> i64 {
        let (a, b) = (2 * low - at, at - 2 * (low + 1));
        a.max(b).max(0)
    };
    let mut best = u64::MAX;
    for row in 0..side {
        for col in 0..side {
            if grid.at(row * side + col) == 0 {
                continue;
            }
            let dr = reach(anchor.0, row as i64);
            let dc = reach(anchor.1, col as i64);
            best = best.min((dr * dr + dc * dc) as u64);
        }
    }
    best
}

pub fn horizon(code: u128, digit: (usize, usize), table: &[(usize, usize)]) -> f64 {
    let point = (digit.0 as f64 / 2.0, digit.1 as f64 / 2.0);
    let mut best = f64::INFINITY;
    for (bit, &other) in table.iter().enumerate().take(9) {
        if code >> bit & 1 == 0 {
            continue;
        }
        if other == digit {
            continue;
        }
        let gap = |at: f64, index: usize| {
            let (low, high) = (index as f64 / 3.0, (index as f64 + 1.0) / 3.0);
            (low - at).max(at - high).max(0.0)
        };
        let (dr, dc) = (gap(point.0, other.0), gap(point.1, other.1));
        best = best.min((dr * dr + dc * dc).sqrt());
    }
    (3.0 * best).min(1.0)
}

pub struct Ripple {
    pub slope: f64,
    pub curve: Vec<f64>,
    pub swing: f64,
    pub drift: f64,
    pub periods: usize,
}

fn sample(mass: &Mass, low: f64, periods: usize, bins: usize) -> Vec<(f64, f64)> {
    let step = (BASE as f64).ln() / bins as f64;
    (0..periods * bins)
        .filter_map(|index| {
            let radius = low * (index as f64 * step).exp();
            let count = mass.at(radius);
            if count == 0 {
                None
            } else {
                Some((radius.ln(), (count as f64).ln()))
            }
        })
        .collect()
}

fn fit(points: &[(f64, f64)]) -> (f64, f64) {
    let n = points.len() as f64;
    let sx: f64 = points.iter().map(|p| p.0).sum();
    let sy: f64 = points.iter().map(|p| p.1).sum();
    let sxx: f64 = points.iter().map(|p| p.0 * p.0).sum();
    let sxy: f64 = points.iter().map(|p| p.0 * p.1).sum();
    let slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    (slope, (sy - slope * sx) / n)
}

fn fold(points: &[(f64, f64)], dimension: f64, bins: usize) -> Vec<f64> {
    let mut sums = vec![0.0; bins];
    let mut hits = vec![0.0; bins];
    for (lx, ly) in points {
        let residual = ly - dimension * lx;
        let phase = (lx / (BASE as f64).ln()).rem_euclid(1.0);
        let bin = ((phase * bins as f64) as usize).min(bins - 1);
        sums[bin] += residual;
        hits[bin] += 1.0;
    }
    let curve: Vec<f64> = sums
        .iter()
        .zip(&hits)
        .map(|(sum, hit)| if *hit == 0.0 { f64::NAN } else { sum / hit })
        .collect();
    let live: Vec<f64> = curve.iter().copied().filter(|v| v.is_finite()).collect();
    let mean = live.iter().sum::<f64>() / live.len().max(1) as f64;
    curve.iter().map(|v| v - mean).collect()
}

/// The whole powers of the base that fit between the window ends, counted by repeated multiplication so no logarithm can round a boundary the wrong way.
pub fn periods(low: f64, high: f64) -> usize {
    let mut count = 0usize;
    let mut edge = low * BASE as f64;
    while edge <= high * (1.0 + 1e-12) {
        count += 1;
        edge *= BASE as f64;
    }
    count
}

pub fn ripple(mass: &Mass, low: f64, high: f64, dimension: f64, bins: usize) -> Ripple {
    let periods = periods(low, high).max(1);
    let points = sample(mass, low, periods, bins);
    let (slope, _) = fit(&points);
    let curve = fold(&points, dimension, bins);
    let swing = spread(&curve);
    let half = periods / 2;
    let drift = if half == 0 {
        f64::NAN
    } else {
        let early = fold(&sample(mass, low, half, bins), dimension, bins);
        let late = fold(
            &sample(
                mass,
                low * (BASE as f64).powi(half as i32),
                periods - half,
                bins,
            ),
            dimension,
            bins,
        );
        early
            .iter()
            .zip(&late)
            .map(|(a, b)| (a - b).abs())
            .fold(
                0.0f64,
                |best, gap| if gap.is_finite() { best.max(gap) } else { best },
            )
    };
    Ripple {
        slope,
        curve,
        swing,
        drift,
        periods,
    }
}

pub fn spread(curve: &[f64]) -> f64 {
    let live: Vec<f64> = curve.iter().copied().filter(|v| v.is_finite()).collect();
    let low = live.iter().copied().fold(f64::INFINITY, f64::min);
    let high = live.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    high - low
}

pub fn distance(left: &[f64], right: &[f64]) -> f64 {
    left.iter()
        .zip(right)
        .map(|(a, b)| (a - b).abs())
        .fold(
            0.0f64,
            |best, gap| if gap.is_finite() { best.max(gap) } else { best },
        )
}

pub fn scaling_error(mass: &Mass, low: f64, high: f64, fill: f64) -> f64 {
    let mut worst: f64 = 0.0;
    let mut radius = low;
    while radius * (BASE as f64) <= high {
        let near = mass.at(radius) as f64;
        let far = mass.at(radius * BASE as f64) as f64;
        if near > 0.0 {
            worst = worst.max((far / (fill * near) - 1.0).abs());
        }
        radius *= 1.05;
    }
    worst
}

#[cfg(test)]
mod tests {
    use super::periods;

    #[test]
    fn the_window_counts_whole_periods_at_every_level() {
        for level in 4..=12u32 {
            let side = 3f64.powi(level as i32);
            assert_eq!(
                periods(27.0, side),
                level as usize - 3,
                "corner level {level}"
            );
            assert_eq!(
                periods(27.0, side / 2.0),
                level as usize - 4,
                "centre level {level}"
            );
        }
        assert_eq!(periods(27.0, 26.0), 0);
        assert_eq!(periods(27.0, 81.0), 1);
        assert_eq!(periods(27.0, 80.9), 0);
    }
}
