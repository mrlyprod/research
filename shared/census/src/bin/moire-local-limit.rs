use census::save;
use mrly_core::error::Result;
use mrly_core::json;
use mrly_math::two::{carpet, Cell2d};

const NAME: &str = "moire-local-limit";
const LOW: usize = 21;
const HIGH: usize = 23;
const PANEL: usize = 410;
const STEPS: usize = 400;
const SCALE: f64 = (LOW * HIGH * LOW * HIGH) as f64;

fn inked(cell: &Cell2d) -> Result<Vec<Vec<bool>>> {
    let n = cell.width();
    let mut rows = vec![vec![false; n]; n];
    for (i, row) in rows.iter_mut().enumerate() {
        for (j, slot) in row.iter_mut().enumerate() {
            *slot = cell.types().get(&[i, j])? != 0;
        }
    }
    Ok(rows)
}

fn pieces(lo: f64, hi: f64) -> Vec<(f64, usize, usize)> {
    let mut cuts = vec![lo, hi];
    for n in [LOW, HIGH] {
        for k in 1..n {
            let x = k as f64 / n as f64;
            if x > lo && x < hi {
                cuts.push(x);
            }
        }
    }
    cuts.sort_by(|a, b| a.total_cmp(b));
    cuts.windows(2)
        .filter(|w| w[1] > w[0])
        .map(|w| {
            let mid = 0.5 * (w[0] + w[1]);
            let a = ((mid * LOW as f64) as usize).min(LOW - 1);
            let b = ((mid * HIGH as f64) as usize).min(HIGH - 1);
            (w[1] - w[0], a, b)
        })
        .collect()
}

fn overlay(low: &[Vec<bool>], high: &[Vec<bool>]) -> Vec<f64> {
    let side = PANEL as f64;
    let spans: Vec<Vec<(f64, usize, usize)>> = (0..PANEL)
        .map(|k| pieces(k as f64 / side, (k + 1) as f64 / side))
        .collect();
    let mut values = Vec::with_capacity(PANEL * PANEL);
    for rows in &spans {
        for cols in &spans {
            let mut area = 0.0;
            for &(ly, ra, rb) in rows {
                for &(lx, ca, cb) in cols {
                    if low[ra][ca] != high[rb][cb] {
                        area += lx * ly;
                    }
                }
            }
            values.push(area * side * side);
        }
    }
    values
}

fn limit(u: f64, v: f64) -> f64 {
    0.5 * (1.0 - (1.0 - 2.0 * u).abs() * (1.0 - 2.0 * v).abs())
}

fn fill(n: usize) -> f64 {
    let q = (n - 1) as f64 / (2 * n) as f64;
    1.0 - q * q
}

fn exact(raw: &[f64]) -> Vec<u32> {
    raw.iter()
        .map(|&v| {
            let n = (v * SCALE).round();
            assert!((v * SCALE - n).abs() < 1e-6);
            n as u32
        })
        .collect()
}

fn octant(whole: &[u32]) -> Vec<u32> {
    let half = PANEL / 2;
    let fold = |k: usize| k.min(PANEL - 1 - k);
    for i in 0..PANEL {
        for j in 0..PANEL {
            let (a, b) = (fold(i), fold(j));
            assert_eq!(whole[i * PANEL + j], whole[a.min(b) * PANEL + a.max(b)]);
        }
    }
    let mut out = Vec::with_capacity(half * (half + 1) / 2);
    for i in 0..half {
        for j in i..half {
            out.push(whole[i * PANEL + j]);
        }
    }
    out
}

fn main() -> Result<()> {
    let low = inked(&carpet(LOW, 1)?)?;
    let high = inked(&carpet(HIGH, 1)?)?;
    assert_eq!(low.iter().flatten().filter(|&&on| on).count(), 341);
    assert_eq!(high.iter().flatten().filter(|&&on| on).count(), 408);

    let raw = overlay(&low, &high);
    let mean = raw.iter().sum::<f64>() / raw.len() as f64;
    let (a, b) = (fill(LOW), fill(HIGH));
    assert!((mean - (a + b - 2.0 * a * b)).abs() < 1e-9);

    let side = PANEL as f64;
    let smooth: Vec<f64> = (0..PANEL * PANEL)
        .map(|k| {
            limit(
                ((k % PANEL) as f64 + 0.5) / side,
                ((k / PANEL) as f64 + 0.5) / side,
            )
        })
        .collect();
    let smean = smooth.iter().sum::<f64>() / smooth.len() as f64;
    assert!((smean - 0.375).abs() < 1e-12);

    let axis: Vec<f64> = (0..PANEL)
        .map(|k| (1.0 - 2.0 * ((k as f64 + 0.5) / side)).abs())
        .collect();
    for (k, &s) in smooth.iter().enumerate() {
        assert_eq!(s, 0.5 * (1.0 - axis[k % PANEL] * axis[k / PANEL]));
    }

    let kappa = 2f64.sqrt() - 1.0;
    assert!((2.0 * kappa - 2f64.sqrt() * (1.0 - kappa)).abs() < 1e-12);
    let arc: Vec<[f64; 2]> = (0..=STEPS)
        .map(|k| {
            let x = kappa + (1.0 - kappa) * k as f64 / STEPS as f64;
            [x, kappa / x]
        })
        .collect();

    save(
        NAME,
        &json!({
            "scale": SCALE as u32,
            "octant": octant(&exact(&raw)),
            "axis": axis,
            "kappa": kappa,
            "arc": arc,
        }),
    )
}
