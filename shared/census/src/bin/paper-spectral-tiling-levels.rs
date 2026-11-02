use census::save;
use mrly_core::error::Result;
use mrly_core::json;

const NAME: &str = "paper-spectral-tiling-levels";
const BASE: u32 = 12;
const DIGITS: [u32; 4] = [0, 3, 6, 9];
const DEPTH: u32 = 5;
const LOW: f64 = -1.25;
const HIGH: f64 = 2.0;
const TERMS: u32 = 16;
const CANVAS: f64 = 1024.0;
const MARGIN: f64 = 0.08;
const GAP: f64 = 4.0;
const SCALE: f64 = 1e9;

type Column = (Vec<f64>, Vec<f64>);

fn lattice() -> Vec<f64> {
    let mut pts = vec![0u64];
    for j in 0..DEPTH {
        let w = u64::from(BASE).pow(j);
        pts = pts
            .iter()
            .flat_map(|&p| DIGITS.iter().map(move |&l| p + w * u64::from(l)))
            .collect();
    }
    pts.sort_unstable();
    pts.into_iter().map(|p| p as f64).collect()
}

fn power(xi: f64) -> f64 {
    let mut v = 1.0;
    let mut t = xi;
    for _ in 0..TERMS {
        t /= f64::from(BASE);
        v *= (std::f64::consts::PI * t).cos() * (8.0 * std::f64::consts::PI * t).cos();
    }
    v * v
}

fn round(q: f64) -> f64 {
    (q * SCALE).round() / SCALE
}

fn stack(bands: &[f64], thin: f64) -> Column {
    let mut kept = Vec::new();
    let mut runs = Vec::new();
    let mut run = 0.0;
    for &q in bands {
        if q < thin {
            run += q;
            continue;
        }
        runs.push(round(run));
        kept.push(round(q));
        run = 0.0;
    }
    runs.push(round(run));
    (kept, runs)
}

fn main() -> Result<()> {
    let base = lattice();
    assert_eq!(
        base.len(),
        4usize.pow(DEPTH),
        "the lattice holds 4^DEPTH points"
    );
    let groups = [
        base.iter().map(|l| 3.0 * l).collect::<Vec<_>>(),
        base.iter().map(|l| 3.0 * l + 0.75).collect::<Vec<_>>(),
    ];
    let side = CANVAS * (1.0 - 2.0 * MARGIN);
    let columns = side.round() as usize;
    let thin = GAP / side;
    let mut heights = Vec::with_capacity(columns);
    let mut stacks: [Vec<Column>; 2] = [Vec::with_capacity(columns), Vec::with_capacity(columns)];
    let mut entries = 0;
    for c in 0..columns {
        let xi = LOW + ((c as f64 + 0.5) / columns as f64) * (HIGH - LOW);
        let parts: Vec<Vec<f64>> = groups
            .iter()
            .map(|g| g.iter().map(|l| power(xi - l)).collect())
            .collect();
        let total: f64 = parts.iter().flatten().sum();
        assert!(
            total <= 1.0 + 1e-9 && total > 1.0 - 1e-2,
            "Psi({xi}) = {total} leaves (1 - 1e-2, 1 + 1e-9]"
        );
        heights.push(total);
        for (g, bands) in parts.iter().enumerate() {
            let column = stack(bands, thin);
            entries += column.0.len() + column.1.len();
            stacks[g].push(column);
        }
    }
    let widest = stacks
        .iter()
        .flatten()
        .map(|(bands, _)| bands.len())
        .max()
        .unwrap_or(0);
    println!("{columns} columns, {entries} entries, {widest} bands at most");
    save(
        NAME,
        &json!({
            "base": BASE,
            "digits": DIGITS,
            "depth": DEPTH,
            "low": LOW,
            "high": HIGH,
            "points": base.len(),
            "columns": columns,
            "gap": GAP,
            "thin": thin,
            "heights": heights,
            "stacks": stacks,
        }),
    )
}
