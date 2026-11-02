use census::save;
use mrly_core::error::Result;
use mrly_core::json;

const NAME: &str = "research-toolpaths";
const DIRS: [usize; 7] = [0, 5, 3, 4, 0, 0, 1];
const FLIPS: [bool; 7] = [false, true, true, false, false, false, true];
const LEVEL: usize = 4;
const STEPS: usize = 2401;
const COUNTS: [usize; 3] = [343, 1371, 686];

fn expand(seq: &[(usize, bool)]) -> Vec<(usize, bool)> {
    let mut out = Vec::with_capacity(seq.len() * 7);
    for &(d, r) in seq {
        for j in 0..7 {
            let i = if r { 6 - j } else { j };
            out.push(((d + DIRS[i]) % 6, FLIPS[i] != r));
        }
    }
    out
}

fn trim(v: f64) -> f64 {
    (v * 1e9).round() / 1e9 + 0.0
}

fn main() -> Result<()> {
    let mut seq = vec![(0usize, false)];
    for _ in 0..LEVEL {
        seq = expand(&seq);
    }
    assert_eq!(seq.len(), STEPS);
    let unit: Vec<(f64, f64)> = (0..6)
        .map(|d| {
            let a = std::f64::consts::PI / 3.0 * d as f64;
            (a.cos(), a.sin())
        })
        .collect();
    let mut pts = vec![(0.0, 0.0)];
    for &(d, _) in &seq {
        let (x, y) = pts[pts.len() - 1];
        pts.push((x + unit[d].0, y + unit[d].1));
    }
    let mut counts = [0usize; 3];
    let mut sharp = Vec::new();
    for i in 1..seq.len() {
        let t = (seq[i].0 + 6 - seq[i - 1].0) % 6;
        let t = t.min(6 - t);
        counts[t] += 1;
        if t == 2 {
            sharp.push(i);
        }
    }
    assert_eq!(counts, COUNTS);
    let path: Vec<[f64; 2]> = pts.iter().map(|&(x, y)| [trim(x), trim(y)]).collect();
    save(NAME, &json!({"path": path, "sharp": sharp}))
}
