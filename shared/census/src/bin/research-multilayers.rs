use census::save;
use mrly_core::error::Result;
use mrly_core::json;

const NAME: &str = "research-multilayers";
const NA: f64 = 2.3;
const NB: f64 = 1.45;
const BASE: usize = 7;
const LEVELS: usize = 12;
const BINS: usize = 128;
const SAMPLES: usize = 48;

type C = (f64, f64);
type M = [C; 4];

fn cm(a: C, b: C) -> C {
    (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0)
}

fn ca(a: C, b: C) -> C {
    (a.0 + b.0, a.1 + b.1)
}

fn mul(p: &M, q: &M) -> M {
    [
        ca(cm(p[0], q[0]), cm(p[1], q[2])),
        ca(cm(p[0], q[1]), cm(p[1], q[3])),
        ca(cm(p[2], q[0]), cm(p[3], q[2])),
        ca(cm(p[2], q[1]), cm(p[3], q[3])),
    ]
}

fn lay(n: f64, ph: f64) -> M {
    let (s, c) = ph.sin_cos();
    [(c, 0.0), (0.0, s / n), (0.0, n * s), (c, 0.0)]
}

fn ratio(m: &M) -> f64 {
    let re = NB * m[0].0 + NB * NB * m[1].0 - m[2].0 - NB * m[3].0;
    let im = NB * m[0].1 + NB * NB * m[1].1 - m[2].1 - NB * m[3].1;
    (re * re + im * im).sqrt() / (2.0 * NB)
}

fn open_levels(digits: &[usize], delta: f64) -> Vec<bool> {
    let mut m = lay(NA, delta);
    let mut out = Vec::with_capacity(LEVELS);
    let mut dead = false;
    let mut span = delta;
    for _ in 0..LEVELS {
        let spacer = lay(NB, span.rem_euclid(std::f64::consts::TAU));
        let mut n = [(1.0, 0.0), (0.0, 0.0), (0.0, 0.0), (1.0, 0.0)];
        for d in 0..BASE {
            n = mul(&n, if digits.contains(&d) { &m } else { &spacer });
        }
        dead |= n.iter().any(|z| !(z.0.abs() < 1e100 && z.1.abs() < 1e100));
        m = n;
        span *= BASE as f64;
        out.push(!dead && ratio(&m) < 1.0);
    }
    out
}

fn brute(digits: &[usize], level: u32, delta: f64) -> f64 {
    let mut m = [(1.0, 0.0), (0.0, 0.0), (0.0, 0.0), (1.0, 0.0)];
    for p in 0..BASE.pow(level) {
        let mut q = p;
        let mut on = true;
        for _ in 0..level {
            on &= digits.contains(&(q % BASE));
            q /= BASE;
        }
        m = mul(&m, &lay(if on { NA } else { NB }, delta));
    }
    ratio(&m)
}

fn recursive(digits: &[usize], level: u32, delta: f64) -> f64 {
    let mut m = lay(NA, delta);
    let mut span = delta;
    for _ in 0..level {
        let spacer = lay(NB, span);
        let mut n = [(1.0, 0.0), (0.0, 0.0), (0.0, 0.0), (1.0, 0.0)];
        for d in 0..BASE {
            n = mul(&n, if digits.contains(&d) { &m } else { &spacer });
        }
        m = n;
        span *= BASE as f64;
    }
    ratio(&m)
}

fn density(digits: &[usize]) -> Vec<Vec<f64>> {
    let width = std::f64::consts::FRAC_PI_2 / BINS as f64;
    let columns: Vec<Vec<f64>> = (0..BINS)
        .map(|bin| {
            let mut column = vec![0.0; LEVELS];
            for s in 0..SAMPLES {
                let delta = width * (bin as f64 + (s as f64 + 0.5) / SAMPLES as f64);
                for (k, open) in open_levels(digits, delta).into_iter().enumerate() {
                    if open {
                        column[k] += 1.0 / SAMPLES as f64;
                    }
                }
            }
            column
        })
        .collect();
    (0..LEVELS)
        .map(|k| columns.iter().map(|column| column[k]).collect())
        .collect()
}

fn main() -> Result<()> {
    let critical = [0, 2, 4, 6];
    let dark = [0, 2, 3, 4, 6];
    for digits in [&critical[..], &dark[..]] {
        for delta in [0.37, 1.21] {
            let (a, b) = (brute(digits, 3, delta), recursive(digits, 3, delta));
            assert!((a - b).abs() <= 1e-9 * a.max(1.0));
        }
    }
    let light = density(&critical);
    let fading = density(&dark);
    assert_eq!(light.len(), LEVELS);
    assert!(light[LEVELS - 1].iter().sum::<f64>() > 1.6 * fading[LEVELS - 1].iter().sum::<f64>());
    save(NAME, &json!({"light": light, "fading": fading}))
}
