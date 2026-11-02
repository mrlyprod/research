use census::save;
use mrly_core::error::Result;
use mrly_core::json;
use std::f64::consts::PI;

const NAME: &str = "paper-first-base-below-a-quarter";
const LOW: usize = 10;
const HIGH: usize = 40;
const SUB: usize = 4;
const CAP: f64 = 6.0e6;
const QUARTER: f64 = 0.25;
const RUNGS: usize = 395;
const GOLD: usize = 163;

// LADDER

fn weights(q: usize, a0: usize, nd: u32, m: usize) -> Vec<f64> {
    let w = q.pow(nd);
    let inv = 1.0 / (q as f64 - 1.0);
    let c = 2.0 * a0 as f64 - q as f64 + 1.0;
    let mass: usize = (0..q).filter(|a| *a != a0).sum();
    let lip = 2.0 * PI * mass as f64 * inv;
    let slack = lip / (2.0 * w as f64 * m as f64) + 1e-12;
    let step = 1.0 / (w as f64 * m as f64);
    let mut g = vec![0.0; w];
    for cell in 0..w.div_ceil(2) {
        let mut top = 0.0f64;
        for r in 0..m {
            let t = ((cell * m + r) as f64 + 0.5) * step;
            let d = (PI * q as f64 * t).sin() / (PI * t).sin();
            let e = d * d - 2.0 * d * (PI * c * t).cos() + 1.0;
            let v = if e > 0.0 { e.sqrt() * inv } else { 0.0 };
            top = top.max(v);
        }
        let u = (top + slack).min(1.0);
        g[cell] = u;
        g[w - 1 - cell] = u;
    }
    g
}

fn sweep(g: &[f64], q: usize, nd: u32, y: &[f64]) -> Vec<f64> {
    let s = q.pow(nd - 1);
    let p = q.pow(nd - 2);
    (0..s)
        .map(|v| {
            let b = (v % p) * q;
            (0..q).map(|c| g[v * q + c] * y[b + c]).sum()
        })
        .collect()
}

fn power(g: &[f64], q: usize, nd: u32, seed: Option<Vec<f64>>, steps: usize) -> Vec<f64> {
    let s = q.pow(nd - 1);
    let mut y = seed.unwrap_or_else(|| vec![1.0; s]);
    for _ in 0..steps {
        let z = sweep(g, q, nd, &y);
        let top = z.iter().copied().fold(0.0f64, f64::max);
        if top <= 0.0 {
            return y;
        }
        y = z.into_iter().map(|x| x / top).collect();
    }
    y
}

fn bound(g: &[f64], q: usize, nd: u32, y: &[f64]) -> f64 {
    let z = sweep(g, q, nd, y);
    let mu = z
        .iter()
        .zip(y.iter())
        .map(|(a, b)| a / b)
        .fold(0.0f64, f64::max);
    mu.ln() / (q as f64).ln()
}

fn lift(y: &[f64], q: usize) -> Vec<f64> {
    (0..y.len() * q).map(|v| y[v / q]).collect()
}

fn exponent(q: usize, a0: usize) -> f64 {
    let mut nd = 3u32;
    let start = weights(q, a0, nd, SUB);
    let mut y = power(&start, q, nd, None, 300);
    let mut e = bound(&start, q, nd, &y);
    while e >= QUARTER && nd < 5 && (q as f64).powi(nd as i32 + 1) <= CAP {
        let g = weights(q, a0, nd + 1, SUB);
        let z = power(&g, q, nd + 1, Some(lift(&y, q)), 40);
        let step = bound(&g, q, nd + 1, &z);
        let gain = e - step;
        nd += 1;
        y = z;
        e = step;
        if e >= QUARTER && gain <= e - QUARTER {
            break;
        }
    }
    e
}

fn ladder() -> Vec<(usize, f64)> {
    let rungs: Vec<(usize, usize)> = (LOW..=HIGH)
        .flat_map(|q| (0..q.div_ceil(2)).map(move |a0| (q, a0)))
        .collect();
    let lanes = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
        .max(1);
    let mut out = vec![0.0f64; rungs.len()];
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..lanes)
            .map(|lane| {
                let rungs = &rungs;
                scope.spawn(move || {
                    rungs
                        .iter()
                        .enumerate()
                        .skip(lane)
                        .step_by(lanes)
                        .map(|(i, (q, a0))| (i, exponent(*q, *a0)))
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        for handle in handles {
            for (i, e) in handle.join().expect("a lane of the ladder") {
                out[i] = e;
            }
        }
    });
    rungs.iter().zip(out).map(|((q, _), e)| (*q, e)).collect()
}

// WRITE

fn main() -> Result<()> {
    let rungs = ladder();
    let below = |q: usize| rungs.iter().filter(|r| r.0 == q && r.1 < QUARTER).count();
    let sets = |q: usize| q.div_ceil(2);
    assert_eq!(rungs.len(), RUNGS, "the ladder has 395 rungs");
    assert!((LOW..21).all(|q| below(q) == 0), "no base under 21 clears");
    assert_eq!(below(21), 1, "base 21 clears at exactly one digit");
    assert!(
        (21..34).all(|q| below(q) < sets(q)),
        "every base from 21 to 33 still misses somewhere"
    );
    assert!(
        (34..=HIGH).all(|q| below(q) == sets(q)),
        "every base from 34 clears at every digit"
    );
    let gold = rungs.iter().filter(|r| r.1 < QUARTER).count();
    assert_eq!(gold, GOLD, "163 rungs sit under the quarter");

    let pairs: Vec<_> = rungs.iter().map(|&(q, e)| json!([q, e])).collect();
    save(NAME, &json!({"rungs": pairs}))
}
