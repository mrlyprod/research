use census::save;
use mrly_core::error::Result;
use mrly_core::json;
use std::f64::consts::PI;

const NAME: &str = "paper-half-interval-mobius";
const GAMMA: f64 = 0.577_215_664_901_532_9;
const LOW: f64 = 2.0;
const HIGH: f64 = 7.0;
const SAMPLES: usize = 1600;
const FIFTH: u64 = 94939;
const QUARTER: u64 = 3789;

// CERTIFICATE

fn rate(base: u64) -> f64 {
    let b = base as f64;
    let fill = (b + 1.0) / 2.0;
    let cp = 2f64.sqrt() - 4.0 / PI;
    let tan = (3.0 * PI / 8.0 + PI / (4.0 * b)).tan();
    let x0 = (b / PI) * ((b + 3.0).ln() + GAMMA + tan.ln()) + cp * (b + 1.0).powi(2) / (8.0 * b);
    (fill + x0 + 0.5 / (PI / (2.0 * b)).sin()) / fill
}

fn margin(base: u64, bar: f64) -> f64 {
    bar - rate(base).ln() / (base as f64).ln()
}

fn odd(x: f64) -> u64 {
    let n = x.round() as u64;
    if n % 2 == 1 {
        n
    } else {
        n + 1
    }
}

// WRITE

fn main() -> Result<()> {
    assert!(margin(FIFTH - 2, 0.2) < 0.0 && margin(FIFTH, 0.2) > 0.0);
    assert!(margin(QUARTER - 2, 0.25) < 0.0 && margin(QUARTER, 0.25) > 0.0);
    assert!((1.3678e-8..1.3679e-8).contains(&margin(FIFTH, 0.2)));
    assert!((7.9625e-6..7.9626e-6).contains(&margin(QUARTER, 0.25)));
    assert!((0.199_387_1..=0.199_387_2).contains(&(0.2 - margin(100_003, 0.2))));

    let mut bases: Vec<u64> = (0..=SAMPLES)
        .map(|i| odd(10f64.powf(LOW + (HIGH - LOW) * i as f64 / SAMPLES as f64)))
        .collect();
    bases.extend([QUARTER - 2, QUARTER, FIFTH - 2, FIFTH]);
    bases.sort_unstable();
    bases.dedup();
    for &b in &bases {
        assert_eq!(margin(b, 0.2) > 0.0, b >= FIFTH);
        assert_eq!(margin(b, 0.25) > 0.0, b >= QUARTER);
    }

    let q25: Vec<f64> = bases.iter().map(|&b| margin(b, 0.25)).collect();
    let q20: Vec<f64> = bases.iter().map(|&b| margin(b, 0.2)).collect();
    save(NAME, &json!({"bases": bases, "q25": q25, "q20": q20}))
}
