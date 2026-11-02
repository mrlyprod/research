use census::save;
use mrly_core::error::Result;
use mrly_core::json;
use mrly_num::factor::totients;
use mrly_num::zeta::{novelty_main, sharp_novelty, smoothed_novelty};
use std::f64::consts::PI;

const NAME: &str = "paper-novelty-meter";
const REACH: usize = 3_000_000;
const LOW: f64 = 8.0;
const HIGH: f64 = 20.5;
const PER_OCTAVE: usize = 16;

fn main() -> Result<()> {
    let phi = totients(REACH);
    assert_eq!(&phi[1..=10], &[1, 1, 2, 2, 4, 2, 6, 4, 6, 4]);
    let mut prefix = vec![0u64; REACH + 1];
    for n in 1..=REACH {
        prefix[n] = prefix[n - 1] + phi[n];
    }
    let total = prefix[REACH] as f64;
    let expect = 3.0 * (REACH as f64).powi(2) / (PI * PI);
    assert!((total - expect).abs() < REACH as f64 * (REACH as f64).ln());

    let main_bump = novelty_main();
    assert!((main_bump - 6.0 / (PI * PI) * 0.575_725_895_994).abs() < 1e-9);

    let samples = ((HIGH - LOW) * PER_OCTAVE as f64).round() as usize + 1;
    assert_eq!(samples, 201);
    let js: Vec<f64> = (0..samples)
        .map(|k| LOW + k as f64 / PER_OCTAVE as f64)
        .collect();
    let smooth: Vec<f64> = js
        .iter()
        .map(|&j| {
            let y = 2.0f64.powf(-j);
            smoothed_novelty(&phi, y, main_bump) / y.powf(1.5)
        })
        .collect();
    let rough: Vec<f64> = js
        .iter()
        .map(|&j| {
            let y = 2.0f64.powf(-j);
            sharp_novelty(&prefix, y) / y
        })
        .collect();
    save(NAME, &json!({"smooth": smooth, "rough": rough}))
}
