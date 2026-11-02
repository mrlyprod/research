use census::save;
use mrly_core::error::Result;
use mrly_core::json;

const NAME: &str = "paper-sparse-mertens-under-grh";
const BASE: i64 = 7;
const LEVEL: u32 = 3;

fn proved() -> f64 {
    let q = BASE as f64;
    let n = ((BASE - 2) as f64 / 2.0).ceil();
    let harmonic = n.ln() + 0.577_215_664_901_532_9 + 1.0 / (2.0 * n);
    let pi = std::f64::consts::PI;
    let phi = (4.0 / pi) * q + (2.0 * q / pi) * harmonic + (1.0 - 2.0 / pi) * (q - 2.0) + 0.727;
    1.0 + phi / q
}

fn main() -> Result<()> {
    let bound = (BASE as f64 * proved()).powi(LEVEL as i32);
    assert!(bound > BASE.pow(LEVEL) as f64);
    save(NAME, &json!({"bound": bound}))
}
