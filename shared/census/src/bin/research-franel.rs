use census::save;
use mrly_core::error::Result;
use mrly_core::json;
use mrly_math::design::elements;
use mrly_num::factor::{gcd, mobius_sieve};
use std::f64::consts::TAU;

const NAME: &str = "research-franel";
const BASE: u64 = 3;
const DIGITS: [u64; 2] = [0, 1];
const ORDER: usize = 81;
const LIT: [usize; 16] = [1, 2, 6, 8, 14, 15, 18, 20, 28, 30, 31, 36, 37, 39, 40, 81];

fn dilated(holds: &[bool], mu: &[i8]) -> Vec<i64> {
    let mut out = vec![0i64; ORDER + 1];
    for d in 1..=ORDER {
        for c in 1..=ORDER / d {
            if holds[d * c] {
                out[d] += mu[c] as i64;
            }
        }
    }
    out
}

fn kernel(d: usize, e: usize) -> f64 {
    let g = gcd(d as u128, e as u128) as f64;
    g * g / (d as f64 * e as f64)
}

fn nodes(dens: &[usize]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for &b in dens {
        for a in 1..=b {
            if gcd(a as u128, b as u128) == 1 {
                out.push((a, b));
            }
        }
    }
    out
}

fn main() -> Result<()> {
    let design: Vec<usize> = elements(BASE, &DIGITS, 5)
        .into_iter()
        .filter(|&n| n as usize <= ORDER)
        .map(|n| n as usize)
        .collect();
    assert_eq!(design.len(), 16);
    let mut holds = vec![false; ORDER + 1];
    for &b in &design {
        holds[b] = true;
    }
    let mu = mobius_sieve(ORDER);
    let x = dilated(&holds, &mu);
    let support: Vec<usize> = (1..=ORDER).filter(|&d| x[d] != 0).collect();
    assert_eq!(support, LIT);
    assert_eq!(x[1], -2);
    assert_eq!(x.iter().map(|v| v * v).sum::<i64>(), 19);

    let mut terms = Vec::new();
    for &d in &support {
        for &e in &support {
            terms.push((d, e, kernel(d, e) * (x[d] * x[e]) as f64));
        }
    }
    assert_eq!(terms.iter().filter(|t| t.2 != 0.0).count(), 256);
    let form: f64 = terms.iter().map(|t| t.2).sum();
    assert!((form - 12.574_087).abs() < 1e-5);

    let farey = nodes(&design);
    assert_eq!(farey.len(), 241);
    for k in 1..=12usize {
        let (mut re, mut im) = (0.0, 0.0);
        for &(a, b) in &farey {
            let angle = TAU * (k * a) as f64 / b as f64;
            re += angle.cos();
            im += angle.sin();
        }
        let divisor: i64 = (1..=k)
            .filter(|d| k % d == 0 && *d <= ORDER)
            .map(|d| d as i64 * x[d])
            .sum();
        assert!((re - divisor as f64).abs() < 1e-9);
        assert!(im.abs() < 1e-9);
    }

    let low = terms.iter().fold(f64::MAX, |a, t| a.min(t.2.abs()));
    let high = terms.iter().fold(0.0f64, |a, t| a.max(t.2.abs()));
    assert!(low > 3e-4 && high == 4.0);

    let lower: Vec<f64> = (1..=ORDER)
        .flat_map(|d| (1..=d).map(move |e| kernel(d, e)))
        .collect();
    let rows: Vec<_> = terms.iter().map(|&(d, e, v)| json!([d, e, v])).collect();
    save(
        NAME,
        &json!({"order": ORDER, "low": low, "high": high, "kernel": lower, "terms": rows}),
    )
}
