use census::save;
use mrly_core::error::Result;
use mrly_core::json;
use mrly_math::design::elements;
use mrly_num::factor::{gcd, mobius_sieve};
use mrly_num::lattice::farey;

const NAME: &str = "paper-franel-converse";
const BASE: u64 = 3;
const DIGITS: [u64; 2] = [0, 1];
const ORDER: usize = 81;

fn nodes_of(dens: &[usize]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for &b in dens {
        for a in 1..=b {
            if gcd(a as u128, b as u128) == 1 {
                out.push((a, b));
            }
        }
    }
    out.sort_by(|p, q| (p.0 * q.1).cmp(&(q.0 * p.1)));
    out
}

fn sawtooth(nodes: &[(usize, usize)]) -> Vec<f64> {
    let m = nodes.len() as f64;
    nodes
        .iter()
        .enumerate()
        .map(|(j, &(a, b))| a as f64 / b as f64 - (j + 1) as f64 / m)
        .collect()
}

fn kernel(dilates: &[i64]) -> f64 {
    let mut total = 0.0;
    for d in 1..dilates.len() {
        if dilates[d] == 0 {
            continue;
        }
        for e in 1..dilates.len() {
            if dilates[e] == 0 {
                continue;
            }
            let g = gcd(d as u128, e as u128) as f64;
            total += g * g / (d as f64 * e as f64) * (dilates[d] * dilates[e]) as f64;
        }
    }
    total
}

fn dilated(holds: &[bool], mu: &[i8]) -> Vec<i64> {
    let q = holds.len() - 1;
    let mut out = vec![0i64; q + 1];
    for d in 1..=q {
        for c in 1..=q / d {
            if holds[d * c] {
                out[d] += mu[c] as i64;
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
    assert_eq!(design[..5], [1, 3, 4, 9, 10]);
    let mut holds = vec![false; ORDER + 1];
    for &b in &design {
        holds[b] = true;
    }
    let thin = nodes_of(&design);
    assert_eq!(thin.len(), 241);
    let full: Vec<(usize, usize)> = farey(ORDER)
        .iter()
        .filter(|n| n.num > 0)
        .map(|n| (n.num as usize, n.den as usize))
        .collect();
    assert_eq!(full.len(), 2020);
    assert_eq!(full, nodes_of(&(1..=ORDER).collect::<Vec<usize>>()));
    assert_eq!(
        thin,
        full.iter()
            .copied()
            .filter(|&(_, b)| holds[b])
            .collect::<Vec<_>>()
    );

    let delta_thin = sawtooth(&thin);
    let delta_full = sawtooth(&full);
    let s2_thin: f64 = delta_thin.iter().map(|v| v * v).sum();
    let s2_full: f64 = delta_full.iter().map(|v| v * v).sum();
    assert!((s2_thin - 0.004_002_104_850).abs() < 1e-9);
    assert!((s2_full - 0.006_524_653_839).abs() < 1e-9);

    let mu = mobius_sieve(ORDER);
    let x_thin = dilated(&holds, &mu);
    assert_eq!(x_thin[1], -2);
    assert_eq!(x_thin.iter().map(|v| v * v).sum::<i64>(), 19);
    let g_thin = kernel(&x_thin);
    assert!((g_thin - 12.574_087).abs() < 1e-5);
    assert!((g_thin - (12.0 * 241.0 * s2_thin + 1.0)).abs() < 1e-8);
    let x_full = dilated(&[true; ORDER + 1], &mu);
    assert_eq!(x_full[1], -4);
    let g_full = kernel(&x_full);
    assert!((g_full - 159.157_609).abs() < 1e-5);
    assert!((g_full - (12.0 * 2020.0 * s2_full + 1.0)).abs() < 1e-8);

    save(NAME, &json!({"thin": delta_thin, "full": delta_full}))
}
