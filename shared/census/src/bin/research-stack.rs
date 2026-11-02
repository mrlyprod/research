use census::save;
use mrly_core::error::Result;
use mrly_core::json;
use mrly_num::gauss::{classes, Ring};
use std::collections::HashSet;

const NAME: &str = "research-stack";
const RING: Ring = Ring::Gaussian;
const BOUND: u64 = 60;

fn nodes(bound: u64) -> Vec<(i64, i64, i64)> {
    let mut seen = HashSet::new();
    let mut out = Vec::new();
    for den in classes(RING, bound) {
        let n = RING.norm(den.0, den.1) as i64;
        for p in 0..n {
            for q in 0..n {
                let g = RING.gaussian_gcd((p, q), den);
                if RING.norm(g.0, g.1) != 1 {
                    continue;
                }
                let x = (p * den.0 + q * den.1).rem_euclid(n);
                let y = (q * den.0 - p * den.1).rem_euclid(n);
                if !seen.insert((x, y, n)) {
                    continue;
                }
                out.push((x, y, n));
            }
        }
    }
    out
}

fn main() -> Result<()> {
    assert_eq!(nodes(50).len(), 672);
    let stack = nodes(BOUND);
    assert_eq!(stack.len(), 880);
    let rows: Vec<_> = stack.iter().map(|&(x, y, n)| json!([x, y, n])).collect();
    save(NAME, &json!({"bound": BOUND, "nodes": rows}))
}
