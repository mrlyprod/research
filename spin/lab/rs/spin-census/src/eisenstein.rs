use crate::gaussian::{factorise, reduce};
use std::collections::HashSet;

pub const HEX_CONSTANT: f64 = 0.638_909_405_445_343_9;
pub const WINDOW: usize = 3;

pub fn loeschian(value: usize, least: &[u32]) -> bool {
    factorise(value, least)
        .iter()
        .all(|(prime, power)| prime % 3 != 2 || power % 2 == 0)
}

pub fn loeschian_prefix(cap: usize, least: &[u32]) -> Vec<u32> {
    let mut out = vec![0u32; cap + 1];
    for norm in 1..=cap {
        out[norm] = out[norm - 1] + u32::from(loeschian(norm, least));
    }
    out
}

pub fn new_rule(scale: usize, least: &[u32]) -> usize {
    let cap = WINDOW * scale * scale;
    let primes: Vec<usize> = factorise(scale, least).iter().map(|(p, _)| *p).collect();
    (1..=cap)
        .filter(|norm| loeschian(*norm, least))
        .filter(|norm| primes.iter().all(|prime| norm % (prime * prime) != 0))
        .count()
}

pub fn mobius_count(scale: usize, prefix: &[u32], least: &[u32], skip: Option<usize>) -> usize {
    let primes: Vec<usize> = factorise(scale, least)
        .iter()
        .map(|(p, _)| *p)
        .filter(|p| Some(*p) != skip)
        .collect();
    let cap = WINDOW * scale * scale;
    let mut total = 0i64;
    for mask in 0..1usize << primes.len() {
        let mut divisor = 1usize;
        let mut sign = 1i64;
        for (index, prime) in primes.iter().enumerate() {
            if mask >> index & 1 == 1 {
                divisor *= prime;
                sign = -sign;
            }
        }
        total += sign * i64::from(prefix[cap / (divisor * divisor)]);
    }
    total as usize
}

pub fn constant_from_primes(cap: usize, least: &[u32]) -> f64 {
    let mut log_b = 0.0f64;
    for prime in (2..=cap).filter(|v| least[*v] as usize == *v && v % 3 == 2) {
        log_b -= (-1.0 / (prime as f64 * prime as f64)).ln_1p();
    }
    (log_b.exp() / (2.0 * 3f64.sqrt())).sqrt()
}

pub fn target(scale: usize, least: &[u32]) -> f64 {
    WINDOW as f64 / 2f64.sqrt() * HEX_CONSTANT * crate::gaussian::jordan(scale, least)
}

pub fn lattice_norms(scale: usize) -> Vec<u64> {
    let cap = (WINDOW * scale * scale) as i64;
    let reach = 2 * scale as i64 + 1;
    let mut seen = vec![false; cap as usize + 1];
    for a in -reach..=reach {
        for b in -reach..=reach {
            let norm = a * a + a * b + b * b;
            if norm >= 1 && norm <= cap {
                seen[norm as usize] = true;
            }
        }
    }
    (1..=cap as usize)
        .filter(|norm| seen[*norm])
        .map(|norm| norm as u64)
        .collect()
}

pub fn union_counts(top: usize) -> Vec<usize> {
    let mut seen: HashSet<(u64, u64)> = HashSet::new();
    let mut out = Vec::new();
    for scale in 1..=top {
        let square = (scale * scale) as u64;
        let mut fresh = 0;
        for norm in lattice_norms(scale) {
            if seen.insert(reduce(norm, square)) {
                fresh += 1;
            }
        }
        out.push(fresh);
    }
    out
}
