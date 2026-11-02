use crate::tile::line_reducible;
use num_bigint::BigUint;

pub fn total(bits: usize) -> BigUint {
    (BigUint::from(1u32) << bits) - BigUint::from(1u32)
}

pub fn tile_total(side: usize) -> BigUint {
    total(side * side)
}

pub fn line_total(side: usize) -> BigUint {
    total(side)
}

pub fn irreducibles(totals: &[BigUint]) -> Vec<BigUint> {
    let mut out: Vec<BigUint> = Vec::new();
    for k in 0..totals.len() {
        let mut value = totals[k].clone();
        for j in 0..k {
            value -= &out[j] * &totals[k - 1 - j];
        }
        out.push(value);
    }
    out
}

pub fn reducible_at(prime: usize, power: usize, plane: bool) -> BigUint {
    let totals: Vec<BigUint> = (1..=power)
        .map(|k| {
            let side = prime.pow(k as u32);
            if plane {
                tile_total(side)
            } else {
                line_total(side)
            }
        })
        .collect();
    let irr = irreducibles(&totals);
    &totals[power - 1] - &irr[power - 1]
}

pub fn line_brute(side: usize) -> usize {
    (1u128..1u128 << side)
        .filter(|mask| line_reducible(*mask, side))
        .count()
}
