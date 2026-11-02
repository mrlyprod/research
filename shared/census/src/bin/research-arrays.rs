use census::save;
use mrly_core::error::Result;
use mrly_core::json;
use std::collections::HashSet;

const NAME: &str = "research-arrays";
const G: [i64; 4] = [0, 1, 2, 4];
const ORDER: u32 = 4;

fn coarray(s: &[i64]) -> HashSet<i64> {
    s.iter()
        .flat_map(|x| s.iter().map(move |y| x - y))
        .collect()
}

fn paired(g: &[i64]) -> Vec<i64> {
    g.iter()
        .copied()
        .filter(|&x| {
            g.iter().any(|&y| {
                y != x
                    && g.iter()
                        .flat_map(|p| g.iter().map(move |q| p - q))
                        .filter(|&d| d == x - y)
                        .count()
                        == 1
            })
        })
        .collect()
}

fn main() -> Result<()> {
    let a = *G.iter().max().unwrap();
    let base = 2 * a + 1;
    assert_eq!(coarray(&G).len() as i64, base);
    let digits: Vec<Vec<i64>> = (0..G.len().pow(ORDER))
        .map(|mut k| {
            (0..ORDER)
                .map(|_| {
                    let d = G[k % G.len()];
                    k /= G.len();
                    d
                })
                .collect()
        })
        .collect();
    let sensors: Vec<i64> = digits
        .iter()
        .map(|d| d.iter().rev().fold(0, |acc, &x| acc * base + x))
        .collect();
    let full = coarray(&sensors);
    assert_eq!(full.len() as i64, base.pow(ORDER));
    let essential: Vec<bool> = sensors
        .iter()
        .map(|&s| {
            let rest: Vec<i64> = sensors.iter().copied().filter(|&x| x != s).collect();
            coarray(&rest) != full
        })
        .collect();
    let u = paired(&G);
    assert_eq!(u, vec![0, 1, 4]);
    for (d, &e) in digits.iter().zip(&essential) {
        assert_eq!(e, d.iter().all(|x| u.contains(x)));
    }
    assert_eq!(essential.iter().filter(|&&e| e).count(), u.len().pow(ORDER));
    let side = a + base * a + 1;
    let cells: Vec<[i64; 3]> = digits
        .iter()
        .zip(&essential)
        .map(|(d, &e)| [d[0] + base * d[2], d[1] + base * d[3], e as i64])
        .collect();
    save(NAME, &json!({"side": side, "cells": cells}))
}
