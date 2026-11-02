use crate::coprime::{bracket, foreign, index};
use mrly_math::bang::baseq::{group_order, representatives};
use mrly_math::bang::factory::{code_to_corners, residue_corners};
use mrly_math::bang::Code;
use std::collections::BTreeMap;
use std::f64::consts::PI;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

const COSET_PRIME: u64 = 3;
const COSET_BASE: usize = 7;
const COSET_TOP: u32 = 7;
const TERMS_BUDGET: u64 = 20_000_000;
const SPLIT: u32 = 2;

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

fn walk(corners: &[Vec<u8>], base: u64, depth: u32, x: [u64; 4], sum: &mut u64) {
    let dim = corners[0].len();
    if depth == 0 {
        if x[..dim].iter().fold(0, |g, &v| gcd(g, v)) == 1 {
            *sum += 1;
        }
        return;
    }
    for corner in corners {
        let mut next = x;
        for (axis, digit) in corner.iter().enumerate() {
            next[axis] = x[axis] * base + u64::from(*digit);
        }
        walk(corners, base, depth - 1, next, sum);
    }
}

pub fn count(corners: &[Vec<u8>], base: u64, level: u32, threads: usize) -> u64 {
    let k = corners.len();
    let dim = corners[0].len();
    let split = SPLIT.min(level);
    let tasks = k.pow(split);
    let next = AtomicUsize::new(0);
    let total = AtomicU64::new(0);
    std::thread::scope(|scope| {
        for _ in 0..threads.max(1) {
            scope.spawn(|| loop {
                let task = next.fetch_add(1, Ordering::Relaxed);
                if task >= tasks {
                    break;
                }
                let mut x = [0u64; 4];
                let mut rest = task;
                for _ in 0..split {
                    let corner = &corners[rest % k];
                    rest /= k;
                    for axis in 0..dim {
                        x[axis] = x[axis] * base + u64::from(corner[axis]);
                    }
                }
                let mut sum = 0;
                walk(corners, base, level - split, x, &mut sum);
                total.fetch_add(sum, Ordering::Relaxed);
            });
        }
    });
    total.into_inner()
}

fn threads() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(8)
}

pub fn terms(base: usize, dimension: usize, code: Code) {
    let corners = code_to_corners(code, dimension, base).expect("the code fits its cells");
    let k = corners.len() as u64;
    let predicted = bracket(&corners, base) * foreign(base, dimension);
    println!(
        "base {base} dim {dimension} code {code} fill {k} index {} bracket {:.6} predicted {:.6}",
        index(&corners, dimension),
        bracket(&corners, base),
        predicted
    );
    let mut level = 0;
    while k.pow(level + 1) <= TERMS_BUDGET {
        level += 1;
        let a = count(&corners, base as u64, level, threads());
        println!(
            "n {level:2} A {a:>10} A/fill^n {:.6} gap {:+.2e}",
            a as f64 / (k as f64).powi(level as i32),
            a as f64 / (k as f64).powi(level as i32) - predicted
        );
    }
}

fn coset_share(corners: &[Vec<u8>], level: u32) -> f64 {
    let p = COSET_PRIME as usize;
    let mut dist = vec![0f64; p * p];
    dist[0] = 1.0;
    for _ in 0..level {
        let mut next = vec![0f64; p * p];
        for (cell, weight) in dist.iter().enumerate() {
            if *weight == 0.0 {
                continue;
            }
            for corner in corners {
                let i = (cell / p + usize::from(corner[0])) % p;
                let j = (cell % p + usize::from(corner[1])) % p;
                next[i * p + j] += weight / corners.len() as f64;
            }
        }
        dist = next;
    }
    dist[0]
}

pub fn limits() {
    let corners: Vec<Vec<u8>> = residue_corners(2, COSET_BASE)
        .into_iter()
        .filter(|v| (usize::from(v[0]) + usize::from(v[1])) % 3 == 1)
        .collect();
    let k = corners.len() as u64;
    let rest = (6.0 / (PI * PI))
        / ((1.0 - 1.0 / (COSET_PRIME * COSET_PRIME) as f64)
            * (1.0 - 1.0 / (COSET_BASE * COSET_BASE) as f64));
    println!(
        "base {COSET_BASE} dim 2 corners with v_1 + v_2 = 1 mod {COSET_PRIME}: fill {k} index {} bracket {:.4} naive delta {:.6}",
        index(&corners, 2),
        bracket(&corners, COSET_BASE),
        rest * (1.0 - 1.0 / (COSET_PRIME * COSET_PRIME) as f64)
    );
    println!("corrected delta_n = (1 - q_n) prod_(p not in {{3, 7}}) (1 - p^-2), q_n the share of level-n points that are 0 mod 3");
    for level in 1..=COSET_TOP {
        let a = count(&corners, COSET_BASE as u64, level, threads());
        let ratio = a as f64 / (k as f64).powi(level as i32);
        let q = coset_share(&corners, level);
        let corrected = (1.0 - q) * rest;
        println!(
            "n {level} n mod 3 {} A {a:>12} A/fill^n {ratio:.6} q_n {q:.6} corrected {corrected:.6} gap {:+.2e}",
            level % 3,
            ratio - corrected
        );
    }
    println!(
        "limits along n mod 3: 1 and 2 at {:.6}, 0 at {:.6}",
        rest,
        rest * (1.0 - 1.0 / COSET_PRIME as f64)
    );
}

pub fn dim4() {
    let reps = representatives(2, 4).expect("the walk stays under the code limit");
    let order = group_order(2, 4).expect("dimension four is under the factorial ceiling");
    println!(
        "base 2 dim 4: {} designs in {} orbits under the {} signed coordinate permutations, least code per orbit",
        1u64 << 16,
        reps.len(),
        order
    );
    let mut rows: Vec<(Code, usize, bool, Vec<u64>)> = Vec::new();
    for (code, _) in reps {
        let corners = code_to_corners(code, 4, 2).expect("the code fits its cells");
        if corners.len() < 2 {
            continue;
        }
        let terms: Vec<u64> = (1..=4).map(|n| count(&corners, 2, n, 1)).collect();
        rows.push((code, corners.len(), index(&corners, 4) == 1, terms));
    }
    let theorem = rows.iter().filter(|r| r.2 && r.1 > 2).count();
    println!(
        "fill >= 2 {}, spanning with fill > 2 {}",
        rows.len(),
        theorem
    );
    let mut groups: BTreeMap<&Vec<u64>, Vec<Code>> = BTreeMap::new();
    for row in &rows {
        groups.entry(&row.3).or_default().push(row.0);
    }
    let mut shared: Vec<(&Vec<u64>, usize)> = groups
        .iter()
        .filter(|(_, codes)| codes.len() > 1)
        .map(|(terms, codes)| (*terms, codes.len()))
        .collect();
    shared.sort_by(|p, q| q.1.cmp(&p.1).then(p.0.cmp(q.0)));
    println!(
        "distinct (A(1)..A(4)) {}, collisions {} covering {} classes",
        groups.len(),
        shared.len(),
        shared.iter().map(|s| s.1).sum::<usize>()
    );
    for (terms, size) in shared.iter().take(4) {
        println!("  {size} classes on {terms:?}");
    }
}
