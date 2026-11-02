use crate::factor::mobius;

// FAMILIES

pub struct Family {
    pub q: u64,
    pub digits: Vec<u64>,
    pub label: String,
    pub lmax: usize,
    pub children: Vec<u64>,
}

pub struct Outcome {
    pub counts: Vec<u64>,
    pub meter: Vec<i64>,
    pub mmax: Vec<u64>,
    pub twisted: Vec<(u64, Vec<i64>)>,
}

fn depth(q: u64, k: usize) -> usize {
    match (q, k) {
        (3, 2) => 24,
        (4, 2) => 22,
        (4, 3) => 14,
        (5, 2) => 21,
        (5, 3) => 13,
        (5, 4) => 11,
        _ => panic!("no depth for q={q} k={k}"),
    }
}

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

pub fn digit_gcd(digits: &[u64]) -> u64 {
    digits.iter().fold(0, |g, &d| gcd(g, d))
}

pub fn make_family(q: u64, digits: Vec<u64>) -> Family {
    let lmax = depth(q, digits.len());
    make_family_depth(q, digits, lmax)
}

pub fn make_family_depth(q: u64, digits: Vec<u64>, lmax: usize) -> Family {
    let label: String = digits.iter().map(|d| d.to_string()).collect();
    let mut children = Vec::new();
    if digit_gcd(&digits) == 1 {
        let top = *digits.iter().max().unwrap();
        let mut a = 2;
        while a * top < q {
            children.push(a);
            a += 1;
        }
    }
    Family {
        q,
        digits,
        label,
        lmax,
        children,
    }
}

pub fn families() -> Vec<Family> {
    let mut out = Vec::new();
    for q in [3u64, 4, 5] {
        for k in 2..=(q as usize - 1) {
            for mask in 0u64..(1 << q) {
                if mask.count_ones() as usize != k {
                    continue;
                }
                let digits: Vec<u64> = (0..q).filter(|d| mask >> d & 1 == 1).collect();
                out.push(make_family(q, digits));
            }
        }
    }
    out
}

// SWEEP

pub fn sweep(q: u64, digits: &[u64], lmax: usize, visit: &mut impl FnMut(u64, usize)) {
    fn rec(
        v: u64,
        len: usize,
        q: u64,
        digits: &[u64],
        lmax: usize,
        visit: &mut impl FnMut(u64, usize),
    ) {
        visit(v, len);
        if len < lmax {
            for &d in digits {
                rec(v * q + d, len + 1, q, digits, lmax, visit);
            }
        }
    }
    for &d in digits {
        if d != 0 {
            rec(d, 1, q, digits, lmax, visit);
        }
    }
}

fn twist(a: u64, v: u64, mu_v: i8) -> i64 {
    match a {
        2 => {
            if v.is_multiple_of(2) {
                0
            } else {
                -(mu_v as i64)
            }
        }
        3 => {
            if v.is_multiple_of(3) {
                0
            } else {
                -(mu_v as i64)
            }
        }
        4 => 0,
        _ => panic!("no twist for a={a}"),
    }
}

// ORDERED SWEEP

pub fn sweep_length(q: u64, digits: &[u64], length: usize, visit: &mut impl FnMut(u64)) {
    fn rec(v: u64, len: usize, q: u64, digits: &[u64], length: usize, visit: &mut impl FnMut(u64)) {
        if len == length {
            visit(v);
            return;
        }
        for &d in digits {
            rec(v * q + d, len + 1, q, digits, length, visit);
        }
    }
    for &d in digits {
        if d != 0 {
            rec(d, 1, q, digits, length, visit);
        }
    }
}

// METER

pub fn run_family(fam: &Family, primes: &[u64]) -> Outcome {
    let l = fam.lmax;
    let has01 = fam.digits.contains(&0) && fam.digits.contains(&1);
    let mut counts = vec![0u64; l + 1];
    let mut meter = vec![0i64; l + 1];
    let mut mmax = vec![0u64; l + 1];
    let mut twisted: Vec<(u64, Vec<i64>)> = fam
        .children
        .iter()
        .map(|&a| (a, vec![0i64; l + 1]))
        .collect();
    let mut run = 0i64;
    let mut peak = 0u64;
    let mut total = 0u64;
    let mut boundary = 1u64;
    for lev in 1..=l {
        let skip = boundary;
        boundary *= fam.q;
        let mut tw_len: Vec<i64> = vec![0; fam.children.len()];
        sweep_length(fam.q, &fam.digits, lev, &mut |v| {
            let mu = mobius(v, primes);
            for (slot, a) in tw_len.iter_mut().zip(fam.children.iter()) {
                *slot += twist(*a, v, mu);
            }
            if has01 && lev > 1 && v == skip {
                return;
            }
            total += 1;
            run += mu as i64;
            peak = peak.max(run.unsigned_abs());
        });
        if has01 {
            total += 1;
            run += mobius(boundary, primes) as i64;
            peak = peak.max(run.unsigned_abs());
        }
        counts[lev] = total;
        meter[lev] = run;
        mmax[lev] = peak;
        for (slot, tl) in twisted.iter_mut().zip(tw_len.iter()) {
            slot.1[lev] = slot.1[lev - 1] + tl;
        }
    }
    Outcome {
        counts,
        meter,
        mmax,
        twisted,
    }
}
