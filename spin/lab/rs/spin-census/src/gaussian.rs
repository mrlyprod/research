use std::collections::HashSet;

pub const LANDAU_RAMANUJAN: f64 = 0.764_223_653_589_220_66;

pub fn least_factors(limit: usize) -> Vec<u32> {
    let mut out = vec![0u32; limit + 1];
    for value in 2..=limit {
        if out[value] == 0 {
            let mut multiple = value;
            while multiple <= limit {
                if out[multiple] == 0 {
                    out[multiple] = value as u32;
                }
                multiple += value;
            }
        }
    }
    out
}

pub fn factorise(mut value: usize, least: &[u32]) -> Vec<(usize, usize)> {
    let mut out: Vec<(usize, usize)> = Vec::new();
    while value > 1 {
        let prime = least[value] as usize;
        let mut power = 0;
        while value.is_multiple_of(prime) {
            value /= prime;
            power += 1;
        }
        out.push((prime, power));
    }
    out
}

pub fn two_square(value: usize, least: &[u32]) -> bool {
    factorise(value, least)
        .iter()
        .all(|(prime, power)| prime % 4 != 3 || power % 2 == 0)
}

pub fn primitive_norm(value: usize, least: &[u32]) -> bool {
    if value == 1 || value == 2 {
        return true;
    }
    let parts = factorise(value, least);
    parts.iter().all(|(prime, power)| match prime % 4 {
        1 => true,
        2 => *power == 1,
        _ => false,
    })
}

pub fn new_disc(scale: usize, least: &[u32]) -> usize {
    let cap = 2 * scale * scale;
    let primes: Vec<usize> = factorise(scale, least).iter().map(|(p, _)| *p).collect();
    (1..=cap)
        .filter(|norm| two_square(*norm, least))
        .filter(|norm| primes.iter().all(|prime| norm % (prime * prime) != 0))
        .count()
}

pub fn two_square_prefix(cap: usize, least: &[u32]) -> Vec<usize> {
    let mut out = vec![0usize; cap + 1];
    for norm in 1..=cap {
        out[norm] = out[norm - 1] + usize::from(two_square(norm, least));
    }
    out
}

pub fn mobius_count(scale: usize, prefix: &[usize], least: &[u32]) -> usize {
    let primes: Vec<usize> = factorise(scale, least).iter().map(|(p, _)| *p).collect();
    let cap = 2 * scale * scale;
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
        total += sign * prefix[cap / (divisor * divisor)] as i64;
    }
    total as usize
}

pub fn jordan(scale: usize, least: &[u32]) -> f64 {
    factorise(scale, least)
        .iter()
        .map(|(prime, _)| 1.0 - 1.0 / (prime * prime) as f64)
        .product()
}

pub fn radical(value: usize, least: &[u32]) -> usize {
    factorise(value, least)
        .iter()
        .map(|(prime, _)| prime)
        .product()
}

pub fn landau_target(scale: usize, least: &[u32]) -> f64 {
    2f64.sqrt() * LANDAU_RAMANUJAN * jordan(scale, least)
}

pub fn scaled(scale: usize, count: usize) -> f64 {
    count as f64 * (scale as f64).ln().sqrt() / (scale * scale) as f64
}

pub fn reduce(mut top: u64, mut bottom: u64) -> (u64, u64) {
    let mut a = top;
    let mut b = bottom;
    while b != 0 {
        let next = a % b;
        a = b;
        b = next;
    }
    top /= a;
    bottom /= a;
    (top, bottom)
}

pub fn disc_norms(scale: usize) -> Vec<u64> {
    let cap = 2 * scale * scale;
    let mut seen = vec![false; cap + 1];
    let reach = (cap as f64).sqrt() as usize + 1;
    for a in 0..=reach {
        for b in 0..=reach {
            let norm = a * a + b * b;
            if norm >= 1 && norm <= cap {
                seen[norm] = true;
            }
        }
    }
    (1..=cap)
        .filter(|norm| seen[*norm])
        .map(|norm| norm as u64)
        .collect()
}

pub fn box_norms(scale: usize) -> Vec<u64> {
    let cap = 2 * scale * scale;
    let mut seen = vec![false; cap + 1];
    for a in 0..=scale {
        for b in 0..=scale {
            let norm = a * a + b * b;
            if norm >= 1 {
                seen[norm] = true;
            }
        }
    }
    (1..=cap)
        .filter(|norm| seen[*norm])
        .map(|norm| norm as u64)
        .collect()
}

pub fn union_counts(top: usize, boxed: bool) -> Vec<usize> {
    let mut seen: HashSet<(u64, u64)> = HashSet::new();
    let mut out = Vec::new();
    for scale in 1..=top {
        let norms = if boxed {
            box_norms(scale)
        } else {
            disc_norms(scale)
        };
        let square = (scale * scale) as u64;
        let mut fresh = 0;
        for norm in norms {
            if seen.insert(reduce(norm, square)) {
                fresh += 1;
            }
        }
        out.push(fresh);
    }
    out
}
