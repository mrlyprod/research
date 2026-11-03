use crate::radii::{self, Live};
use mrly_num::factor::gcd;
use std::collections::HashMap;

const GASKET: [(u64, u64); 3] = [(0, 0), (1, 0), (0, 1)];
const FREE: [(u64, u64); 9] = [
    (0, 0),
    (1, 0),
    (2, 0),
    (0, 1),
    (1, 1),
    (2, 1),
    (0, 2),
    (1, 2),
    (2, 2),
];
const THETA: [i128; 4] = [-2, 0, -1, 1];
const MATCH: f64 = 1e-9;
const PRIME_BITS: f64 = 61.0;
const RATIO_LENGTH: usize = 20;
const BOUND_LENGTH: usize = 40;
const CLOSED_LENGTH: usize = 60;
const CROSS_STATES: usize = 20;

#[derive(Clone, Copy, PartialEq)]
pub enum Digits {
    Gasket,
    Free,
}

impl Digits {
    fn name(self) -> &'static str {
        match self {
            Digits::Gasket => "A",
            Digits::Free => "B",
        }
    }

    fn alphabet(self) -> &'static [(u64, u64)] {
        match self {
            Digits::Gasket => &GASKET,
            Digits::Free => &FREE,
        }
    }
}

fn encode(c: [u64; 4]) -> i64 {
    (c[0] | c[1] << 8 | c[2] << 16 | c[3] << 24) as i64
}

pub fn build(s: u64, t: u64, digits: Digits) -> (Live, Live) {
    let mut index = HashMap::new();
    let mut order: Vec<[u64; 4]> = vec![[0; 4]];
    let mut edges: Vec<Vec<usize>> = Vec::new();
    index.insert([0u64; 4], 0usize);
    let mut i = 0;
    while i < order.len() {
        let c = order[i];
        let mut out = Vec::new();
        for &(d1, d2) in digits.alphabet() {
            let x = [s * d1 + c[0], s * d2 + c[1], t * d1 + c[2], t * d2 + c[3]];
            if !GASKET.contains(&(x[0] % 3, x[1] % 3)) || !GASKET.contains(&(x[2] % 3, x[3] % 3)) {
                continue;
            }
            let next = [x[0] / 3, x[1] / 3, x[2] / 3, x[3] / 3];
            let j = *index.entry(next).or_insert_with(|| {
                order.push(next);
                order.len() - 1
            });
            out.push(j);
        }
        edges.push(out);
        i += 1;
    }
    let mut back = vec![Vec::new(); order.len()];
    for (i, out) in edges.iter().enumerate() {
        for &j in out {
            back[j].push(i);
        }
    }
    let mut seen = vec![false; order.len()];
    seen[0] = true;
    let mut stack = vec![0usize];
    while let Some(x) = stack.pop() {
        for &y in &back[x] {
            if !seen[y] {
                seen[y] = true;
                stack.push(y);
            }
        }
    }
    let keep: Vec<usize> = (0..order.len()).filter(|&i| seen[i]).collect();
    let mut place = vec![usize::MAX; order.len()];
    for (k, &i) in keep.iter().enumerate() {
        place[i] = k;
    }
    let core = Live {
        states: keep.iter().map(|&i| encode(order[i])).collect(),
        edges: keep
            .iter()
            .map(|&i| {
                edges[i]
                    .iter()
                    .filter(|&&j| seen[j])
                    .map(|&j| place[j])
                    .collect()
            })
            .collect(),
    };
    let reach = Live {
        states: order.iter().map(|&c| encode(c)).collect(),
        edges,
    };
    (reach, core)
}

fn component_rho(local: &[Vec<usize>]) -> (f64, f64) {
    let m = local.len();
    let mut v = vec![1.0f64; m];
    let (mut lo, mut hi) = (0.0f64, f64::INFINITY);
    for _ in 0..200_000 {
        let mut w: Vec<f64> = v.clone();
        for (i, out) in local.iter().enumerate() {
            for &j in out {
                w[i] += v[j];
            }
        }
        lo = f64::INFINITY;
        hi = 0.0;
        for i in 0..m {
            let r = w[i] / v[i];
            lo = lo.min(r);
            hi = hi.max(r);
        }
        let norm: f64 = w.iter().sum();
        for x in w.iter_mut() {
            *x /= norm;
        }
        v = w;
        if hi - lo < 1e-13 {
            break;
        }
    }
    ((lo + hi) / 2.0 - 1.0, hi - lo)
}

pub fn radius(edges: &[Vec<usize>]) -> (f64, f64) {
    let mut best = (0.0f64, 0.0f64);
    for comp in radii::components(edges) {
        let place: HashMap<usize, usize> = comp.iter().enumerate().map(|(k, &v)| (v, k)).collect();
        let local: Vec<Vec<usize>> = comp
            .iter()
            .map(|&v| {
                edges[v]
                    .iter()
                    .filter_map(|w| place.get(w).copied())
                    .collect()
            })
            .collect();
        if local.iter().all(|out| out.is_empty()) {
            continue;
        }
        let (rho, width) = component_rho(&local);
        if rho > best.0 {
            best = (rho, width);
        }
    }
    best
}

fn walks(live: &Live, length: usize) -> Vec<Vec<u128>> {
    let n = live.edges.len();
    let mut rows = vec![vec![1u128; n]];
    for w in 1..=length {
        let prev = &rows[w - 1];
        let next: Vec<u128> = live
            .edges
            .iter()
            .map(|out| out.iter().map(|&j| prev[j]).sum())
            .collect();
        rows.push(next);
    }
    rows
}

fn closed(live: &Live, length: usize) -> u128 {
    let n = live.edges.len();
    let mut cur = vec![0u128; n];
    cur[0] = 1;
    for _ in 0..length {
        cur = live
            .edges
            .iter()
            .map(|out| out.iter().map(|&j| cur[j]).sum())
            .collect();
    }
    cur[0]
}

struct Big(Vec<u64>);

impl Big {
    fn mul_add(&mut self, m: u64, a: u64) {
        let mut carry = a as u128;
        for limb in self.0.iter_mut() {
            let x = *limb as u128 * m as u128 + carry;
            *limb = x as u64;
            carry = x >> 64;
        }
        if carry > 0 {
            self.0.push(carry as u64);
        }
    }

    fn below(&self, other: &Big) -> bool {
        let top = |v: &[u64]| v.iter().rposition(|&x| x != 0).map_or(0, |i| i + 1);
        let (a, b) = (top(&self.0), top(&other.0));
        if a != b {
            return a < b;
        }
        for i in (0..a).rev() {
            if self.0[i] != other.0[i] {
                return self.0[i] < other.0[i];
            }
        }
        false
    }
}

fn nonnegative(residues: &[u64], ps: &[u64]) -> bool {
    let m = ps.len();
    let mut digits = vec![0u64; m];
    for i in 0..m {
        let p = ps[i];
        let mut t = residues[i] % p;
        for j in 0..i {
            let inv = radii::powmod(ps[j] % p, p - 2, p);
            t = radii::mulmod(radii::sub(t, digits[j] % p, p), inv, p);
        }
        digits[i] = t;
    }
    let mut value = Big(vec![digits[m - 1]]);
    for i in (0..m - 1).rev() {
        value.mul_add(ps[i], digits[i]);
    }
    let mut total = Big(vec![1]);
    for &p in ps {
        total.mul_add(p, 0);
    }
    value.mul_add(2, 0);
    value.below(&total)
}

fn shift_mod(c: &[u64], a: u64, p: u64) -> Vec<u64> {
    let mut c = c.to_vec();
    let n = c.len() - 1;
    for k in 0..n {
        for j in (k..n).rev() {
            c[j] = (c[j] + radii::mulmod(a, c[j + 1], p)) % p;
        }
    }
    c
}

fn primes_for(n: usize) -> usize {
    let bits = n as f64 * (1.0 + 90f64.sqrt()).log2() + 8.0;
    (bits / PRIME_BITS).ceil() as usize + 1
}

#[derive(Clone, Copy)]
pub struct Cert {
    pub at3: bool,
    pub le2: bool,
    pub eq2: bool,
    pub theta: bool,
}

pub fn certify(live: &Live, ps: &[u64]) -> Cert {
    let n = live.edges.len();
    let ps = &ps[..primes_for(n).max(2)];
    let polys: Vec<Vec<u64>> = ps.iter().map(|&p| radii::charpoly_mod(live, p)).collect();
    let shifted = |a: u64| {
        let qs: Vec<Vec<u64>> = polys
            .iter()
            .zip(ps)
            .map(|(c, &p)| shift_mod(c, a, p))
            .collect();
        let zero = qs.iter().all(|q| q[0] == 0);
        let nonneg = (0..=n).all(|j| {
            let residues: Vec<u64> = qs.iter().map(|q| q[j]).collect();
            nonnegative(&residues, ps)
        });
        (zero, nonneg)
    };
    let (zero3, nonneg3) = shifted(3);
    let (zero2, nonneg2) = shifted(2);
    let theta = n >= 3
        && polys
            .iter()
            .zip(ps)
            .all(|(c, &p)| radii::remainder_mod(c, &THETA, p).iter().all(|&x| x == 0));
    Cert {
        at3: zero3 && nonneg3,
        le2: nonneg2,
        eq2: nonneg2 && zero2,
        theta,
    }
}

fn shift_exact(c: &[i128], a: i128) -> Vec<i128> {
    let mut c = c.to_vec();
    let n = c.len() - 1;
    for k in 0..n {
        for j in (k..n).rev() {
            c[j] += a * c[j + 1];
        }
    }
    c
}

fn cross_check(live: &Live, ps: &[u64], cert: Cert) -> Option<bool> {
    let n = live.edges.len();
    if !(2..=CROSS_STATES).contains(&n) {
        return None;
    }
    let exact = radii::charpoly_exact(live, ps)?;
    let verdict = |a: i128| {
        let q = shift_exact(&exact, a);
        (q[0] == 0, q.iter().all(|&x| x >= 0))
    };
    let (zero3, nonneg3) = verdict(3);
    let (zero2, nonneg2) = verdict(2);
    let theta = radii::divides(&THETA, &exact);
    Some(
        cert.at3 == (zero3 && nonneg3)
            && cert.le2 == nonneg2
            && cert.eq2 == (nonneg2 && zero2)
            && cert.theta == theta,
    )
}

pub struct Row {
    pub s: u64,
    pub t: u64,
    pub live: usize,
    pub reach: usize,
    pub rho: f64,
    pub width: f64,
    pub reach_rho: f64,
    pub cert: Cert,
    pub k: u32,
    pub ratio: f64,
    pub holds: bool,
    pub agrees: Option<bool>,
}

fn is_shift(s: u64, t: u64) -> bool {
    s == 1 && t == 3u64.pow(radii::v3(t))
}

fn k_of(s: u64, t: u64) -> u32 {
    let (s1, t1) = (s / 3u64.pow(radii::v3(s)), t / 3u64.pow(radii::v3(t)));
    radii::v3(s * t) + radii::v3(s1.abs_diff(t1))
}

pub fn row(s: u64, t: u64, digits: Digits, ps: &[u64]) -> Row {
    let (reach, live) = build(s, t, digits);
    let (rho, width) = radius(&live.edges);
    let (reach_rho, _) = radius(&reach.edges);
    let cert = certify(&live, ps);
    let agrees = cross_check(&live, ps, cert);
    let k = k_of(s, t);
    let p = walks(&live, BOUND_LENGTH);
    let holds = is_shift(s, t)
        || p.iter()
            .enumerate()
            .all(|(w, counts)| counts.iter().all(|&c| (c << k) <= (3u128.pow(k) << w)));
    let top = *p[RATIO_LENGTH].iter().max().unwrap() as f64;
    let ratio = top * 2f64.powi(k as i32) / (3f64.powi(k as i32) * 2f64.powi(RATIO_LENGTH as i32));
    Row {
        s,
        t,
        live: live.edges.len(),
        reach: reach.edges.len(),
        rho,
        width,
        reach_rho,
        cert,
        k,
        ratio,
        holds,
        agrees,
    }
}

fn parallel<T: Send>(
    pairs: &[(u64, u64)],
    threads: usize,
    f: impl Fn(u64, u64) -> T + Sync,
) -> Vec<T> {
    let mut tagged: Vec<(usize, T)> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|lane| {
                let f = &f;
                scope.spawn(move || {
                    pairs
                        .iter()
                        .enumerate()
                        .skip(lane)
                        .step_by(threads)
                        .map(|(i, &(s, t))| (i, f(s, t)))
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap())
            .collect()
    });
    tagged.sort_by_key(|(i, _)| *i);
    tagged.into_iter().map(|(_, v)| v).collect()
}

fn theta() -> f64 {
    let (mut lo, mut hi) = (1.0f64, 2.0f64);
    for _ in 0..200 {
        let mid = (lo + hi) / 2.0;
        if mid * mid * mid - mid * mid - 2.0 > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    (lo + hi) / 2.0
}

fn list(rows: &[&Row]) -> String {
    rows.iter()
        .map(|r| format!("({},{})", r.s, r.t))
        .collect::<Vec<_>>()
        .join(" ")
}

fn describe(r: &Row) -> String {
    format!(
        "({},{}) live {} reachable {} rho {:.12} width {:.1e} K {} ratio {:.7}{}{}{}{}",
        r.s,
        r.t,
        r.live,
        r.reach,
        r.rho,
        r.width,
        r.k,
        r.ratio,
        if r.cert.at3 { " exactly-3" } else { "" },
        if r.cert.eq2 { " exactly-2" } else { "" },
        if r.cert.le2 && !r.cert.eq2 {
            " below-2"
        } else {
            ""
        },
        if r.cert.theta { " theta-divides" } else { "" },
    )
}

fn summarize(digits: Digits, rows: &[Row]) -> Vec<(u64, u64)> {
    let name = digits.name();
    let theta = theta();
    let live_max = rows.iter().map(|r| r.live).max().unwrap();
    let reach_max = rows.iter().map(|r| r.reach).max().unwrap();
    let largest: Vec<&Row> = rows.iter().filter(|r| r.live == live_max).collect();
    let widest: Vec<&Row> = rows.iter().filter(|r| r.reach == reach_max).collect();
    let worst = rows.iter().map(|r| r.width).fold(0.0f64, f64::max);
    println!(
        "{} live states (reachable from the start and returning to it) max {} at {}, reachable states max {} at {}, worst Perron bracket {:.1e}",
        name,
        live_max,
        list(&largest),
        reach_max,
        list(&widest),
        worst
    );
    let differ: Vec<&Row> = rows
        .iter()
        .filter(|r| (r.rho - r.reach_rho).abs() > MATCH)
        .collect();
    println!(
        "{} reachable states that never return carry cycles above the live radius on {} pairs {}",
        name,
        differ.len(),
        list(&differ)
    );
    let at3: Vec<&Row> = rows.iter().filter(|r| r.cert.at3).collect();
    let eq2: Vec<&Row> = rows.iter().filter(|r| r.cert.eq2).collect();
    let below: Vec<&Row> = rows.iter().filter(|r| r.cert.le2 && !r.cert.eq2).collect();
    let open: Vec<&Row> = rows.iter().filter(|r| !r.cert.le2 && !r.cert.at3).collect();
    println!(
        "{} radius exactly 3 (charpoly root 3, shift by 3 nonnegative): {} {}",
        name,
        at3.len(),
        list(&at3)
    );
    println!(
        "{} radius exactly 2 (charpoly root 2, shift by 2 nonnegative): {} {}",
        name,
        eq2.len(),
        list(&eq2)
    );
    println!(
        "{} radius below 2 (shift by 2 nonnegative, 2 no root): {}, neither certified, so in (2,3) or above: {} {}",
        name,
        below.len(),
        open.len(),
        list(&open)
    );
    let top = below.iter().map(|r| r.rho).fold(0.0f64, f64::max);
    let attainers: Vec<&Row> = below
        .iter()
        .filter(|r| (r.rho - top).abs() < MATCH)
        .copied()
        .collect();
    println!(
        "{} largest radius below 2: {:.12} at {}, x^3 - x^2 - 2 divides their charpolys {}",
        name,
        top,
        list(&attainers),
        attainers.iter().all(|r| r.cert.theta)
    );
    let reach: Vec<&Row> = below
        .iter()
        .filter(|r| r.rho > theta - MATCH)
        .copied()
        .collect();
    let above: Vec<&Row> = reach
        .iter()
        .filter(|r| r.rho > theta + MATCH)
        .copied()
        .collect();
    let at: Vec<&Row> = reach
        .iter()
        .filter(|r| (r.rho - theta).abs() <= MATCH && r.cert.theta)
        .copied()
        .collect();
    println!(
        "{} below 2 and at or above theta {:.14}: {}, strictly above {} (theta-divisible among them {}), exactly at theta with x^3 - x^2 - 2 dividing {}",
        name,
        theta,
        reach.len(),
        above.len(),
        above.iter().filter(|r| r.cert.theta).count(),
        at.len()
    );
    let checked: Vec<&Row> = rows.iter().filter(|r| r.agrees.is_some()).collect();
    println!(
        "{} modular certificates against the exact i128 charpoly shifted by 2 and 3 on the {} automata with 2 to {} live states: agree {}",
        name,
        checked.len(),
        CROSS_STATES,
        checked.iter().all(|r| r.agrees == Some(true))
    );
    let off: Vec<&Row> = rows.iter().filter(|r| !is_shift(r.s, r.t)).collect();
    let violations: Vec<&Row> = off.iter().filter(|r| !r.holds).copied().collect();
    println!(
        "{} bound P_w <= (3/2)^K 2^w at every live state to w = {} on {} non-shift pairs: violations {} {}",
        name,
        BOUND_LENGTH,
        off.len(),
        violations.len(),
        list(&violations)
    );
    let mut ranked: Vec<&Row> = off.clone();
    ranked.sort_by(|x, y| y.ratio.partial_cmp(&x.ratio).unwrap());
    println!(
        "{} closest ratios P_{} / ((3/2)^K 2^{}) over the live states:",
        name, RATIO_LENGTH, RATIO_LENGTH
    );
    for r in ranked.iter().take(8) {
        println!("  ({},{}) K {} ratio {:.7}", r.s, r.t, r.k, r.ratio);
    }
    println!(
        "{} rows at radius 3, at 2, and below 2 at or above theta:",
        name
    );
    for r in rows
        .iter()
        .filter(|r| r.cert.at3 || r.cert.eq2 || (r.cert.le2 && r.rho > theta - MATCH))
    {
        println!("  {}", describe(r));
    }
    eq2.iter().map(|r| (r.s, r.t)).collect()
}

pub fn report(height: u64, threads: usize) {
    let pairs: Vec<(u64, u64)> = (2..=height)
        .flat_map(|t| (1..t).map(move |s| (s, t)))
        .filter(|&(s, t)| gcd(s as u128, t as u128) == 1)
        .collect();
    println!(
        "multiplier pairs 1 <= s < t <= {} coprime unordered {}; A reads z over the gasket digits, B over all of {{0,1,2}}^2, both count z with sz and tz in the gasket",
        height,
        pairs.len()
    );
    let ps = radii::primes(primes_for(1000));
    println!(
        "certificates: charpoly modulo primes near 2^61, Taylor shift by 2 or 3 reconstructed by Garner under the Hadamard bound (1 + sqrt 90)^n, radius <= a iff every coefficient of p(x + a) is nonnegative, theta divisibility by modular remainder; match tolerance {:.0e}",
        MATCH
    );
    let mut sets = Vec::new();
    for digits in [Digits::Gasket, Digits::Free] {
        let rows = parallel(&pairs, threads, |s, t| row(s, t, digits, &ps));
        sets.push(summarize(digits, &rows));
    }
    println!("A and B attain 2 on the same pairs: {}", sets[0] == sets[1]);
    let (_, b413) = build(4, 13, Digits::Free);
    println!(
        "B(4,13) closed paths at n = {}: {}, theta^{} = {:.3e}",
        CLOSED_LENGTH,
        closed(&b413, CLOSED_LENGTH),
        CLOSED_LENGTH,
        theta().powi(CLOSED_LENGTH as i32)
    );
    let (_, a14) = build(1, 4, Digits::Gasket);
    let poly = radii::charpoly_exact(&a14, &ps).unwrap();
    let g = walks(&a14, BOUND_LENGTH);
    let start: Vec<u128> = g.iter().map(|row| row[0]).collect();
    let recurrence = (2..=BOUND_LENGTH).all(|m| start[m] == start[m - 1] + 2 * start[m - 2]);
    let closed_form = (0..=BOUND_LENGTH).all(|w| {
        let sign: i128 = if w % 2 == 0 { -1 } else { 1 };
        start[w] as i128 * 3 == (1i128 << (w + 2)) + sign
    });
    println!(
        "A(1,4) live states {} charpoly low to high {:?}, G_w from the start {:?}, G_m = G_(m-1) + 2 G_(m-2) to m = {}: {}, G_w = (2^(w+2) + (-1)^(w+1))/3: {}",
        a14.edges.len(),
        poly,
        &start[..8],
        BOUND_LENGTH,
        recurrence,
        closed_form
    );
    for (s, t) in [(2u64, 9u64), (1, 18)] {
        for digits in [Digits::Gasket, Digits::Free] {
            let (_, live) = build(s, t, digits);
            let (rho, _) = radius(&live.edges);
            println!(
                "{}({},{}) live states {} radius {:.12}, the same automaton as ({},{}) with the carries swapped",
                digits.name(),
                s,
                t,
                live.edges.len(),
                rho,
                t,
                s
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cert(s: u64, t: u64, digits: Digits) -> Cert {
        let (_, live) = build(s, t, digits);
        certify(&live, &radii::primes(primes_for(live.edges.len()).max(2)))
    }

    #[test]
    fn shift_pair_has_radius_three() {
        let c = cert(1, 3, Digits::Free);
        assert!(c.at3 && !c.le2);
    }

    #[test]
    fn one_four_has_radius_exactly_two() {
        let c = cert(1, 4, Digits::Gasket);
        assert!(c.eq2 && !c.at3 && !c.theta);
    }

    #[test]
    fn four_thirteen_sits_at_theta_for_gasket_digits_only() {
        assert!(cert(4, 13, Digits::Gasket).theta);
        assert!(!cert(4, 13, Digits::Free).theta);
    }

    #[test]
    fn one_four_walks_follow_the_recurrence() {
        let (_, live) = build(1, 4, Digits::Gasket);
        let g: Vec<u128> = walks(&live, 10).iter().map(|row| row[0]).collect();
        assert_eq!(g[..6], [1, 3, 5, 11, 21, 43]);
    }
}
