use num_bigint::{BigInt, Sign};

use crate::{bigpow, gcd, ratio_f64};

// DIGIT SUMS

fn sum_dist(digits: &[u64], r: usize) -> Vec<u128> {
    let top = *digits.iter().max().unwrap() as usize;
    let mut out = vec![1u128];
    for _ in 0..r {
        let mut cur = vec![0u128; out.len() + top];
        for (i, &c) in out.iter().enumerate() {
            if c == 0 {
                continue;
            }
            for &f in digits {
                cur[i + f as usize] += c;
            }
        }
        out = cur;
    }
    out
}

fn digit_energy(digits: &[u64]) -> u128 {
    sum_dist(digits, 2).iter().map(|c| c * c).sum()
}

// TRANSFER

pub struct Transfer {
    states: Vec<(usize, usize)>,
    m: Vec<Vec<u128>>,
}

fn ordered(a: &[u128], q: usize, c: usize, cp: usize, d: usize, dp: usize) -> u128 {
    let mut s = 0u128;
    for rho in 0..q {
        let i = d * q + rho;
        let j = dp * q + rho;
        if i < c || j < cp {
            continue;
        }
        let (i, j) = (i - c, j - cp);
        if i >= a.len() || j >= a.len() {
            continue;
        }
        s += a[i] * a[j];
    }
    s
}

pub fn transfer(q: u64, digits: &[u64], p: usize) -> Transfer {
    assert!(p >= 2 && p.is_multiple_of(2));
    let r = p / 2;
    let a = sum_dist(digits, r);
    let mut states = Vec::new();
    for c in 0..r {
        for cp in c..r {
            states.push((c, cp));
        }
    }
    let n = states.len();
    let mut m = vec![vec![0u128; n]; n];
    for (i, &(c, cp)) in states.iter().enumerate() {
        for (j, &(d, dp)) in states.iter().enumerate() {
            let mut v = ordered(&a, q as usize, c, cp, d, dp);
            if c != cp {
                v += ordered(&a, q as usize, cp, c, d, dp);
            }
            m[i][j] = v;
        }
    }
    Transfer { states, m }
}

impl Transfer {
    fn n(&self) -> usize {
        self.states.len()
    }

    fn step(&self, w: &[BigInt]) -> Vec<BigInt> {
        let n = self.n();
        let mut out = vec![BigInt::from(0); n];
        for (i, wi) in w.iter().enumerate().take(n) {
            if *wi == BigInt::from(0) {
                continue;
            }
            for (j, slot) in out.iter_mut().enumerate() {
                if self.m[i][j] != 0 {
                    *slot += wi * BigInt::from(self.m[i][j]);
                }
            }
        }
        out
    }

    pub fn energies(&self, lmax: usize) -> (Vec<BigInt>, Vec<BigInt>) {
        let n = self.n();
        let mut w = vec![BigInt::from(0); n];
        w[0] = BigInt::from(1);
        let mut emod = Vec::new();
        let mut eint = Vec::new();
        for _ in 0..=lmax {
            let mut sm = BigInt::from(0);
            let mut si = BigInt::from(0);
            for (j, &(d, dp)) in self.states.iter().enumerate() {
                if d == dp {
                    si += &w[j];
                    sm += &w[j];
                } else {
                    sm += &w[j] * 2;
                }
            }
            emod.push(sm);
            eint.push(si);
            w = self.step(&w);
        }
        (emod, eint)
    }

    fn max_row_sum(&self) -> u128 {
        self.m.iter().map(|r| r.iter().sum::<u128>()).max().unwrap()
    }

    fn charpoly(&self) -> Vec<BigInt> {
        let n = self.n();
        let a: Vec<Vec<BigInt>> = self
            .m
            .iter()
            .map(|r| r.iter().map(|&v| BigInt::from(v)).collect())
            .collect();
        let ident = |n: usize| -> Vec<Vec<BigInt>> {
            (0..n)
                .map(|i| (0..n).map(|j| BigInt::from((i == j) as u8)).collect())
                .collect()
        };
        let mul = |x: &Vec<Vec<BigInt>>, y: &Vec<Vec<BigInt>>| -> Vec<Vec<BigInt>> {
            let mut z = vec![vec![BigInt::from(0); n]; n];
            for i in 0..n {
                for l in 0..n {
                    if x[i][l] == BigInt::from(0) {
                        continue;
                    }
                    for j in 0..n {
                        z[i][j] += &x[i][l] * &y[l][j];
                    }
                }
            }
            z
        };
        let mut c = vec![BigInt::from(0); n + 1];
        c[n] = BigInt::from(1);
        let mut mk = ident(n);
        for step in 1..=n {
            if step > 1 {
                let am = mul(&a, &mk);
                mk = am;
                for (i, row) in mk.iter_mut().enumerate().take(n) {
                    row[i] += &c[n - step + 1];
                }
            }
            let am = mul(&a, &mk);
            let mut tr = BigInt::from(0);
            for (i, row) in am.iter().enumerate().take(n) {
                tr += &row[i];
            }
            let kk = BigInt::from(step as u64);
            assert!(&tr % &kk == BigInt::from(0), "leverrier division not exact");
            c[n - step] = -(&tr / &kk);
        }
        c
    }

    fn reachable(&self) -> Vec<usize> {
        let n = self.n();
        let mut seen = vec![false; n];
        let mut stack = vec![0usize];
        seen[0] = true;
        while let Some(i) = stack.pop() {
            for (j, reached) in seen.iter_mut().enumerate() {
                if self.m[i][j] != 0 && !*reached {
                    *reached = true;
                    stack.push(j);
                }
            }
        }
        (0..n).filter(|&i| seen[i]).collect()
    }

    fn power(&self, idx: &[usize]) -> Vec<f64> {
        let n = idx.len();
        let mut v = vec![1.0f64; n];
        for _ in 0..4000 {
            let mut nv = vec![0.0f64; n];
            for a in 0..n {
                let mut s = v[a];
                for b in 0..n {
                    s += self.m[idx[a]][idx[b]] as f64 * v[b];
                }
                nv[a] = s;
            }
            let mx = nv.iter().cloned().fold(0.0f64, f64::max);
            for x in nv.iter_mut() {
                *x /= mx;
            }
            v = nv;
        }
        v
    }

    pub fn bracket(&self) -> (Frac, Frac) {
        let reach = self.reachable();
        let v0 = self.power(&reach);
        let idx: Vec<usize> = reach
            .iter()
            .zip(v0.iter())
            .filter(|(_, &x)| x > 1e-9)
            .map(|(&i, _)| i)
            .collect();
        let v = self.power(&idx);
        let n = idx.len();
        let vv: Vec<BigInt> = v
            .iter()
            .map(|&x| {
                let s = (x * 2f64.powi(60)).round() as u128;
                BigInt::from(s.max(1))
            })
            .collect();
        let mut lo: Option<Frac> = None;
        let mut hi: Option<Frac> = None;
        for a in 0..n {
            let mut mv = BigInt::from(0);
            for b in 0..n {
                mv += BigInt::from(self.m[idx[a]][idx[b]]) * &vv[b];
            }
            let f = Frac {
                num: mv,
                den: vv[a].clone(),
            };
            lo = Some(match lo {
                None => f.clone(),
                Some(l) => {
                    if f.lt(&l) {
                        f.clone()
                    } else {
                        l
                    }
                }
            });
            hi = Some(match hi {
                None => f,
                Some(h) => {
                    if h.lt(&f) {
                        f
                    } else {
                        h
                    }
                }
            });
        }
        (lo.unwrap(), hi.unwrap())
    }
}

// FRACTIONS

#[derive(Clone)]
pub struct Frac {
    pub num: BigInt,
    pub den: BigInt,
}

impl Frac {
    fn int(v: i64) -> Frac {
        Frac {
            num: BigInt::from(v),
            den: BigInt::from(1),
        }
    }

    fn lt(&self, o: &Frac) -> bool {
        &self.num * &o.den < &o.num * &self.den
    }

    fn mul(&self, o: &Frac) -> Frac {
        Frac {
            num: &self.num * &o.num,
            den: &self.den * &o.den,
        }
    }

    fn add(&self, o: &Frac) -> Frac {
        Frac {
            num: &self.num * &o.den + &o.num * &self.den,
            den: &self.den * &o.den,
        }
    }

    fn sub(&self, o: &Frac) -> Frac {
        Frac {
            num: &self.num * &o.den - &o.num * &self.den,
            den: &self.den * &o.den,
        }
    }

    fn abs(&self) -> Frac {
        Frac {
            num: if self.num.sign() == Sign::Minus {
                -&self.num
            } else {
                self.num.clone()
            },
            den: self.den.clone(),
        }
    }

    fn sign(&self) -> Sign {
        self.num.sign()
    }

    fn to_f64(&self) -> f64 {
        ratio_f64(&self.num, &self.den)
    }

    fn floor_digits(&self, digits: usize) -> String {
        let scaled = &self.num * bigpow(&BigInt::from(10), digits) / &self.den;
        place(&scaled, digits)
    }

    fn ceil_digits(&self, digits: usize) -> String {
        let ten = bigpow(&BigInt::from(10), digits);
        let scaled = (&self.num * &ten + &self.den - BigInt::from(1)) / &self.den;
        place(&scaled, digits)
    }
}

fn eval_frac(c: &[BigInt], x: &Frac) -> Frac {
    let mut acc = Frac::int(0);
    for coef in c.iter().rev() {
        acc = acc.mul(x).add(&Frac {
            num: coef.clone(),
            den: BigInt::from(1),
        });
    }
    acc
}

fn place(scaled: &BigInt, digits: usize) -> String {
    let s = scaled.to_string();
    if digits == 0 {
        return s;
    }
    let s = if s.len() <= digits {
        format!("{}{}", "0".repeat(digits + 1 - s.len()), s)
    } else {
        s
    };
    let (a, b) = s.split_at(s.len() - digits);
    format!("{a}.{b}")
}

fn frac_band(lo: &Frac, hi: &Frac, digits: usize) -> String {
    let l = lo.floor_digits(digits);
    let h = hi.ceil_digits(digits);
    if l == h {
        l
    } else {
        format!("{l}..{h}")
    }
}

pub(crate) fn band(lo: f64, hi: f64, digits: usize) -> String {
    let scale = 10f64.powi(digits as i32);
    let l = (lo * scale).floor() / scale;
    let h = (hi * scale).ceil() / scale;
    if (l - h).abs() < 0.5 / scale {
        format!("{l:.digits$}")
    } else {
        format!("{l:.digits$}..{h:.digits$}")
    }
}

// FACTORS

fn eval_big(c: &[BigInt], x: &BigInt) -> BigInt {
    let mut acc = BigInt::from(0);
    for coef in c.iter().rev() {
        acc = acc * x + coef;
    }
    acc
}

fn eval_i128(c: &[i128], x: i128) -> Option<i128> {
    let mut acc: i128 = 0;
    for &coef in c.iter().rev() {
        acc = acc.checked_mul(x)?.checked_add(coef)?;
    }
    Some(acc)
}

fn to_i128(c: &[BigInt]) -> Option<Vec<i128>> {
    let mut out = Vec::new();
    for x in c {
        let s = x.to_string();
        out.push(s.parse::<i128>().ok()?);
    }
    Some(out)
}

fn divide_root(c: &[BigInt], r: &BigInt) -> Vec<BigInt> {
    let n = c.len() - 1;
    let mut qv = vec![BigInt::from(0); n];
    let mut carry = BigInt::from(0);
    for i in (0..n).rev() {
        carry = &carry * r + &c[i + 1];
        qv[i] = carry.clone();
    }
    assert!(&carry * r + &c[0] == BigInt::from(0));
    qv
}

fn quad_rem_i128(c: &[i128], b: i128, cc: i128) -> Option<(i128, i128)> {
    let n = c.len() - 1;
    let mut q = vec![0i128; n - 1];
    for i in (0..=n - 2).rev() {
        let q1 = if i + 1 < q.len() { q[i + 1] } else { 0 };
        let q2 = if i + 2 < q.len() { q[i + 2] } else { 0 };
        q[i] = c[i + 2]
            .checked_sub(b.checked_mul(q1)?)?
            .checked_sub(cc.checked_mul(q2)?)?;
    }
    let q1 = if 1 < q.len() { q[1] } else { 0 };
    let r1 = c[1]
        .checked_sub(b.checked_mul(q[0])?)?
        .checked_sub(cc.checked_mul(q1)?)?;
    let r0 = c[0].checked_sub(cc.checked_mul(q[0])?)?;
    Some((r1, r0))
}

fn divide_quad(c: &[BigInt], b: &BigInt, cc: &BigInt) -> Vec<BigInt> {
    let n = c.len() - 1;
    let mut q = vec![BigInt::from(0); n - 1];
    for i in (0..=n - 2).rev() {
        let q1 = if i + 1 < q.len() {
            q[i + 1].clone()
        } else {
            BigInt::from(0)
        };
        let q2 = if i + 2 < q.len() {
            q[i + 2].clone()
        } else {
            BigInt::from(0)
        };
        q[i] = &c[i + 2] - b * q1 - cc * q2;
    }
    let q1 = if 1 < q.len() {
        q[1].clone()
    } else {
        BigInt::from(0)
    };
    assert!(&c[1] - b * &q[0] - cc * q1 == BigInt::from(0));
    assert!(&c[0] - cc * &q[0] == BigInt::from(0));
    q
}

fn divisors(n: u128) -> Vec<u128> {
    let mut out = vec![1u128];
    let mut m = n;
    let mut f = 2u128;
    while f * f <= m {
        if m.is_multiple_of(f) {
            let mut e = 0;
            while m.is_multiple_of(f) {
                m /= f;
                e += 1;
            }
            let base = out.clone();
            let mut pw = 1u128;
            for _ in 0..e {
                pw *= f;
                out.extend(base.iter().map(|d| d * pw));
            }
        }
        f += 1;
    }
    if m > 1 {
        let base = out.clone();
        out.extend(base.iter().map(|d| d * m));
    }
    out.sort();
    out
}

pub struct Factors {
    pub linear: Vec<BigInt>,
    pub quads: Vec<Vec<BigInt>>,
    pub rest: Vec<BigInt>,
    pub scanned: bool,
    pub witness: Option<u64>,
}

pub fn factor_bounded(c: &[BigInt], bound: u128) -> Factors {
    let mut cur: Vec<BigInt> = c.to_vec();
    let mut linear = Vec::new();
    let mut quads = Vec::new();
    let small = to_i128(c);
    if small.is_none() && bound > 100000 {
        return Factors {
            linear,
            quads,
            rest: cur,
            scanned: false,
            witness: None,
        };
    }
    let bound_i = bound as i128;
    let mut r: i128 = -bound_i;
    while r <= bound_i {
        let rb = BigInt::from(r);
        let is_root = match small.as_ref().and_then(|s| eval_i128(s, r)) {
            Some(v) => v == 0,
            None => eval_big(c, &rb) == BigInt::from(0),
        };
        if is_root {
            while cur.len() > 1 && eval_big(&cur, &rb) == BigInt::from(0) {
                cur = divide_root(&cur, &rb);
                linear.push(rb.clone());
            }
        }
        r += 1;
    }
    let mut again = cur.len() >= 4;
    while again {
        again = false;
        let ci = match to_i128(&cur) {
            Some(v) => v,
            None => break,
        };
        let a0 = ci[0].unsigned_abs();
        if a0 == 0 || a0 > 1_000_000_000_000_000 {
            break;
        }
        let divs = divisors(a0);
        let bb = bound_i.checked_mul(bound_i).unwrap_or(i128::MAX);
        'search: for &d in &divs {
            if d as i128 > bb {
                break;
            }
            for cc in [d as i128, -(d as i128)] {
                let mut b = -2 * bound_i;
                while b <= 2 * bound_i {
                    let hit = match quad_rem_i128(&ci, b, cc) {
                        Some((r1, r0)) => r1 == 0 && r0 == 0,
                        None => false,
                    };
                    if hit {
                        let bq = BigInt::from(b);
                        let cq = BigInt::from(cc);
                        cur = divide_quad(&cur, &bq, &cq);
                        quads.push(vec![cq, bq, BigInt::from(1)]);
                        again = cur.len() >= 4;
                        break 'search;
                    }
                    b += 1;
                }
            }
        }
    }
    let witness = if cur.len() >= 5 {
        irreducibility_witness(&cur)
    } else {
        None
    };
    Factors {
        linear,
        quads,
        rest: cur,
        scanned: true,
        witness,
    }
}

impl Factors {
    fn all(&self) -> Vec<(Vec<BigInt>, bool)> {
        let mut out: Vec<(Vec<BigInt>, bool)> = Vec::new();
        for r in &self.linear {
            out.push((vec![-r, BigInt::from(1)], true));
        }
        for q in &self.quads {
            out.push((q.clone(), true));
        }
        if self.rest.len() >= 2 {
            let certified = self.scanned && (self.rest.len() <= 4 || self.witness.is_some());
            out.push((self.rest.clone(), certified));
        }
        out
    }
}

fn nonzero_on(g: &[BigInt], lo: &Frac, hi: &Frac) -> bool {
    let at_lo = eval_frac(g, lo).abs();
    let mut slope = Frac::int(0);
    let mut pw = Frac::int(1);
    for (i, coef) in g.iter().enumerate().skip(1) {
        let term = Frac {
            num: BigInt::from(i as u64) * coef,
            den: BigInt::from(1),
        }
        .abs()
        .mul(&pw);
        slope = slope.add(&term);
        pw = pw.mul(hi);
    }
    let width = hi.sub(lo);
    width.mul(&slope).lt(&at_lo)
}

pub fn locate(factors: &[(Vec<BigInt>, bool)], lo: &Frac, hi: &Frac) -> Option<usize> {
    let mut found = None;
    for (i, (g, _)) in factors.iter().enumerate() {
        let sl = eval_frac(g, lo).sign();
        let sh = eval_frac(g, hi).sign();
        let changes = sl == Sign::NoSign || sh == Sign::NoSign || sl != sh;
        if changes {
            if found.is_some() {
                return None;
            }
            found = Some(i);
        } else if !nonzero_on(g, lo, hi) {
            return None;
        }
    }
    found
}

// IRREDUCIBILITY MOD P

fn pm_trim(a: &mut Vec<u64>) {
    while a.len() > 1 && *a.last().unwrap() == 0 {
        a.pop();
    }
}

fn pm_rem(a: &[u64], m: &[u64], p: u64) -> Vec<u64> {
    let mut a = a.to_vec();
    pm_trim(&mut a);
    let dm = m.len() - 1;
    if dm == 0 {
        return vec![0];
    }
    let inv = modinv(m[dm], p);
    while a.len() > dm {
        let lead = a[a.len() - 1];
        if lead != 0 {
            let f = lead * inv % p;
            let shift = a.len() - 1 - dm;
            for i in 0..=dm {
                a[shift + i] = (a[shift + i] + p - f * m[i] % p) % p;
            }
        }
        a.pop();
        if a.is_empty() {
            a.push(0);
        }
    }
    pm_trim(&mut a);
    a
}

fn modinv(a: u64, p: u64) -> u64 {
    let mut r = 1u64;
    let mut b = a % p;
    let mut e = p - 2;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % p;
        }
        b = b * b % p;
        e >>= 1;
    }
    r
}

fn pm_mulmod(a: &[u64], b: &[u64], m: &[u64], p: u64) -> Vec<u64> {
    let mut c = vec![0u64; a.len() + b.len() - 1];
    for (i, &x) in a.iter().enumerate() {
        if x == 0 {
            continue;
        }
        for (j, &y) in b.iter().enumerate() {
            c[i + j] = (c[i + j] + x * y) % p;
        }
    }
    pm_rem(&c, m, p)
}

fn pm_powmod(base: &[u64], mut e: u64, m: &[u64], p: u64) -> Vec<u64> {
    let mut r = vec![1u64];
    let mut b = pm_rem(base, m, p);
    while e > 0 {
        if e & 1 == 1 {
            r = pm_mulmod(&r, &b, m, p);
        }
        b = pm_mulmod(&b, &b, m, p);
        e >>= 1;
    }
    r
}

fn pm_gcd(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
    let mut a = a.to_vec();
    let mut b = b.to_vec();
    pm_trim(&mut a);
    pm_trim(&mut b);
    while !(b.len() == 1 && b[0] == 0) {
        let r = pm_rem(&a, &b, p);
        a = b;
        b = r;
    }
    a
}

fn pm_sub(a: &[u64], b: &[u64], p: u64) -> Vec<u64> {
    let n = a.len().max(b.len());
    let mut c = vec![0u64; n];
    for i in 0..n {
        let x = if i < a.len() { a[i] } else { 0 };
        let y = if i < b.len() { b[i] } else { 0 };
        c[i] = (x + p - y) % p;
    }
    pm_trim(&mut c);
    c
}

fn reduce_mod(c: &[BigInt], p: u64) -> Vec<u64> {
    let pb = BigInt::from(p);
    c.iter()
        .map(|x| {
            let m = ((x % &pb) + &pb) % &pb;
            m.to_string().parse::<u64>().unwrap()
        })
        .collect()
}

fn prime_factors(mut d: u64) -> Vec<u64> {
    let mut out = Vec::new();
    let mut f = 2;
    while f * f <= d {
        if d.is_multiple_of(f) {
            out.push(f);
            while d.is_multiple_of(f) {
                d /= f;
            }
        }
        f += 1;
    }
    if d > 1 {
        out.push(d);
    }
    out
}

pub fn irreducible_mod(c: &[BigInt], p: u64) -> bool {
    let m = reduce_mod(c, p);
    let d = m.len() - 1;
    if d < 1 || m[d] == 0 {
        return false;
    }
    if d == 1 {
        return true;
    }
    let deriv: Vec<u64> = (1..=d).map(|i| m[i] * (i as u64 % p) % p).collect();
    let g = pm_gcd(&m, &deriv, p);
    if g.len() > 1 {
        return false;
    }
    let x = vec![0u64, 1];
    let mut powers = vec![x.clone()];
    for _ in 0..d {
        let last = powers.last().unwrap().clone();
        powers.push(pm_powmod(&last, p, &m, p));
    }
    let top = pm_sub(&powers[d], &x, p);
    if !(top.len() == 1 && top[0] == 0) {
        return false;
    }
    for l in prime_factors(d as u64) {
        let i = d / l as usize;
        let diff = pm_sub(&powers[i], &x, p);
        let g = pm_gcd(&m, &diff, p);
        if g.len() > 1 {
            return false;
        }
    }
    true
}

const PRIMES: [u64; 30] = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
    101, 103, 107, 109, 113,
];

pub fn irreducibility_witness(c: &[BigInt]) -> Option<u64> {
    if c.len() <= 2 {
        return Some(0);
    }
    PRIMES.iter().cloned().find(|&p| irreducible_mod(c, p))
}

// RENDER

fn poly_string(c: &[BigInt]) -> String {
    let mut parts = Vec::new();
    for i in (0..c.len()).rev() {
        if c[i] == BigInt::from(0) {
            continue;
        }
        let sign = if c[i].sign() == Sign::Minus { "-" } else { "+" };
        let mag = if c[i].sign() == Sign::Minus {
            -&c[i]
        } else {
            c[i].clone()
        };
        let term = match i {
            0 => format!("{mag}"),
            1 => {
                if mag == BigInt::from(1) {
                    "x".to_string()
                } else {
                    format!("{mag} x")
                }
            }
            _ => {
                if mag == BigInt::from(1) {
                    format!("x^{i}")
                } else {
                    format!("{mag} x^{i}")
                }
            }
        };
        if parts.is_empty() {
            parts.push(if sign == "-" {
                format!("-{term}")
            } else {
                term
            });
        } else {
            parts.push(format!("{sign} {term}"));
        }
    }
    parts.join(" ")
}

// STUDY

pub struct Family {
    pub q: u64,
    pub digits: Vec<u64>,
    pub label: &'static str,
}

pub struct Moment {
    pub p: usize,
    pub states: usize,
    pub charpoly: Vec<BigInt>,
    pub factors: Factors,
    pub lo: Frac,
    pub hi: Frac,
    pub located: Option<usize>,
    pub emod: Vec<BigInt>,
    pub eint: Vec<BigInt>,
}

pub fn moment(fam: &Family, p: usize, lmax: usize) -> Moment {
    let t = transfer(fam.q, &fam.digits, p);
    let (emod, eint) = t.energies(lmax);
    let charpoly = t.charpoly();
    let n = t.n();
    for l in 0..=(lmax - n) {
        let mut sm = BigInt::from(0);
        let mut si = BigInt::from(0);
        for i in 0..=n {
            sm += &charpoly[i] * &emod[l + i];
            si += &charpoly[i] * &eint[l + i];
        }
        assert!(sm == BigInt::from(0), "recurrence fails for E_mod");
        assert!(si == BigInt::from(0), "recurrence fails for E_int");
    }
    let factors = factor_bounded(&charpoly, t.max_row_sum());
    let (lo, hi) = t.bracket();
    let located = if factors.scanned {
        locate(&factors.all(), &lo, &hi)
    } else {
        None
    };
    Moment {
        p,
        states: n,
        charpoly,
        factors,
        lo,
        hi,
        located,
        emod,
        eint,
    }
}

fn assert_bounds(fam: &Family, m: &Moment) {
    let q = BigInt::from(fam.q);
    let k = BigInt::from(fam.digits.len() as u64);
    let r = m.p / 2;
    for l in 0..m.emod.len() {
        let ql = bigpow(&q, l);
        let kl = bigpow(&k, l);
        let s = &ql * &m.emod[l];
        let hoelder = &ql * bigpow(&kl, r);
        let zero = bigpow(&kl, m.p);
        let upper = &ql * bigpow(&kl, m.p - 1);
        assert!(s >= hoelder, "hoelder floor fails");
        assert!(s >= zero, "zero frequency floor fails");
        assert!(s <= upper, "trivial ceiling fails");
        assert!(m.eint[l] <= m.emod[l], "integer energy exceeds modular");
    }
}

pub(crate) fn alpha(fam: &Family) -> f64 {
    (fam.digits.len() as f64).ln() / (fam.q as f64).ln()
}

pub(crate) fn theta_band(fam: &Family, m: &Moment) -> (f64, f64) {
    let lq = (fam.q as f64).ln();
    let lo = 1.0 + m.lo.to_f64().ln() / lq - 1e-12;
    let hi = 1.0 + m.hi.to_f64().ln() / lq + 1e-12;
    (lo, hi)
}

fn minpoly_string(m: &Moment) -> String {
    if !m.factors.scanned {
        return "divides charpoly, scan skipped".to_string();
    }
    let all = m.factors.all();
    match m.located {
        Some(i) => {
            let (g, certified) = &all[i];
            if *certified {
                match (g.len(), m.factors.witness) {
                    (l, Some(p)) if l >= 5 => format!("{} [irreducible mod {p}]", poly_string(g)),
                    _ => poly_string(g),
                }
            } else {
                format!("divides {} [not certified]", poly_string(g))
            }
        }
        None => "not located".to_string(),
    }
}

fn factors_string(m: &Moment) -> String {
    let mut parts = Vec::new();
    for r in &m.factors.linear {
        parts.push(format!("(x - {r})").replace("- -", "+ "));
    }
    for q in &m.factors.quads {
        parts.push(format!("({})", poly_string(q)));
    }
    if m.factors.rest.len() >= 2 {
        parts.push(format!("({})", poly_string(&m.factors.rest)));
    }
    parts.join(" ")
}

pub fn algebra_row(fam: &Family, m: &Moment) -> String {
    let carries = 2 * *fam.digits.iter().max().unwrap() >= fam.q;
    format!(
        "| {} | {} | {} | {} | {} | {} | {} | {} | {} |",
        fam.q,
        fam.label,
        fam.digits.len(),
        m.p,
        if carries { "yes" } else { "no" },
        m.states,
        poly_string(&m.charpoly),
        factors_string(m),
        minpoly_string(m)
    )
}

pub fn growth_row(fam: &Family, m: &Moment) -> String {
    let q = fam.q as f64;
    let k = fam.digits.len() as f64;
    let (tlo, thi) = theta_band(fam, m);
    let al = alpha(fam);
    let lam_lo = Frac {
        num: &m.lo.num * BigInt::from(fam.q),
        den: m.lo.den.clone(),
    };
    let lam_hi = Frac {
        num: &m.hi.num * BigInt::from(fam.q),
        den: m.hi.den.clone(),
    };
    let r = (m.p / 2) as i32;
    let floor = (k.powi(m.p as i32)).max(q * k.powi(r));
    let ceil = q * k.powi(m.p as i32 - 1);
    let root_lo = (lam_lo.to_f64()).powf(1.0 / m.p as f64) / k;
    let root_hi = (lam_hi.to_f64()).powf(1.0 / m.p as f64) / k;
    format!(
        "| {} | {} | {} | {} | {} | {} | {} | {} | {} |",
        fam.q,
        fam.label,
        m.p,
        frac_band(&lam_lo, &lam_hi, 9),
        floor,
        ceil,
        band(tlo, thi, 9),
        band(al - 1e-12, al + 1e-12, 6),
        band(root_lo - 1e-12, root_hi + 1e-12, 6)
    )
}

pub fn typeii_row(fam: &Family, m: &Moment) -> String {
    let (tlo, thi) = theta_band(fam, m);
    let al = alpha(fam);
    let nu_lo = tlo - 4.0 * al - 1e-12;
    let nu_hi = thi - 4.0 * al + 1e-12;
    let t2_lo = tlo / 4.0 + 0.25;
    let t2_hi = thi / 4.0 + 0.25;
    let miss_lo = t2_lo - al - 1e-12;
    let miss_hi = t2_hi - al + 1e-12;
    let need_lo = 1.0 - al / 2.0 - 1e-12;
    let need_hi = 1.0 - al / 2.0 + 1e-12;
    let eta_lo = (1.0 + 3.0 * al - thi) / 2.0 - 1e-12;
    let eta_hi = (1.0 + 3.0 * al - tlo) / 2.0 + 1e-12;
    format!(
        "| {} | {} | {} | {} | {} | {} | {} | {} | {} |",
        fam.q,
        fam.label,
        band(al - 1e-12, al + 1e-12, 6),
        band(tlo, thi, 6),
        band(nu_lo, nu_hi, 6),
        band(t2_lo, t2_hi, 6),
        band(miss_lo, miss_hi, 6),
        band(need_lo, need_hi, 6),
        band(eta_lo, eta_hi, 6)
    )
}

pub fn values_row(fam: &Family, m: &Moment, l: usize) -> String {
    let s = bigpow(&BigInt::from(fam.q), l) * &m.emod[l];
    format!(
        "| {} | {} | {} | {} | {} | {} |",
        fam.q, fam.label, l, m.emod[l], m.eint[l], s
    )
}

pub fn ratio_row(fam: &Family, m: &Moment) -> String {
    let lmax = m.emod.len() - 1;
    let q = BigInt::from(fam.q);
    let k4 = bigpow(&BigInt::from(fam.digits.len() as u64), m.p);
    let s = |l: usize| bigpow(&q, l) * &m.emod[l];
    let prim = |l: usize| s(l) - &k4 * s(l - 1);
    let growth = ratio_f64(&m.emod[lmax], &m.emod[lmax - 1]);
    let prim_growth = {
        let a = prim(lmax);
        let b = prim(lmax - 1);
        if b == BigInt::from(0) {
            0.0
        } else {
            ratio_f64(&a, &b)
        }
    };
    let mi = ratio_f64(&m.emod[lmax], &m.eint[lmax]);
    let rho = m.lo.to_f64();
    let zero_share = {
        let z = bigpow(&BigInt::from(fam.digits.len() as u64), m.p * lmax);
        ratio_f64(&z, &s(lmax))
    };
    format!(
        "| {} | {} | {} | {} | {:.3e} | {:.6} | {:.6} | {:.3e} |",
        fam.q,
        fam.label,
        m.p,
        lmax,
        (growth - rho).abs(),
        mi,
        prim_growth / (fam.q as f64 * rho),
        zero_share
    )
}

// ARCS

pub(crate) fn transform_table(q: u64, digits: &[u64], l: usize) -> Vec<(f64, f64)> {
    let n = (q as usize).pow(l as u32);
    let mut g = vec![(0.0f64, 0.0f64); n];
    for (b, slot) in g.iter_mut().enumerate() {
        let y = b as f64 / n as f64;
        let mut re = 0.0;
        let mut im = 0.0;
        for &f in digits {
            let t = 2.0 * std::f64::consts::PI * f as f64 * y;
            re += t.cos();
            im += t.sin();
        }
        *slot = (re, im);
    }
    g
}

pub struct Arcs {
    pub l: usize,
    pub dmax: u64,
    pub points: usize,
    pub major: usize,
    pub parseval_rel: f64,
    pub fourth_rel: f64,
    pub major_l2: f64,
    pub major_l4: f64,
    pub zero_l4: f64,
}

pub fn arcs(fam: &Family, l: usize, s4_exact: &BigInt) -> Arcs {
    let q = fam.q;
    let n = (q as usize).pow(l as u32);
    let g = transform_table(q, &fam.digits, l);
    let mut sq = vec![0.0f64; n];
    for (a, slot) in sq.iter_mut().enumerate() {
        let mut b = a;
        let mut v = 1.0f64;
        for _ in 0..l {
            let (re, im) = g[b];
            v *= re * re + im * im;
            b = (b * q as usize) % n;
        }
        *slot = v;
    }
    let x = n as f64;
    let dmax = x.powf(0.4).floor() as u64;
    let qq = x.powf(0.6);
    let mut major = vec![false; n];
    for d in 1..=dmax {
        for lnum in 0..d {
            if gcd(lnum, d) != 1 {
                continue;
            }
            let c = lnum as f64 * x / d as f64;
            let h = x / (d as f64 * qq);
            let lo = (c - h).ceil() as i64;
            let hi = (c + h).floor() as i64;
            for a in lo..=hi {
                let idx = a.rem_euclid(n as i64) as usize;
                major[idx] = true;
            }
        }
    }
    let total2: f64 = sq.iter().sum();
    let total4: f64 = sq.iter().map(|v| v * v).sum();
    let m2: f64 = (0..n).filter(|&a| major[a]).map(|a| sq[a]).sum();
    let m4: f64 = (0..n).filter(|&a| major[a]).map(|a| sq[a] * sq[a]).sum();
    let k = fam.digits.len() as f64;
    let parseval = x * k.powi(l as i32);
    let s4 = ratio_f64(s4_exact, &BigInt::from(1));
    Arcs {
        l,
        dmax,
        points: n,
        major: major.iter().filter(|&&b| b).count(),
        parseval_rel: (total2 - parseval).abs() / parseval,
        fourth_rel: (total4 - s4).abs() / s4,
        major_l2: m2 / total2,
        major_l4: m4 / total4,
        zero_l4: sq[0] * sq[0] / total4,
    }
}

pub fn arcs_row(fam: &Family, a: &Arcs) -> String {
    let share = (fam.digits.len() as f64 / fam.q as f64).powi(a.l as i32);
    format!(
        "| {} | {} | {} | {} | {} | {} | {:.1e} | {:.1e} | {:.4} | {:.4} | {:.4} | {:.4} |",
        fam.q,
        fam.label,
        a.l,
        a.dmax,
        a.points,
        a.major,
        a.parseval_rel,
        a.fourth_rel,
        1.0 - a.major_l2,
        share,
        a.major_l4,
        a.zero_l4
    )
}

// FAMILIES

fn ex(q: u64, e: u64) -> Vec<u64> {
    (0..q).filter(|&f| f != e).collect()
}

pub fn families() -> Vec<Family> {
    vec![
        Family {
            q: 3,
            digits: vec![0, 1],
            label: "01",
        },
        Family {
            q: 3,
            digits: vec![0, 2],
            label: "02",
        },
        Family {
            q: 3,
            digits: vec![1, 2],
            label: "12",
        },
        Family {
            q: 4,
            digits: vec![0, 1, 2],
            label: "012",
        },
        Family {
            q: 5,
            digits: vec![0, 1, 2, 3],
            label: "0123",
        },
        Family {
            q: 5,
            digits: vec![0, 2, 4],
            label: "024",
        },
        Family {
            q: 10,
            digits: ex(10, 0),
            label: "ex0",
        },
        Family {
            q: 10,
            digits: ex(10, 1),
            label: "ex1",
        },
        Family {
            q: 10,
            digits: ex(10, 3),
            label: "ex3",
        },
        Family {
            q: 10,
            digits: ex(10, 5),
            label: "ex5",
        },
        Family {
            q: 10,
            digits: ex(10, 7),
            label: "ex7",
        },
        Family {
            q: 10,
            digits: ex(10, 9),
            label: "ex9",
        },
        Family {
            q: 100,
            digits: ex(100, 37),
            label: "ex37",
        },
        Family {
            q: 100,
            digits: (0..50).collect(),
            label: "0to49",
        },
        Family {
            q: 100,
            digits: vec![0, 1],
            label: "01",
        },
    ]
}

fn find<'a>(fams: &'a [Family], q: u64, label: &str) -> &'a Family {
    fams.iter().find(|f| f.q == q && f.label == label).unwrap()
}

pub fn run() {
    let lmax = 60;
    let fams = families();
    let fourth: Vec<Moment> = fams.iter().map(|f| moment(f, 4, lmax)).collect();
    for (f, m) in fams.iter().zip(fourth.iter()) {
        assert_bounds(f, m);
    }
    for (a, b) in [((3, "02"), (3, "01")), ((3, "12"), (3, "01"))] {
        let ia = fams
            .iter()
            .position(|f| f.q == a.0 && f.label == a.1)
            .unwrap();
        let ib = fams
            .iter()
            .position(|f| f.q == b.0 && f.label == b.1)
            .unwrap();
        assert!(fourth[ia].emod == fourth[ib].emod && fourth[ia].eint == fourth[ib].eint);
    }
    {
        let aux = Family {
            q: 5,
            digits: vec![0, 1, 2],
            label: "012",
        };
        let ma = moment(&aux, 4, lmax);
        let i = fams
            .iter()
            .position(|f| f.q == 5 && f.label == "024")
            .unwrap();
        assert!(fourth[i].emod == ma.emod && fourth[i].eint == ma.eint);
    }
    for (f, m) in fams.iter().zip(fourth.iter()) {
        if 2 * *f.digits.iter().max().unwrap() < f.q {
            let e = BigInt::from(digit_energy(&f.digits));
            for l in 0..=lmax {
                assert!(m.emod[l] == bigpow(&e, l) && m.eint[l] == bigpow(&e, l));
            }
        }
    }
    println!("riesz algebra");
    println!(
        "| q | F | k | p | carries | states | charpoly | factors | minimal polynomial of rho |"
    );
    for (f, m) in fams.iter().zip(fourth.iter()) {
        println!("{}", algebra_row(f, m));
    }
    println!("riesz growth");
    println!("| q | F | p | Lambda(p) | floor | ceiling | theta_p | alpha | Lambda^(1/p)/k |");
    for (f, m) in fams.iter().zip(fourth.iter()) {
        println!("{}", growth_row(f, m));
    }
    println!("riesz values");
    println!("| q | F | L | E_mod | E_int | S_4 |");
    for (f, m) in fams.iter().zip(fourth.iter()) {
        for l in [1usize, 2, 3] {
            println!("{}", values_row(f, m, l));
        }
    }
    println!("riesz ratios");
    println!("| q | F | p | L | growth-rho | E_mod/E_int | prim growth/Lambda | zero share |");
    for (f, m) in fams.iter().zip(fourth.iter()) {
        println!("{}", ratio_row(f, m));
    }
    println!("riesz type II");
    println!("| q | F | alpha | theta_4 | nu_4 | L4 exponent | miss | delta needed | eta_4 |");
    for (f, m) in fams.iter().zip(fourth.iter()) {
        println!("{}", typeii_row(f, m));
    }
    println!("riesz higher moments");
    println!(
        "| q | F | k | p | carries | states | charpoly | factors | minimal polynomial of rho |"
    );
    let higher: Vec<(&Family, Moment)> = [
        (3, "01", 6),
        (3, "01", 8),
        (3, "01", 10),
        (3, "12", 6),
        (4, "012", 6),
        (5, "0123", 6),
        (10, "ex7", 6),
        (100, "01", 6),
    ]
    .iter()
    .map(|&(q, label, p)| {
        let f = find(&fams, q, label);
        let m = moment(f, p, lmax);
        assert_bounds(f, &m);
        (f, m)
    })
    .collect();
    for (f, m) in higher.iter() {
        println!("{}", algebra_row(f, m));
    }
    println!("| q | F | p | Lambda(p) | floor | ceiling | theta_p | alpha | Lambda^(1/p)/k |");
    for (f, m) in higher.iter() {
        println!("{}", growth_row(f, m));
    }
    println!("| q | F | p | L | growth-rho | E_mod/E_int | prim growth/Lambda | zero share |");
    for (f, m) in higher.iter() {
        println!("{}", ratio_row(f, m));
    }
    println!("riesz arcs");
    println!("| q | F | L | D | points | major | parseval rel | fourth rel | minor l2 | (k/q)^L | major l4 | zero l4 |");
    for (q, label, l) in [
        (3, "01", 8),
        (3, "01", 12),
        (4, "012", 6),
        (4, "012", 9),
        (5, "0123", 6),
        (5, "0123", 8),
        (10, "ex7", 4),
        (10, "ex7", 6),
        (100, "01", 3),
        (100, "0to49", 2),
        (100, "0to49", 3),
    ] {
        let i = fams
            .iter()
            .position(|f| f.q == q && f.label == label)
            .unwrap();
        let s4 = bigpow(&BigInt::from(q), l) * &fourth[i].emod[l];
        let a = arcs(&fams[i], l, &s4);
        println!("{}", arcs_row(&fams[i], &a));
    }
}

// TESTS

#[cfg(test)]
mod tests {
    use super::*;

    fn roots_string(r: &[BigInt]) -> String {
        r.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    }

    fn strings(q: u64, digits: &[u64], l: usize) -> Vec<u128> {
        let mut out = vec![0u128];
        for _ in 0..l {
            let mut next = Vec::with_capacity(out.len() * digits.len());
            for &v in &out {
                for &f in digits {
                    next.push(v * q as u128 + f as u128);
                }
            }
            out = next;
        }
        out
    }

    fn histogram_energy(q: u64, digits: &[u64], l: usize, r: usize) -> (u128, u128) {
        let n = (q as u128).pow(l as u32);
        let vals = strings(q, digits, l);
        let mut hist = vec![1u128];
        for _ in 0..r {
            let mut next = vec![0u128; hist.len() + n as usize];
            for (s, &c) in hist.iter().enumerate() {
                if c == 0 {
                    continue;
                }
                for &v in &vals {
                    next[s + v as usize] += c;
                }
            }
            hist = next;
        }
        let eint: u128 = hist.iter().map(|c| c * c).sum();
        let mut modhist = vec![0u128; n as usize];
        for (s, &c) in hist.iter().enumerate() {
            modhist[s % n as usize] += c;
        }
        let emod: u128 = modhist.iter().map(|c| c * c).sum();
        (emod, eint)
    }

    fn direct_fourth(q: u64, digits: &[u64], l: usize, grid: usize) -> f64 {
        let n = (q as usize).pow(l as u32);
        let mut total = 0.0;
        for a in 0..grid {
            let mut v = 1.0f64;
            for j in 0..l {
                let y = a as f64 * (q as f64).powi(j as i32) / grid as f64;
                let mut re = 0.0;
                let mut im = 0.0;
                for &f in digits {
                    let t = 2.0 * std::f64::consts::PI * f as f64 * y;
                    re += t.cos();
                    im += t.sin();
                }
                v *= re * re + im * im;
            }
            total += v * v;
        }
        let _ = n;
        total
    }

    #[test]
    fn energy_matches_histogram() {
        for (q, digits, l) in [
            (3u64, vec![0u64, 1], 7usize),
            (3, vec![1, 2], 7),
            (3, vec![0, 2], 6),
            (4, vec![0, 1, 2], 5),
            (5, vec![0, 1, 2, 3], 4),
            (5, vec![0, 2, 4], 4),
            (10, ex(10, 7), 3),
            (10, ex(10, 0), 3),
            (100, ex(100, 37), 1),
            (100, (0..50).collect(), 1),
            (100, vec![0, 1], 2),
        ] {
            let t = transfer(q, &digits, 4);
            let (emod, eint) = t.energies(l);
            for ll in 0..=l {
                let (hm, hi) = histogram_energy(q, &digits, ll, 2);
                assert_eq!(emod[ll], BigInt::from(hm), "E_mod q={q} L={ll}");
                assert_eq!(eint[ll], BigInt::from(hi), "E_int q={q} L={ll}");
            }
        }
    }

    #[test]
    fn higher_energy_matches_histogram() {
        for (q, digits, p, l) in [
            (3u64, vec![0u64, 1], 6usize, 5usize),
            (3, vec![0, 1], 8, 4),
            (3, vec![1, 2], 6, 4),
            (4, vec![0, 1, 2], 6, 3),
            (5, vec![0, 1, 2, 3], 6, 3),
        ] {
            let t = transfer(q, &digits, p);
            let (emod, eint) = t.energies(l);
            for ll in 0..=l {
                let (hm, hi) = histogram_energy(q, &digits, ll, p / 2);
                assert_eq!(emod[ll], BigInt::from(hm));
                assert_eq!(eint[ll], BigInt::from(hi));
            }
        }
    }

    #[test]
    fn fourth_moment_matches_direct_sum() {
        for (q, digits, l) in [
            (3u64, vec![0u64, 1], 6usize),
            (3, vec![1, 2], 6),
            (4, vec![0, 1, 2], 5),
            (5, vec![0, 1, 2, 3], 4),
            (10, ex(10, 7), 4),
            (100, vec![0, 1], 2),
        ] {
            let t = transfer(q, &digits, 4);
            let (emod, eint) = t.energies(l);
            let n = (q as usize).pow(l as u32);
            let s4 = ratio_f64(&(bigpow(&BigInt::from(q), l) * &emod[l]), &BigInt::from(1));
            let grid = direct_fourth(q, &digits, l, n);
            assert!((grid - s4).abs() <= 1e-9 * s4, "S_4 q={q} L={l}");
            let e_int = ratio_f64(&eint[l], &BigInt::from(1));
            let integral = direct_fourth(q, &digits, l, 2 * n) / (2 * n) as f64;
            assert!(
                (integral - e_int).abs() <= 1e-9 * e_int,
                "E_int q={q} L={l}"
            );
        }
    }

    #[test]
    fn charpoly_and_scaling_pins() {
        let f = Family {
            q: 3,
            digits: vec![0, 1],
            label: "01",
        };
        let m = moment(&f, 4, 20);
        assert_eq!(poly_string(&m.charpoly), "x^3 - 8 x^2 + 13 x - 6");
        assert_eq!(roots_string(&m.factors.linear), "1,1,6");
        assert_eq!(minpoly_string(&m), "x - 6");
        assert_eq!(digit_energy(&[0, 1]), 6);
        assert_eq!(digit_energy(&(0..50).collect::<Vec<u64>>()), 83350);
        let g = Family {
            q: 3,
            digits: vec![1, 2],
            label: "12",
        };
        let mg = moment(&g, 4, 20);
        assert_eq!(mg.emod, m.emod);
        assert_eq!(minpoly_string(&mg), "x - 6");
        let h = Family {
            q: 3,
            digits: vec![0, 2],
            label: "02",
        };
        let mh = moment(&h, 4, 20);
        assert_eq!(poly_string(&mh.charpoly), "x^3 - 15 x^2 + 74 x - 120");
        assert_eq!(minpoly_string(&mh), "x - 6");
    }

    #[test]
    fn irreducibility_tester_pins() {
        let poly = |c: &[i64]| -> Vec<BigInt> { c.iter().map(|&x| BigInt::from(x)).collect() };
        assert!(irreducible_mod(&poly(&[1, 0, 1]), 3));
        assert!(!irreducible_mod(&poly(&[1, 0, 1]), 5));
        assert!(irreducible_mod(&poly(&[-2, 0, 1]), 3));
        assert!(!irreducible_mod(&poly(&[-2, 0, 1]), 7));
        assert!(irreducible_mod(&poly(&[-2, 0, 0, 1]), 7));
        assert!(!irreducible_mod(&poly(&[-1, 0, 0, 1]), 7));
        for p in [2u64, 3, 5, 7, 11, 13] {
            assert!(!irreducible_mod(&poly(&[1, 0, 0, 0, 1]), p));
        }
        assert!(irreducible_mod(&poly(&[1, 1, 1, 1, 1]), 3));
        assert_eq!(irreducibility_witness(&poly(&[1, 1, 1, 1, 1])), Some(2));
        let f = factor_bounded(&poly(&[-6, 11, -6, 1]), 10);
        assert!(f.scanned);
        assert_eq!(roots_string(&f.linear), "1,2,3");
        assert_eq!(f.rest.len(), 1);
        let g = factor_bounded(&poly(&[2, 0, 3, 0, 1]), 5);
        assert!(g.linear.is_empty());
        assert_eq!(g.quads.len(), 1);
        assert_eq!(poly_string(&g.rest), "x^2 + 2");
        assert_eq!(divisors(12), vec![1, 2, 3, 4, 6, 12]);
    }

    #[test]
    fn rows_pinned() {
        let fams = families();
        let row = |q: u64, label: &str, p: usize| {
            let f = find(&fams, q, label);
            let m = moment(f, p, 30);
            (
                algebra_row(f, &m),
                values_row(f, &m, 3),
                typeii_row(f, &m),
                m,
            )
        };
        let (a, v, t, m) = row(3, "01", 4);
        assert_eq!(a, "| 3 | 01 | 2 | 4 | no | 3 | x^3 - 8 x^2 + 13 x - 6 | (x - 1) (x - 1) (x - 6) | x - 6 |");
        assert_eq!(v, "| 3 | 01 | 3 | 216 | 216 | 5832 |");
        assert_eq!(t, "| 3 | 01 | 0.630929..0.630930 | 2.630929..2.630930 | 0.107210..0.107211 | 0.907732..0.907733 | 0.276802..0.276803 | 0.684535..0.684536 | 0.130929..0.130930 |");
        assert_eq!(frac_band(&m.lo, &m.hi, 9), "6.000000000");
        let (a, _, _, _) = row(3, "02", 4);
        assert_eq!(a, "| 3 | 02 | 2 | 4 | yes | 3 | x^3 - 15 x^2 + 74 x - 120 | (x - 4) (x - 5) (x - 6) | x - 6 |");
        let (a, _, _, _) = row(5, "024", 4);
        assert_eq!(a, "| 5 | 024 | 3 | 4 | yes | 3 | x^3 - 42 x^2 + 563 x - 2394 | (x - 9) (x - 14) (x - 19) | x - 19 |");
        let (a, v, _, m) = row(4, "012", 4);
        assert_eq!(a, "| 4 | 012 | 3 | 4 | yes | 3 | x^3 - 27 x^2 + 136 x - 176 | (x - 4) (x^2 - 23 x + 44) | x^2 - 23 x + 44 |");
        assert_eq!(v, "| 4 | 012 | 3 | 9173 | 8203 | 587072 |");
        let rho = (23.0 + 353f64.sqrt()) / 2.0;
        assert!(m.lo.to_f64() <= rho + 1e-9 && rho - 1e-9 <= m.hi.to_f64());
        let (a, v, _, m) = row(5, "0123", 4);
        assert_eq!(a, "| 5 | 0123 | 4 | 4 | yes | 3 | x^3 - 64 x^2 + 659 x - 1476 | (x - 9) (x^2 - 55 x + 164) | x^2 - 55 x + 164 |");
        assert_eq!(v, "| 5 | 0123 | 3 | 139752 | 116864 | 17469000 |");
        let rho = (55.0 + 2369f64.sqrt()) / 2.0;
        assert!(m.lo.to_f64() <= rho + 1e-9 && rho - 1e-9 <= m.hi.to_f64());
        let (a, v, _, m) = row(10, "ex7", 4);
        assert_eq!(a, "| 10 | ex7 | 9 | 4 | yes | 3 | x^3 - 731 x^2 + 49474 x - 434304 | (x - 64) (x^2 - 667 x + 6786) | x^2 - 667 x + 6786 |");
        assert_eq!(v, "| 10 | ex7 | 3 | 283307409 | 188677161 | 283307409000 |");
        let rho = (667.0 + 417745f64.sqrt()) / 2.0;
        assert!(m.lo.to_f64() <= rho + 1e-8 && rho - 1e-8 <= m.hi.to_f64());
        let (a, v, _, m) = row(100, "ex37", 4);
        assert_eq!(a, "| 100 | ex37 | 99 | 4 | yes | 3 | x^3 - 970301 x^2 + 9322916414 x - 925656819304 | (x - 9604) (x^2 - 960697 x + 96382426) | x^2 - 960697 x + 96382426 |");
        assert_eq!(v, "| 100 | ex37 | 3 | 886386992219168729 | 588561774108546671 | 886386992219168729000000 |");
        let rho = (960697.0 + (960697f64 * 960697.0 - 4.0 * 96382426.0).sqrt()) / 2.0;
        assert!(m.lo.to_f64() <= rho + 1e-4 && rho - 1e-4 <= m.hi.to_f64());
        let (a, _, _, _) = row(100, "0to49", 4);
        assert_eq!(a, "| 100 | 0to49 | 50 | 4 | no | 3 | x^3 - 83350 x^2 | (x - 0) (x - 0) (x - 83350) | x - 83350 |");
        for (q, label, p, want) in [
            (3, "01", 6, "x^2 - 26 x + 90"),
            (3, "01", 8, "x^2 - 99 x + 1134"),
            (3, "01", 10, "x^3 - 392 x^2 + 17469 x - 96228"),
            (4, "012", 6, "x^2 - 197 x + 2604"),
            (5, "0123", 6, "x^3 - 858 x^2 + 31475 x - 62700"),
            (10, "ex7", 6, "x^3 - 53697 x^2 + 29368547 x - 460410771"),
            (100, "01", 6, "x - 20"),
        ] {
            let (_, _, _, m) = row(q, label, p);
            assert_eq!(minpoly_string(&m), want, "{q} {label} {p}");
        }
        let (_, _, _, m) = row(3, "01", 6);
        let rho = 13.0 + 79f64.sqrt();
        assert!(m.lo.to_f64() <= rho + 1e-9 && rho - 1e-9 <= m.hi.to_f64());
    }

    #[test]
    fn arcs_pinned() {
        let fams = families();
        let f = find(&fams, 3, "01");
        let m = moment(f, 4, 12);
        let s4 = bigpow(&BigInt::from(3), 12) * &m.emod[12];
        let a = arcs(f, 12, &s4);
        assert_eq!((a.dmax, a.points, a.major), (195, 531441, 46525));
        assert!(a.parseval_rel < 1e-10 && a.fourth_rel < 1e-10);
        assert!((a.major_l2 - 0.2005).abs() < 5e-4);
        assert!((a.major_l4 - 0.6462).abs() < 5e-4);
        assert!((a.zero_l4 - (16f64 / 18.0).powi(12)).abs() < 1e-9);
    }

    #[test]
    fn bracket_contains_growth() {
        for f in families() {
            let m = moment(&f, 4, 40);
            let g = ratio_f64(&m.emod[40], &m.emod[39]);
            let rho = m.lo.to_f64();
            assert!(
                (g - rho).abs() <= 1e-3 * rho,
                "growth off at {} {}",
                f.q,
                f.label
            );
            let width = m.hi.to_f64() - rho;
            assert!(
                width <= 1e-9 * rho,
                "bracket too wide at {} {}",
                f.q,
                f.label
            );
        }
    }
}
