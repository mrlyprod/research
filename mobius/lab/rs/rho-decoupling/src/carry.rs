use num_bigint::BigInt;
use std::collections::BTreeMap;

use crate::{bigpow, pow_checked, ratio_f64};

// MODULAR ARITHMETIC

const MR_LIMIT: u128 = 3_317_044_064_679_887_385_961_981;
const MR_BASES: [u128; 13] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41];

fn addmod(a: u128, b: u128, m: u128) -> u128 {
    let s = a + b;
    if s >= m {
        s - m
    } else {
        s
    }
}

fn mulmod(a: u128, b: u128, m: u128) -> u128 {
    if m < 1u128 << 64 {
        return (a % m) * (b % m) % m;
    }
    let mut r = 0u128;
    let mut x = a % m;
    let mut y = b % m;
    while y > 0 {
        if y & 1 == 1 {
            r = addmod(r, x, m);
        }
        x = addmod(x, x, m);
        y >>= 1;
    }
    r
}

fn powmod(a: u128, mut e: u128, m: u128) -> u128 {
    let mut r = 1u128 % m;
    let mut b = a % m;
    while e > 0 {
        if e & 1 == 1 {
            r = mulmod(r, b, m);
        }
        b = mulmod(b, b, m);
        e >>= 1;
    }
    r
}

fn gcd_u128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

// PRIMALITY

fn is_prime(n: u128) -> bool {
    if n < 2 {
        return false;
    }
    for &p in MR_BASES.iter() {
        if n.is_multiple_of(p) {
            return n == p;
        }
    }
    let mut d = n - 1;
    let mut r = 0;
    while d.is_multiple_of(2) {
        d /= 2;
        r += 1;
    }
    'outer: for &a in MR_BASES.iter() {
        let mut x = powmod(a, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..r {
            x = mulmod(x, x, n);
            if x == n - 1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

// FACTORISATION

fn absdiff(a: u128, b: u128) -> u128 {
    a.abs_diff(b)
}

fn brent(n: u128, c: u128) -> u128 {
    let step = |x: u128| addmod(mulmod(x, x, n), c % n, n);
    let mut y = 2u128;
    let mut g = 1u128;
    let mut r = 1u128;
    let mut q = 1u128;
    let mut x = 0u128;
    let mut ys = 0u128;
    while g == 1 {
        x = y;
        for _ in 0..r {
            y = step(y);
        }
        let mut k = 0u128;
        while k < r && g == 1 {
            ys = y;
            let lim = if r - k < 128 { r - k } else { 128 };
            for _ in 0..lim {
                y = step(y);
                q = mulmod(q, absdiff(x, y), n);
            }
            g = gcd_u128(q, n);
            k += lim;
        }
        r *= 2;
    }
    if g == n {
        g = 1;
        while g == 1 {
            ys = step(ys);
            let d = absdiff(x, ys);
            if d == 0 {
                return n;
            }
            g = gcd_u128(d, n);
        }
    }
    g
}

fn split(n: u128) -> u128 {
    let mut c = 1u128;
    loop {
        let d = brent(n, c);
        if d > 1 && d < n {
            return d;
        }
        c += 1;
    }
}

fn factor_into(n: u128, out: &mut BTreeMap<u128, u32>, unknown: &mut bool) {
    if n == 1 {
        return;
    }
    if is_prime(n) {
        if n >= MR_LIMIT {
            *unknown = true;
        }
        *out.entry(n).or_insert(0) += 1;
        return;
    }
    let d = split(n);
    factor_into(d, out, unknown);
    factor_into(n / d, out, unknown);
}

fn factorise(mut n: u128) -> (BTreeMap<u128, u32>, bool) {
    assert!(n < 1u128 << 126, "modulus outside the mulmod envelope");
    let mut out = BTreeMap::new();
    let mut unknown = false;
    let mut p = 2u128;
    while p <= 100_000 && p * p <= n {
        while n.is_multiple_of(p) {
            *out.entry(p).or_insert(0) += 1;
            n /= p;
        }
        p += 1;
    }
    factor_into(n, &mut out, &mut unknown);
    (out, unknown)
}

fn quot_map(m: &BTreeMap<u128, u32>, g: u64) -> BTreeMap<u128, u32> {
    let mut out = m.clone();
    let mut r = g as u128;
    let mut p = 2u128;
    while r > 1 {
        while r.is_multiple_of(p) {
            let c = out.get_mut(&p).expect("g fails to divide q^t - 1");
            *c -= 1;
            if *c == 0 {
                out.remove(&p);
            }
            r /= p;
        }
        p += 1;
    }
    out
}

fn mu_map(m: &BTreeMap<u128, u32>) -> i8 {
    let mut s = 1i8;
    for e in m.values() {
        if *e >= 2 {
            return 0;
        }
        s = -s;
    }
    s
}

fn fac_string(m: &BTreeMap<u128, u32>) -> String {
    m.iter()
        .map(|(p, e)| {
            if *e == 1 {
                p.to_string()
            } else {
                format!("{p}^{e}")
            }
        })
        .collect::<Vec<String>>()
        .join(" * ")
}

// CYCLOTOMIC SPLIT

fn cyclotomic(q: u64, tmax: usize) -> Vec<BigInt> {
    let mut phi: Vec<BigInt> = vec![BigInt::from(1); tmax + 1];
    for d in 1..=tmax {
        let mut v = bigpow(&BigInt::from(q), d) - BigInt::from(1);
        for (e, p) in phi.iter().enumerate().take(d).skip(1) {
            if d % e == 0 {
                v /= p;
            }
        }
        phi[d] = v;
    }
    phi
}

fn to_u128(x: &BigInt) -> u128 {
    x.to_string()
        .parse::<u128>()
        .expect("cyclotomic value overflows u128")
}

// CARRY DP

fn poly_powers(digits: &[u64], nmax: usize) -> Vec<Vec<u128>> {
    let top = *digits.iter().max().unwrap() as usize;
    let mut out: Vec<Vec<u128>> = vec![vec![1u128]];
    for n in 1..=nmax {
        let prev = &out[n - 1];
        let mut cur = vec![0u128; prev.len() + top];
        for (i, &c) in prev.iter().enumerate() {
            if c == 0 {
                continue;
            }
            for &f in digits {
                cur[i + f as usize] += c;
            }
        }
        out.push(cur);
    }
    out
}

fn count_div(q: u64, polys: &[Vec<u128>], l: usize, t: usize, g: u64) -> u128 {
    let s = l / t;
    let u = l - s * t;
    let dig = (q - 1) / g;
    let cmax = s + 1;
    let mut state = vec![0u128; cmax + 1];
    let mut next = vec![0u128; cmax + 1];
    let mut td = vec![0u64; t];
    let mut total: u128 = 0;
    for j in 0..=(s as u64 + 1) * g {
        let mut carry = 0u64;
        for slot in td.iter_mut() {
            let p = dig * j + carry;
            *slot = p % q;
            carry = p / q;
        }
        let high = carry as usize;
        if high > cmax {
            continue;
        }
        state.iter_mut().for_each(|x| *x = 0);
        state[0] = 1;
        for (c, &want) in td.iter().enumerate().take(t) {
            let sc = if c < u { s + 1 } else { s };
            let poly = &polys[sc];
            next.iter_mut().for_each(|x| *x = 0);
            for (cr, &v) in state.iter().enumerate().take(cmax + 1) {
                if v == 0 {
                    continue;
                }
                let mut m = ((want + q - (cr as u64) % q) % q) as usize;
                while m < poly.len() {
                    let w = poly[m];
                    if w != 0 {
                        let nc = (((cr + m) as u64) - want) / q;
                        next[nc as usize] += v * w;
                    }
                    m += q as usize;
                }
            }
            std::mem::swap(&mut state, &mut next);
        }
        total += state[high];
    }
    total
}

// STUDY

struct Ctx {
    q: u64,
    digits: Vec<u64>,
    label: &'static str,
    gs: Vec<u64>,
    polys: Vec<Vec<u128>>,
    qt: Vec<BigInt>,
    fac: Vec<(BTreeMap<u128, u32>, bool)>,
}

impl Ctx {
    fn new(q: u64, digits: Vec<u64>, label: &'static str, lmax: usize) -> Ctx {
        let gs: Vec<u64> = (1..=q - 1).filter(|g| (q - 1).is_multiple_of(*g)).collect();
        let polys = poly_powers(&digits, lmax);
        let qt: Vec<BigInt> = (0..=lmax)
            .map(|t| bigpow(&BigInt::from(q), t) - BigInt::from(1))
            .collect();
        let phi = cyclotomic(q, lmax);
        let split: Vec<(BTreeMap<u128, u32>, bool)> = (0..=lmax)
            .map(|d| {
                if d == 0 {
                    (BTreeMap::new(), false)
                } else {
                    factorise(to_u128(&phi[d]))
                }
            })
            .collect();
        let mut fac: Vec<(BTreeMap<u128, u32>, bool)> = Vec::with_capacity(lmax + 1);
        for t in 0..=lmax {
            let mut m: BTreeMap<u128, u32> = BTreeMap::new();
            let mut unk = false;
            for (d, (primes, unknown)) in split.iter().enumerate().take(t + 1).skip(1) {
                if t % d == 0 {
                    unk |= unknown;
                    for (p, e) in primes {
                        *m.entry(*p).or_insert(0) += e;
                    }
                }
            }
            fac.push((m, unk));
        }
        Ctx {
            q,
            digits,
            label,
            gs,
            polys,
            qt,
            fac,
        }
    }

    fn k(&self) -> u64 {
        self.digits.len() as u64
    }
}

struct Term {
    t: usize,
    e: BigInt,
    mu: i8,
    n: u128,
    et: BigInt,
    aet: BigInt,
}

struct Res {
    l: usize,
    ratio: f64,
    absnorm: f64,
    live: usize,
    nterms: usize,
    unknown: usize,
    top: Vec<Term>,
}

fn build(c: &Ctx, l: usize) -> (Vec<Term>, usize) {
    let kl = BigInt::from(pow_checked(c.k(), l));
    let two = BigInt::from(2);
    let zero = BigInt::from(0);
    let mut out: Vec<Term> = Vec::new();
    let mut unknown = 0usize;
    for t in 1..=l {
        for &g in &c.gs {
            let e = &c.qt[t] / g;
            if e < two {
                continue;
            }
            if c.fac[t].1 {
                unknown += 1;
                continue;
            }
            let mu = mu_map(&quot_map(&c.fac[t].0, g));
            if mu == 0 {
                continue;
            }
            let n = count_div(c.q, &c.polys, l, t, g);
            let et = &e * BigInt::from(n) - &kl;
            let aet = if et < zero { -&et } else { et.clone() };
            out.push(Term {
                t,
                e,
                mu,
                n,
                et,
                aet,
            });
        }
    }
    (out, unknown)
}

fn study(c: &Ctx, l: usize) -> Res {
    let (mut terms, unknown) = build(c, l);
    let kl = BigInt::from(pow_checked(c.k(), l));
    let zero = BigInt::from(0);
    let mut den = BigInt::from(1);
    for t in terms.iter() {
        den *= &t.e;
    }
    let mut sig = BigInt::from(0);
    let mut abs = BigInt::from(0);
    for t in terms.iter() {
        let w = &den / &t.e;
        sig += BigInt::from(t.mu) * &t.et * &w;
        abs += &t.aet * &w;
    }
    let ratio = if abs == zero {
        0.0
    } else {
        let neg = sig < zero;
        let mag = if neg { -&sig } else { sig.clone() };
        let v = ratio_f64(&mag, &abs);
        if neg {
            -v
        } else {
            v
        }
    };
    let absnorm = ratio_f64(&abs, &(&den * &kl));
    terms.sort_by(|a, b| (&b.aet * &a.e).cmp(&(&a.aet * &b.e)));
    let nterms = terms.len();
    let live = if terms.is_empty() || terms[0].aet == zero {
        0
    } else {
        let (bm, be) = (terms[0].aet.clone(), terms[0].e.clone());
        terms
            .iter()
            .filter(|t| BigInt::from(10) * &t.aet * &be >= &bm * &t.e)
            .count()
    };
    let top: Vec<Term> = terms
        .into_iter()
        .take(4)
        .map(|t| Term {
            t: t.t,
            e: t.e,
            mu: t.mu,
            n: t.n,
            et: t.et,
            aet: t.aet,
        })
        .collect();
    Res {
        l,
        ratio,
        absnorm,
        live,
        nterms,
        unknown,
        top,
    }
}

// RENDER

fn depth_row(c: &Ctx, r: &Res) -> String {
    format!(
        "| {} | {} | {} | {:+.3} | {:.1e} | {} | {} | {} |",
        c.q, c.label, r.l, r.ratio, r.absnorm, r.live, r.nterms, r.unknown
    )
}

fn top_row(c: &Ctx, t: &Term) -> String {
    format!(
        "| {} | {} | {} | {} | {:+} | {} | {} |",
        c.q, c.label, t.t, t.e, t.mu, t.n, t.et
    )
}

fn mu_row(c: &Ctx, g: u64, lmax: usize) -> String {
    let two = BigInt::from(2);
    let mut plus: Vec<String> = Vec::new();
    let mut minus: Vec<String> = Vec::new();
    let mut zeros = 0usize;
    for t in 1..=lmax {
        if &c.qt[t] / g < two {
            continue;
        }
        match mu_map(&quot_map(&c.fac[t].0, g)) {
            1 => plus.push(t.to_string()),
            -1 => minus.push(t.to_string()),
            _ => zeros += 1,
        }
    }
    format!(
        "| {} | Q_t/{} | {} | {} | {} |",
        c.q,
        g,
        plus.join(","),
        minus.join(","),
        zeros
    )
}

fn fac_row(name: &str, m: &BTreeMap<u128, u32>) -> String {
    format!("| {} | {} |", name, fac_string(m))
}

fn seq_row(c: &Ctx, res: &[Res], lo: usize, hi: usize) -> String {
    let vals: Vec<String> = res
        .iter()
        .filter(|r| r.l >= lo && r.l <= hi)
        .map(|r| format!("{:+.2}", r.ratio))
        .collect();
    format!(
        "| {} | {} | {}..{} | {} |",
        c.q,
        c.label,
        lo,
        hi,
        vals.join(" ")
    )
}

fn ex7_digits() -> Vec<u64> {
    (0..10).filter(|&f| f != 7).collect()
}

fn contexts(lmax: usize) -> Vec<Ctx> {
    vec![
        Ctx::new(3, vec![0, 1], "01", lmax),
        Ctx::new(10, ex7_digits(), "ex7", lmax),
    ]
}

fn fac_rows(cs: &[Ctx]) -> Vec<String> {
    let mut out = Vec::new();
    for &t in &[7usize, 37, 39] {
        out.push(fac_row(&format!("3^{t} - 1"), &cs[0].fac[t].0));
    }
    for &t in &[19usize, 23, 31, 37] {
        out.push(fac_row(&format!("R_{t}"), &quot_map(&cs[1].fac[t].0, 9)));
    }
    out
}

pub fn run() {
    let lmax = 40;
    let cs = contexts(lmax);
    let res: Vec<Vec<Res>> = cs
        .iter()
        .map(|c| (3..=lmax).map(|l| study(c, l)).collect())
        .collect();
    println!("carry depth");
    println!("| q | F | L | ratio | abs/k^L | live | terms | unknown |");
    for (i, c) in cs.iter().enumerate() {
        for &l in &[10usize, 20, 30, 40] {
            println!("{}", depth_row(c, &res[i][l - 3]));
        }
    }
    println!("carry top L=40");
    println!("| q | F | t | e | mu | N | eT |");
    for (i, c) in cs.iter().enumerate() {
        for t in res[i][lmax - 3].top.iter() {
            println!("{}", top_row(c, t));
        }
    }
    println!("carry mobius");
    println!("| q | e | mu=+1 | mu=-1 | zeros |");
    for c in cs.iter() {
        for &g in &c.gs {
            println!("{}", mu_row(c, g, lmax));
        }
    }
    println!("carry factors");
    println!("| n | factorisation |");
    for row in fac_rows(&cs) {
        println!("{row}");
    }
    println!("carry ratio sequence");
    println!("| q | F | L | Sigma_L/Abs_L |");
    for (i, c) in cs.iter().enumerate() {
        for &(lo, hi) in &[(3usize, 12usize), (13, 22), (23, 32), (33, 40)] {
            println!("{}", seq_row(c, &res[i], lo, hi));
        }
    }
}

// TESTS

#[cfg(test)]
mod tests {
    use super::*;
    use crate::n_div;

    fn brute(q: u64, digits: &[u64], l: usize, e: u128) -> u128 {
        fn rec(v: u128, len: usize, q: u64, digits: &[u64], l: usize, e: u128, hits: &mut u128) {
            if len == l {
                if v.is_multiple_of(e) {
                    *hits += 1;
                }
                return;
            }
            for &f in digits {
                rec(v * q as u128 + f as u128, len + 1, q, digits, l, e, hits);
            }
        }
        let mut hits = 0;
        rec(0, 0, q, digits, l, e, &mut hits);
        hits
    }

    #[test]
    fn carry_matches_brute() {
        for (q, digits, l) in [
            (3u64, vec![0u64, 1], 8usize),
            (10, ex7_digits(), 5),
            (5, vec![0, 2, 4], 6),
        ] {
            let polys = poly_powers(&digits, l);
            let gs: Vec<u64> = (1..=q - 1).filter(|g| (q - 1) % g == 0).collect();
            for t in 1..=l {
                let qt = (q as u128).pow(t as u32) - 1;
                for &g in &gs {
                    let e = qt / g as u128;
                    if e < 2 {
                        continue;
                    }
                    assert_eq!(
                        count_div(q, &polys, l, t, g),
                        brute(q, &digits, l, e),
                        "q={q} l={l} t={t} g={g}"
                    );
                }
            }
        }
    }

    #[test]
    fn carry_matches_residue_dp() {
        let mut checks = 0;
        for (q, digits) in [(3u64, vec![0u64, 1]), (10, ex7_digits())] {
            let polys = poly_powers(&digits, 12);
            let gs: Vec<u64> = (1..=q - 1).filter(|g| (q - 1) % g == 0).collect();
            for &l in &[8usize, 12] {
                for t in 1..=l {
                    let qt = (q as u128).pow(t as u32) - 1;
                    for &g in &gs {
                        let e = qt / g as u128;
                        if !(2..=30000).contains(&e) {
                            continue;
                        }
                        assert_eq!(
                            count_div(q, &polys, l, t, g),
                            n_div(q, &digits, l, e as u64),
                            "q={q} l={l} t={t} g={g}"
                        );
                        checks += 1;
                    }
                }
            }
        }
        assert!(
            checks >= 40,
            "residue cross-check too thin at {checks} cells"
        );
    }

    #[test]
    fn mobius_is_certified() {
        for c in contexts(40).iter() {
            for t in 1..=40 {
                assert!(!c.fac[t].1, "uncertified cofactor at q={} t={t}", c.q);
            }
        }
    }

    #[test]
    fn mobius_rows_pinned() {
        let cs = contexts(40);
        let mut got: Vec<String> = Vec::new();
        for c in cs.iter() {
            for &g in &c.gs {
                got.push(mu_row(c, g, 40));
            }
        }
        got.extend(fac_rows(&cs));
        let want = [
            "| 3 | Q_t/1 | 3,7,13,21,27,29,31 | 1,9,11,17,19,23,33,37 | 25 |",
            "| 3 | Q_t/2 | 9,11,17,19,23,33,37 | 3,7,13,21,27,29,31 | 25 |",
            "| 10 | Q_t/1 |  |  | 40 |",
            "| 10 | Q_t/3 | 2,13,19,20,23,25,29,31,32,35,37,38,40 | 1,4,5,7,8,10,11,14,16,17,26,28,34 | 14 |",
            "| 10 | Q_t/9 | 3,4,5,7,8,10,11,14,15,16,17,24,26,28,33,34,39 | 2,6,12,13,19,20,21,23,25,29,30,31,32,35,37,38,40 | 5 |",
            "| 3^7 - 1 | 2 * 1093 |",
            "| 3^37 - 1 | 2 * 13097927 * 17189128703 |",
            "| 3^39 - 1 | 2 * 13^2 * 313 * 6553 * 7333 * 797161 |",
            "| R_19 | 1111111111111111111 |",
            "| R_23 | 11111111111111111111111 |",
            "| R_31 | 2791 * 6943319 * 57336415063790604359 |",
            "| R_37 | 2028119 * 247629013 * 2212394296770203368013 |",
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn carry_rows_pinned() {
        let cs = contexts(40);
        let res: Vec<Vec<Res>> = cs
            .iter()
            .map(|c| (3..=40).map(|l| study(c, l)).collect())
            .collect();
        let mut got: Vec<String> = Vec::new();
        for (i, c) in cs.iter().enumerate() {
            for &l in &[10usize, 20, 30, 40] {
                got.push(depth_row(c, &res[i][l - 3]));
            }
        }
        for (i, c) in cs.iter().enumerate() {
            for t in res[i][37].top.iter() {
                got.push(top_row(c, t));
            }
        }
        for (i, c) in cs.iter().enumerate() {
            for &(lo, hi) in &[(3usize, 12usize), (13, 22), (23, 32), (33, 40)] {
                got.push(seq_row(c, &res[i], lo, hi));
            }
        }
        let want = [
            "| 3 | 01 | 10 | -0.211 | 4.7e-2 | 4 | 7 | 0 |",
            "| 3 | 01 | 20 | -0.123 | 4.1e-3 | 6 | 15 | 0 |",
            "| 3 | 01 | 30 | +0.069 | 1.0e-3 | 6 | 23 | 0 |",
            "| 3 | 01 | 40 | -0.498 | 2.1e-4 | 6 | 29 | 0 |",
            "| 10 | ex7 | 10 | +0.812 | 1.1e-5 | 5 | 15 | 0 |",
            "| 10 | ex7 | 20 | -0.495 | 3.5e-8 | 2 | 31 | 0 |",
            "| 10 | ex7 | 30 | -0.127 | 4.7e-10 | 2 | 44 | 0 |",
            "| 10 | ex7 | 40 | -0.192 | 3.9e-12 | 4 | 60 | 0 |",
            "| 3 | 01 | 7 | 1093 | -1 | 1059181242 | 58173469730 |",
            "| 3 | 01 | 9 | 19682 | -1 | 107075926 | 1007956747756 |",
            "| 3 | 01 | 7 | 2186 | +1 | 460789966 | -92224762100 |",
            "| 3 | 01 | 9 | 9841 | +1 | 148454776 | 361431822840 |",
            "| 10 | ex7 | 5 | 33333 | -1 | 4434309484921670553282063418321580 | 8646548121236467809716529928539 |",
            "| 10 | ex7 | 8 | 33333333 | -1 | 4434111044140794648321659572234 | -5129421023116418959440572469821679 |",
            "| 10 | ex7 | 10 | 3333333333 | -1 | 44419676099006184777353657556 | 256757567534800575426920054498816547 |",
            "| 10 | ex7 | 7 | 3333333 | -1 | 44342590913069608659343530190385 | -207818310865474807662686276694396 |",
            "| 3 | 01 | 3..12 | -0.33 -0.64 -1.00 -1.00 -0.81 -0.83 -0.33 -0.21 -0.09 -0.09 |",
            "| 3 | 01 | 13..22 | -0.13 -0.24 -0.36 -0.41 -0.40 -0.33 -0.26 -0.12 +0.02 +0.14 |",
            "| 3 | 01 | 23..32 | +0.19 +0.18 +0.19 +0.17 +0.12 +0.10 +0.09 +0.07 +0.07 +0.06 |",
            "| 3 | 01 | 33..40 | +0.02 -0.05 -0.17 -0.23 -0.31 -0.38 -0.44 -0.50 |",
            "| 10 | ex7 | 3..12 | +0.76 +0.12 +0.02 +0.15 -0.31 +0.57 +0.53 +0.81 +0.92 +0.13 |",
            "| 10 | ex7 | 13..22 | +0.57 +0.82 +0.20 -0.82 -0.79 -0.73 -0.24 -0.49 -0.86 -0.91 |",
            "| 10 | ex7 | 23..32 | -0.95 -0.31 +0.68 +0.63 +0.58 -0.96 -0.99 -0.13 -0.10 -0.27 |",
            "| 10 | ex7 | 33..40 | -0.31 -0.45 -0.60 -0.64 -0.65 -0.59 -0.23 -0.19 |",
        ];
        assert_eq!(got, want);
    }
}
