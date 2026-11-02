use num_bigint::BigInt;

mod carry;
mod large;
mod menergy;
mod riesz;
mod signed;
mod vaughan;

// EXACT DP

fn pow_checked(k: u64, l: usize) -> u128 {
    let mut p: u128 = 1;
    for _ in 0..l {
        p = p.checked_mul(k as u128).expect("k^L overflows u128");
    }
    p
}

fn residue_counts(q: u64, digits: &[u64], l: usize, d: u64) -> Vec<u128> {
    let d = d as usize;
    let mut state = vec![0u128; d];
    state[0] = 1;
    let mut next = vec![0u128; d];
    for _ in 0..l {
        next.iter_mut().for_each(|x| *x = 0);
        for (r, &count) in state.iter().enumerate() {
            if count == 0 {
                continue;
            }
            let base = (r as u64 * q) % d as u64;
            for &f in digits {
                let idx = ((base + f) % d as u64) as usize;
                next[idx] += count;
            }
        }
        std::mem::swap(&mut state, &mut next);
    }
    state
}

fn n_div(q: u64, digits: &[u64], l: usize, d: u64) -> u128 {
    residue_counts(q, digits, l, d)[0]
}

// ARITHMETIC HELPERS

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

fn diff_gcd(digits: &[u64]) -> u64 {
    let mut g = 0;
    for w in digits.windows(2) {
        g = gcd(g, w[1] - w[0]);
    }
    g.max(1)
}

fn mu_sieve(n: usize) -> Vec<i8> {
    let mut mu = vec![1i8; n + 1];
    let mut primes: Vec<usize> = Vec::new();
    let mut composite = vec![false; n + 1];
    for i in 2..=n {
        if !composite[i] {
            primes.push(i);
            mu[i] = -1;
        }
        for &p in &primes {
            let ip = i * p;
            if ip > n {
                break;
            }
            composite[ip] = true;
            if i % p == 0 {
                mu[ip] = 0;
                break;
            }
            mu[ip] = -mu[i];
        }
    }
    mu[0] = 0;
    mu
}

fn mult_order(q: u64, d: u64) -> u64 {
    let mut x = q % d;
    let mut t = 1;
    while x != 1 {
        x = (x * q) % d;
        t += 1;
        if t > d {
            return 0;
        }
    }
    t
}

// EXACT READOUT

fn bigpow(base: &BigInt, l: usize) -> BigInt {
    let mut p = BigInt::from(1);
    for _ in 0..l {
        p *= base;
    }
    p
}

fn ratio_f64(num: &BigInt, den: &BigInt) -> f64 {
    let scaled = (num * bigpow(&BigInt::from(10), 40)) / den;
    let v: f64 = scaled.to_string().parse().unwrap();
    v / 1e40
}

struct Frac {
    num: BigInt,
    den: BigInt,
}

impl Frac {
    fn zero() -> Frac {
        Frac {
            num: BigInt::from(0),
            den: BigInt::from(1),
        }
    }

    fn add(&mut self, num: BigInt, den: BigInt) {
        self.num = &self.num * &den + num * &self.den;
        self.den = &self.den * den;
    }

    fn to_f64(&self) -> f64 {
        let neg = self.num < BigInt::from(0);
        let mag = if neg { -&self.num } else { self.num.clone() };
        let v = ratio_f64(&mag, &self.den);
        if neg {
            -v
        } else {
            v
        }
    }
}

fn gamma_emp(digits: &[u64], d: u64) -> f64 {
    let k = digits.len() as f64;
    let mut best = 0.0f64;
    for a in 1..d {
        let mut re = 0.0;
        let mut im = 0.0;
        for &f in digits {
            let t = 2.0 * std::f64::consts::PI * (a as f64) * (f as f64) / (d as f64);
            re += t.cos();
            im += t.sin();
        }
        best = best.max(re.hypot(im) / k);
    }
    best
}

fn assert_lemma_a(errnum: &BigInt, kl: &BigInt, k: u64, d: u64, l: usize) {
    let kd2 = BigInt::from(k * k * d * d);
    let kd2m8 = &kd2 - BigInt::from(8);
    let lhs = errnum * bigpow(&kd2, l);
    let rhs = BigInt::from(d) * kl * bigpow(&kd2m8, l);
    assert!(lhs <= rhs, "lemma A fails at k={k} d={d} l={l}");
}

// CENSUS

fn census(q: u64, digits: &[u64], label: &str, ls: &[usize], dmax: u64, mu: &[i8]) {
    let k = digits.len() as u64;
    let delta = diff_gcd(digits);
    let mut lastcop: Vec<(usize, f64)> = Vec::new();
    for &l in ls {
        let kl_u = pow_checked(k, l);
        let kl = BigInt::from(kl_u);
        let mut worst = (0.0f64, 0u64);
        let mut worstcop = (0.0f64, 0u64);
        let mut slack = 0.0f64;
        let mut s1: i128 = 0;
        let mut t = Frac::zero();
        let mut ts = Frac::zero();
        for d in 2..=dmax {
            if gcd(d, q) != 1 {
                continue;
            }
            let n = n_div(q, digits, l, d);
            let signed = BigInt::from(d) * BigInt::from(n) - &kl;
            let errnum = if signed < BigInt::from(0) {
                -&signed
            } else {
                signed.clone()
            };
            let norm = ratio_f64(&errnum, &kl);
            if norm > worst.0 {
                worst = (norm, d);
            }
            if gcd(d, delta) == 1 {
                assert_lemma_a(&errnum, &kl, k, d, l);
                if norm > worstcop.0 {
                    worstcop = (norm, d);
                }
                let kd2 = (k * k * d * d) as f64;
                let bound = d as f64 * ((kd2 - 8.0) / kd2).powi(l as i32);
                if bound > 0.0 {
                    slack = slack.max(norm / bound);
                }
            }
            if mu[d as usize] != 0 {
                s1 += mu[d as usize] as i128 * n as i128;
                t.add(
                    BigInt::from(mu[d as usize]) * &signed,
                    BigInt::from(d) * &kl,
                );
                ts.add(errnum, BigInt::from(d) * &kl);
            }
        }
        let rate = worstcop.0.powf(1.0 / l as f64);
        let gam = gamma_emp(digits, worstcop.1);
        let ord = mult_order(q, worstcop.1);
        println!(
            "row q={q} F={label} L={l} D={dmax} worstcop={:.4e} at d={} ord={ord} rate={rate:.4} gam={gam:.4} slack={slack:.3}",
            worstcop.0, worstcop.1
        );
        if delta > 1 {
            println!(
                "wallrow q={q} F={label} L={l} D={dmax} delta={delta} worst={:.4e} at d={}",
                worst.0, worst.1
            );
        }
        let tv = t.to_f64();
        let tsv = ts.to_f64();
        let ratio = if tsv > 0.0 { tv.abs() / tsv } else { 0.0 };
        println!("typeI q={q} F={label} L={l} D={dmax} S1={s1} relT={tv:.4e} reltriv={tsv:.4e} ratio={ratio:.3e}");
        lastcop.push((l, worstcop.0));
    }
    if lastcop.len() >= 2 {
        let (l1, w1) = lastcop[lastcop.len() - 2];
        let (l2, w2) = lastcop[lastcop.len() - 1];
        if w1 > 0.0 && w2 > 0.0 {
            let per = (w2 / w1).powf(1.0 / (l2 - l1) as f64);
            println!("decay q={q} F={label} D={dmax} L={l1}..{l2} factor={per:.4}");
        }
    }
    if gcd(7, q) == 1 && gcd(7, delta) == 1 {
        let l = *ls.last().unwrap();
        let kl = BigInt::from(pow_checked(k, l));
        let n = n_div(q, digits, l, 7);
        let signed = BigInt::from(7) * BigInt::from(n) - &kl;
        let errnum = if signed < BigInt::from(0) {
            -signed
        } else {
            signed
        };
        let rate = (ratio_f64(&errnum, &kl) / 7.0).powf(1.0 / l as f64);
        println!(
            "decouple q={q} k={k} F={label} L={l} d=7 rate={rate:.4} gam={:.4}",
            gamma_emp(digits, 7)
        );
    }
}

fn pinned(q: u64, digits: &[u64], label: &str, l: usize, d: u64) {
    let k = digits.len() as u64;
    let kl = BigInt::from(pow_checked(k, l));
    let n = n_div(q, digits, l, d);
    let signed = BigInt::from(d) * BigInt::from(n) - &kl;
    let errnum = if signed < BigInt::from(0) {
        -signed
    } else {
        signed
    };
    let rate = (ratio_f64(&errnum, &kl) / d as f64).powf(1.0 / l as f64);
    let mean = rate * (d as f64).powf(1.0 / l as f64);
    let ord = mult_order(q, d);
    println!("pinned q={q} F={label} L={l} d={d} ord={ord} rate={rate:.4} mean={mean:.4}");
}

fn main() {
    let mu = mu_sieve(500);
    let ex = |q: u64, e: u64| -> Vec<u64> { (0..q).filter(|&f| f != e).collect() };
    census(3, &[0, 1], "01", &[8, 16, 24, 32], 200, &mu);
    census(3, &[0, 2], "02", &[8, 16, 24, 32], 200, &mu);
    census(3, &[1, 2], "12", &[8, 16, 24, 32], 200, &mu);
    census(4, &[0, 1, 2], "012", &[8, 16, 24], 200, &mu);
    census(5, &[0, 1, 2, 3], "0123", &[8, 16, 24], 200, &mu);
    census(5, &[0, 2, 4], "024", &[8, 16, 20], 200, &mu);
    census(10, &ex(10, 7), "ex7", &[6, 12, 18, 24], 300, &mu);
    census(100, &ex(100, 37), "ex37", &[4, 8, 12, 16], 500, &mu);
    census(
        100,
        &(0..50).collect::<Vec<u64>>(),
        "0to49",
        &[4, 8, 12, 16],
        500,
        &mu,
    );
    census(100, &[0, 1], "01", &[16, 32, 64, 96], 500, &mu);
    pinned(100, &ex(100, 37), "ex37", 12, 101);
    pinned(100, &ex(100, 37), "ex37", 12, 9999);
    pinned(100, &ex(100, 37), "ex37", 12, 3367);
    pinned(100, &ex(100, 37), "ex37", 12, 999999);
    carry::run();
    riesz::run();
    large::run();
    menergy::run();
    signed::run();
    vaughan::run();
}

// TESTS

#[cfg(test)]
mod tests {
    use super::*;

    fn brute(q: u64, digits: &[u64], l: usize, d: u64) -> u128 {
        fn rec(v: u128, len: usize, q: u64, digits: &[u64], l: usize, d: u64, hits: &mut u128) {
            if len == l {
                if v.is_multiple_of(d as u128) {
                    *hits += 1;
                }
                return;
            }
            for &f in digits {
                rec(v * q as u128 + f as u128, len + 1, q, digits, l, d, hits);
            }
        }
        let mut hits = 0;
        rec(0, 0, q, digits, l, d, &mut hits);
        hits
    }

    #[test]
    fn dp_matches_brute() {
        for d in 1..=40 {
            assert_eq!(n_div(3, &[0, 1], 8, d), brute(3, &[0, 1], 8, d));
        }
        for d in 1..=30 {
            assert_eq!(n_div(5, &[1, 3, 4], 6, d), brute(5, &[1, 3, 4], 6, d));
        }
        let f10: Vec<u64> = (0..10).filter(|&f| f != 7).collect();
        for d in 1..=30 {
            assert_eq!(n_div(10, &f10, 4, d), brute(10, &f10, 4, d));
        }
        for d in 1..=25 {
            assert_eq!(n_div(100, &[0, 17], 3, d), brute(100, &[0, 17], 3, d));
        }
    }

    #[test]
    fn residues_sum_to_kl() {
        for (q, digits, l, d) in [
            (3u64, vec![0u64, 1], 12usize, 35u64),
            (10, vec![1, 4, 9], 8, 77),
            (100, vec![0, 1], 20, 99),
        ] {
            let total: u128 = residue_counts(q, &digits, l, d).iter().sum();
            assert_eq!(total, pow_checked(digits.len() as u64, l));
        }
    }

    #[test]
    fn crt_reduction_exact() {
        let f = [0u64, 1, 3, 4];
        let low = residue_counts(6, &f, 4, 5);
        let split = low[0] + low[1];
        assert_eq!(n_div(6, &f, 5, 10), split);
    }

    #[test]
    fn mu_pins() {
        let mu = mu_sieve(100);
        assert_eq!(&mu[1..13], &[1, -1, -1, 0, -1, 1, -1, 0, 0, 1, -1, 0]);
        let m100: i32 = (1..=100).map(|d| mu[d] as i32).sum();
        assert_eq!(m100, 1);
    }

    #[test]
    fn lemma_a_exact_small() {
        let l = 16;
        let kl = BigInt::from(pow_checked(2, l));
        for d in 2..=60u64 {
            if gcd(d, 3) != 1 {
                continue;
            }
            let n = n_div(3, &[0, 1], l, d);
            let signed = BigInt::from(d) * BigInt::from(n) - &kl;
            let errnum = if signed < BigInt::from(0) {
                -signed
            } else {
                signed
            };
            assert_lemma_a(&errnum, &kl, 2, d, l);
        }
    }

    #[test]
    fn helpers_pin() {
        assert_eq!(pow_checked(3, 4), 81);
        assert_eq!(diff_gcd(&[0, 2, 4]), 2);
        assert_eq!(diff_gcd(&[0, 1]), 1);
        assert_eq!(mult_order(10, 7), 6);
        assert_eq!(mult_order(3, 8), 2);
    }
}
