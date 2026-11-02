use std::collections::HashMap;
use std::env;
use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::thread;
use std::time::Instant;

// DESIGNS

const DESIGNS: usize = 57;
const COIN: usize = 57;
const STREAMS: usize = 58;
const LAGS: usize = 8;
const GOLDEN: u64 = 0x9e37_79b9_7f4a_7c15;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Design {
    b: u64,
    mask: u32,
}

impl Design {
    fn all() -> Vec<Design> {
        (2..=6u64)
            .flat_map(|b| (1..1u32 << (b - 1)).map(move |mask| Design { b, mask }))
            .collect()
    }

    fn find(name: &str) -> Design {
        Design::all()
            .into_iter()
            .find(|d| d.name() == name)
            .unwrap()
    }

    fn index(self) -> usize {
        offset(self.b) + self.mask as usize - 1
    }

    fn has(self, d: u64) -> bool {
        d > 0 && (self.mask >> (d - 1)) & 1 == 1
    }

    fn sign(self, d: u64) -> i64 {
        if self.has(d) {
            -1
        } else {
            1
        }
    }

    fn parity(self, mut v: u128) -> u64 {
        let b = self.b as u128;
        let mut p = 0;
        while v > 0 {
            p ^= self.has((v % b) as u64) as u64;
            v /= b;
        }
        p
    }

    fn name(self) -> String {
        let digits: String = (1..self.b)
            .filter(|&d| self.has(d))
            .map(|d| char::from(b'0' + d as u8))
            .collect();
        format!("{}:{}", self.b, digits)
    }

    fn mirror(self) -> Design {
        let top = self.sign(self.b - 1);
        let mask = (1..self.b)
            .filter(|&d| top * self.sign(self.b - 1 - d) == -1)
            .map(|d| 1 << (d - 1))
            .sum();
        Design { b: self.b, mask }
    }

    fn twin(self) -> Option<Design> {
        if self.b % 2 == 0 {
            return None;
        }
        let odd: u32 = (1..self.b)
            .filter(|d| d % 2 == 1)
            .map(|d| 1 << (d - 1))
            .sum();
        let mask = self.mask ^ odd;
        (mask != 0).then_some(Design { b: self.b, mask })
    }
}

fn offset(b: u64) -> usize {
    (2..b).map(|c| (1usize << (c - 1)) - 1).sum()
}

fn label(i: usize) -> String {
    if i == COIN {
        "coin".into()
    } else {
        Design::all()[i].name()
    }
}

fn mix(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

fn coin(n: u64, salt: u64) -> u64 {
    mix(salt.wrapping_add((n + 1).wrapping_mul(GOLDEN))) >> 63
}

fn tags() -> Vec<&'static str> {
    let all = Design::all();
    let seqs: Vec<Vec<u64>> = all
        .iter()
        .map(|d| (0..1u128 << 12).map(|n| d.parity(n)).collect())
        .collect();
    (0..DESIGNS)
        .map(|i| {
            if (0..i).any(|j| seqs[j] == seqs[i]) {
                "copy"
            } else if seqs[i].iter().enumerate().all(|(n, &p)| p == n as u64 % 2) {
                "folded"
            } else if all[i].twin().is_some() {
                "twin"
            } else {
                ""
            }
        })
        .collect()
}

// COMPLEX

#[derive(Clone, Copy, Debug, PartialEq)]
struct Cx {
    re: f64,
    im: f64,
}

impl Cx {
    const ZERO: Cx = Cx { re: 0.0, im: 0.0 };
    const ONE: Cx = Cx { re: 1.0, im: 0.0 };

    fn real(x: f64) -> Cx {
        Cx { re: x, im: 0.0 }
    }

    fn turn(r: u64, q: u64) -> Cx {
        let t = std::f64::consts::TAU * (r % q) as f64 / q as f64;
        Cx {
            re: t.cos(),
            im: t.sin(),
        }
    }

    fn add(self, o: Cx) -> Cx {
        Cx {
            re: self.re + o.re,
            im: self.im + o.im,
        }
    }

    fn sub(self, o: Cx) -> Cx {
        Cx {
            re: self.re - o.re,
            im: self.im - o.im,
        }
    }

    fn mul(self, o: Cx) -> Cx {
        Cx {
            re: self.re * o.re - self.im * o.im,
            im: self.re * o.im + self.im * o.re,
        }
    }

    fn scale(self, s: f64) -> Cx {
        Cx {
            re: self.re * s,
            im: self.im * s,
        }
    }

    fn div(self, o: Cx) -> Cx {
        let d = o.re * o.re + o.im * o.im;
        Cx {
            re: (self.re * o.re + self.im * o.im) / d,
            im: (self.im * o.re - self.re * o.im) / d,
        }
    }

    fn pow(self, n: u32) -> Cx {
        (0..n).fold(Cx::ONE, |a, _| a.mul(self))
    }

    fn abs(self) -> f64 {
        self.re.hypot(self.im)
    }

    fn turns(self) -> f64 {
        self.im.atan2(self.re) / std::f64::consts::TAU
    }
}

// ARITHMETIC

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn inverse(a: u64, m: u64) -> u64 {
    let (mut r0, mut r1) = ((a % m) as i128, m as i128);
    let (mut s0, mut s1) = (1i128, 0i128);
    while r1 != 0 {
        let q = r0 / r1;
        (r0, r1) = (r1, r0 - q * r1);
        (s0, s1) = (s1, s0 - q * s1);
    }
    assert!(m == 1 || r0 == 1);
    s0.rem_euclid(m as i128) as u64
}

fn gauss(k: u32, a: u64, q: u64) -> Cx {
    (0..q).fold(Cx::ZERO, |s, x| {
        s.add(Cx::turn(
            (a as u128 * (x as u128).pow(k) % q as u128) as u64,
            q,
        ))
    })
}

// POINTS

struct Point {
    num: u64,
    den: u64,
    tile: Vec<Cx>,
    sigma: Cx,
    eta: f64,
}

impl Point {
    fn all(d: Design) -> Vec<Point> {
        let b = d.b;
        let big = b * b;
        let m = big - 1;
        (0..m)
            .map(|c| {
                let tile: Vec<Cx> = (0..big)
                    .map(|t| Cx::turn(c * t % m, m).scale((d.sign(t % b) * d.sign(t / b)) as f64))
                    .collect();
                let sigma = tile.iter().fold(Cx::ZERO, |s, &y| s.add(y));
                let eta = if sigma.abs() < 1e-9 {
                    f64::INFINITY
                } else {
                    1.0 - sigma.abs().ln() / (big as f64).ln()
                };
                let g = gcd(c, m);
                Point {
                    num: c / g,
                    den: m / g,
                    tile,
                    sigma,
                    eta,
                }
            })
            .collect()
    }

    fn riesz(&self) -> Vec<Cx> {
        let big = self.tile.len() as f64;
        self.tile
            .iter()
            .map(|y| y.scale(big).div(self.sigma))
            .collect()
    }

    fn weight(&self, k: u32) -> Cx {
        gauss(k, (self.den - self.num % self.den) % self.den, self.den).scale(1.0 / self.den as f64)
    }

    fn name(&self) -> String {
        format!("{}/{}", self.num, self.den)
    }
}

// ADIC

fn units(p: u64, k: u32) -> (u32, Vec<u64>) {
    let c = match p {
        2 => 3,
        3 => 2,
        _ => 1,
    };
    let m = p.pow(c);
    let mut s: Vec<u64> = (1..m)
        .filter(|u| u % p != 0)
        .map(|u| u.pow(k) % m)
        .collect();
    s.sort_unstable();
    s.dedup();
    (c, s)
}

fn expect(p: u64, t: u32, e: u32, r: &[Cx], k: u32) -> Cx {
    let (c, allowed) = units(p, k);
    let big = p.pow(t);
    let digits = (e + c + t - 1) / t + 1;
    let span = p.pow(t * digits - e);
    let (pc, pe) = (p.pow(c), p.pow(e));
    let mut sum = Cx::ZERO;
    let mut count = 0u64;
    for w in (0..span).filter(|w| w % p != 0 && allowed.contains(&(w % pc))) {
        let mut v = pe * w;
        let mut prod = Cx::ONE;
        for _ in 0..digits {
            prod = prod.mul(r[(v % big) as usize]);
            v /= big;
        }
        sum = sum.add(prod);
        count += 1;
    }
    sum.scale(1.0 / count as f64)
}

fn prime_power(p: u64, t: u32, r: &[Cx], k: u32) -> Cx {
    let big = p.pow(t) as f64;
    let mut sum = Cx::ZERO;
    for j in 0..t {
        let (a, e) = (k * j / t, k * j % t);
        let weight = (1.0 - 1.0 / p as f64) / (p as f64).powi(j as i32);
        sum = sum.add(r[0].pow(a).mul(expect(p, t, e, r, k)).scale(weight));
    }
    sum.div(Cx::ONE.sub(r[0].pow(k).scale(1.0 / big)))
}

struct Six<'a> {
    r: &'a [Cx],
    memo: HashMap<(u32, u32, u64, u32, u64), Cx>,
    inv9: Vec<u64>,
    inv4: Vec<u64>,
}

impl Six<'_> {
    fn new(r: &[Cx]) -> Six<'_> {
        Six {
            r,
            memo: HashMap::new(),
            inv9: (0..64).map(|n| inverse(9, 1u64 << n)).collect(),
            inv4: (0..41).map(|n| inverse(4, 3u64.pow(n))).collect(),
        }
    }

    fn class(&mut self, left: u32, s2: u32, t2: u64, s3: u32, t3: u64) -> Cx {
        if left == 0 || (s2 == 0 && s3 == 0) {
            return Cx::ONE;
        }
        if let Some(&v) = self.memo.get(&(left, s2, t2, s3, t3)) {
            return v;
        }
        let (l2, l3) = (s2.min(2), s3.min(2));
        let (m2, m3) = (1u64 << l2, 3u64.pow(l3));
        let (q2, q3) = (1u128 << s2, 3u128.pow(s3));
        let (n2, n3) = (s2 - l2, s3 - l3);
        let mut sum = Cx::ZERO;
        let mut count = 0;
        for d in (0..36u64).filter(|d| d % m2 == t2 % m2 && d % m3 == t3 % m3) {
            let u2 = if n2 > 0 {
                let diff = (t2 as i128 - d as i128).rem_euclid(q2 as i128) as u128 / 4;
                (diff * self.inv9[n2 as usize] as u128 % (1u128 << n2)) as u64
            } else {
                0
            };
            let u3 = if n3 > 0 {
                let diff = (t3 as i128 - d as i128).rem_euclid(q3 as i128) as u128 / 9;
                (diff * self.inv4[n3 as usize] as u128 % 3u128.pow(n3)) as u64
            } else {
                0
            };
            let next = if left == u32::MAX { left } else { left - 1 };
            let inner = self.class(next, n2, u2, n3, u3);
            sum = sum.add(self.r[d as usize].mul(inner));
            count += 1;
        }
        let v = sum.scale(1.0 / count as f64);
        self.memo.insert((left, s2, t2, s3, t3), v);
        v
    }
}

fn coset(g: u64, i: u32, allowed: &[u64], m: u64) -> Vec<u64> {
    let gi = (0..i).fold(1u64, |a, _| a * g % m);
    let mut s: Vec<u64> = allowed.iter().map(|a| a * gi % m).collect();
    s.sort_unstable();
    s.dedup();
    s
}

fn six_adic(r: &[Cx], k: u32, emax: u32) -> (Cx, Vec<Cx>) {
    let (c2, a2) = units(2, k);
    let (c3, a3) = units(3, k);
    let mut six = Six::new(r);
    let denom = Cx::ONE.sub(r[0].pow(k).scale(1.0 / 36.0)).scale(3.0);
    let mut terms = vec![Cx::ZERO; emax as usize + 1];
    for j in 0..2u32 {
        let (f, a) = (k * j % 2, k * j / 2);
        let lead = r[0].pow(a).scale(6f64.powi(-(j as i32)));
        for e in 0..=emax {
            let pairs = if e == 0 {
                vec![(0, 0)]
            } else {
                vec![(e, 0), (0, e)]
            };
            for (e2, e3) in pairs {
                let weight = 2f64.powi(-(e2 as i32)) * 3f64.powi(-(e3 as i32));
                let (v2, v3) = (f + k * e2, f + k * e3);
                let s2 = coset(3, f + k * e3, &a2, 1 << c2);
                let s3 = coset(2, f + k * e2, &a3, 3u64.pow(c3));
                let mut avg = Cx::ZERO;
                for &x2 in &s2 {
                    for &x3 in &s3 {
                        let t2 = ((1u128 << v2) * x2 as u128 % (1u128 << (v2 + c2))) as u64;
                        let t3 = (3u128.pow(v3) * x3 as u128 % 3u128.pow(v3 + c3)) as u64;
                        avg = avg.add(six.class(u32::MAX, v2 + c2, t2, v3 + c3, t3));
                    }
                }
                let part = avg.scale(weight / (s2.len() * s3.len()) as f64);
                terms[e as usize] = terms[e as usize].add(lead.mul(part).div(denom));
            }
        }
    }
    (terms.iter().fold(Cx::ZERO, |s, &t| s.add(t)), terms)
}

#[cfg(test)]
fn six_local(r: &[Cx], depth: u32) -> Cx {
    let (c2, a2) = units(2, 2);
    let (c3, a3) = units(3, 2);
    let mut six = Six::new(r);
    let mut total = r[0].pow(depth).scale(6f64.powi(-(depth as i32)));
    for j in 0..depth {
        let lead = r[0].pow(j).scale(6f64.powi(-(j as i32)) / 3.0);
        for e in 0..=24u32 {
            let pairs = if e == 0 {
                vec![(0, 0)]
            } else if e <= 18 {
                vec![(e, 0), (0, e)]
            } else {
                vec![(e, 0)]
            };
            for (e2, e3) in pairs {
                let weight = 2f64.powi(-(e2 as i32)) * 3f64.powi(-(e3 as i32));
                let (v2, v3) = (2 * e2, 2 * e3);
                let s2 = coset(3, 2 * e3, &a2, 1 << c2);
                let s3 = coset(2, 2 * e2, &a3, 3u64.pow(c3));
                let mut avg = Cx::ZERO;
                for &x2 in &s2 {
                    for &x3 in &s3 {
                        let t2 = ((1u128 << v2) * x2 as u128 % (1u128 << (v2 + c2))) as u64;
                        let t3 = (3u128.pow(v3) * x3 as u128 % 3u128.pow(v3 + c3)) as u64;
                        avg = avg.add(six.class(depth - j, v2 + c2, t2, v3 + c3, t3));
                    }
                }
                total = total.add(lead.mul(avg).scale(weight / (s2.len() * s3.len()) as f64));
            }
        }
    }
    total
}

fn adic(b: u64, p: &Point, k: u32, emax: u32) -> (Cx, Vec<Cx>) {
    let r = p.riesz();
    match b {
        2 => (prime_power(2, 2, &r, k), vec![]),
        3 => (prime_power(3, 2, &r, k), vec![]),
        4 => (prime_power(2, 4, &r, k), vec![]),
        5 => (prime_power(5, 2, &r, k), vec![]),
        _ => six_adic(&r, k, emax),
    }
}

fn local(p: &Point, k: u32, digits: u32) -> Cx {
    let big = p.tile.len() as u64;
    let q = big.pow(digits);
    let mut sum = Cx::ZERO;
    for n in 0..q {
        let mut v = ((n as u128).pow(k) % q as u128) as u64;
        let mut prod = Cx::ONE;
        for _ in 0..digits {
            prod = prod.mul(p.tile[(v % big) as usize]);
            v /= big;
        }
        sum = sum.add(prod);
    }
    sum.div(p.sigma.pow(digits))
}

fn emax(k: u32) -> u32 {
    if k == 2 {
        14
    } else {
        10
    }
}

#[derive(Clone)]
struct Arc {
    num: u64,
    den: u64,
    eta: f64,
    weight: Cx,
    adic: Cx,
}

impl Arc {
    fn amp(&self) -> Cx {
        self.weight.mul(self.adic)
    }

    fn name(&self) -> String {
        format!("{}/{}", self.num, self.den)
    }
}

fn reach(k: u32) -> f64 {
    if k == 2 {
        0.5
    } else {
        0.25
    }
}

fn arcs_of(d: Design, k: u32) -> Vec<Arc> {
    let mut out: Vec<Arc> = Point::all(d)
        .into_iter()
        .filter(|p| p.eta < reach(k) - 1e-9)
        .map(|p| {
            let weight = p.weight(k);
            let adic = if weight.abs() < 1e-12 {
                Cx::ZERO
            } else {
                adic(d.b, &p, k, emax(k)).0
            };
            Arc {
                num: p.num,
                den: p.den,
                eta: p.eta,
                weight,
                adic,
            }
        })
        .collect();
    out.sort_by(|x, y| {
        x.eta
            .total_cmp(&y.eta)
            .then(x.den.cmp(&y.den))
            .then(x.num.cmp(&y.num))
    });
    out
}

fn all_arcs(k: u32) -> Vec<Vec<Arc>> {
    let file = data_dir().join(format!("arcs-{}-{}.txt", k, emax(k)));
    if let Ok(text) = fs::read_to_string(&file) {
        let mut out = vec![Vec::new(); DESIGNS];
        for line in text.lines() {
            let v: Vec<&str> = line.split_whitespace().collect();
            let f = |j: usize| v[j].parse::<f64>().unwrap();
            out[v[0].parse::<usize>().unwrap()].push(Arc {
                num: v[1].parse().unwrap(),
                den: v[2].parse().unwrap(),
                eta: f(3),
                weight: Cx { re: f(4), im: f(5) },
                adic: Cx { re: f(6), im: f(7) },
            });
        }
        return out;
    }
    let clock = Instant::now();
    let all = Design::all();
    let next = AtomicUsize::new(0);
    let found = Mutex::new(vec![Vec::new(); DESIGNS]);
    thread::scope(|sc| {
        for _ in 0..workers() {
            sc.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= DESIGNS {
                    break;
                }
                let j = DESIGNS - 1 - i;
                let a = arcs_of(all[j], k);
                found.lock().unwrap()[j] = a;
            });
        }
    });
    let out = found.into_inner().unwrap();
    let mut text = String::new();
    for (i, arcs) in out.iter().enumerate() {
        for a in arcs {
            text.push_str(&format!(
                "{} {} {} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}\n",
                i, a.num, a.den, a.eta, a.weight.re, a.weight.im, a.adic.re, a.adic.im
            ));
        }
    }
    fs::create_dir_all(data_dir()).unwrap();
    fs::write(&file, text).unwrap();
    eprintln!("arcs along n^{}: {:.1} s", k, clock.elapsed().as_secs_f64());
    out
}

// WEIGHTS

const MOM: usize = 44;

fn moments(d: Design, num: u64, den: u64, levels: usize) -> Vec<[Cx; MOM]> {
    let b = d.b;
    let bf = b as f64;
    let mut binom = [[0.0f64; MOM]; MOM];
    for i in 0..MOM {
        binom[i][0] = 1.0;
        for l in 1..=i {
            binom[i][l] = binom[i - 1][l - 1] + if l < i { binom[i - 1][l] } else { 0.0 };
        }
    }
    let inv: Vec<f64> = (0..MOM).map(|l| bf.powi(-(l as i32))).collect();
    let mut mu = vec![[Cx::ZERO; MOM]; levels + 1];
    for (i, m) in mu[0].iter_mut().enumerate() {
        *m = Cx::real((-0.5f64).powi(i as i32));
    }
    let mut shift = num % den;
    for j in 1..=levels {
        let mut next = [Cx::ZERO; MOM];
        for dg in 0..b {
            let ph = Cx::turn(dg * shift % den, den).scale(d.sign(dg) as f64);
            let a = (2.0 * dg as f64 + 1.0 - bf) / (2.0 * bf);
            let pa: Vec<f64> = (0..MOM).map(|m| a.powi(m as i32)).collect();
            for (i, slot) in next.iter_mut().enumerate() {
                let mut acc = Cx::ZERO;
                for l in 0..=i {
                    acc = acc.add(mu[j - 1][l].scale(binom[i][l] * inv[l] * pa[i - l]));
                }
                *slot = slot.add(ph.mul(acc));
            }
        }
        mu[j] = next;
        shift = shift * b % den;
    }
    mu
}

fn digits(mut m: u128, b: u128) -> Vec<u64> {
    let mut out = Vec::new();
    while m > 0 {
        out.push((m % b) as u64);
        m /= b;
    }
    out
}

fn weighted(d: Design, num: u64, den: u64, k: u32, tops: &[u128]) -> Vec<Cx> {
    let b = d.b as u128;
    let s = 1.0 / k as f64 - 1.0;
    let levels = tops.iter().map(|&m| digits(m, b).len()).max().unwrap_or(1);
    let mu = moments(d, num, den, levels);
    let mut gen = [0.0f64; MOM];
    gen[0] = 1.0;
    for i in 1..MOM {
        gen[i] = gen[i - 1] * (s - (i - 1) as f64) / i as f64;
    }
    let mut bj = vec![1u64; levels + 1];
    for j in 1..=levels {
        bj[j] = bj[j - 1] * d.b % den;
    }
    let block = |c: u128, j: usize| -> Cx {
        let cf = c as f64 + 0.5;
        let mut acc = Cx::ZERO;
        let mut pw = 1.0;
        for i in 0..MOM {
            acc = acc.add(mu[j][i].scale(gen[i] * pw));
            pw /= cf;
        }
        let ph = Cx::turn((c % den as u128) as u64 * bj[j] % den * num % den, den);
        let sign = 1.0 - 2.0 * d.parity(c) as f64;
        acc.mul(ph)
            .scale(sign / k as f64 * (d.b as f64).powf(j as f64 * s) * cf.powf(s))
    };
    let mut low = vec![Cx::ZERO; levels + 1];
    for j in 1..=levels {
        let mut add = Cx::ZERO;
        for c in 1..d.b as u128 {
            add = add.add(block(c, j - 1));
        }
        low[j] = low[j - 1].add(add);
    }
    tops.iter()
        .map(|&m| {
            let dg = digits(m, b);
            let mut out = Cx::ZERO;
            let mut prefix = 0u128;
            for j in (0..dg.len()).rev() {
                for x in 0..dg[j] as u128 {
                    let c = prefix * b + x;
                    out = out.add(if c == 0 { low[j] } else { block(c, j) });
                }
                prefix = prefix * b + dg[j] as u128;
            }
            out
        })
        .collect()
}

fn model(d: Design, k: u32, arcs: &[Arc], ns: &[u64]) -> Vec<(f64, Vec<f64>)> {
    let tops: Vec<u128> = ns.iter().map(|&n| (n as u128).pow(k)).collect();
    let parts: Vec<Vec<f64>> = arcs
        .iter()
        .map(|a| {
            weighted(d, a.num, a.den, k, &tops)
                .into_iter()
                .map(|u| a.amp().mul(u).re)
                .collect()
        })
        .collect();
    ns.iter()
        .enumerate()
        .map(|(i, _)| {
            let each: Vec<f64> = parts.iter().map(|p| p[i]).collect();
            (1.0 + each.iter().sum::<f64>(), each)
        })
        .collect()
}

// WALKER

const W: usize = 8;
const BASES: [u64; 3] = [3, 5, 6];
const LIMB: [u64; 3] = [531_441, 390_625, 279_936];
const EVEN: u128 = 0x5555_5555_5555_5555_5555_5555_5555_5555;

struct Tables {
    t: [Vec<u32>; 3],
    four: [u64; 8],
}

fn table(b: u64, size: u64) -> Vec<u32> {
    let mut digit = vec![0u32; b as usize];
    for mask in 1..1u32 << (b - 1) {
        let d = Design { b, mask };
        for r in 1..b {
            if d.has(r) {
                digit[r as usize] |= 1 << (mask - 1);
            }
        }
    }
    let mut t = vec![0u32; size as usize];
    for c in 1..size as usize {
        t[c] = t[c / b as usize] ^ digit[c % b as usize];
    }
    t
}

impl Tables {
    fn new() -> Tables {
        let mut four = [0u64; 8];
        for (c, f) in four.iter_mut().enumerate() {
            for mask in 1..8u64 {
                *f |= ((mask & c as u64).count_ones() as u64 & 1) << (mask - 1);
            }
        }
        Tables {
            t: [table(3, LIMB[0]), table(5, LIMB[1]), table(6, LIMB[2])],
            four,
        }
    }
}

#[derive(Clone)]
struct Walker {
    k: usize,
    wide: [u128; 4],
    limbs: [[[u64; W]; 4]; 3],
    len: [[usize; 4]; 3],
}

impl Walker {
    fn new(n: u64, k: usize) -> Walker {
        let mut d: Vec<u128> = (0..=k as u128)
            .map(|i| (n as u128 + i).pow(k as u32))
            .collect();
        for i in 1..=k {
            for j in (i..=k).rev() {
                d[j] -= d[j - 1];
            }
        }
        let mut w = Walker {
            k,
            wide: [0; 4],
            limbs: [[[0; W]; 4]; 3],
            len: [[0; 4]; 3],
        };
        for (i, &v0) in d.iter().enumerate() {
            w.wide[i] = v0;
            for (bi, &p) in LIMB.iter().enumerate() {
                let mut v = v0;
                for l in 0..W {
                    w.limbs[bi][i][l] = (v % p as u128) as u64;
                    v /= p as u128;
                    if w.limbs[bi][i][l] > 0 {
                        w.len[bi][i] = l + 1;
                    }
                }
                assert_eq!(v, 0);
            }
        }
        w
    }

    fn step(&mut self) {
        for i in 0..self.k {
            self.wide[i] += self.wide[i + 1];
            for bi in 0..3 {
                let p = LIMB[bi];
                let lb = self.len[bi][i + 1];
                let (lo, hi) = self.limbs[bi].split_at_mut(i + 1);
                let (a, b) = (&mut lo[i], &hi[0]);
                let mut carry = 0;
                let mut l = 0;
                while l < lb || carry > 0 {
                    let s = a[l] + if l < lb { b[l] } else { 0 } + carry;
                    if s >= p {
                        a[l] = s - p;
                        carry = 1;
                    } else {
                        a[l] = s;
                        carry = 0;
                    }
                    l += 1;
                }
                if l > self.len[bi][i] {
                    self.len[bi][i] = l;
                }
            }
        }
    }

    fn word(&self, t: &Tables) -> u64 {
        let v = self.wide[0];
        let (lo, hi) = (v & EVEN, (v >> 1) & EVEN);
        let c = ((lo & !hi).count_ones() & 1)
            | ((hi & !lo).count_ones() & 1) << 1
            | ((hi & lo).count_ones() & 1) << 2;
        let mut w = (v.count_ones() & 1) as u64 | t.four[c as usize] << offset(4);
        for (bi, &b) in BASES.iter().enumerate() {
            let mut a = 0u32;
            for &limb in &self.limbs[bi][0][..self.len[bi][0]] {
                a ^= t.t[bi][limb as usize];
            }
            w |= (a as u64) << offset(b);
        }
        w
    }
}

// CENSUS

struct Slices {
    planes: [u64; 16],
    adds: u32,
}

impl Slices {
    fn new() -> Slices {
        Slices {
            planes: [0; 16],
            adds: 0,
        }
    }

    fn add(&mut self, w: u64, total: &mut [u64; 64]) {
        let mut carry = w;
        for p in self.planes.iter_mut() {
            let t = *p & carry;
            *p ^= carry;
            carry = t;
            if carry == 0 {
                break;
            }
        }
        self.adds += 1;
        if self.adds == u16::MAX as u32 {
            self.flush(total);
        }
    }

    fn flush(&mut self, total: &mut [u64; 64]) {
        for (k, p) in self.planes.iter().enumerate() {
            for (i, t) in total.iter_mut().enumerate() {
                *t += ((p >> i) & 1) << k;
            }
        }
        self.planes = [0; 16];
        self.adds = 0;
    }
}

#[derive(Clone)]
struct Block {
    count: u64,
    ones: [[u64; 64]; LAGS],
    grams: Vec<[u64; 256]>,
}

impl Block {
    fn new(picks: usize) -> Block {
        Block {
            count: 0,
            ones: [[0; 64]; LAGS],
            grams: vec![[0; 256]; picks],
        }
    }

    fn merge(&mut self, o: &Block) {
        self.count += o.count;
        for (a, b) in self.ones.iter_mut().zip(&o.ones) {
            for (x, y) in a.iter_mut().zip(b) {
                *x += y;
            }
        }
        for (a, b) in self.grams.iter_mut().zip(&o.grams) {
            for (x, y) in a.iter_mut().zip(b) {
                *x += y;
            }
        }
    }

    fn sum(&self, lag: usize, s: usize) -> i64 {
        self.count as i64 - 2 * self.ones[lag][s] as i64
    }
}

fn workers() -> usize {
    thread::available_parallelism()
        .map_or(4, |n| n.get())
        .min(8)
}

fn chunk(s: u64, e: u64, k: usize, salt: u64, picks: &[usize], t: &Tables) -> Block {
    let mut out = Block::new(picks.len());
    let mut walker = Walker::new(s, k);
    let mut ring = [0u64; 8];
    let mut slices: Vec<Slices> = (0..LAGS).map(|_| Slices::new()).collect();
    let mut regs = vec![0usize; picks.len()];
    for n in s..e + 7 {
        let w = walker.word(t) | coin(n, salt) << COIN;
        walker.step();
        if n < e {
            slices[0].add(w, &mut out.ones[0]);
        }
        for h in 1..LAGS as u64 {
            if n >= s + h && n - h < e {
                slices[h as usize].add(w ^ ring[((n - h) & 7) as usize], &mut out.ones[h as usize]);
            }
        }
        ring[(n & 7) as usize] = w;
        for (r, &p) in regs.iter_mut().zip(picks) {
            *r = (*r >> 1) | (((w >> p) & 1) as usize) << 7;
        }
        if n >= s + 7 {
            for (g, &r) in out.grams.iter_mut().zip(&regs) {
                g[r] += 1;
            }
        }
    }
    for (h, sl) in slices.iter_mut().enumerate() {
        sl.flush(&mut out.ones[h]);
    }
    out.count = e - s;
    out
}

fn tally(k: usize, top: u32, salt: u64, picks: &[usize]) -> Vec<Block> {
    let t = Tables::new();
    let size: u64 = 1 << 16;
    let mut jobs: Vec<(u64, u64, usize)> = vec![(0, 1, 0)];
    for b in 1..=top.min(16) {
        jobs.push((1 << (b - 1), 1 << b, b as usize));
    }
    let mut s = size;
    while s < 1u64 << top {
        jobs.push((s, s + size, 64 - s.leading_zeros() as usize));
        s += size;
    }
    let total = Mutex::new(vec![Block::new(picks.len()); top as usize + 1]);
    let next = AtomicUsize::new(0);
    thread::scope(|sc| {
        for _ in 0..workers() {
            sc.spawn(|| {
                let mut mine = vec![Block::new(picks.len()); top as usize + 1];
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= jobs.len() {
                        break;
                    }
                    let (s0, e0, blk) = jobs[i];
                    mine[blk].merge(&chunk(s0, e0, k, salt, picks, &t));
                }
                let mut all = total.lock().unwrap();
                for (a, m) in all.iter_mut().zip(&mine) {
                    a.merge(m);
                }
            });
        }
    });
    let blocks = total.into_inner().unwrap();
    let mut cum = Vec::new();
    let mut run = Block::new(picks.len());
    for b in &blocks {
        run.merge(b);
        cum.push(run.clone());
    }
    cum
}

fn data_dir() -> PathBuf {
    let here = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let repo = here.ancestors().nth(5).unwrap().to_path_buf();
    repo.join("data").join(here.strip_prefix(&repo).unwrap())
}

fn cached(k: usize, top: u32, picks: &[usize]) -> Vec<Block> {
    let file = data_dir().join(format!("power-{}-{}.txt", k, top));
    if let Ok(text) = fs::read_to_string(&file) {
        let nums: Vec<u64> = text
            .split_whitespace()
            .map(|x| x.parse().unwrap())
            .collect();
        let per = 1 + LAGS * 64 + picks.len() * 256;
        if nums.len() == per * (top as usize + 1) {
            return nums
                .chunks(per)
                .map(|c| {
                    let mut b = Block::new(picks.len());
                    b.count = c[0];
                    for h in 0..LAGS {
                        b.ones[h].copy_from_slice(&c[1 + 64 * h..1 + 64 * (h + 1)]);
                    }
                    for (g, gram) in b.grams.iter_mut().enumerate() {
                        let at = 1 + LAGS * 64 + 256 * g;
                        gram.copy_from_slice(&c[at..at + 256]);
                    }
                    b
                })
                .collect();
        }
    }
    let clock = Instant::now();
    let cum = tally(k, top, k as u64, picks);
    eprintln!(
        "tally along n^{} to 2^{}: {:.1} s",
        k,
        top,
        clock.elapsed().as_secs_f64()
    );
    fs::create_dir_all(data_dir()).unwrap();
    let mut text = String::new();
    for b in &cum {
        text.push_str(&b.count.to_string());
        for h in 0..LAGS {
            for x in b.ones[h] {
                text.push(' ');
                text.push_str(&x.to_string());
            }
        }
        for g in &b.grams {
            for x in g {
                text.push(' ');
                text.push_str(&x.to_string());
            }
        }
        text.push('\n');
    }
    fs::write(&file, text).unwrap();
    cum
}

fn walsh(c: &[u64; 256]) -> [i64; 256] {
    let mut t = c.map(|x| x as i64);
    let mut h = 1;
    while h < 256 {
        for i in 0..256 {
            if i & h == 0 {
                let (a, b) = (t[i], t[i | h]);
                t[i] = a + b;
                t[i | h] = a - b;
            }
        }
        h <<= 1;
    }
    t
}

fn chi(c: &[u64; 256]) -> f64 {
    let n: u64 = c.iter().sum();
    let t = walsh(c);
    (1..256).map(|s| (t[s] as f64).powi(2)).sum::<f64>() / n as f64
}

fn longer(d: Design) -> Vec<(u32, u64, u64, f64)> {
    let b = d.b;
    let phi2 = |x: f64| -> f64 {
        let (mut re, mut im) = (0.0f64, 0.0f64);
        for t in 0..b {
            let a = std::f64::consts::TAU * t as f64 * x;
            re += d.sign(t) as f64 * a.cos();
            im += d.sign(t) as f64 * a.sin();
        }
        (re * re + im * im) / (b * b) as f64
    };
    let mut out = Vec::new();
    let mut k = 3u32;
    while b.pow(k) <= 1 << 14 {
        let m = b.pow(k) - 1;
        let g: Vec<f64> = (0..m)
            .map(|c| phi2(c as f64 / m as f64).max(1e-300).ln())
            .collect();
        'c: for c in 0..m {
            let mut x = c;
            let mut s = 0.0;
            for j in 0..k {
                if x < c || (j > 0 && x == c) {
                    continue 'c;
                }
                s += g[x as usize];
                x = x * b % m;
            }
            let eta = -s / (2.0 * k as f64 * (b as f64).ln());
            if eta < 0.25 {
                out.push((k, c / gcd(c, m), m / gcd(c, m), eta));
            }
        }
        k += 1;
    }
    out
}

// VERBS

fn checkpoints(k: u32, top: u32) -> Vec<u32> {
    let lo = if k == 2 { 12 } else { 10 };
    (lo..=top).filter(|j| (top - j) % 2 == 0).collect()
}

fn coefficients() {
    let all = Design::all();
    let tag = tags();
    for k in [2u32, 3] {
        let table = all_arcs(k);
        println!("# COEFFICIENTS ALONG n^{}", k);
        println!();
        println!(
            "the model S(N) = 1 + Re sum_gamma A(gamma) U(gamma; N^{k}), U(gamma; M) = sum_(1<=m<M) x(m) e(m gamma) m^(1/{k} - 1)/{k}; gamma = c/q runs over the points c/(b^2 - 1), each fixed by B = b^2, with eta(gamma) = 1 - log_B abs(sigma) < {reach}, sigma = sum_(t<B) x(t) e(t gamma); A = (G(-c, q)/q) K, G the degree-{k} Gauss sum; K = lim_J E_(n in Z_B) prod_(i<J) y(digit_i(n^{k}))/phi(0), y(t) = x(t) e(t gamma), phi(0) = sigma/B, in closed form for B = p^t and by the 2-adic and 3-adic valuations of n for B = 36, the e-series cut at e = {e}; listed: the arcs with eta < {lim}, best first; copies omitted",
            k = k,
            reach = reach(k),
            e = emax(k),
            lim = if k == 2 { "1/4" } else { "1/4 and a nonzero Gauss sum" }
        );
        println!();
        println!(
            "{:<8} {:<6} {:>8} {:>8} {:>7} {:>9} {:>8} {:>9}",
            "design", "class", "gamma", "eta", "|G|/q", "|K|", "arg/2pi", "|A|"
        );
        for (i, d) in all.iter().enumerate() {
            if tag[i] == "copy" {
                continue;
            }
            let shown: Vec<&Arc> = table[i]
                .iter()
                .filter(|a| a.eta < 0.25 && (k == 2 || a.weight.abs() > 1e-12))
                .collect();
            if shown.is_empty() {
                let best = table[i]
                    .iter()
                    .find(|a| a.weight.abs() > 1e-12)
                    .map(|a| format!("{} {:.4}", a.name(), a.eta))
                    .unwrap_or("-".into());
                println!(
                    "{:<8} {:<6} none; best with G != 0 below {}: {}",
                    d.name(),
                    tag[i],
                    reach(k),
                    best
                );
                continue;
            }
            for (j, a) in shown.iter().enumerate() {
                println!(
                    "{:<8} {:<6} {:>8} {:>8.4} {:>7.4} {:>9.4} {:>8.4} {:>9.4}",
                    if j == 0 { d.name() } else { String::new() },
                    if j == 0 { tag[i] } else { "" },
                    a.name(),
                    a.eta,
                    a.weight.abs(),
                    a.adic.abs(),
                    a.adic.turns(),
                    a.amp().abs()
                );
            }
        }
        println!();
        if k == 2 {
            let mut pairs = Vec::new();
            for (i, d) in all.iter().enumerate() {
                let m = d.mirror();
                let j = m.index();
                if m.mask <= d.mask || tag[i] == "copy" || tag[j] == "copy" {
                    continue;
                }
                let top = |t: &[Arc]| {
                    t.iter()
                        .filter(|a| a.eta < 0.25)
                        .fold(0f64, |x, a| x.max(a.amp().abs()))
                };
                let (x, y) = (top(&table[i]), top(&table[j]));
                if x.max(y) > 0.0 {
                    pairs.push(format!("{}~{} {:.4}/{:.4}", d.name(), m.name(), x, y));
                }
            }
            println!(
                "mirror pairs with an arc below 1/4, largest |A| below 1/4 of each: {}",
                pairs.join(", ")
            );
            println!();
        }
    }
    println!("# LONGER ORBITS");
    println!();
    println!("orbits of x -> b x of least period 3 and more, searched while b^k <= 2^14 (period 14 at b = 2, 5 at b = 5 and 6), with eta < 1/4 in floating point, one line per design that has any: period, denominator, number of points, least eta, largest quadratic Gauss weight; the model above keeps the periods 1 and 2, and orbits shadowing the best one continue past the search");
    for (i, d) in all.iter().enumerate() {
        if tag[i] == "copy" || tag[i] == "folded" {
            continue;
        }
        let found = longer(*d);
        if found.is_empty() {
            continue;
        }
        let mut groups: Vec<(u32, u64, usize, f64, f64)> = Vec::new();
        for (k, c, q, e) in found {
            let w = gauss(2, (q - c % q) % q, q).abs() / q as f64;
            match groups.iter_mut().find(|g| g.0 == k && g.1 == q) {
                Some(g) => {
                    g.2 += 1;
                    g.3 = g.3.min(e);
                    g.4 = g.4.max(w);
                }
                None => groups.push((k, q, 1, e, w)),
            }
        }
        let text: Vec<String> = groups
            .iter()
            .map(|g| {
                format!(
                    "period {} q {}: {} points, eta {:.4}, abs(G)/q <= {:.4}",
                    g.0,
                    g.1,
                    g.2 * g.0 as usize,
                    g.3,
                    g.4
                )
            })
            .collect();
        println!("{:<8} {}", d.name(), text.join("; "));
    }
    println!();
    println!("# CHECKS");
    println!();
    println!("K_J = sum_(n mod B^J) prod_(i<J) y(digit_i(n^2 mod B^J)) / sigma^J, the tree cut at depth J, against the limit K and its e-terms");
    for (name, c, js) in [
        ("2:1", 1usize, vec![2u32, 4, 6, 8, 10]),
        ("6:5", 0, vec![1, 2, 3, 4, 5]),
        ("6:12345", 0, vec![1, 2, 3, 4, 5]),
        ("6:14", 14, vec![1, 2, 3, 4]),
        ("6:23", 28, vec![1, 2, 3, 4]),
    ] {
        let d = Design::find(name);
        let pts = Point::all(d);
        let p = &pts[c];
        let (lim, terms) = adic(d.b, p, 2, emax(2));
        let row: Vec<String> = js
            .iter()
            .map(|&j| {
                let v = local(p, 2, j);
                format!("J={} {:.5}{:+.5}i", j, v.re, v.im)
            })
            .collect();
        println!(
            "{} at {}: {}; limit {:.5}{:+.5}i",
            name,
            p.name(),
            row.join(", "),
            lim.re,
            lim.im
        );
        if !terms.is_empty() {
            let t: Vec<String> = terms.iter().map(|t| format!("{:.1e}", t.abs())).collect();
            println!("  e-terms, abs: {}", t.join(" "));
        }
    }
    let tm = Design::find("2:1");
    let p = &Point::all(tm)[1];
    let w = Cx::turn(1, 3);
    let k0 = adic(2, p, 2, 0).0;
    let closed = Cx::real(4.0 / 3.0).mul(Cx::ONE.sub(w));
    let amp = p.weight(2).mul(k0);
    println!(
        "2:1 at 1/3: K = {:.12}{:+.12}i against (4/3)(1 - e(1/3)) = {:.12}{:+.12}i; A = G(-1,3)/3 . K = {:.12}{:+.12}i against (4/3) e(-1/3) = {:.12}{:+.12}i",
        k0.re, k0.im, closed.re, closed.im, amp.re, amp.im, Cx::turn(2, 3).scale(4.0 / 3.0).re, Cx::turn(2, 3).scale(4.0 / 3.0).im
    );
    let tm = Design::find("2:1");
    let ms: Vec<u64> = vec![1 << 12, 1 << 16, 1 << 20, 1 << 24];
    let us = weighted(
        tm,
        1,
        3,
        2,
        &ms.iter().map(|&m| m as u128).collect::<Vec<_>>(),
    );
    for j in [1u32, 4, 8] {
        let q = 4u64.pow(j);
        let mut rho = vec![0f64; q as usize];
        for n in 0..q {
            rho[(n * n % q) as usize] += 1.0;
        }
        let kj = local(&Point::all(tm)[1], 2, j);
        let mut acc = Cx::ZERO;
        let mut row = Vec::new();
        let mut at = 0;
        for m in 1..*ms.last().unwrap() {
            let y = Cx::turn(m % 3, 3).scale(1.0 - 2.0 * (m.count_ones() & 1) as f64);
            acc = acc.add(y.scale(rho[(m % q) as usize] / (2.0 * (m as f64).sqrt())));
            if m + 1 == ms[at] {
                let gap = acc.sub(kj.mul(us[at]));
                row.push(format!(
                    "M=2^{} gap {:.3} ({:.1e} of the tree)",
                    63 - ms[at].leading_zeros(),
                    gap.abs(),
                    gap.abs() / acc.abs()
                ));
                at += 1;
            }
        }
        println!(
            "2:1 at 1/3, tree cut at depth J={} summed over m < M against K_J U: {}",
            j,
            row.join(", ")
        );
    }
    println!();
    println!("# VANISHING");
    println!();
    for k in [2u32, 3] {
        let table = all_arcs(k);
        let mut zero = Vec::new();
        for (i, d) in all.iter().enumerate() {
            if tag[i] == "copy" {
                continue;
            }
            for a in table[i]
                .iter()
                .filter(|a| a.eta < 0.25 && a.weight.abs() < 1e-12)
            {
                zero.push(format!("{} {} ({:.4})", d.name(), a.name(), a.eta));
            }
        }
        println!(
            "along n^{}: arcs with eta < 1/4 and G(-c, q) = 0: {}",
            k,
            zero.join(", ")
        );
    }
    println!();
}

fn compare(k: u32, top: u32, cut: f64) {
    let all = Design::all();
    let tag = tags();
    let picks: Vec<usize> = if k == 2 {
        vec![]
    } else {
        [
            "2:1", "4:13", "5:124", "6:5", "6:12345", "6:14", "6:135", "6:345",
        ]
        .iter()
        .map(|n| Design::find(n).index())
        .chain([COIN])
        .collect()
    };
    let cum = cached(k as usize, top, &picks);
    let ks = checkpoints(k, top);
    let ns: Vec<u64> = ks.iter().map(|&j| 1u64 << j).collect();
    println!(
        "# {}",
        if k == 2 {
            "MODEL ALONG THE SQUARES"
        } else {
            "MODEL ALONG THE CUBES"
        }
    );
    println!();
    println!(
        "observed z = S(N)/sqrt(N), S(N) = sum_(n<N) x(n^{k}), against the model 1 + Re sum_gamma A(gamma) U(gamma; N^{k}) over sqrt(N), every arc with eta < {cut}; per design two rows, obs and model, at N = 2^j; lit: abs(z) > 4 at some checkpoint",
        k = k,
        cut = cut
    );
    println!();
    print!("{:<8} {:<6} {:>7} {:<5}", "design", "class", "eta", "row");
    for j in &ks {
        print!(" {:>7}", format!("2^{}", j));
    }
    println!();
    let mut resid = Vec::new();
    let mut table = Vec::new();
    let table_arcs = all_arcs(k);
    for (i, d) in all.iter().enumerate() {
        if tag[i] == "copy" {
            continue;
        }
        let arcs: Vec<Arc> = table_arcs[i]
            .iter()
            .filter(|a| a.eta < cut)
            .cloned()
            .collect();
        let eta = table_arcs[i].first().map_or(f64::INFINITY, |a| a.eta);
        let pred = model(*d, k, &arcs, &ns);
        let obs: Vec<f64> = ks
            .iter()
            .map(|&j| cum[j as usize].sum(0, i) as f64 / ((1u64 << j) as f64).sqrt())
            .collect();
        let mz: Vec<f64> = pred
            .iter()
            .zip(&ns)
            .map(|(p, &n)| p.0 / (n as f64).sqrt())
            .collect();
        for (o, m) in obs.iter().zip(&mz) {
            if tag[i] != "folded" {
                resid.push(o - m);
            }
        }
        for (row, v) in [("obs", &obs), ("model", &mz)] {
            print!(
                "{:<8} {:<6} {:>7} {:<5}",
                if row == "obs" {
                    d.name()
                } else {
                    String::new()
                },
                if row == "obs" { tag[i] } else { "" },
                if row == "obs" {
                    format!("{:.4}", eta)
                } else {
                    String::new()
                },
                row
            );
            for x in v.iter() {
                print!(" {:>7.2}", x);
            }
            println!();
        }
        table.push((d.name(), tag[i], eta, obs, mz));
    }
    let c = cum[top as usize].sum(0, COIN) as f64 / ((1u64 << top) as f64).sqrt();
    println!(
        "{:<8} {:<6} {:>7} {:<5} z at 2^{} = {:.2}",
        "coin", "", "", "obs", top, c
    );
    println!();
    let n = resid.len() as f64;
    let rms = (resid.iter().map(|r| r * r).sum::<f64>() / n).sqrt();
    let worst = resid.iter().fold(0f64, |m, r| m.max(r.abs()));
    println!(
        "residual obs - model over {} cells of the non-folded designs: rms {:.3}, largest {:.2}",
        resid.len(),
        rms,
        worst
    );
    let split = |want: bool| {
        let r: Vec<f64> = table
            .iter()
            .filter(|r| r.1 != "folded" && r.4.iter().any(|z| z.abs() > 4.0) == want)
            .flat_map(|r| r.3.iter().zip(&r.4).map(|(o, m)| o - m).collect::<Vec<_>>())
            .collect();
        (
            (r.iter().map(|x| x * x).sum::<f64>() / r.len().max(1) as f64).sqrt(),
            r.len(),
        )
    };
    let (lit_rms, lit_n) = split(true);
    let (dark_rms, dark_n) = split(false);
    println!("residual rms over the designs the model lights {:.3} ({} cells), over the rest {:.3} ({} cells)", lit_rms, lit_n, dark_rms, dark_n);
    let ratios: Vec<String> = table
        .iter()
        .filter(|r| r.4.last().map_or(false, |z| z.abs() > 8.0))
        .map(|r| format!("{} {:.2}", r.0, r.3.last().unwrap() / r.4.last().unwrap()))
        .collect();
    println!(
        "obs/model at 2^{} where abs(model) > 8: {}",
        top,
        ratios.join(", ")
    );
    for th in [0.5, 1.0, 1.5, 2.0] {
        let cells: Vec<(f64, f64, String)> = table
            .iter()
            .filter(|r| r.1 != "folded")
            .flat_map(|r| {
                r.3.iter()
                    .zip(&r.4)
                    .filter(|(_, m)| m.abs() > th)
                    .map(|(o, m)| (*o, *m, r.0.clone()))
                    .collect::<Vec<_>>()
            })
            .collect();
        let mut names: Vec<&String> = cells.iter().map(|c| &c.2).collect();
        names.sort();
        names.dedup();
        let rel = (cells
            .iter()
            .map(|(o, m, _)| (o / m - 1.0).powi(2))
            .sum::<f64>()
            / cells.len().max(1) as f64)
            .sqrt();
        println!(
            "cells with abs(model) > {}: {} in {} designs, obs of the same sign {}, within a factor 2 {}, rms of obs/model - 1 {:.3}",
            th,
            cells.len(),
            names.len(),
            cells.iter().filter(|(o, m, _)| o * m > 0.0).count(),
            cells.iter().filter(|(o, m, _)| o / m > 0.5 && o / m < 2.0).count(),
            rel
        );
    }
    let near: Vec<String> = table
        .iter()
        .filter(|r| r.1 != "folded")
        .flat_map(|r| {
            let mo = r.3.iter().fold(0f64, |m, z| m.max(z.abs()));
            let mm = r.4.iter().fold(0f64, |m, z| m.max(z.abs()));
            let mut v = Vec::new();
            if (mo - 4.0).abs() < 0.1 {
                v.push(format!("{} obs {:.4}", r.0, mo));
            }
            if (mm - 4.0).abs() < 0.1 {
                v.push(format!("{} model {:.4}", r.0, mm));
            }
            v
        })
        .collect();
    println!(
        "designs within 0.1 of the lit cut 4: {}",
        if near.is_empty() {
            "none".to_string()
        } else {
            near.join(", ")
        }
    );
    if k == 3 {
        for (i, _) in all.iter().enumerate() {
            if tag[i] == "copy" || tag[i] == "folded" {
                continue;
            }
            for h in 1..LAGS {
                let zs: Vec<f64> = ks
                    .iter()
                    .map(|&j| cum[j as usize].sum(h, i) as f64 / ((1u64 << j) as f64).sqrt())
                    .collect();
                if zs.iter().any(|z| z.abs() > 3.0) {
                    println!(
                        "lag row {} h={}: {}",
                        label(i),
                        h,
                        zs.iter()
                            .map(|z| format!("{:.2}", z))
                            .collect::<Vec<_>>()
                            .join(" ")
                    );
                }
            }
        }
    }
    let mut cells = [[0usize; 2]; 2];
    let mut odd = Vec::new();
    for (name, t, eta, obs, mz) in &table {
        if *t == "folded" {
            continue;
        }
        let lo = obs.iter().any(|z| z.abs() > 4.0) as usize;
        let lm = mz.iter().any(|z| z.abs() > 4.0) as usize;
        cells[lo][lm] += 1;
        if lo != lm {
            let mo = obs.iter().fold(0f64, |m, z| m.max(z.abs()));
            let mm = mz.iter().fold(0f64, |m, z| m.max(z.abs()));
            odd.push(format!(
                "{} (eta {:.4}, max obs {:.2}, max model {:.2})",
                name, eta, mo, mm
            ));
        }
    }
    println!(
        "lit by both {}, by neither {}, observed only {}, model only {}: {}",
        cells[1][1],
        cells[0][0],
        cells[1][0],
        cells[0][1],
        if odd.is_empty() {
            "none".to_string()
        } else {
            odd.join("; ")
        }
    );
    let mut big: Vec<(f64, String)> = table
        .iter()
        .filter(|r| r.1 != "folded")
        .map(|r| {
            let m =
                r.3.iter()
                    .zip(&r.4)
                    .fold(0f64, |m, (o, p)| m.max((o - p).abs()));
            (m, r.0.clone())
        })
        .collect();
    big.sort_by(|a, b| b.0.total_cmp(&a.0));
    println!(
        "largest residuals by design: {}",
        big.iter()
            .take(6)
            .map(|b| format!("{} {:.2}", b.1, b.0))
            .collect::<Vec<_>>()
            .join(", ")
    );
    let lagmax = (1..LAGS)
        .flat_map(|h| {
            (0..STREAMS)
                .filter(|&s| s == COIN || (tag[s] != "copy" && tag[s] != "folded"))
                .map(move |s| (h, s))
        })
        .map(|(h, s)| {
            (
                cum[top as usize].sum(h, s) as f64 / ((1u64 << top) as f64).sqrt(),
                h,
                s,
            )
        })
        .fold(
            (0f64, 0, 0),
            |m, x| if x.0.abs() > m.0.abs() { x } else { m },
        );
    println!(
        "lags h <= 7 at 2^{}: largest abs(z) over the non-folded designs and the coin {:.2} at {} h={}",
        top,
        lagmax.0,
        label(lagmax.2),
        lagmax.1
    );
    if !picks.is_empty() {
        println!();
        println!(
            "8-blocks along n^{} at 2^{}: X2 (iid mean 255, sd 38.7) and the lags",
            k, top
        );
        for (g, &p) in picks.iter().enumerate() {
            let lags: Vec<String> = (1..LAGS)
                .map(|h| {
                    format!(
                        "{:.2}",
                        cum[top as usize].sum(h, p) as f64 / ((1u64 << top) as f64).sqrt()
                    )
                })
                .collect();
            println!(
                "{:<8} X2 {:>8.1}  zB {:>6.2}  zC1..7 {}",
                label(p),
                chi(&cum[top as usize].grams[g]),
                cum[top as usize].sum(0, p) as f64 / ((1u64 << top) as f64).sqrt(),
                lags.join(" ")
            );
        }
    }
    println!();
}

fn terms(top: u32, qmax: u64) {
    let ks: Vec<u32> = vec![12, 18, 24, 30, top];
    let ns: Vec<u64> = ks.iter().map(|&j| 1u64 << j).collect();
    let tops: Vec<u128> = ns.iter().map(|&n| (n as u128) * (n as u128)).collect();
    let cum = cached(2, top, &[]);
    let table_arcs = all_arcs(2);
    println!("# TERMS");
    println!();
    println!(
        "the q-terms of the model along the squares: T_q(N) = sum_(a mod q, (a,q)=1) (G(-a,q)/q) U(a/q; N^2) over sqrt(N), and the running sum over q, against the trees of every arc below 1/4 and the census, at N = 2^j"
    );
    for name in ["2:1", "6:5", "6:12345"] {
        let d = Design::find(name);
        let i = d.index();
        println!();
        print!("{:<8} {:>5}", name, "q");
        for j in &ks {
            print!(" {:>9} {:>9}", format!("T 2^{}", j), "sum");
        }
        println!();
        let mut run = vec![0.0f64; ks.len()];
        for q in 1..=qmax {
            let mut tq = vec![0.0f64; ks.len()];
            for a in (0..q).filter(|&a| gcd(a, q) == 1) {
                let g = gauss(2, (q - a) % q, q).scale(1.0 / q as f64);
                if g.abs() < 1e-12 {
                    continue;
                }
                for (x, u) in tq.iter_mut().zip(weighted(d, a, q, 2, &tops)) {
                    *x += g.mul(u).re;
                }
            }
            if tq.iter().all(|x| x.abs() < 1e-9) {
                continue;
            }
            print!("{:<8} {:>5}", "", q);
            for (c, (x, r)) in tq.iter().zip(run.iter_mut()).enumerate() {
                *r += x;
                let s = (ns[c] as f64).sqrt();
                print!(" {:>9.3} {:>9.3}", x / s, *r / s);
            }
            println!();
        }
        let arcs: Vec<Arc> = table_arcs[i]
            .iter()
            .filter(|a| a.eta < 0.25)
            .cloned()
            .collect();
        let pred = model(d, 2, &arcs, &ns);
        print!("{:<8} {:>5}", "", "trees");
        for (c, p) in pred.iter().enumerate() {
            print!(" {:>9} {:>9.3}", "", (p.0 - 1.0) / (ns[c] as f64).sqrt());
        }
        println!();
        print!("{:<8} {:>5}", "", "obs");
        for (c, &j) in ks.iter().enumerate() {
            print!(
                " {:>9} {:>9.3}",
                "",
                (cum[j as usize].sum(0, i) - 1) as f64 / (ns[c] as f64).sqrt()
            );
        }
        println!();
        let lead: Vec<String> = arcs
            .iter()
            .take(4)
            .map(|a| {
                format!(
                    "{} eta {:.4} A {:.4}{:+.4}i",
                    a.name(),
                    a.eta,
                    a.amp().re,
                    a.amp().im
                )
            })
            .collect();
        println!("{:<8} arcs: {}", "", lead.join("; "));
    }
    println!();
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let clock = Instant::now();
    let num = |i: usize, d: u32| args.get(i).and_then(|s| s.parse().ok()).unwrap_or(d);
    match args.get(1).map(String::as_str) {
        Some("coefficients") => coefficients(),
        Some("model") => compare(
            2,
            num(2, 32),
            args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0.25),
        ),
        Some("terms") => terms(num(2, 32), num(3, 36) as u64),
        Some("cubes") => compare(
            3,
            num(2, 32),
            args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0.25),
        ),
        _ => {
            eprintln!("verbs: coefficients, model [k], terms [k q], cubes [k]");
            std::process::exit(2);
        }
    }
    eprintln!("{:.2} s", clock.elapsed().as_secs_f64());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn walker_reads_the_digits() {
        let t = Tables::new();
        let all = Design::all();
        for k in [2usize, 3] {
            for start in [0u64, 977, 1 << 21, (1 << 33) + 12345] {
                let mut w = Walker::new(start, k);
                for n in start..start + 300 {
                    let v = (n as u128).pow(k as u32);
                    let word = w.word(&t);
                    for (i, d) in all.iter().enumerate() {
                        assert_eq!((word >> i) & 1, d.parity(v), "{} at {}^{}", d.name(), n, k);
                    }
                    w.step();
                }
            }
        }
    }

    #[test]
    fn slices_count_the_bits() {
        let mut total = [0u64; 64];
        let mut s = Slices::new();
        let mut want = [0u64; 64];
        for i in 0..200_000u64 {
            let w = mix(i);
            s.add(w, &mut total);
            for (b, x) in want.iter_mut().enumerate() {
                *x += (w >> b) & 1;
            }
        }
        s.flush(&mut total);
        assert_eq!(total, want);
    }

    #[test]
    fn weighted_sums_match_the_direct_sum() {
        for (name, num, den, k) in [
            ("2:1", 1u64, 3u64, 2u32),
            ("6:5", 0, 1, 2),
            ("6:135", 3, 7, 3),
            ("5:14", 2, 9, 2),
        ] {
            let d = Design::find(name);
            let tops = [1u128, 2, 37, 1000, 77_777, 300_001];
            let got = weighted(d, num, den, k, &tops);
            for (&m, g) in tops.iter().zip(got) {
                let mut want = Cx::ZERO;
                for x in 1..m as u64 {
                    let s = 1.0 - 2.0 * d.parity(x as u128) as f64;
                    want = want.add(
                        Cx::turn(x % den * num % den, den)
                            .scale(s * (x as f64).powf(1.0 / k as f64 - 1.0) / k as f64),
                    );
                }
                assert!(
                    g.sub(want).abs() < 1e-9 * (1.0 + want.abs()),
                    "{} {} {:?} {:?}",
                    name,
                    m,
                    g,
                    want
                );
            }
        }
    }

    #[test]
    fn closed_forms_match_the_local_sums() {
        for (name, c, j) in [
            ("2:1", 1usize, 10u32),
            ("3:2", 2, 6),
            ("4:13", 6, 5),
            ("5:12", 6, 4),
            ("5:14", 8, 4),
        ] {
            let d = Design::find(name);
            let pts = Point::all(d);
            let p = &pts[c];
            let lim = adic(d.b, p, 2, emax(2)).0;
            let near = local(p, 2, j);
            let ratio = ((d.b * d.b) as f64).sqrt() / p.sigma.abs();
            assert!(
                lim.sub(near).abs() < 3.0 * ratio.powf(j as f64) * lim.abs().max(1.0),
                "{} {:?} {:?}",
                name,
                lim,
                near
            );
        }
    }

    #[test]
    fn six_classes_match_the_enumeration() {
        let d = Design::find("6:12345");
        let pts = Point::all(d);
        for c in [0usize, 7, 15] {
            let r = pts[c].riesz();
            let mut six = Six::new(&r);
            for (s2, t2, s3, t3) in [
                (3u32, 1u64, 1u32, 1u64),
                (5, 4, 1, 1),
                (3, 1, 3, 9),
                (7, 16, 1, 1),
                (3, 1, 5, 81),
            ] {
                let digits = (s2.max(s3) + 1) / 2;
                let q = 36u64.pow(digits);
                let mut sum = Cx::ZERO;
                let mut n = 0;
                for v in 0..q {
                    if v % (1 << s2) != t2 || v % 3u64.pow(s3) != t3 {
                        continue;
                    }
                    let mut x = v;
                    let mut prod = Cx::ONE;
                    for _ in 0..digits {
                        prod = prod.mul(r[(x % 36) as usize]);
                        x /= 36;
                    }
                    sum = sum.add(prod);
                    n += 1;
                }
                let want = sum.scale(1.0 / n as f64);
                let got = six.class(u32::MAX, s2, t2, s3, t3);
                assert!(
                    got.sub(want).abs() < 1e-9 * (1.0 + want.abs()),
                    "{} {:?} {:?}",
                    c,
                    got,
                    want
                );
            }
        }
    }

    #[test]
    fn six_valuations_match_the_local_sums() {
        for (name, c) in [
            ("6:5", 0usize),
            ("6:12345", 0),
            ("6:135", 15),
            ("6:14", 14),
            ("6:23", 28),
        ] {
            let d = Design::find(name);
            let pts = Point::all(d);
            let p = &pts[c];
            for depth in 1..=if c == 0 { 4u32 } else { 3 } {
                let want = local(p, 2, depth);
                let got = six_local(&p.riesz(), depth);
                assert!(
                    got.sub(want).abs() < 1e-6 * (1.0 + want.abs()),
                    "{} {} {:?} {:?}",
                    name,
                    depth,
                    got,
                    want
                );
            }
        }
    }

    #[test]
    fn local_sums_count_the_square_roots() {
        let q = 36u64;
        for m in 0..q {
            let rho = (0..q).filter(|n| n * n % q == m).count() as f64;
            let mut s = Cx::ZERO;
            for den in (1..=q).filter(|d| q % d == 0) {
                for a in (0..den).filter(|&a| gcd(a, den) == 1) {
                    s = s.add(
                        gauss(2, (den - a) % den, den)
                            .scale(1.0 / den as f64)
                            .mul(Cx::turn(a * m % den, den)),
                    );
                }
            }
            assert!((s.re - rho).abs() < 1e-9 && s.im.abs() < 1e-9);
        }
    }
}
