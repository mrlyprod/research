use std::env;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;
use std::thread;
use std::time::Instant;

// DESIGNS

const DESIGNS: usize = 57;
const RS: usize = 57;
const COIN: usize = 58;
const STREAMS: usize = 59;
const GOLDEN: u64 = 0x9e37_79b9_7f4a_7c15;
const SALT_SQ: u64 = 1;
const SALT_PR: u64 = 2;

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

    fn parity(self, mut v: u64) -> u64 {
        let mut p = 0;
        while v > 0 {
            p ^= self.has(v % self.b) as u64;
            v /= self.b;
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

fn rudin(v: u64) -> u64 {
    ((v & (v >> 1)).count_ones() & 1) as u64
}

fn label(i: usize) -> String {
    match i {
        RS => "rudin".into(),
        COIN => "coin".into(),
        _ => Design::all()[i].name(),
    }
}

fn mix(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

fn coin(i: u64, salt: u64) -> u64 {
    mix(salt.wrapping_add((i + 1).wrapping_mul(GOLDEN))) >> 63
}

// READER

const P3: u64 = 1_594_323;
const P4: u64 = 1 << 20;
const P5: u64 = 1_953_125;
const P6: u64 = 1_679_616;

fn offset(b: u64) -> usize {
    (2..b).map(|c| (1usize << (c - 1)) - 1).sum()
}

fn build(b: u64, size: u64) -> Vec<u32> {
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

fn fold<const B: u64, T: Copy + Into<u32>>(t: &[T], mut v: u64) -> u64 {
    let mut a = 0u32;
    while v > 0 {
        a ^= t[(v % B) as usize].into();
        v /= B;
    }
    a as u64
}

struct Reader {
    t3: Vec<u8>,
    t4: Vec<u8>,
    t5: Vec<u16>,
    t6: Vec<u32>,
}

impl Reader {
    fn new() -> Reader {
        Reader {
            t3: build(3, P3).into_iter().map(|x| x as u8).collect(),
            t4: build(4, P4).into_iter().map(|x| x as u8).collect(),
            t5: build(5, P5).into_iter().map(|x| x as u16).collect(),
            t6: build(6, P6),
        }
    }

    fn word(&self, v: u64) -> u64 {
        (v.count_ones() & 1) as u64
            | fold::<P3, u8>(&self.t3, v) << offset(3)
            | fold::<P4, u8>(&self.t4, v) << offset(4)
            | fold::<P5, u16>(&self.t5, v) << offset(5)
            | fold::<P6, u32>(&self.t6, v) << offset(6)
            | rudin(v) << RS
    }
}

// TALLY

type Grid = Vec<[u64; 256]>;
const CHUNK: u32 = 16;

fn workers() -> usize {
    thread::available_parallelism()
        .map_or(4, |n| n.get())
        .min(8)
}

fn chunks(top: u32) -> Vec<(u64, u64, usize)> {
    let mut out = vec![(0, 1, 0)];
    for k in 1..=top.min(CHUNK) {
        out.push((1 << (k - 1), 1 << k, k as usize));
    }
    let mut s = 1u64 << CHUNK;
    while s < 1u64 << top {
        out.push((s, s + (1 << CHUNK), 64 - s.leading_zeros() as usize));
        s += 1 << CHUNK;
    }
    out
}

fn tally<F: Fn(u64) -> u64 + Sync>(top: u32, word: F) -> Grid {
    let jobs = chunks(top);
    let cells = (top as usize + 1) * STREAMS;
    let total = Mutex::new(vec![[0u64; 256]; cells]);
    let next = AtomicUsize::new(0);
    thread::scope(|s| {
        for _ in 0..workers() {
            s.spawn(|| {
                let mut mine = vec![[0u64; 256]; cells];
                let mut buf: Vec<u64> = Vec::with_capacity((1 << CHUNK) + 8);
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= jobs.len() {
                        break;
                    }
                    let (s0, s1, blk) = jobs[i];
                    buf.clear();
                    buf.extend((s0..s1 + 7).map(&word));
                    let len = (s1 - s0) as usize;
                    for d in 0..STREAMS {
                        let cnt = &mut mine[blk * STREAMS + d];
                        let mut w = 0usize;
                        for (j, x) in buf[..8].iter().enumerate() {
                            w |= ((x >> d) & 1) as usize * (1 << j);
                        }
                        cnt[w] += 1;
                        for x in &buf[8..len + 7] {
                            w = (w >> 1) | (((x >> d) & 1) as usize) << 7;
                            cnt[w] += 1;
                        }
                    }
                }
                let mut all = total.lock().unwrap();
                for (a, m) in all.iter_mut().zip(&mine) {
                    for (x, y) in a.iter_mut().zip(m) {
                        *x += y;
                    }
                }
            });
        }
    });
    total.into_inner().unwrap()
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

struct Stats {
    n: f64,
    t: [i64; 256],
    seen: usize,
}

impl Stats {
    fn at(grid: &Grid, k: u32, d: usize) -> Stats {
        let mut c = [0u64; 256];
        for blk in 0..=k as usize {
            for (x, y) in c.iter_mut().zip(&grid[blk * STREAMS + d]) {
                *x += y;
            }
        }
        Stats {
            n: (1u64 << k) as f64,
            t: walsh(&c),
            seen: c.iter().filter(|&&x| x > 0).count(),
        }
    }

    fn z(&self, s: usize) -> f64 {
        self.t[s] as f64 / self.n.sqrt()
    }

    fn bias(&self) -> f64 {
        self.z(1)
    }

    fn lag(&self, h: usize) -> f64 {
        self.z(1 | 1 << h)
    }

    fn chi(&self) -> f64 {
        (1..256).map(|s| (self.t[s] as f64).powi(2)).sum::<f64>() / self.n
    }

    fn order(&self, lo: u32, hi: u32) -> f64 {
        (1..256usize)
            .filter(|s| (lo..=hi).contains(&s.count_ones()))
            .map(|s| (self.t[s] as f64).powi(2))
            .sum::<f64>()
            / self.n
    }
}

fn null() -> (f64, f64) {
    let mut mean = 0.0;
    let mut var = 0.0;
    for s in (1..256usize).filter(|s| s & 1 == 1) {
        let m = (8 - (63 - s.leading_zeros() as usize)) as f64;
        mean += m;
        var += 2.0 * m * m;
    }
    (mean, var.sqrt())
}

// PRIMES

fn primes(count: usize) -> Vec<u32> {
    let n = count as f64;
    let limit = (n * (n.ln() + n.ln().ln()) + 1000.0) as u64;
    assert!(limit < 1 << 32);
    let root = (limit as f64).sqrt() as u64 + 2;
    let mut small = vec![true; root as usize + 1];
    let mut base = Vec::new();
    for i in 2..=root as usize {
        if small[i] {
            if i > 2 {
                base.push(i as u64);
            }
            let mut j = i * i;
            while j <= root as usize {
                small[j] = false;
                j += i;
            }
        }
    }
    let seg: u64 = 1 << 20;
    let segments = limit / (2 * seg) + 1;
    let next = AtomicUsize::new(0);
    let found = Mutex::new(Vec::new());
    thread::scope(|s| {
        for _ in 0..workers() {
            s.spawn(|| {
                let mut sieve = vec![true; seg as usize];
                loop {
                    let k = next.fetch_add(1, Ordering::Relaxed) as u64;
                    if k >= segments {
                        break;
                    }
                    let lo = 2 * seg * k;
                    sieve.fill(true);
                    for &p in &base {
                        if p * p >= lo + 2 * seg {
                            break;
                        }
                        let mut m = (lo + p - 1) / p * p;
                        m = m.max(p * p);
                        if m % 2 == 0 {
                            m += p;
                        }
                        while m < lo + 2 * seg {
                            sieve[((m - lo) / 2) as usize] = false;
                            m += 2 * p;
                        }
                    }
                    let mut out = Vec::new();
                    if k == 0 {
                        out.push(2u32);
                        sieve[0] = false;
                    }
                    for (i, &ok) in sieve.iter().enumerate() {
                        if ok {
                            out.push((lo + 2 * i as u64 + 1) as u32);
                        }
                    }
                    found.lock().unwrap().push((k, out));
                }
            });
        }
    });
    let mut parts = found.into_inner().unwrap();
    parts.sort_by_key(|p| p.0);
    let mut all: Vec<u32> = parts.into_iter().flat_map(|p| p.1).collect();
    assert!(all.len() >= count + 7, "sieve too short");
    all.truncate(count + 7);
    all
}

// SPECTRUM

#[derive(Clone, Copy, Debug)]
struct Iv {
    lo: f64,
    hi: f64,
}

impl Iv {
    fn pt(x: f64) -> Iv {
        Iv { lo: x, hi: x }
    }

    fn add(self, o: Iv) -> Iv {
        Iv {
            lo: (self.lo + o.lo).next_down(),
            hi: (self.hi + o.hi).next_up(),
        }
    }

    fn neg(self) -> Iv {
        Iv {
            lo: -self.hi,
            hi: -self.lo,
        }
    }

    fn mul(self, o: Iv) -> Iv {
        let c = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        Iv {
            lo: c.iter().copied().fold(f64::INFINITY, f64::min).next_down(),
            hi: c
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up(),
        }
    }

    fn div(self, o: Iv) -> Iv {
        assert!(o.lo > 0.0);
        let c = [
            self.lo / o.lo,
            self.lo / o.hi,
            self.hi / o.lo,
            self.hi / o.hi,
        ];
        Iv {
            lo: c.iter().copied().fold(f64::INFINITY, f64::min).next_down(),
            hi: c
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max)
                .next_up(),
        }
    }
}

fn two_pi() -> Iv {
    Iv {
        lo: 2.0 * std::f64::consts::PI,
        hi: 2.0 * std::f64::consts::PI.next_up(),
    }
}

fn taylor(x: Iv, odd: bool) -> Iv {
    let x2 = x.mul(x);
    let mut term = if odd { x } else { Iv::pt(1.0) };
    let mut sum = term;
    let mut k = if odd { 1.0 } else { 0.0 };
    for _ in 0..15 {
        term = term.mul(x2).div(Iv::pt((k + 1.0) * (k + 2.0))).neg();
        k += 2.0;
        sum = sum.add(term);
    }
    sum.add(Iv {
        lo: -1e-30,
        hi: 1e-30,
    })
}

fn atanh_series(y: Iv) -> Iv {
    let y2 = y.mul(y);
    let mut term = y;
    let mut sum = y;
    for k in 1..40 {
        term = term.mul(y2);
        sum = sum.add(term.div(Iv::pt((2 * k + 1) as f64)));
    }
    let tail = y.hi.abs().max(y.lo.abs()).powi(81) * 2.0;
    sum.add(Iv {
        lo: -tail - 1e-300,
        hi: tail + 1e-300,
    })
    .mul(Iv::pt(2.0))
}

fn ln_iv(x: Iv) -> Iv {
    assert!(x.lo > 0.0);
    let ln2 = atanh_series(Iv::pt(1.0).div(Iv::pt(3.0)));
    let end = |v: f64| {
        let e = v.log2().floor() as i32;
        let m = v / 2f64.powi(e);
        let y = Iv::pt(m).add(Iv::pt(-1.0)).div(Iv::pt(m).add(Iv::pt(1.0)));
        ln2.mul(Iv::pt(e as f64)).add(atanh_series(y))
    };
    Iv {
        lo: end(x.lo).lo,
        hi: end(x.hi).hi,
    }
}

fn cos2pi(p: i128, q: i128) -> Iv {
    let mut p = p.rem_euclid(q);
    let mut q = q;
    if 2 * p > q {
        p = q - p;
    }
    let mut flip = false;
    if 4 * p > q {
        p = q - 2 * p;
        q *= 2;
        flip = true;
    }
    assert!(q < 1 << 52);
    let v = if 8 * p > q {
        taylor(
            two_pi()
                .mul(Iv::pt((q - 4 * p) as f64))
                .div(Iv::pt((4 * q) as f64)),
            true,
        )
    } else {
        taylor(two_pi().mul(Iv::pt(p as f64)).div(Iv::pt(q as f64)), false)
    };
    if flip {
        v.neg()
    } else {
        v
    }
}

struct Tile {
    b: u64,
    a: Vec<i64>,
}

impl Tile {
    fn of(d: Design) -> Tile {
        let b = d.b;
        let a = (0..b)
            .map(|s| (0..b - s).map(|r| d.sign(r) * d.sign(r + s)).sum())
            .collect();
        Tile { b, a }
    }

    fn phi2(&self, x: f64) -> f64 {
        let mut s = self.b as f64;
        for (k, &a) in self.a.iter().enumerate().skip(1) {
            s += 2.0 * a as f64 * (2.0 * std::f64::consts::PI * k as f64 * x).cos();
        }
        (s / (self.b * self.b) as f64).clamp(0.0, 1.0)
    }

    fn phi2_iv(&self, p: i128, q: i128) -> Iv {
        let mut s = Iv::pt(self.b as f64);
        for (k, &a) in self.a.iter().enumerate().skip(1) {
            if a != 0 {
                s = s.add(Iv::pt(2.0 * a as f64).mul(cos2pi(k as i128 * p, q)));
            }
        }
        let r = s.div(Iv::pt((self.b * self.b) as f64));
        Iv {
            lo: r.lo.max(0.0),
            hi: r.hi.min(1.0),
        }
    }

    fn slope(&self) -> f64 {
        self.a
            .iter()
            .enumerate()
            .skip(1)
            .map(|(k, &a)| k as f64 * a.abs() as f64)
            .sum::<f64>()
    }
}

// ORBITS

#[derive(Clone, Copy, Debug)]
struct Orbit {
    k: u32,
    a: u64,
    q: u64,
    eta: f64,
}

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn squarefree(mut q: u64) -> bool {
    let mut p = 2;
    while p * p <= q {
        if q % (p * p) == 0 {
            return false;
        }
        if q % p == 0 {
            q /= p;
        }
        p += 1;
    }
    true
}

const REACH: u64 = 1 << 20;

fn periodic(t: &Tile) -> Vec<Orbit> {
    let b = t.b;
    let mut out = Vec::new();
    let mut k = 1u32;
    while b.pow(k) <= REACH {
        let m = b.pow(k) - 1;
        let g: Vec<f64> = (0..m).map(|c| t.phi2(c as f64 / m as f64).ln()).collect();
        'a: for a in 0..m {
            let mut c = a;
            let mut s = 0.0;
            for j in 0..k {
                if c < a || (j > 0 && c == a) {
                    continue 'a;
                }
                s += g[c as usize];
                c = c * b % m;
            }
            let g0 = gcd(a, m);
            out.push(Orbit {
                k,
                a: a / g0,
                q: m / g0,
                eta: -s / (2.0 * k as f64 * (b as f64).ln()),
            });
        }
        k += 1;
    }
    out
}

fn best<F: Fn(&Orbit) -> bool>(orbits: &[Orbit], keep: F) -> Orbit {
    *orbits
        .iter()
        .filter(|o| keep(o))
        .min_by(|x, y| x.eta.total_cmp(&y.eta).then(x.q.cmp(&y.q)))
        .unwrap()
}

fn point(t: &Tile, num: u64, den: u64) -> f64 {
    let mut seen = vec![usize::MAX; den as usize];
    let mut path = Vec::new();
    let mut x = num % den;
    while seen[x as usize] == usize::MAX {
        seen[x as usize] = path.len();
        path.push(x);
        x = x * t.b % den;
    }
    let vals: Vec<f64> = path
        .iter()
        .map(|&y| t.phi2(y as f64 / den as f64))
        .collect();
    if vals.iter().any(|&v| v < 1e-20) {
        return f64::INFINITY;
    }
    let cyc = &vals[seen[x as usize]..];
    -cyc.iter().map(|v| v.ln()).sum::<f64>() / (2.0 * cyc.len() as f64 * (t.b as f64).ln())
}

fn lag(t: &Tile, h: u64) -> f64 {
    (0..2 * h)
        .map(|k| point(t, k, 2 * h))
        .fold(f64::INFINITY, f64::min)
}

struct Spectrum {
    top: Orbit,
    sq: Orbit,
    pr: Orbit,
    lags: [f64; 8],
}

fn spectrum(d: Design) -> Spectrum {
    let t = Tile::of(d);
    let orbits = periodic(&t);
    let mut lags = [f64::INFINITY; 8];
    for (h, e) in lags.iter_mut().enumerate().skip(1) {
        *e = lag(&t, h as u64);
    }
    Spectrum {
        top: best(&orbits, |_| true),
        sq: best(&orbits, |o| o.q % 4 != 2),
        pr: best(&orbits, |o| squarefree(o.q)),
        lags,
    }
}

// ENERGY

fn energy(d: Design) -> [[i64; 2]; 2] {
    let b = d.b as usize;
    let al: Vec<i64> = (0..=b)
        .map(|s| {
            (0..b.saturating_sub(s))
                .map(|r| d.sign(r as u64) * d.sign((r + s) as u64))
                .sum()
        })
        .collect();
    let t11 = al[0] * al[0] + 2 * (1..b).map(|t| al[t] * al[t]).sum::<i64>();
    let t12 = 2 * (0..b).map(|s| al[s] * al[b - s]).sum::<i64>();
    let t21 = (0..b - 1)
        .map(|s| al[s] * al[s + 1] + al[b - s] * al[b - s - 1])
        .sum::<i64>()
        + al[1] * al[0];
    let t22 = (0..b - 1)
        .map(|s| al[s] * al[b - s - 1] + al[b - s] * al[s + 1])
        .sum::<i64>()
        + al[b - 1] * al[0];
    [[t11, t12], [t21, t22]]
}

fn dim2(d: Design) -> f64 {
    let t = energy(d);
    let tr = (t[0][0] + t[1][1]) as f64;
    let det = (t[0][0] * t[1][1] - t[0][1] * t[1][0]) as f64;
    let disc = tr * tr - 4.0 * det;
    let rho = if disc >= 0.0 {
        ((tr + disc.sqrt()) / 2.0)
            .abs()
            .max(((tr - disc.sqrt()) / 2.0).abs())
    } else {
        det.sqrt()
    };
    3.0 - rho.ln() / (d.b as f64).ln()
}

fn ranks(v: &[f64]) -> Vec<f64> {
    let mut idx: Vec<usize> = (0..v.len()).collect();
    idx.sort_by(|&a, &b| v[a].total_cmp(&v[b]));
    let mut out = vec![0.0; v.len()];
    let mut i = 0;
    while i < idx.len() {
        let mut j = i;
        while j + 1 < idx.len() && v[idx[j + 1]] == v[idx[i]] {
            j += 1;
        }
        for &k in &idx[i..=j] {
            out[k] = (i + j) as f64 / 2.0;
        }
        i = j + 1;
    }
    out
}

fn spearman(x: &[f64], y: &[f64]) -> f64 {
    let (rx, ry) = (ranks(x), ranks(y));
    let n = x.len() as f64;
    let (mx, my) = (rx.iter().sum::<f64>() / n, ry.iter().sum::<f64>() / n);
    let cov: f64 = rx.iter().zip(&ry).map(|(a, b)| (a - mx) * (b - my)).sum();
    let vx: f64 = rx.iter().map(|a| (a - mx).powi(2)).sum();
    let vy: f64 = ry.iter().map(|b| (b - my).powi(2)).sum();
    cov / (vx * vy).sqrt()
}

// CERTIFICATE

fn evaluate(b: usize, wt: &[f64], pol: &[usize], eta: &mut [f64], h: &mut [f64]) {
    let n = wt.len();
    let next = |i: usize| (i * b + pol[i]) % n;
    let mut state = vec![0u8; n];
    let mut stack = Vec::new();
    for s in 0..n {
        if state[s] != 0 {
            continue;
        }
        stack.clear();
        let mut i = s;
        while state[i] == 0 {
            state[i] = 1;
            stack.push(i);
            i = next(i);
        }
        if state[i] == 1 {
            let pos = stack.iter().position(|&x| x == i).unwrap();
            let len = (stack.len() - pos) as f64;
            let mu = stack[pos..].iter().map(|&c| wt[c]).sum::<f64>() / len;
            h[i] = 0.0;
            eta[i] = mu;
            state[i] = 2;
            for &c in stack[pos + 1..].iter().rev() {
                h[c] = wt[c] - mu + h[next(c)];
                eta[c] = mu;
                state[c] = 2;
            }
            stack.truncate(pos);
        }
        while let Some(c) = stack.pop() {
            let j = next(c);
            eta[c] = eta[j];
            h[c] = wt[c] - eta[j] + h[j];
            state[c] = 2;
        }
    }
}

fn howard(b: usize, wt: &[f64]) -> (Vec<f64>, Vec<f64>, Vec<usize>, usize) {
    let n = wt.len();
    let succ = |i: usize, d: usize| (i * b + d) % n;
    let mut pol: Vec<usize> = (0..n)
        .map(|i| {
            (0..b)
                .max_by(|&x, &y| wt[succ(i, x)].total_cmp(&wt[succ(i, y)]))
                .unwrap()
        })
        .collect();
    let mut eta = vec![0.0; n];
    let mut h = vec![0.0; n];
    let mut rounds = 0;
    loop {
        evaluate(b, wt, &pol, &mut eta, &mut h);
        rounds += 1;
        let mut changed = false;
        for i in 0..n {
            let cur = succ(i, pol[i]);
            let (mut be, mut bh, mut pick) = (eta[cur], h[cur], pol[i]);
            for d in 0..b {
                let j = succ(i, d);
                if eta[j] > be + 1e-12 || (eta[j] >= be - 1e-12 && h[j] > bh + 1e-9) {
                    be = eta[j];
                    bh = h[j];
                    pick = d;
                }
            }
            if pick != pol[i] {
                pol[i] = pick;
                changed = true;
            }
        }
        if !changed || rounds >= 1000 {
            break;
        }
    }
    (eta, h, pol, rounds)
}

fn depth(b: u64) -> u32 {
    let mut r = 1;
    while b.pow(r + 1) <= 1 << 17 {
        r += 1;
    }
    r
}

struct Bracket {
    lo: f64,
    hi: f64,
    orbit: (u32, u64, u64),
    cycle: usize,
    rounds: usize,
}

fn certify(t: &Tile, k: u32, a: u64) -> f64 {
    let b = t.b;
    let m = b.pow(k) - 1;
    let mut prod = Iv::pt(1.0);
    let mut c = a;
    for _ in 0..k {
        prod = prod.mul(t.phi2_iv(c as i128, m as i128));
        c = c * b % m;
    }
    if prod.lo <= 0.0 {
        return f64::INFINITY;
    }
    ln_iv(Iv::pt(prod.lo))
        .neg()
        .div(ln_iv(Iv::pt(b as f64)).mul(Iv::pt(2.0 * k as f64)))
        .hi
}

fn bracket(d: Design, guess: Orbit) -> Bracket {
    let t = Tile::of(d);
    let b = t.b;
    let r = depth(b);
    let n = b.pow(r) as usize;
    let lip = two_pi()
        .mul(Iv::pt(t.slope()))
        .div(Iv::pt((b * b) as f64 * n as f64))
        .hi;
    let g2: Vec<f64> = (0..n)
        .map(|i| {
            let c = t.phi2_iv(2 * i as i128 + 1, 2 * n as i128).hi;
            (c + lip).next_up().min(1.0)
        })
        .collect();
    let wt: Vec<f64> = g2.iter().map(|g| g.ln()).collect();
    let (eta, h, pol, rounds) = howard(b as usize, &wt);
    let hmax = h.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let v: Vec<f64> = h.iter().map(|x| (x - hmax).exp().max(1e-300)).collect();
    let mut c = 0.0f64;
    for i in 0..n {
        let mx = (0..b as usize)
            .map(|dd| v[(i * b as usize + dd) % n])
            .fold(0.0, f64::max);
        let u = (g2[i] * mx).next_up();
        c = c.max((u / v[i]).next_up());
    }
    let lo = ln_iv(Iv::pt(c))
        .neg()
        .div(ln_iv(Iv::pt(b as f64)).mul(Iv::pt(2.0)))
        .lo
        .max(0.0);
    let start = (0..n).max_by(|&x, &y| eta[x].total_cmp(&eta[y])).unwrap();
    let mut seen = vec![usize::MAX; n];
    let mut walk = Vec::new();
    let mut i = start;
    while seen[i] == usize::MAX {
        seen[i] = walk.len();
        walk.push(i);
        i = (i * b as usize + pol[i]) % n;
    }
    let cyc: Vec<usize> = walk[seen[i]..].to_vec();
    let mut orbit = (guess.k, guess.a, guess.q);
    let mut hi = certify(&t, guess.k, guess.a * ((b.pow(guess.k) - 1) / guess.q));
    let k = cyc.len() as u32;
    if (k as f64) * (b as f64).log2() < 46.0 {
        let mut digits: Vec<u64> = Vec::new();
        let mut x = cyc[0];
        for _ in 0..r {
            digits.push((x % b as usize) as u64);
            x /= b as usize;
        }
        digits.reverse();
        for &c in &cyc {
            digits.push(pol[c] as u64);
        }
        let word = digits[..k as usize]
            .iter()
            .fold(0u64, |acc, &x| acc * b + x);
        let alt = certify(&t, k, word);
        if alt < hi {
            hi = alt;
            let m = b.pow(k) - 1;
            orbit = (k, word / gcd(word, m), m / gcd(word, m));
        }
    }
    Bracket {
        lo,
        hi,
        orbit,
        cycle: cyc.len(),
        rounds,
    }
}

// DIAGONAL

struct Corr {
    b: u64,
    m: Vec<[[f64; 2]; 2]>,
    c1: f64,
}

impl Corr {
    fn of(d: Design) -> Corr {
        let b = d.b;
        let mut al: Vec<f64> = (0..=b)
            .map(|s| {
                (0..b.saturating_sub(s))
                    .map(|r| (d.sign(r) * d.sign(r + s)) as f64)
                    .sum()
            })
            .collect();
        al.push(0.0);
        let bf = b as f64;
        let m = (0..b as usize)
            .map(|s| {
                [
                    [al[s] / bf, al[b as usize - s] / bf],
                    [al[s + 1] / bf, al[b as usize - s - 1] / bf],
                ]
            })
            .collect();
        Corr {
            b,
            m,
            c1: al[1] / (bf - al[b as usize - 1]),
        }
    }

    fn row(&self, w: [f64; 2], s: usize) -> [f64; 2] {
        let m = &self.m[s];
        [
            w[0] * m[0][0] + w[1] * m[1][0],
            w[0] * m[0][1] + w[1] * m[1][1],
        ]
    }

    fn tail(&self, mut w: [f64; 2], mut q: u64) -> f64 {
        while q > 0 {
            w = self.row(w, (q % self.b) as usize);
            q /= self.b;
        }
        w[0] + w[1] * self.c1
    }

    fn at(&self, g: u64) -> f64 {
        self.tail([1.0, 0.0], g)
    }

    fn progression(&self, n: u64, a: u64, c: u64) -> f64 {
        let b = self.b;
        let mut nd = Vec::new();
        let mut x = n;
        while x > 0 {
            nd.push(x % b);
            x /= b;
        }
        let cap = (a + 2) as usize;
        let mut cur = vec![[[0.0f64; 2]; 3]; cap];
        cur[0][1] = [1.0, 0.0];
        let mut cq = c;
        for &nj in &nd {
            let cj = cq % b;
            cq /= b;
            let mut nxt = vec![[[0.0f64; 2]; 3]; cap];
            for (carry, states) in cur.iter().enumerate() {
                for (cmp, w) in states.iter().enumerate() {
                    if w[0] == 0.0 && w[1] == 0.0 {
                        continue;
                    }
                    for dgt in 0..b {
                        let t = a * dgt + carry as u64 + cj;
                        let cmp2 = if dgt < nj {
                            0
                        } else if dgt > nj {
                            2
                        } else {
                            cmp
                        };
                        let v = self.row(*w, (t % b) as usize);
                        let slot = &mut nxt[(t / b) as usize][cmp2];
                        slot[0] += v[0];
                        slot[1] += v[1];
                    }
                }
            }
            cur = nxt;
        }
        cur.iter()
            .enumerate()
            .map(|(carry, states)| self.tail(states[0], carry as u64 + cq))
            .sum()
    }
}

fn diagonal(top: u32, ptop: u32) {
    let all = Design::all();
    let reader = Reader::new();
    let sq = tally(top, |n| reader.word(n * n));
    let ps = primes(1 << ptop);
    let pr = tally(ptop, |i| reader.word(ps[i as usize] as u64));
    let n = 1u64 << top;
    let np = 1usize << ptop;
    println!("# DIAGONAL");
    println!();
    println!(
        "squares: model z = sqrt(N)^-1 sum_(n<N) C_x(2hn + h^2), C_x the design's autocorrelation by its 2 x 2 digit matrices, summed exactly by a digit walk over n < N = 2^{} with the carry of 2hn + h^2; primes: model z = sqrt(N)^-1 sum_(i<N) x(r) x((r + g) mod b) C_x(floor((r + g)/b)), r = p_i mod b, g = p_(i+h) - p_i, N = 2^{}; each cell model/observed",
        top, ptop
    );
    println!();
    let mut lit = 0;
    let mut worst: f64 = 0.0;
    let mut worst_at = String::new();
    let mut rms = 0.0;
    let mut cells = 0;
    let mut plit = 0;
    let mut pworst: f64 = 0.0;
    let mut pworst_at = String::new();
    let mut prms = 0.0;
    for (i, d) in all.iter().enumerate() {
        let corr = Corr::of(*d);
        let s = Stats::at(&sq, top, i);
        let p = Stats::at(&pr, ptop, i);
        print!("{:<7}", d.name());
        for h in 1..8u64 {
            let model = corr.progression(n, 2 * h, h * h) / (n as f64).sqrt();
            let seen = s.lag(h as usize);
            let gap = model - seen;
            cells += 1;
            rms += gap * gap;
            if model.abs() > 4.0 {
                lit += 1;
            }
            if gap.abs() > worst.abs() {
                worst = gap;
                worst_at = format!("{} h={}", d.name(), h);
            }
            print!(" {:>7.2}/{:<7.2}", model, seen);
        }
        let b = d.b;
        print!("  |");
        for h in 1..4usize {
            let mut sum = 0.0;
            for j in 0..np {
                let r = ps[j] as u64 % b;
                let g = (ps[j + h] - ps[j]) as u64;
                let sign = (d.sign(r) * d.sign((r + g) % b)) as f64;
                sum += sign * corr.at((r + g) / b);
            }
            let model = sum / (np as f64).sqrt();
            let seen = p.lag(h);
            let gap = model - seen;
            prms += gap * gap;
            if model.abs() > 4.0 {
                plit += 1;
            }
            if gap.abs() > pworst.abs() {
                pworst = gap;
                pworst_at = format!("{} h={}", d.name(), h);
            }
            print!(" {:>8.1}/{:<8.1}", model, seen);
        }
        println!();
    }
    println!();
    println!(
        "squares: {} cells, {} with |model| > 4, rms of model - observed {:.3}, largest {:.2} at {}",
        cells,
        lit,
        (rms / cells as f64).sqrt(),
        worst,
        worst_at
    );
    println!(
        "primes: {} cells, {} with |model| > 4, rms of model - observed {:.3}, largest {:.2} at {}",
        3 * DESIGNS,
        plit,
        (prms / (3 * DESIGNS) as f64).sqrt(),
        pworst,
        pworst_at
    );
}

// CENSUS

fn sequences(all: &[Design], len: u64) -> Vec<Vec<u64>> {
    all.iter()
        .map(|d| (0..len).map(|n| d.parity(n)).collect())
        .collect()
}

fn classes(all: &[Design]) -> Vec<String> {
    classify(all, true)
}

fn classify(all: &[Design], show: bool) -> Vec<String> {
    let len = 1 << 16;
    let seqs = sequences(all, len);
    let mut tags = vec![String::new(); DESIGNS];
    let mut distinct = 0;
    let mut folded = 0;
    let mut pairs = 0;
    for i in 0..DESIGNS {
        if let Some(j) = (0..i).find(|&j| seqs[j] == seqs[i]) {
            tags[i] = format!("copy {}", all[j].name());
            continue;
        }
        distinct += 1;
        if (0..len).all(|n| seqs[i][n as usize] == n % 2) {
            tags[i] = "folded".into();
            folded += 1;
            continue;
        }
        if let Some(tw) = all[i].twin() {
            let j = all.iter().position(|&e| e == tw).unwrap();
            assert!((0..len).all(|n| seqs[i][n as usize] ^ (n % 2) == seqs[j][n as usize]));
            tags[i] = format!("twin {}", tw.name());
            pairs += 1;
        }
    }
    if !show {
        return tags;
    }
    println!("# CLASSES");
    println!();
    println!(
        "designs {}  distinct sequences on n < 2^16 {}  folded (n mod 2) {}  designs in twin pairs {}",
        DESIGNS, distinct, folded, pairs
    );
    for (i, tag) in tags.iter().enumerate() {
        if !tag.is_empty() {
            println!("  {:<7} {}", all[i].name(), tag);
        }
    }
    let mirrors: Vec<String> = all
        .iter()
        .filter(|d| d.mirror() != **d && d.mirror().mask > d.mask)
        .map(|d| format!("{}~{}", d.name(), d.mirror().name()))
        .collect();
    let selfs = all.iter().filter(|d| d.mirror() == **d).count();
    println!(
        "mirror pairs, v*(d) = v(b-1) v(b-1-d), same tile autocorrelation hence the same spectral measure: {} pairs, {} self-mirror: {}",
        mirrors.len(),
        selfs,
        mirrors.join(" ")
    );
    println!();
    tags
}

fn seen_along_n(reader: &Reader) -> Vec<usize> {
    let len = 1u64 << 20;
    let words: Vec<u64> = (0..len + 7)
        .map(|n| reader.word(n) | coin(n, 0) << COIN)
        .collect();
    (0..STREAMS)
        .map(|d| {
            let mut hit = [false; 256];
            let mut w = 0usize;
            for (j, x) in words.iter().enumerate() {
                w = (w >> 1) | (((x >> d) & 1) as usize) << 7;
                if j >= 7 {
                    hit[w] = true;
                }
            }
            hit.iter().filter(|&&x| x).count()
        })
        .collect()
}

fn status(d: usize, tags: &[String]) -> &'static str {
    if d == RS {
        return "proved";
    }
    if d >= DESIGNS {
        return "";
    }
    let b = Design::all()[d].b;
    if tags[d] == "folded" || tags[d] == "copy 3:1" {
        "periodic"
    } else if b % 2 == 0 {
        "proved"
    } else {
        "open"
    }
}

fn fmt(e: f64) -> String {
    if e.is_finite() {
        format!("{:.4}", e)
    } else {
        "-".into()
    }
}

fn census(top: u32, ptop: u32) {
    let all = Design::all();
    assert_eq!(all.len(), DESIGNS);
    let tags = classes(&all);
    let clock = Instant::now();
    let reader = Reader::new();
    let along = seen_along_n(&reader);
    let sq = tally(top, |n| reader.word(n * n) | coin(n, SALT_SQ) << COIN);
    let tsq = clock.elapsed().as_secs_f64();
    let ps = primes(1 << ptop);
    let pr = tally(ptop, |i| {
        reader.word(ps[i as usize] as u64) | coin(i, SALT_PR) << COIN
    });
    let tpr = clock.elapsed().as_secs_f64() - tsq;
    let (mean, sd) = null();
    println!("# SQUARES");
    println!();
    println!(
        "y(n) = x(n^2), 8-blocks at every start n < N = 2^{}; z = sqrt(N) . mean, so z_B = S(N)/sqrt(N); X2 the overlapping 8-block chi-square, iid mean {:.0} and sd {:.1}; X2 = X2(1) + X2(2) + X2(3+) by the order of the Walsh coefficients; seen counts the 8-blocks met along n^2 and, in the last column, along n < 2^20; tally {:.1} s",
        top, mean, sd, tsq
    );
    println!();
    println!(
        "{:<7} {:<6} {:<8} {:>8} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>9} {:>9} {:>9} {:>9} {:>4} {:>4}",
        "design", "class", "normal", "zB", "zC1", "zC2", "zC3", "zC4", "zC5", "zC6", "zC7", "X2", "X2(1)", "X2(2)", "X2(3+)", "sq", "n"
    );
    for d in 0..STREAMS {
        let s = Stats::at(&sq, top, d);
        let tag = if d < DESIGNS { tags[d].as_str() } else { "" };
        print!(
            "{:<7} {:<6} {:<8} {:>8.2}",
            label(d),
            tag.split(' ').next().unwrap(),
            status(d, &tags),
            s.bias()
        );
        for h in 1..8 {
            print!(" {:>7.2}", s.lag(h));
        }
        println!(
            " {:>9.1} {:>9.1} {:>9.1} {:>9.1} {:>4} {:>4}",
            s.chi(),
            s.order(1, 1),
            s.order(2, 2),
            s.order(3, 8),
            s.seen,
            along[d]
        );
    }
    println!();
    println!("# GROWTH ALONG THE SQUARES");
    println!();
    let ks: Vec<u32> = (12..=top).filter(|k| (top - k) % 2 == 0).collect();
    print!("{:<7} {:<5}", "design", "stat");
    for k in &ks {
        print!(" {:>8}", format!("2^{}", k));
    }
    println!();
    for d in 0..STREAMS {
        for stat in ["zB", "X2"] {
            print!("{:<7} {:<5}", label(d), stat);
            for &k in &ks {
                let s = Stats::at(&sq, k, d);
                let v = if stat == "zB" { s.bias() } else { s.chi() };
                print!(" {:>8.1}", v);
            }
            println!();
        }
    }
    println!();
    println!("# PRIMES");
    println!();
    println!(
        "y(i) = x(p_i) along the consecutive primes, i < N = 2^{}, p_N = {}; same columns; sieve and tally {:.1} s",
        ptop,
        ps[(1usize << ptop) - 1],
        tpr
    );
    println!();
    println!(
        "{:<7} {:>8} {:>8} {:>8} {:>8} {:>10} {:>10} {:>10} {:>10}",
        "design", "zB", "zC1", "zC2", "zC3", "X2", "X2(1)", "X2(2)", "X2(3+)"
    );
    for d in 0..STREAMS {
        let s = Stats::at(&pr, ptop, d);
        println!(
            "{:<7} {:>8.1} {:>8.1} {:>8.1} {:>8.1} {:>10.0} {:>10.0} {:>10.0} {:>10.0}",
            label(d),
            s.bias(),
            s.lag(1),
            s.lag(2),
            s.lag(3),
            s.chi(),
            s.order(1, 1),
            s.order(2, 2),
            s.order(3, 8)
        );
    }
    println!();
    println!("# GROWTH ALONG THE PRIMES");
    println!();
    let kp: Vec<u32> = (12..=ptop).filter(|k| (ptop - k) % 2 == 0).collect();
    print!("{:<7} {:<5}", "design", "stat");
    for k in &kp {
        print!(" {:>9}", format!("2^{}", k));
    }
    println!();
    for d in 0..STREAMS {
        print!("{:<7} {:<5}", label(d), "zB");
        for &k in &kp {
            print!(" {:>9.1}", Stats::at(&pr, k, d).bias());
        }
        println!();
    }
    println!();
    law(&all, &sq, top, &pr, ptop, &tags);
}

fn slope(grid: &Grid, d: usize, ks: &[u32], stat: usize) -> f64 {
    let pts: Vec<(f64, f64)> = ks
        .iter()
        .map(|&k| {
            let s = Stats::at(grid, k, d);
            (k as f64, s.z(stat).abs().max(1e-9).log2())
        })
        .collect();
    let n = pts.len() as f64;
    let mx = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
    pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum::<f64>()
        / pts.iter().map(|p| (p.0 - mx).powi(2)).sum::<f64>()
}

fn law(all: &[Design], sq: &Grid, top: u32, pr: &Grid, ptop: u32, tags: &[String]) {
    let (mean, sd) = null();
    let floor = mean + 4.0 * sd;
    println!("# LAW");
    println!();
    println!(
        "eta: best periodic orbit a/q of the Fourier products, period <= log_b 2^20; D2 = 3 - log_b rho(T), T the 2 x 2 energy recursion of the autocorrelation; sq: best orbit with q != 2 mod 4, pred = 1/2 - 2 eta_sq; pr: best orbit with q squarefree, pred = 1/2 - eta_pr; slope: least-squares slope of log2 |z_B| over the last 9 checkpoints; lit: lags h <= 7 with |z_Ch| > 4 at N = 2^{}; e<1/4: lags whose best orbit through k/(2h) has eta < 1/4",
        top
    );
    println!();
    println!(
        "{:<7} {:<6} {:>7} {:>8} {:>7} {:>7} {:>7} {:>7} {:>7} {:>8} {:>7} {:>7} {:>7} {:>7}  {:<10} {}",
        "design", "class", "eta", "orbit", "D2", "eta_sq", "pred", "zB", "slope", "X2", "eta_pr", "pred", "zB pr", "slope", "lit", "e<1/4"
    );
    let kq: Vec<u32> = (top - 8..=top).collect();
    let kr: Vec<u32> = (ptop - 8..=ptop).collect();
    let mut rows = Vec::new();
    let mut lit_ok = 0;
    let mut lit_off = Vec::new();
    let mut quiet = Vec::new();
    for (i, d) in all.iter().enumerate() {
        let sp = spectrum(*d);
        let s = Stats::at(sq, top, i);
        let p = Stats::at(pr, ptop, i);
        let tag = tags[i].split(' ').next().unwrap();
        let lit: Vec<usize> = (1..8).filter(|&h| s.lag(h).abs() > 4.0).collect();
        let want: Vec<usize> = (1..8).filter(|&h| sp.lags[h] < 0.25).collect();
        let show = |v: &[usize]| {
            v.iter()
                .map(|h| h.to_string())
                .collect::<Vec<_>>()
                .join(",")
        };
        println!(
            "{:<7} {:<6} {:>7} {:>8} {:>7.4} {:>7} {:>7.3} {:>7.1} {:>7.3} {:>8.0} {:>7} {:>7.3} {:>7.1} {:>7.3}  {:<10} {}",
            d.name(),
            tag,
            fmt(sp.top.eta),
            format!("{}/{}", sp.top.a, sp.top.q),
            dim2(*d),
            fmt(sp.sq.eta),
            0.5 - 2.0 * sp.sq.eta,
            s.bias(),
            slope(sq, i, &kq, 1),
            s.chi(),
            fmt(sp.pr.eta),
            0.5 - sp.pr.eta,
            p.bias(),
            slope(pr, i, &kr, 1),
            show(&lit),
            show(&want)
        );
        if tag == "copy" || tag == "folded" {
            continue;
        }
        for &h in &lit {
            if want.contains(&h) {
                lit_ok += 1;
            } else {
                lit_off.push(format!("{} h={} z={:.1}", d.name(), h, s.lag(h)));
            }
        }
        for &h in &want {
            if !lit.contains(&h) {
                quiet.push(format!("{} h={}", d.name(), h));
            }
        }
        rows.push((d.name(), sp.top.eta, dim2(*d), s.chi(), p.bias().abs()));
    }
    for d in [RS, COIN] {
        let s = Stats::at(sq, top, d);
        let p = Stats::at(pr, ptop, d);
        println!(
            "{:<7} {:<6} {:>7} {:>8} {:>7} {:>7} {:>7} {:>7.1} {:>7.3} {:>8.0} {:>7} {:>7} {:>7.1} {:>7.3}",
            label(d),
            "",
            "",
            "",
            "",
            "",
            "",
            s.bias(),
            slope(sq, d, &kq, 1),
            s.chi(),
            "",
            "",
            p.bias(),
            slope(pr, d, &kr, 1)
        );
    }
    println!();
    let eta: Vec<f64> = rows.iter().map(|r| r.1).collect();
    let d2: Vec<f64> = rows.iter().map(|r| r.2).collect();
    let chi: Vec<f64> = rows.iter().map(|r| r.3).collect();
    let zp: Vec<f64> = rows.iter().map(|r| r.4).collect();
    let above: Vec<&(String, f64, f64, f64, f64)> = rows.iter().filter(|r| r.1 > 0.25).collect();
    let below: Vec<&(String, f64, f64, f64, f64)> = rows.iter().filter(|r| r.1 < 0.25).collect();
    let (top_above, top_name) = above
        .iter()
        .map(|r| (r.3, r.0.clone()))
        .fold((0.0, String::new()), |a, b| if b.0 > a.0 { b } else { a });
    println!(
        "distinct non-folded designs {}; X2 floor mean + 4 sd = {:.1}",
        rows.len(),
        floor
    );
    println!(
        "eta > 1/4: {} designs, {} above the floor, largest X2 {:.1} at {}",
        above.len(),
        above.iter().filter(|r| r.3 > floor).count(),
        top_above,
        top_name
    );
    let mut crossed = Vec::new();
    for (i, d) in all.iter().enumerate() {
        let tag = tags[i].split(' ').next().unwrap();
        if tag == "copy" || tag == "folded" || spectrum(*d).top.eta <= 0.25 {
            continue;
        }
        let over: Vec<(u32, f64)> = (12..=top)
            .map(|k| (k, Stats::at(sq, k, i).chi()))
            .filter(|c| c.1 > floor)
            .collect();
        if let Some(peak) = over.iter().max_by(|a, b| a.1.total_cmp(&b.1)) {
            crossed.push(format!(
                "{} above at {} of the {} checkpoints 2^12..2^{}, peak {:.1} at 2^{}, last 2^{}",
                d.name(),
                over.len(),
                top - 11,
                top,
                peak.1,
                peak.0,
                over.last().unwrap().0
            ));
        }
    }
    println!(
        "eta > 1/4 above the floor at an earlier checkpoint: {}",
        crossed.join("; ")
    );
    let near: Vec<String> = rows
        .iter()
        .filter(|r| (r.3 - floor).abs() <= sd)
        .map(|r| format!("{} {:.1}", r.0, r.3))
        .collect();
    println!("within 1 sd of the floor at 2^{}: {}", top, near.join(", "));
    println!(
        "eta < 1/4: {} designs, {} above the floor; quiet: {}",
        below.len(),
        below.iter().filter(|r| r.3 > floor).count(),
        below
            .iter()
            .filter(|r| r.3 <= floor)
            .map(|r| format!("{} ({:.4}, X2 {:.0})", r.0, r.1, r.3))
            .collect::<Vec<_>>()
            .join(", ")
    );
    println!(
        "spearman over the {} designs: eta vs X2 along the squares {:.3}; D2 vs X2 {:.3}; eta vs |z_B| along the primes {:.3}; eta vs D2 {:.3}",
        rows.len(),
        spearman(&eta, &chi),
        spearman(&d2, &chi),
        spearman(&eta, &zp),
        spearman(&eta, &d2)
    );
    let mut pairs = Vec::new();
    for i in 0..rows.len() {
        for j in 0..rows.len() {
            let (a, b) = (&rows[i], &rows[j]);
            if i != j && (a.1 - b.1).abs() < 1e-4 && a.3 > 10.0 * b.3 {
                pairs.push(format!(
                    "{} / {} (eta {:.4}, X2 {:.0} / {:.0})",
                    a.0, b.0, a.1, a.3, b.3
                ));
            }
        }
    }
    println!(
        "equal eta within 1e-4, X2 apart by 10x: {}",
        pairs.join("; ")
    );
    let mut far: Vec<(f64, String)> = Vec::new();
    for (i, d) in all.iter().enumerate() {
        let m = d.mirror();
        let j = all.iter().position(|&e| e == m).unwrap();
        if m.mask > d.mask && !tags[i].starts_with("copy") && !tags[j].starts_with("copy") {
            let (a, b) = (Stats::at(sq, top, i), Stats::at(sq, top, j));
            let ratio = a.chi().max(b.chi()) / a.chi().min(b.chi());
            far.push((
                ratio,
                format!(
                    "{}~{} X2 {:.0}/{:.0} zB {:.1}/{:.1}",
                    d.name(),
                    m.name(),
                    a.chi(),
                    b.chi(),
                    a.bias(),
                    b.bias()
                ),
            ));
        }
    }
    far.sort_by(|x, y| y.0.total_cmp(&x.0));
    println!(
        "mirror pairs by X2 ratio: {}",
        far.iter()
            .take(6)
            .map(|f| format!("{} ({:.1}x)", f.1, f.0))
            .collect::<Vec<_>>()
            .join("; ")
    );
    let mut inverted = 0;
    let mut total = 0;
    for a in &rows {
        for b in &rows {
            if a.1 < b.1 - 0.01 && a.3.max(b.3) > floor {
                total += 1;
                if b.3 > 2.0 * a.3 {
                    inverted += 1;
                }
            }
        }
    }
    println!(
        "pairs with eta apart by > 0.01 and one X2 above the floor: {}, ranked backwards by a factor 2 in X2: {}",
        total, inverted
    );
    let mut twins = Vec::new();
    for (i, d) in all.iter().enumerate() {
        let sp = spectrum(*d);
        if tags[i].starts_with("twin") && sp.top.eta < 0.25 {
            let (k, m) = (12..=top)
                .map(|k| (k, Stats::at(sq, k, i).bias().abs()))
                .fold((0, 0.0), |a, b| if b.1 > a.1 { b } else { a });
            twins.push(format!(
                "{} orbit {}/{} G {} max {:.2} at 2^{}",
                d.name(),
                sp.top.a,
                sp.top.q,
                if sp.top.q % 4 == 2 { "zero" } else { "nonzero" },
                m,
                k
            ));
        }
    }
    println!(
        "twins with eta < 1/4, max abs(z_B) over 2^12..2^{}: {}",
        top,
        twins.join("; ")
    );
    let open: Vec<usize> = (0..DESIGNS)
        .filter(|&d| status(d, tags) == "open")
        .collect();
    let seen_all = open.iter().all(|&d| Stats::at(sq, top, d).seen == 256);
    let (omax, oname) = open
        .iter()
        .map(|&d| (Stats::at(sq, top, d).chi(), all[d].name()))
        .fold((0.0, String::new()), |a, b| if b.0 > a.0 { b } else { a });
    let high: Vec<usize> = open
        .iter()
        .copied()
        .filter(|&d| spectrum(all[d]).top.eta > 0.25)
        .collect();
    let (hmax, hname) = high
        .iter()
        .map(|&d| (Stats::at(sq, top, d).chi(), all[d].name()))
        .fold((0.0, String::new()), |a, b| if b.0 > a.0 { b } else { a });
    println!(
        "open class, odd base and not folded: {} designs, every 8-block met along n^2 below 2^{}: {}, largest X2 {:.0} at {}; of them {} with eta > 1/4, largest X2 {:.1} at {}; proved class, even base: {} designs",
        open.len(),
        top,
        seen_all,
        omax,
        oname,
        high.len(),
        hmax,
        hname,
        (0..DESIGNS).filter(|&d| status(d, tags) == "proved").count()
    );
    println!(
        "lit lags (|z| > 4) on an orbit with eta < 1/4 through k/(2h): {}; lit elsewhere: {}; predicted and quiet: {}",
        lit_ok,
        if lit_off.is_empty() { "none".to_string() } else { lit_off.join(", ") },
        quiet.len()
    );
    println!();
}

// SPEED

fn decay(grid: &Grid, d: usize, set: usize, ks: &[u32]) -> f64 {
    let pts: Vec<(f64, f64)> = ks
        .iter()
        .map(|&k| {
            let st = Stats::at(grid, k, d);
            (k as f64, ((st.t[set] as f64).abs().max(1.0) / st.n).log2())
        })
        .collect();
    let n = pts.len() as f64;
    let mx = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
    pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum::<f64>()
        / pts.iter().map(|p| (p.0 - mx).powi(2)).sum::<f64>()
}

fn speed(top: u32) {
    let all = Design::all();
    let tags = classify(&all, false);
    let reader = Reader::new();
    let sq = tally(top, |n| reader.word(n * n) | coin(n, SALT_SQ) << COIN);
    let wide: Vec<u32> = (18..=top).collect();
    let late: Vec<u32> = (top - 6..=top).collect();
    let floor = 8.0;
    println!("# SPEED");
    println!();
    println!(
        "the decay exponent of the mean m(N) = T_S/N of a statistic along n^2: least-squares slope of log2 abs(m(2^k)) over every k in 18..{} (fit) and {}..{} (late), against -2 eta_sq for the bias and -2 e_h for the lag h; only statistics with abs(z) >= {} at 2^{}, so the noise moves log2 abs(z) by less than 0.2; the noise alone has -0.5",
        top, top - 6, top, floor, top
    );
    println!();
    println!(
        "{:<7} {:>7} {:<5} {:>8} {:>8} {:>8} {:>8} {:>8}",
        "design", "eta", "stat", "z", "fit", "late", "pred", "gap"
    );
    let mut gaps: Vec<(f64, String)> = Vec::new();
    let mut late_gaps: Vec<f64> = Vec::new();
    for (i, d) in all.iter().enumerate() {
        let tag = tags[i].split(' ').next().unwrap();
        let sp = spectrum(*d);
        if tag == "copy" || tag == "folded" || sp.top.eta >= 0.25 {
            continue;
        }
        let st = Stats::at(&sq, top, i);
        let mut cells: Vec<(String, usize, f64)> = Vec::new();
        if st.bias().abs() >= floor {
            cells.push(("B".into(), 1, -2.0 * sp.sq.eta));
        }
        for h in 1..8 {
            if st.lag(h).abs() >= floor {
                cells.push((format!("C{}", h), 1 | 1 << h, -2.0 * sp.lags[h]));
            }
        }
        for (name, set, pred) in cells {
            let fit = decay(&sq, i, set, &wide);
            let lat = decay(&sq, i, set, &late);
            println!(
                "{:<7} {:>7.4} {:<5} {:>8.2} {:>8.3} {:>8.3} {:>8.3} {:>8.3}",
                d.name(),
                sp.top.eta,
                name,
                st.z(set),
                fit,
                lat,
                pred,
                fit - pred
            );
            gaps.push((fit - pred, format!("{} {}", d.name(), name)));
            late_gaps.push(lat - pred);
        }
    }
    for d in [RS, COIN] {
        let st = Stats::at(&sq, top, d);
        println!(
            "{:<7} {:>7} {:<5} {:>8.2} {:>8.3} {:>8.3} {:>8} {:>8}",
            label(d),
            "",
            "B",
            st.bias(),
            decay(&sq, d, 1, &wide),
            decay(&sq, d, 1, &late),
            "-0.500",
            ""
        );
    }
    println!();
    let n = gaps.len() as f64;
    let mean = gaps.iter().map(|g| g.0).sum::<f64>() / n;
    let mabs = gaps.iter().map(|g| g.0.abs()).sum::<f64>() / n;
    let worst = gaps
        .iter()
        .max_by(|a, b| a.0.abs().total_cmp(&b.0.abs()))
        .unwrap();
    let near = gaps.iter().filter(|g| g.0.abs() <= 0.05).count();
    let lmabs = late_gaps.iter().map(|g| g.abs()).sum::<f64>() / n;
    let mut names: Vec<&str> = gaps
        .iter()
        .map(|g| g.1.split(' ').next().unwrap())
        .collect();
    names.dedup();
    println!(
        "{} statistics fitted in {} designs: mean of fit - pred {:+.3}, mean abs {:.3}, within 0.05: {}, largest {:+.3} at {}; over the late window mean abs {:.3}",
        gaps.len(),
        names.len(),
        mean,
        mabs,
        near,
        worst.0,
        worst.1,
        lmabs
    );
}

// THUE

fn thue_job(s: u64, e: u64) -> [i64; 3] {
    let mut sq = (s as u128) * (s as u128);
    let mut n = s;
    let mut hist = 0u64;
    let mut out = [0i64; 3];
    while n < e + 6 {
        let bit = (((sq as u64) ^ ((sq >> 64) as u64)).count_ones() & 1) as u64;
        hist = (hist << 1) | bit;
        if n < e {
            out[0] += 1 - 2 * bit as i64;
        }
        if n >= s + 3 && n - 3 < e {
            out[1] += 1 - 2 * ((bit ^ (hist >> 3)) & 1) as i64;
        }
        if n >= s + 6 && n - 6 < e {
            out[2] += 1 - 2 * ((bit ^ (hist >> 6)) & 1) as i64;
        }
        sq += 2 * n as u128 + 1;
        n += 1;
    }
    out
}

fn thue_sums(top: u32) -> Vec<[i64; 3]> {
    let chunk: u64 = 1 << 24;
    let blocks = top as usize + 1;
    let mut sums = vec![[0i64; 3]; blocks];
    for n in 0..chunk.min(1 << top) {
        let blk = if n == 0 {
            0
        } else {
            64 - n.leading_zeros() as usize
        };
        let one = thue_job(n, n + 1);
        for t in 0..3 {
            sums[blk][t] += one[t];
        }
    }
    let starts: Vec<u64> = (1..(1u64 << top) / chunk).map(|j| j * chunk).collect();
    let next = AtomicUsize::new(0);
    let total = Mutex::new(vec![[0i64; 3]; blocks]);
    thread::scope(|sc| {
        for _ in 0..workers() {
            sc.spawn(|| {
                let mut mine = vec![[0i64; 3]; blocks];
                loop {
                    let j = next.fetch_add(1, Ordering::Relaxed);
                    if j >= starts.len() {
                        break;
                    }
                    let s0 = starts[j];
                    let blk = 64 - s0.leading_zeros() as usize;
                    let r = thue_job(s0, s0 + chunk);
                    for t in 0..3 {
                        mine[blk][t] += r[t];
                    }
                }
                let mut all = total.lock().unwrap();
                for (a, m) in all.iter_mut().zip(&mine) {
                    for t in 0..3 {
                        a[t] += m[t];
                    }
                }
            });
        }
    });
    for (a, m) in sums.iter_mut().zip(total.into_inner().unwrap()) {
        for t in 0..3 {
            a[t] += m[t];
        }
    }
    let mut cum = vec![[0i64; 3]; blocks];
    let mut run = [0i64; 3];
    for k in 0..blocks {
        for t in 0..3 {
            run[t] += sums[k][t];
        }
        cum[k] = run;
    }
    cum
}

fn thue(top: u32) {
    let cum = thue_sums(top);
    let blocks = top as usize + 1;
    println!("# THUE");
    println!();
    println!(
        "Thue-Morse along n^2 alone, n^2 in 128 bits: T(N) = sum_(n<N) t(n^2), and the lags 3, 6, at every N = 2^k up to 2^{}; z = T/sqrt(N); decay = least-squares slope of log2 abs(T/N) over a window of k, against -log2(4/3) = -0.4150",
        top
    );
    println!();
    println!(
        "{:>4} {:>16} {:>9} {:>9} {:>9}",
        "k", "S(N)", "z_B", "z_C3", "z_C6"
    );
    for k in 12..blocks {
        let sn = (1u64 << k) as f64;
        println!(
            "{:>4} {:>16} {:>9.2} {:>9.2} {:>9.2}",
            k,
            cum[k][0],
            cum[k][0] as f64 / sn.sqrt(),
            cum[k][1] as f64 / sn.sqrt(),
            cum[k][2] as f64 / sn.sqrt()
        );
    }
    let fit = |t: usize, lo: u32, hi: u32| {
        let pts: Vec<(f64, f64)> = (lo..=hi)
            .map(|k| {
                (
                    k as f64,
                    ((cum[k as usize][t] as f64).abs().max(1.0) / (1u64 << k) as f64).log2(),
                )
            })
            .collect();
        let n = pts.len() as f64;
        let mx = pts.iter().map(|p| p.0).sum::<f64>() / n;
        let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
        pts.iter().map(|p| (p.0 - mx) * (p.1 - my)).sum::<f64>()
            / pts.iter().map(|p| (p.0 - mx).powi(2)).sum::<f64>()
    };
    println!();
    for (lo, hi) in [(18, 30.min(top)), (18, top), (30.min(top), top)] {
        if hi > lo {
            println!(
                "decay over 2^{}..2^{}: bias {:.4}, C3 {:.4}, C6 {:.4}",
                lo,
                hi,
                fit(0, lo, hi),
                fit(1, lo, hi),
                fit(2, lo, hi)
            );
        }
    }
}

// ETA

fn eta() {
    println!("# ETA");
    println!();
    println!(
        "eta = -lim (1/L) log_b sup_alpha prod_(j<L) |phi(b^j alpha)|, phi(alpha) = b^-1 sum_(d<b) (-1)^[d in D] e(d alpha); lo from a max-plus certificate on the b-adic cylinders of depth r (b^r <= 2^17), sup of |phi|^2 on a cylinder bounded by the centre value plus a Lipschitz term, every float rounded outward; hi from the certified product along the orbit a/q; both printed outward, lo floored and hi ceiled at 6 places"
    );
    println!();
    println!(
        "{:<7} {:>9} {:>9} {:>9} {:>4} {:>16} {:>6} {:>6} {:>8} {:>16}",
        "design", "eta lo", "eta hi", "width", "k", "orbit a/q", "cycle", "rounds", "D2", "T"
    );
    let mut widest = 0.0f64;
    let floor6 = |x: f64| (((x - 1e-15) * 1e6).floor() / 1e6).max(0.0);
    let ceil6 = |x: f64| ((x + 1e-15) * 1e6).ceil() / 1e6;
    for d in Design::all() {
        let sp = spectrum(d);
        let br = bracket(d, sp.top);
        let (lo, hi) = (floor6(br.lo), ceil6(br.hi));
        widest = widest.max(hi - lo);
        let t = energy(d);
        println!(
            "{:<7} {:>9.6} {:>9.6} {:>9.2e} {:>4} {:>16} {:>6} {:>6} {:>8.5} {:>16}",
            d.name(),
            lo,
            hi,
            hi - lo,
            br.orbit.0,
            format!("{}/{}", br.orbit.1, br.orbit.2),
            br.cycle,
            br.rounds,
            dim2(d),
            format!("{} {} {} {}", t[0][0], t[0][1], t[1][0], t[1][1])
        );
    }
    println!();
    println!("widest bracket {:.2e}", widest);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let clock = Instant::now();
    let num = |i: usize, d: u32| args.get(i).and_then(|s| s.parse().ok()).unwrap_or(d);
    match args.get(1).map(String::as_str) {
        Some("census") => census(num(2, 28), num(3, 24)),
        Some("eta") => eta(),
        Some("diagonal") => diagonal(num(2, 24), num(3, 22)),
        Some("speed") => speed(num(2, 30)),
        Some("thue") => thue(num(2, 34)),
        _ => {
            eprintln!("verbs: census [k p], eta, diagonal [k p], speed [k], thue [k]");
            std::process::exit(2);
        }
    }
    eprintln!("{:.2} s", clock.elapsed().as_secs_f64());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reader_matches_the_digits() {
        let all = Design::all();
        let reader = Reader::new();
        let mut values: Vec<u64> = (0..5000).collect();
        values.extend(
            (0..2000u64)
                .map(|i| (1u64 << 30) - 1 - 977 * i)
                .map(|n| n * n),
        );
        values.extend((0..2000u64).map(|i| mix(i) >> 4));
        for v in values {
            let w = reader.word(v);
            for (i, d) in all.iter().enumerate() {
                assert_eq!((w >> i) & 1, d.parity(v), "{} at {}", d.name(), v);
            }
            assert_eq!((w >> RS) & 1, rudin(v));
        }
    }

    #[test]
    fn walsh_chi_square_is_pearson() {
        let mut c = [0u64; 256];
        for (i, x) in c.iter_mut().enumerate() {
            *x = mix(i as u64) % 1000;
        }
        let n: u64 = c.iter().sum();
        let s = Stats {
            n: n as f64,
            t: walsh(&c),
            seen: 256,
        };
        let e = n as f64 / 256.0;
        let pearson: f64 = c.iter().map(|&x| (x as f64 - e).powi(2) / e).sum();
        assert!((s.chi() - pearson).abs() < 1e-6 * pearson);
        assert_eq!(null().0, 255.0);
    }

    #[test]
    fn tally_counts_every_block_once() {
        let grid = tally(17, |n| mix(n));
        for d in [0, 30, STREAMS - 1] {
            let mut c = [0u64; 256];
            for n in 0..1u64 << 17 {
                let w = (0..8).fold(0usize, |w, i| w | (((mix(n + i) >> d) & 1) as usize) << i);
                c[w] += 1;
            }
            let mut got = [0u64; 256];
            for blk in 0..=17 {
                for (x, y) in got.iter_mut().zip(&grid[blk * STREAMS + d]) {
                    *x += y;
                }
            }
            assert_eq!(got, c);
        }
    }

    #[test]
    fn intervals_hold_the_logarithms() {
        for x in [1e-9, 0.01, 0.4330127, 0.75, 1.0, 1.5, 2.0, 6.0, 1e6] {
            let iv = ln_iv(Iv::pt(x));
            assert!(
                iv.lo <= x.ln() + 1e-15 && x.ln() - 1e-15 <= iv.hi && iv.hi - iv.lo < 1e-12,
                "{} {:?}",
                x,
                iv
            );
        }
    }

    #[test]
    fn intervals_hold_the_cosines() {
        for q in [3i128, 7, 12, 1000, 59049, 1 << 20] {
            for p in [0, 1, q / 8, q / 4, q / 3, q / 2, q - 1, 5 * q + 7] {
                let iv = cos2pi(p, q);
                let x = (2.0 * std::f64::consts::PI * (p.rem_euclid(q)) as f64 / q as f64).cos();
                assert!(iv.lo <= x + 1e-15 && x - 1e-15 <= iv.hi && iv.hi - iv.lo < 1e-14);
            }
        }
    }

    #[test]
    fn digit_walk_sums_the_autocorrelation() {
        let tm = Corr::of(Design { b: 2, mask: 1 });
        let t3 = Corr::of(Design { b: 3, mask: 2 });
        assert!((tm.at(1) + 1.0 / 3.0).abs() < 1e-15 && (tm.at(3) - 1.0 / 3.0).abs() < 1e-15);
        assert!(t3.at(1).abs() < 1e-15 && (t3.at(2) + 1.0 / 3.0).abs() < 1e-15);
        for d in [
            Design { b: 2, mask: 1 },
            Design { b: 5, mask: 6 },
            Design { b: 6, mask: 6 },
        ] {
            let corr = Corr::of(d);
            let len = 1u64 << 16;
            for g in [1u64, 2, 7, 30, 1000] {
                let emp = (0..len)
                    .map(|m| 1.0 - 2.0 * (d.parity(m) ^ d.parity(m + g)) as f64)
                    .sum::<f64>()
                    / len as f64;
                assert!(
                    (emp - corr.at(g)).abs() < 0.02,
                    "{} {} {} {}",
                    d.name(),
                    g,
                    emp,
                    corr.at(g)
                );
            }
            for (n, a, c) in [(1000u64, 6u64, 9u64), (4097, 14, 49), (777, 2, 1)] {
                let direct: f64 = (0..n).map(|k| corr.at(a * k + c)).sum();
                assert!((corr.progression(n, a, c) - direct).abs() < 1e-9 * n as f64);
            }
        }
    }

    #[test]
    fn energy_recursion_matches_the_blocks() {
        for d in [
            Design { b: 2, mask: 1 },
            Design { b: 5, mask: 6 },
            Design { b: 6, mask: 6 },
        ] {
            let t = energy(d);
            let (mut e, mut r) = (1i128, 0i128);
            for l in 1..=4u32 {
                (e, r) = (
                    t[0][0] as i128 * e + t[0][1] as i128 * r,
                    t[1][0] as i128 * e + t[1][1] as i128 * r,
                );
                let len = d.b.pow(l) as i64;
                let x: Vec<i128> = (0..len)
                    .map(|m| 1 - 2 * d.parity(m as u64) as i128)
                    .collect();
                let a: Vec<i128> = (-len..=len)
                    .map(|h| {
                        (0..len)
                            .filter(|&m| (0..len).contains(&(m + h)))
                            .map(|m| x[m as usize] * x[(m + h) as usize])
                            .sum()
                    })
                    .collect();
                let de: i128 = a.iter().map(|v| v * v).sum();
                let dr: i128 = a.windows(2).map(|w| w[0] * w[1]).sum();
                assert_eq!((e, r), (de, dr), "{} level {}", d.name(), l);
            }
        }
        assert!(
            (dim2(Design { b: 2, mask: 1 }) - (3.0 - (1.0 + 17f64.sqrt()).log2())).abs() < 1e-12
        );
    }

    #[test]
    fn mirrors_share_the_spectrum() {
        for d in Design::all() {
            let m = d.mirror();
            assert_eq!(m.mirror(), d);
            assert_eq!(Tile::of(d).a, Tile::of(m).a, "{}", d.name());
            let (cd, cm) = (Corr::of(d), Corr::of(m));
            for g in 0..500 {
                assert!((cd.at(g) - cm.at(g)).abs() < 1e-12);
            }
            for len in [5u32, 6] {
                for n in 0..d.b.pow(len) {
                    let top = if d.sign(d.b - 1) == -1 {
                        len as u64 % 2
                    } else {
                        0
                    };
                    assert_eq!(m.parity(n), d.parity(d.b.pow(len) - 1 - n) ^ top);
                }
            }
        }
    }

    #[test]
    fn wide_squares_match_the_tally() {
        let reader = Reader::new();
        let grid = tally(22, |n| reader.word(n * n));
        let cum = thue_sums(22);
        for k in [12u32, 17, 22] {
            let st = Stats::at(&grid, k, 0);
            assert_eq!(
                [st.t[1], st.t[1 | 1 << 3], st.t[1 | 1 << 6]],
                cum[k as usize]
            );
        }
    }

    #[test]
    fn brackets_hold_the_closed_forms() {
        let tm = Design { b: 2, mask: 1 };
        let t3 = Design { b: 3, mask: 2 };
        for (d, exact) in [
            (tm, (2.0 / 3f64.sqrt()).log2()),
            (t3, (3.0 / 5f64.sqrt()).ln() / 3f64.ln()),
        ] {
            let br = bracket(d, spectrum(d).top);
            assert!(
                br.lo <= exact && exact <= br.hi && br.hi - br.lo < 1e-3,
                "{} {} {}",
                d.name(),
                br.lo,
                br.hi
            );
        }
    }
}
