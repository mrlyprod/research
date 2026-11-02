use std::env;
use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Mutex, OnceLock};
use std::thread;
use std::time::Instant;

// DESIGNS

const DESIGNS: usize = 57;
const RS: usize = 57;
const COIN: usize = 58;
const STREAMS: usize = 59;
const GOLDEN: u64 = 0x9e37_79b9_7f4a_7c15;
const SALT: u64 = 3;

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

    fn tile(self) -> Vec<f64> {
        (0..self.b).map(|d| self.sign(d) as f64).collect()
    }

    fn name(self) -> String {
        let digits: String = (1..self.b)
            .filter(|&d| self.has(d))
            .map(|d| char::from(b'0' + d as u8))
            .collect();
        format!("{}:{}", self.b, digits)
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

fn find(name: &str) -> usize {
    (0..STREAMS).find(|&i| label(i) == name).unwrap()
}

fn mix(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

fn coin(n: u64) -> u64 {
    mix(SALT.wrapping_add(n.wrapping_mul(GOLDEN))) >> 63
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
            | coin(v) << COIN
    }
}

fn sign(w: u64, i: usize) -> i64 {
    1 - 2 * ((w >> i) & 1) as i64
}

// CENSUS

const ZS: [u64; 10] = [7, 11, 13, 17, 31, 47, 73, 127, 257, 1021];
const BINS: usize = ZS.len();
const SPAN: u64 = 1 << 20;
const MARK: u64 = 1_339_484_198;

fn workers() -> usize {
    thread::available_parallelism()
        .map_or(4, |n| n.get())
        .min(8)
}

fn odd_primes(limit: u64) -> Vec<u64> {
    let mut is = vec![true; limit as usize + 1];
    let mut out = Vec::new();
    for i in 2..=limit as usize {
        if is[i] {
            if i > 2 {
                out.push(i as u64);
            }
            let mut j = i * i;
            while j <= limit as usize {
                is[j] = false;
                j += i;
            }
        }
    }
    out
}

fn level(p: u64) -> u8 {
    ZS.iter().filter(|&&z| p > z).count() as u8
}

fn ratio(z: u64) -> f64 {
    let mut r = 1.0;
    for p in 2..=z {
        if (2..p).take_while(|d| d * d <= p).all(|d| p % d != 0) {
            r *= p as f64 / (p - 1) as f64;
        }
    }
    r
}

fn ratios() -> &'static [f64; BINS] {
    static R: OnceLock<[f64; BINS]> = OnceLock::new();
    R.get_or_init(|| ZS.map(ratio))
}

fn checkpoints(top: u32) -> Vec<u64> {
    let mut xs: Vec<u64> = (64..=16 * top)
        .map(|c| 2f64.powf(c as f64 / 16.0).round() as u64)
        .collect();
    if MARK < 1 << top {
        xs.push(MARK);
    }
    xs.sort();
    xs.dedup();
    xs
}

#[derive(Clone)]
struct Acc {
    primes: u64,
    theta: f64,
    lam2: f64,
    pi: [i64; STREAMS],
    psi: [f64; STREAMS],
    mc: [i64; BINS],
    mu: [f64; BINS],
    m: [[i64; STREAMS]; BINS],
    u: [[f64; STREAMS]; BINS],
}

impl Acc {
    fn new() -> Acc {
        Acc {
            primes: 0,
            theta: 0.0,
            lam2: 0.0,
            pi: [0; STREAMS],
            psi: [0.0; STREAMS],
            mc: [0; BINS],
            mu: [0.0; BINS],
            m: [[0; STREAMS]; BINS],
            u: [[0.0; STREAMS]; BINS],
        }
    }

    fn add(&mut self, o: &Acc) {
        self.primes += o.primes;
        self.theta += o.theta;
        self.lam2 += o.lam2;
        for s in 0..STREAMS {
            self.pi[s] += o.pi[s];
            self.psi[s] += o.psi[s];
        }
        for b in 0..BINS {
            self.mc[b] += o.mc[b];
            self.mu[b] += o.mu[b];
            for s in 0..STREAMS {
                self.m[b][s] += o.m[b][s];
                self.u[b][s] += o.u[b][s];
            }
        }
    }

    fn sigma(&self) -> f64 {
        self.lam2.sqrt()
    }

    fn model(&self, lvl: usize, s: usize) -> f64 {
        ratios()[lvl] * (lvl..BINS).map(|b| self.m[b][s] as f64).sum::<f64>()
    }

    fn line(&self, x: u64) -> String {
        let mut out = format!("{} {} {} {}", x, self.primes, self.theta, self.lam2);
        for s in 0..STREAMS {
            out += &format!(" {} {}", self.pi[s], self.psi[s]);
        }
        for b in 0..BINS {
            out += &format!(" {} {}", self.mc[b], self.mu[b]);
            for s in 0..STREAMS {
                out += &format!(" {} {}", self.m[b][s], self.u[b][s]);
            }
        }
        out
    }

    fn parse(line: &str) -> (u64, Acc) {
        let mut t = line.split_whitespace();
        let mut int = || t.next().unwrap().parse::<f64>().unwrap();
        let x = int() as u64;
        let mut a = Acc::new();
        a.primes = int() as u64;
        a.theta = int();
        a.lam2 = int();
        for s in 0..STREAMS {
            a.pi[s] = int() as i64;
            a.psi[s] = int();
        }
        for b in 0..BINS {
            a.mc[b] = int() as i64;
            a.mu[b] = int();
            for s in 0..STREAMS {
                a.m[b][s] = int() as i64;
                a.u[b][s] = int();
            }
        }
        (x, a)
    }
}

struct Scratch {
    comp: Vec<u8>,
    low: Vec<u8>,
    mc: Vec<u32>,
    ms: Vec<f64>,
    pc: Vec<u32>,
    ps: Vec<f64>,
}

impl Scratch {
    fn new() -> Scratch {
        let half = (SPAN / 2 + 1) as usize;
        Scratch {
            comp: vec![0; half],
            low: vec![0; half],
            mc: vec![0; BINS * 2048],
            ms: vec![0.0; BINS * 2048],
            pc: vec![0; 2048],
            ps: vec![0.0; 2048],
        }
    }
}

fn spread(
    cnt: &mut [u32],
    sum: &mut [f64],
    ones: &mut [i64; STREAMS],
    weight: &mut [f64; STREAMS],
) {
    for k in 0..8 {
        for byte in 0..256 {
            let c = cnt[k * 256 + byte];
            if c == 0 {
                continue;
            }
            let s = sum[k * 256 + byte];
            for bit in 0..8 {
                let st = 8 * k + bit;
                if (byte >> bit) & 1 == 1 && st < STREAMS {
                    ones[st] += c as i64;
                    weight[st] += s;
                }
            }
            cnt[k * 256 + byte] = 0;
            sum[k * 256 + byte] = 0.0;
        }
    }
}

fn segment(
    reader: &Reader,
    base: &[u64],
    small: &[(u64, u8)],
    lo: u64,
    hi: u64,
    sc: &mut Scratch,
) -> Acc {
    let mut acc = Acc::new();
    let first = lo.max(11) | 1;
    if first >= hi {
        return acc;
    }
    let len = (hi - first).div_ceil(2) as usize;
    let comp = &mut sc.comp[..len];
    let low = &mut sc.low[..len];
    comp.fill(0);
    for &p in base {
        if p * p >= hi {
            break;
        }
        let mut start = (p * p).max(first.div_ceil(p) * p);
        if start % 2 == 0 {
            start += p;
        }
        let mut i = ((start - first) / 2) as usize;
        while i < len {
            comp[i] = 1;
            i += p as usize;
        }
    }
    low.fill(BINS as u8);
    for &(p, lv) in small {
        let mut start = first.div_ceil(p) * p;
        if start.is_multiple_of(2) {
            start += p;
        }
        let mut i = ((start - first) / 2) as usize;
        while i < len {
            low[i] = lv;
            i += p as usize;
        }
    }
    for i in 0..len {
        let lv = low[i];
        if lv == 0 {
            continue;
        }
        let n = first + 2 * i as u64;
        let w = reader.word(n);
        let lnn = (n as f64).ln();
        let inv = 1.0 / lnn;
        let bin = (lv - 1) as usize;
        acc.mc[bin] += 1;
        acc.mu[bin] += inv;
        let at = bin * 2048;
        for k in 0..8 {
            let idx = at + k * 256 + ((w >> (8 * k)) & 255) as usize;
            sc.mc[idx] += 1;
            sc.ms[idx] += inv;
        }
        if comp[i] == 0 {
            acc.primes += 1;
            acc.theta += lnn;
            acc.lam2 += lnn * lnn;
            for k in 0..8 {
                let idx = k * 256 + ((w >> (8 * k)) & 255) as usize;
                sc.pc[idx] += 1;
                sc.ps[idx] += lnn;
            }
        }
    }
    let (mut ones, mut weight) = ([0i64; STREAMS], [0f64; STREAMS]);
    spread(&mut sc.pc, &mut sc.ps, &mut ones, &mut weight);
    for s in 0..STREAMS {
        acc.pi[s] = acc.primes as i64 - 2 * ones[s];
        acc.psi[s] = acc.theta - 2.0 * weight[s];
    }
    for b in 0..BINS {
        let (mut ones, mut weight) = ([0i64; STREAMS], [0f64; STREAMS]);
        spread(
            &mut sc.mc[b * 2048..(b + 1) * 2048],
            &mut sc.ms[b * 2048..(b + 1) * 2048],
            &mut ones,
            &mut weight,
        );
        for s in 0..STREAMS {
            acc.m[b][s] = acc.mc[b] - 2 * ones[s];
            acc.u[b][s] = acc.mu[b] - 2.0 * weight[s];
        }
    }
    acc
}

fn extras(reader: &Reader, top: u64) -> Vec<(u64, Acc)> {
    let mut out = Vec::new();
    let mut one = Acc::new();
    let w = reader.word(1);
    one.mc[BINS - 1] = 1;
    for s in 0..STREAMS {
        one.m[BINS - 1][s] = sign(w, s);
    }
    out.push((1, one));
    let mut ps = vec![2u64];
    ps.extend(odd_primes(65_536));
    for &p in &ps {
        let l = (p as f64).ln();
        let mut n = p;
        let mut k = 1;
        while n < top {
            let w = reader.word(n);
            let mut a = Acc::new();
            a.lam2 = l * l;
            for s in 0..STREAMS {
                a.psi[s] = sign(w, s) as f64 * l;
            }
            if k == 1 && p < 11 {
                a.primes = 1;
                a.theta = l;
                for s in 0..STREAMS {
                    a.pi[s] = sign(w, s);
                }
            }
            if k > 1 || p < 11 {
                out.push((n, a));
            }
            n *= p;
            k += 1;
        }
    }
    out
}

fn census(top: u32) -> Vec<(u64, Acc)> {
    let clock = Instant::now();
    let xs = checkpoints(top);
    let reader = Reader::new();
    let base = odd_primes(1 << 16);
    let mut small: Vec<(u64, u8)> = odd_primes(1021)
        .into_iter()
        .map(|p| (p, if p <= 7 { 0 } else { level(p) }))
        .collect();
    small.reverse();
    let mut tasks = Vec::new();
    let mut lo = 0;
    for (j, &x) in xs.iter().enumerate() {
        while lo < x {
            let hi = (lo + SPAN).min(x);
            tasks.push((j, lo, hi));
            lo = hi;
        }
    }
    let next = AtomicUsize::new(0);
    let found = Mutex::new(Vec::new());
    thread::scope(|s| {
        for _ in 0..workers() {
            s.spawn(|| {
                let mut sc = Scratch::new();
                loop {
                    let t = next.fetch_add(1, Ordering::Relaxed);
                    if t >= tasks.len() {
                        break;
                    }
                    let (_, lo, hi) = tasks[t];
                    let acc = segment(&reader, &base, &small, lo, hi, &mut sc);
                    found.lock().unwrap().push((t, acc));
                }
            });
        }
    });
    let mut found = found.into_inner().unwrap();
    found.sort_by_key(|f| f.0);
    let mut part = vec![Acc::new(); xs.len()];
    for (t, acc) in &found {
        part[tasks[*t].0].add(acc);
    }
    for (n, acc) in extras(&reader, 1 << top) {
        part[xs.partition_point(|&x| x <= n)].add(&acc);
    }
    let mut rows = Vec::with_capacity(xs.len());
    let mut run = Acc::new();
    for (j, &x) in xs.iter().enumerate() {
        run.add(&part[j]);
        rows.push((x, run.clone()));
    }
    eprintln!(
        "census 2^{}: {} checkpoints, {} tasks, {:.1} s",
        top,
        xs.len(),
        tasks.len(),
        clock.elapsed().as_secs_f64()
    );
    rows
}

fn data_dir() -> PathBuf {
    let here = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let repo = here.ancestors().nth(5).unwrap().to_path_buf();
    repo.join("data").join(here.strip_prefix(&repo).unwrap())
}

fn cached(top: u32) -> Vec<(u64, Acc)> {
    let file = data_dir().join(format!("census-{}.txt", top));
    if let Ok(text) = fs::read_to_string(&file) {
        return text.lines().map(Acc::parse).collect();
    }
    let rows = census(top);
    fs::create_dir_all(data_dir()).unwrap();
    let text: String = rows.iter().map(|(x, a)| a.line(*x) + "\n").collect();
    fs::write(&file, text).unwrap();
    rows
}

fn at(rows: &[(u64, Acc)], x: u64) -> &Acc {
    &rows.iter().find(|r| r.0 == x).unwrap().1
}

// FOURIER

#[derive(Clone, Copy, Debug, PartialEq)]
struct Cx {
    re: f64,
    im: f64,
}

impl Cx {
    const ZERO: Cx = Cx { re: 0.0, im: 0.0 };
    const ONE: Cx = Cx { re: 1.0, im: 0.0 };

    fn turn(r: u64, q: u64) -> Cx {
        let t = 2.0 * std::f64::consts::PI * r as f64 / q as f64;
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
}

fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn factor(mut q: u64) -> Vec<(u64, u32)> {
    let mut out = Vec::new();
    let mut p = 2;
    while p * p <= q {
        if q.is_multiple_of(p) {
            let mut e = 0;
            while q.is_multiple_of(p) {
                q /= p;
                e += 1;
            }
            out.push((p, e));
        }
        p += 1;
    }
    if q > 1 {
        out.push((q, 1));
    }
    out
}

fn mobius(q: u64) -> i64 {
    let f = factor(q);
    if f.iter().any(|&(_, e)| e > 1) {
        0
    } else if f.len().is_multiple_of(2) {
        1
    } else {
        -1
    }
}

fn totient(q: u64) -> u64 {
    factor(q).iter().fold(q, |t, &(p, _)| t / p * (p - 1))
}

fn digits(mut x: u64, b: u64) -> Vec<u64> {
    let mut out = Vec::new();
    while x > 0 {
        out.push(x % b);
        x /= b;
    }
    out
}

fn sums(tiles: &[Vec<f64>], b: usize, a: u64, q: u64, tab: &[Cx], dig: &[u64], out: &mut [Cx]) {
    let l = dig.len();
    let mut es = [[Cx::ZERO; 6]; 40];
    let mut r = a % q;
    for row in es.iter_mut().take(l) {
        let mut idx = 0u64;
        for e in row.iter_mut().take(b) {
            *e = tab[idx as usize];
            idx += r;
            if idx >= q {
                idx -= q;
            }
        }
        r = r * b as u64 % q;
    }
    for (t, tile) in tiles.iter().enumerate() {
        let mut pows = [Cx::ONE; 41];
        for j in 0..l {
            let mut p = Cx::ZERO;
            for d in 0..b {
                p = p.add(es[j][d].scale(tile[d]));
            }
            pows[j + 1] = pows[j].mul(p);
        }
        let mut s = Cx::ZERO;
        let mut f = Cx::ONE;
        for j in (0..l).rev() {
            let xj = dig[j] as usize;
            let mut part = Cx::ZERO;
            for d in 0..xj {
                part = part.add(es[j][d].scale(tile[d]));
            }
            s = s.add(f.mul(part).mul(pows[j]));
            f = f.mul(es[j][xj].scale(tile[xj]));
        }
        out[t] = s;
    }
}

fn rudin_sum(a: u64, q: u64, tab: &[Cx], dig: &[u64]) -> Cx {
    let l = dig.len();
    let mut e = [Cx::ONE; 40];
    let mut r = a % q;
    for v in e.iter_mut().take(l) {
        *v = tab[r as usize];
        r = r * 2 % q;
    }
    let mut a0 = [Cx::ONE; 41];
    let mut a1 = [Cx::ONE; 41];
    for j in 1..=l {
        let t = e[j - 1].mul(a1[j - 1]);
        a0[j] = a0[j - 1].add(t);
        a1[j] = a0[j - 1].add(t.scale(-1.0));
    }
    let mut s = Cx::ZERO;
    let mut f = Cx::ONE;
    let mut above = 0;
    for j in (0..l).rev() {
        if dig[j] == 1 {
            s = s.add(f.mul(a0[j]));
            f = f.mul(e[j]).scale(if above == 1 { -1.0 } else { 1.0 });
            above = 1;
        } else {
            above = 0;
        }
    }
    s
}

fn terms(q: u64, xs: &[u64]) -> Vec<Vec<f64>> {
    let all = Design::all();
    let mut out = vec![vec![0.0; COIN]; xs.len()];
    let mu = mobius(q);
    if mu == 0 {
        return out;
    }
    let w = mu as f64 / totient(q) as f64;
    let tab: Vec<Cx> = (0..q).map(|r| Cx::turn(r, q)).collect();
    let units: Vec<u64> = (0..q).filter(|&a| 2 * a <= q && gcd(a, q) == 1).collect();
    for b in 2..=6u64 {
        let ids: Vec<usize> = (0..DESIGNS).filter(|&i| all[i].b == b).collect();
        let tiles: Vec<Vec<f64>> = ids.iter().map(|&i| all[i].tile()).collect();
        let mut s = vec![Cx::ZERO; ids.len()];
        for (k, &x) in xs.iter().enumerate() {
            let dig = digits(x, b);
            for &a in &units {
                let m = if q <= 2 || 2 * a == q { w } else { 2.0 * w };
                sums(&tiles, b as usize, a, q, &tab, &dig, &mut s);
                for (j, &i) in ids.iter().enumerate() {
                    out[k][i] += m * s[j].re;
                }
                if b == 2 {
                    out[k][RS] += m * rudin_sum(a, q, &tab, &dig).re;
                }
            }
        }
    }
    out
}

fn squarefree(limit: u64) -> Vec<u64> {
    (1..=limit).filter(|&q| mobius(q) != 0).collect()
}

fn parallel_terms(qs: &[u64], xs: &[u64]) -> Vec<Vec<Vec<f64>>> {
    let next = AtomicUsize::new(0);
    let found = Mutex::new(vec![Vec::new(); qs.len()]);
    thread::scope(|s| {
        for _ in 0..workers() {
            s.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= qs.len() {
                    break;
                }
                let t = terms(qs[i], xs);
                found.lock().unwrap()[i] = t;
            });
        }
    });
    found.into_inner().unwrap()
}

fn model(qmax: u64, xs: &[u64]) -> Vec<Vec<f64>> {
    let qs = squarefree(qmax);
    let next = AtomicUsize::new(0);
    let parts = Mutex::new(Vec::new());
    thread::scope(|s| {
        for _ in 0..workers() {
            s.spawn(|| {
                let mut mine = vec![vec![0.0; COIN]; xs.len()];
                loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= qs.len() {
                        break;
                    }
                    for (m, v) in mine.iter_mut().zip(terms(qs[i], xs)) {
                        for (a, b) in m.iter_mut().zip(v) {
                            *a += b;
                        }
                    }
                }
                parts.lock().unwrap().push(mine);
            });
        }
    });
    let mut parts = parts.into_inner().unwrap();
    parts.sort_by(|p, q| p[0][0].total_cmp(&q[0][0]));
    let mut total = vec![vec![0.0; COIN]; xs.len()];
    for p in parts {
        for (t, v) in total.iter_mut().zip(p) {
            for (a, b) in t.iter_mut().zip(v) {
                *a += b;
            }
        }
    }
    total
}

fn grid(top: u32) -> Vec<u64> {
    let mut xs: Vec<u64> = (32..=4 * top)
        .map(|c| 2f64.powf(c as f64 / 4.0).round() as u64)
        .collect();
    if MARK < 1 << top {
        xs.push(MARK);
    }
    xs.sort();
    xs
}

fn flat(xs: &[u64], m: &[f64]) -> Vec<f64> {
    let w = |u: f64| 2f64.powf(u) / (u * u * std::f64::consts::LN_2);
    let mut out = Vec::with_capacity(xs.len());
    let mut integral = 0.0;
    for i in 0..xs.len() {
        if i > 0 {
            let (u0, u1) = ((xs[i - 1] as f64).log2(), (xs[i] as f64).log2());
            let (m0, m1) = (m[i - 1] / xs[i - 1] as f64, m[i] / xs[i] as f64);
            let steps = 32;
            let h = (u1 - u0) / steps as f64;
            for j in 0..=steps {
                let c = if j == 0 || j == steps {
                    1.0
                } else if j % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                let sfrac = j as f64 / steps as f64;
                integral += c * h / 3.0 * (m0 * (1.0 - sfrac) + m1 * sfrac) * w(u0 + j as f64 * h);
            }
        }
        out.push(m[i] / (xs[i] as f64).ln() + integral + m[0] / (xs[0] as f64).ln());
    }
    out
}

// VERBS

const Q: u64 = 1024;
const FROM: u64 = 1 << 20;

fn summary(top: u32) {
    let rows = census(top);
    fs::create_dir_all(data_dir()).unwrap();
    let text: String = rows.iter().map(|(x, a)| a.line(*x) + "\n").collect();
    fs::write(data_dir().join(format!("census-{}.txt", top)), text).unwrap();
    let (x, a) = rows.last().unwrap();
    println!("# CENSUS");
    println!();
    println!(
        "n < x = {}: {} checkpoints round(2^(c/16)) and {}; pi(x) = {}, theta(x) = {:.1}, sum Lambda^2 = {:.1}, sigma = {:.1}",
        x,
        rows.len(),
        MARK,
        a.primes,
        a.theta,
        a.lam2,
        a.sigma()
    );
    if let Some((_, m)) = rows.iter().find(|r| r.0 == MARK) {
        println!();
        println!(
            "along the first {} primes (p < {}), z_B = sum x(p)/sqrt(pi):",
            m.primes, MARK
        );
        for name in ["2:1", "6:2", "6:3", "6:4", "6:23", "6:1", "6:5"] {
            let s = find(name);
            print!(" {} {:.1}", name, m.pi[s] as f64 / (m.primes as f64).sqrt());
        }
        println!();
    }
}

fn pick(qs: &[u64], t: &[Vec<Vec<f64>>], xi: usize, d: usize, keep: impl Fn(u64) -> bool) -> f64 {
    qs.iter()
        .zip(t)
        .filter(|(q, _)| keep(**q))
        .map(|(_, v)| v[xi][d])
        .sum()
}

fn thue(top: u32) {
    let rows = cached(top);
    let tm = 0;
    let theta = 3f64.log2() / 2.0;
    let ks: Vec<u32> = (12..=top).collect();
    let xs: Vec<u64> = ks.iter().map(|&k| 1u64 << k).collect();
    let qs: Vec<u64> = squarefree(210)
        .into_iter()
        .filter(|q| 210 % q == 0)
        .collect();
    let t = parallel_terms(&qs, &xs);
    let full = model(Q, &xs);
    let mut small = (0f64, 0f64);
    println!("# THUE-MORSE ALONG THE PRIMES");
    println!();
    println!(
        "Psi(x) = sum_(n<x) t(n) Lambda(n) at x = 2^k, sigma = sqrt(sum_(n<x) Lambda(n)^2), z = Psi/sigma; T_q = (mu(q)/phi(q)) sum_((a,q)=1) S(a/q; x) in units of sigma; M_Q = sum_(q<=Q) T_q, Q = {}; r_X = (Psi - X)/sigma for X = T_3 alone, the q | 6 part, the q | 30 part and M_Q; G = -M_Q/x^(log_4 3); flat: r of sum_(p<x) t(p) against the model weighted by 1/log n, in units of sqrt(pi(x))",
        Q
    );
    println!();
    println!(
        "{:>3} {:>12} {:>8} {:>8} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7}",
        "k", "Psi", "z", "T3", "T6", "T5", "T10", "rest", "r_T3", "r_6", "r_30", "r_Q", "G_obs", "G", "flat"
    );
    let grid = grid(top);
    let gm = model(Q, &grid);
    let fl = flat(&grid, &gm.iter().map(|v| v[tm]).collect::<Vec<f64>>());
    let mut worst = (0f64, 0f64);
    let mut domain = (0f64, 0f64);
    for (i, &k) in ks.iter().enumerate() {
        let a = at(&rows, xs[i]);
        let sg = a.sigma();
        let psi = a.psi[tm];
        let q = |qq: u64| pick(&qs, &t, i, tm, |x| x == qq);
        let m6 = pick(&qs, &t, i, tm, |x| 6 % x == 0);
        let m30 = pick(&qs, &t, i, tm, |x| 30 % x == 0);
        small.0 = small.0.max((q(15) + q(30)).abs() / sg);
        small.1 = small.1.max(pick(&qs, &t, i, tm, |x| x % 7 == 0).abs() / sg);
        let mq = full[i][tm];
        let m = (k / 4) as i32;
        let five = -5f64.powi(m) * [1.0, 0.5, 0.0, 2.5][(k % 4) as usize];
        worst.0 = worst
            .0
            .max((m6 + 3f64.powi(k as i32 / 2)).abs() / 3f64.powi(k as i32 / 2));
        worst.1 = worst.1.max((q(5) + q(10) - five).abs() / 5f64.powi(m));
        let gi = grid.iter().position(|&g| g == xs[i]).unwrap();
        let flat_r = (a.pi[tm] as f64 - fl[gi]) / (a.primes as f64).sqrt();
        if (Q * Q) as f64 <= xs[i] as f64 {
            domain.0 = domain.0.max((psi - mq).abs() / sg);
        } else {
            domain.1 = domain.1.max((psi - mq).abs() / sg);
        }
        println!(
            "{:>3} {:>12.0} {:>8.2} {:>8.2} {:>7.2} {:>7.2} {:>7.2} {:>7.2} {:>7.2} {:>7.2} {:>7.2} {:>7.2} {:>7.4} {:>7.4} {:>7.2}",
            k,
            psi,
            psi / sg,
            q(3) / sg,
            q(6) / sg,
            q(5) / sg,
            q(10) / sg,
            (mq - m30) / sg,
            (psi - q(3)) / sg,
            (psi - m6) / sg,
            (psi - m30) / sg,
            (psi - mq) / sg,
            -psi / (xs[i] as f64).powf(theta),
            -mq / (xs[i] as f64).powf(theta),
            flat_r
        );
    }
    println!();
    println!(
        "exact at 2^k, k = 12..{}: the q | 6 part is -3^floor(k/2) to relative {:.1e}; T_5 + T_10 is -5^floor(k/4) (1, 1/2, 0, 5/2)[k mod 4] to relative {:.1e}",
        top, worst.0, worst.1
    );
    println!(
        "over k = 12..{}: |T_15 + T_30| <= {:.1e} sigma, and the eight terms q | 210 with 7 | q sum to at most {:.4} sigma",
        top, small.0, small.1
    );
    println!(
        "max |r_Q| over the k with Q <= sqrt(x), k = 20..{}: {:.2}; over k = 12..19, where Q > sqrt(x): {:.2}",
        top, domain.0, domain.1
    );
    let top_a = 1usize << 22;
    let mut a_small = vec![0i64; top_a + 1];
    let (mut run, mut lo4, mut hi4, mut off) = (0i64, i64::MAX, i64::MIN, 0usize);
    for n in 0..=4 * top_a {
        if n % 4 == 0 && n > 0 {
            let d = 3 * run - 9 * a_small[n / 4];
            lo4 = lo4.min(d);
            hi4 = hi4.max(d);
            off += (d != 0) as usize;
        }
        if n <= top_a {
            a_small[n] = run;
        }
        if n % 2 == 1 && n % 3 != 0 {
            run += 1 - 2 * (n.count_ones() & 1) as i64;
        }
    }
    println!(
        "the q | 6 part A(x) = 3 sum_(n<x, (n,6)=1) t(n) has A(4x) - 3 A(x) in [{}, {}] for x = 1..2^22, nonzero at {} of the {} x; only T_3 scales exactly",
        lo4, hi4, off, top_a
    );
    let fine: Vec<u64> = (0..2048)
        .map(|i| 2f64.powf(30.0 + i as f64 / 1024.0).round() as u64)
        .collect();
    let tf = parallel_terms(&[1, 2, 3, 6], &fine);
    let gs: Vec<f64> = (0..fine.len())
        .map(|i| -(0..4).map(|j| tf[j][i][tm]).sum::<f64>() / (fine[i] as f64).powf(theta))
        .collect();
    let lo = gs.iter().copied().fold(f64::INFINITY, f64::min);
    let hi = gs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let mean = gs.iter().sum::<f64>() / gs.len() as f64;
    let (ilo, ihi) = (
        gs.iter().position(|&g| g == lo).unwrap(),
        gs.iter().position(|&g| g == hi).unwrap(),
    );
    let centre = 30.0 + ilo as f64 / 1024.0;
    let near: Vec<u64> = (0..4096)
        .map(|i| {
            2f64.powf(centre + (i as f64 - 2048.0) / 2048.0 / 1024.0)
                .round() as u64
        })
        .collect();
    let tn = parallel_terms(&[1, 2, 3, 6], &near);
    let inf = (0..near.len())
        .map(|i| -(0..4).map(|j| tn[j][i][tm]).sum::<f64>() / (near[i] as f64).powf(theta))
        .fold(f64::INFINITY, f64::min);
    println!(
        "G_6(u) = -(q | 6 part)/x^(log_4 3) over one period u = log2 x in [30, 32), 2048 points: min {:.4} at u = {:.3}, max {:.4} at u = {:.3}, mean {:.4}; refined on 4096 points within 2^-10 of the minimum: {:.5}; G_6 = 1 at even u, 1/sqrt 3 = {:.4} at odd u",
        lo,
        centre,
        hi,
        30.0 + ihi as f64 / 1024.0,
        mean,
        inf,
        1.0 / 3f64.sqrt()
    );
    let tmd = Design::all()[tm];
    let big = terms(65535, &[1u64 << top]);
    let sgt = at(&rows, 1u64 << top).sigma();
    println!(
        "orbits on q = 2^m - 1: g = {:.3} at q = 255, {:.3} at q = 65535 (squarefree); T_65535 at 2^{} is {:.2} sigma, {:.2} sqrt(x)",
        growth(tmd, 255),
        growth(tmd, 65535),
        top,
        big[0][tm] / sgt,
        big[0][tm] / ((1u64 << top) as f64).sqrt()
    );
    for name in ["2:1", "6:2", "6:3", "6:4", "6:23", "6:1", "6:5"] {
        let d = find(name);
        let gi = grid.iter().position(|&g| g == MARK).unwrap();
        let a = at(&rows, MARK);
        let fm = flat(&grid, &gm.iter().map(|v| v[d]).collect::<Vec<f64>>());
        print!(
            " {} {:.1}/{:.1}",
            name,
            a.pi[d] as f64 / (a.primes as f64).sqrt(),
            fm[gi] / (a.primes as f64).sqrt()
        );
    }
    println!("  (z_B along the first 2^26 primes, observed/model)");
}

fn law(top: u32) {
    let rows = cached(top);
    let x = 1u64 << top;
    let a = at(&rows, x);
    let sg = a.sigma();
    let grid = grid(top);
    let gm = model(Q, &grid);
    let qs = squarefree(Q);
    let tq = parallel_terms(&qs, &[x]);
    let late: Vec<usize> = (0..grid.len())
        .filter(|&i| grid[i] >= FROM && grid[i] != MARK)
        .collect();
    let h: f64 = qs.iter().map(|&q| 1.0 / totient(q) as f64).sum();
    println!("# LAW");
    println!();
    println!(
        "x = 2^{}, sigma = {:.1}, sqrt(pi) = {:.1}; z = Psi/sigma, model = M_Q/sigma, Q = {}, r = (Psi - M_Q)/sigma; max and rms of r over the {} quarter-octave checkpoints from 2^20; top: the largest T_q; unit: the q | rad(b) part, the restriction to the residues prime to b; flat: sum_(p<x) x(p)/sqrt(pi) against the model weighted by 1/log n",
        top,
        sg,
        (a.primes as f64).sqrt(),
        Q,
        late.len()
    );
    println!();
    println!(
        "{:<8} {:>9} {:>9} {:>6} {:>5} {:>5} {:>13} {:>9} {:>9} {:>9} {:>6}",
        "design", "z", "model", "r", "max", "rms", "top", "unit", "z_flat", "flat", "r_flat"
    );
    let all = Design::all();
    let mut summary = Vec::new();
    for s in 0..COIN {
        let m: Vec<f64> = gm.iter().map(|v| v[s]).collect();
        let fl = flat(&grid, &m);
        let r: Vec<f64> = late
            .iter()
            .map(|&i| {
                let b = at(&rows, grid[i]);
                (b.psi[s] - m[i]) / b.sigma()
            })
            .collect();
        let rf: Vec<f64> = late
            .iter()
            .map(|&i| {
                let b = at(&rows, grid[i]);
                (b.pi[s] as f64 - fl[i]) / (b.primes as f64).sqrt()
            })
            .collect();
        let max = r.iter().fold(0f64, |m, v| m.max(v.abs()));
        let maxf = rf.iter().fold(0f64, |m, v| m.max(v.abs()));
        let rms = (r.iter().map(|v| v * v).sum::<f64>() / r.len() as f64).sqrt();
        let last = grid.iter().position(|&g| g == x).unwrap();
        let (tq_best, tv) = qs
            .iter()
            .zip(&tq)
            .map(|(q, v)| (*q, v[0][s]))
            .max_by(|p, q| p.1.abs().total_cmp(&q.1.abs()))
            .unwrap();
        let rad = if s < DESIGNS {
            factor(all[s].b).iter().map(|f| f.0).product::<u64>()
        } else {
            2
        };
        let unit: f64 = qs
            .iter()
            .zip(&tq)
            .filter(|(q, _)| rad % **q == 0)
            .map(|(_, v)| v[0][s])
            .sum();
        summary.push((s, max, maxf, rms, r[r.len() - 1]));
        println!(
            "{:<8} {:>9.2} {:>9.2} {:>6.2} {:>5.2} {:>5.2} {:>13} {:>9.2} {:>9.2} {:>9.2} {:>6.2}",
            label(s),
            a.psi[s] / sg,
            m[last] / sg,
            (a.psi[s] - m[last]) / sg,
            max,
            rms,
            format!("{}:{:.2}", tq_best, tv / sg),
            unit / sg,
            a.pi[s] as f64 / (a.primes as f64).sqrt(),
            fl[last] / (a.primes as f64).sqrt(),
            rf[rf.len() - 1]
        );
        let _ = maxf;
    }
    println!();
    let distinct: Vec<&(usize, f64, f64, f64, f64)> = summary
        .iter()
        .filter(|w| w.0 < DESIGNS && !["4:12", "5:13"].contains(&label(w.0).as_str()))
        .collect();
    for cut in [1.0, 2.0, 2.5, 3.0, 4.0] {
        println!(
            "of the 55 distinct designs, max |r| <= {}: {}, max |r_flat| <= {}: {}",
            cut,
            distinct.iter().filter(|w| w.1 <= cut).count(),
            cut,
            distinct.iter().filter(|w| w.2 <= cut).count()
        );
    }
    let at_x: Vec<f64> = distinct.iter().map(|w| w.4).collect();
    let pooled = (distinct.iter().map(|w| w.3 * w.3).sum::<f64>() / distinct.len() as f64).sqrt();
    let rms_x = (at_x.iter().map(|v| v * v).sum::<f64>() / at_x.len() as f64).sqrt();
    let iid = (1.0 - a.theta * h / a.lam2).sqrt();
    let sd = |b: &Acc| (1.0 - b.theta * h / b.lam2).sqrt();
    let floor = (late
        .iter()
        .map(|&i| sd(at(&rows, grid[i])).powi(2))
        .sum::<f64>()
        / late.len() as f64)
        .sqrt();
    let twins: Vec<(usize, usize)> = distinct
        .iter()
        .filter_map(|w| {
            let d = all[w.0];
            if d.b.is_multiple_of(2) {
                return None;
            }
            let odd: u32 = (1..d.b).filter(|e| e % 2 == 1).map(|e| 1 << (e - 1)).sum();
            let t = (0..DESIGNS).find(|&i| {
                all[i]
                    == Design {
                        b: d.b,
                        mask: d.mask ^ odd,
                    }
            })?;
            (w.0 < t && distinct.iter().any(|v| v.0 == t)).then_some((w.0, t))
        })
        .collect();
    let lone: Vec<f64> = distinct
        .iter()
        .filter(|w| !twins.iter().any(|p| p.1 == w.0))
        .map(|w| w.3)
        .collect();
    let pooled_lone = (lone.iter().map(|v| v * v).sum::<f64>() / lone.len() as f64).sqrt();
    println!(
        "iid signs give the residual sd sqrt(1 - theta H/sum Lambda^2), H = sum_(q<=Q) mu^2(q)/phi(q) = {:.3}: {:.3} at 2^20 and {:.3} at 2^{}, pooled over the checkpoints {:.3}",
        h,
        sd(at(&rows, FROM)),
        iid,
        top,
        floor
    );
    println!(
        "r over the 55: pooled rms over the checkpoints {:.3} against the pooled floor {:.3}, ratio {:.3}; at 2^{} alone rms {:.3} against {:.3}, chi2/55 = {:.2}",
        pooled,
        floor,
        pooled / floor,
        top,
        rms_x,
        iid,
        rms_x * rms_x / (iid * iid)
    );
    println!(
        "{} sign-twin pairs along the odd primes ({}), so {} independent streams, pooled rms {:.3}; rudin: rms {:.3}, max {:.2}",
        twins.len(),
        twins
            .iter()
            .map(|p| format!("{}/{}", label(p.0), label(p.1)))
            .collect::<Vec<_>>()
            .join(", "),
        lone.len(),
        pooled_lone,
        summary[RS].3,
        summary[RS].1
    );
    let coin: Vec<f64> = late
        .iter()
        .map(|&i| {
            let b = at(&rows, grid[i]);
            b.psi[COIN] / b.sigma()
        })
        .collect();
    println!(
        "coin, no model: z = Psi/sigma at 2^{} {:.2}, rms over the checkpoints {:.3}, max {:.2}",
        top,
        a.psi[COIN] / sg,
        (coin.iter().map(|v| v * v).sum::<f64>() / coin.len() as f64).sqrt(),
        coin.iter().fold(0f64, |m, v| m.max(v.abs()))
    );
    let sieve: Vec<f64> = (0..COIN + 1)
        .map(|s| (a.psi[s] - a.model(BINS - 1, s)) / sg)
        .collect();
    let cut7 = parallel_terms(
        &squarefree(210)
            .into_iter()
            .filter(|q| 210 % q == 0)
            .collect::<Vec<u64>>(),
        &grid,
    );
    let mut dev = 0f64;
    for (i, &g) in grid.iter().enumerate() {
        let b = at(&rows, g);
        for s in 0..COIN {
            let f: f64 = cut7.iter().map(|v| v[i][s]).sum();
            dev = dev.max((f - b.model(0, s)).abs() / b.sigma());
        }
    }
    let density = ratios()[BINS - 1] * a.mc[BINS - 1] as f64 / x as f64;
    println!(
        "sieve model M_z = (P(z)/phi(P(z))) sum_(n<x, (n,P(z))=1) x(n): at z = 7 it equals the q | 210 part of the expansion to {:.1e} sigma at each of the {} quarter-octave checkpoints; at z = 1021, r = {:.2} for the folded class 3:1, {:.2} for Thue-Morse, {:.2} for the coin; iid signs give sd {:.3} there; its density (P/phi(P)) #(n < x, (n,P)=1)/x is {:.5}",
        dev,
        grid.len(),
        sieve[find("3:1")],
        sieve[0],
        sieve[COIN],
        (1.0 - ratios()[BINS - 1] * a.theta / a.lam2).sqrt(),
        density
    );
}

fn lucas(n: u32) -> (f64, f64) {
    let (mut f, mut g) = (0f64, 1f64);
    for _ in 0..n {
        (f, g) = (g, f + g);
    }
    (f, 2.0 * g - f)
}

fn exact() {
    let all = Design::all();
    let ks: Vec<u32> = (1..=12).collect();
    let xs: Vec<u64> = ks.iter().map(|&k| 6u64.pow(k)).collect();
    let qs = [1u64, 2, 3, 6, 5, 10, 15, 30];
    let t = parallel_terms(&qs, &xs);
    let part = |d: usize, k: usize, keep: &dyn Fn(u64) -> bool| -> f64 {
        qs.iter()
            .zip(&t)
            .filter(|(q, _)| keep(**q))
            .map(|(_, v)| v[k][d])
            .sum()
    };
    println!("# EXACT");
    println!();
    println!("at x = 6^k, k = 1..12, for the 31 designs of base 6, P(z) = sum_(d<6) v(d) z^d, zeta = e(1/5): the q | 6 part against 3 (v(1) + v(5)) P(1)^(k-1), the q in (5, 10, 15, 30) part against -(3/4) sum_(a=1..4) (v(1) zeta^a + v(5)) P(zeta^a)^(k-1); deviation relative to max(1, |prediction|)");
    let (mut du, mut df) = (0f64, 0f64);
    for d in (0..DESIGNS).filter(|&d| all[d].b == 6) {
        let v = all[d].tile();
        for (i, &k) in ks.iter().enumerate() {
            let p1: f64 = v.iter().sum();
            let unit = 3.0 * (v[1] + v[5]) * p1.powi(k as i32 - 1);
            du = du.max((part(d, i, &|q| 6 % q == 0) - unit).abs() / unit.abs().max(1.0));
            let mut fifth = Cx::ZERO;
            for a in 1..5u64 {
                let mut pz = Cx::ZERO;
                for (dd, &vd) in v.iter().enumerate() {
                    pz = pz.add(Cx::turn(a * dd as u64 % 5, 5).scale(vd));
                }
                let mut pk = Cx::ONE;
                for _ in 1..k {
                    pk = pk.mul(pz);
                }
                fifth = fifth.add(Cx::turn(a, 5).scale(v[1]).add(Cx::ONE.scale(v[5])).mul(pk));
            }
            let pred = -0.75 * fifth.re;
            df = df.max((part(d, i, &|q| q % 5 == 0) - pred).abs() / pred.abs().max(1.0));
        }
    }
    println!(
        "unit part: max deviation {:.1e}; fifths part: max deviation {:.1e}",
        du, df
    );
    let (d23, d14) = (find("6:23"), find("6:14"));
    let mut dl = 0f64;
    for (i, &k) in ks.iter().enumerate() {
        let (_, l) = lucas(3 * k - 1);
        let (f, _) = lucas(3 * k - 2);
        dl = dl.max((part(d23, i, &|q| q % 5 == 0) + 0.75 * l).abs() / l);
        dl = dl.max((part(d14, i, &|q| q % 5 == 0) + 3.75 * f).abs() / f);
    }
    println!("6:23: the fifths part is -(3/4) L_(3k-1), 6:14: -(15/4) F_(3k-2), Lucas and Fibonacci: max relative deviation {:.1e}", dl);
    println!();
    let eights: Vec<u64> = (1..=6).map(|m| 8u64.pow(m)).collect();
    let t7 = terms(7, &eights);
    println!(
        "Thue-Morse T_7 at 8^m, m = 1..6, its orbit product -i sqrt 7 turning a quarter each step: {}",
        t7.iter().map(|v| format!("{:.3}", v[0])).collect::<Vec<_>>().join(", ")
    );
    println!();
    println!(
        "{:>3} {:>16} {:>16} {:>16} {:>16}",
        "k", "6:23 fifths", "-(3/4)L", "6:14 fifths", "-(15/4)F"
    );
    for (i, &k) in ks.iter().enumerate().filter(|(_, k)| **k <= 6) {
        let (_, l) = lucas(3 * k - 1);
        let (f, _) = lucas(3 * k - 2);
        println!(
            "{:>3} {:>16.3} {:>16.3} {:>16.3} {:>16.3}",
            k,
            part(d23, i, &|q| q % 5 == 0),
            -0.75 * l,
            part(d14, i, &|q| q % 5 == 0),
            -3.75 * f
        );
    }
}

fn growth(d: Design, q: u64) -> f64 {
    let v = d.tile();
    let mut best = f64::NEG_INFINITY;
    for a in (0..q).filter(|&a| gcd(a, q) == 1) {
        let mut seen = vec![usize::MAX; q as usize];
        let mut path = Vec::new();
        let mut c = a;
        while seen[c as usize] == usize::MAX {
            seen[c as usize] = path.len();
            path.push(c);
            c = c * d.b % q;
        }
        let cycle = &path[seen[c as usize]..];
        let mut sum = 0.0;
        for &c in cycle {
            let mut p = Cx::ZERO;
            for (k, &vk) in v.iter().enumerate() {
                p = p.add(Cx::turn(c * k as u64 % q, q).scale(vk));
            }
            sum += p.re.hypot(p.im).max(1e-300).ln();
        }
        best = best.max(sum / (cycle.len() as f64 * (d.b as f64).ln()));
    }
    best
}

const LADDER: [u64; 9] = [6, 30, 64, 128, 256, 512, 1024, 2048, 4096];

fn ladder(top: u32, qmax: u64) {
    let rows = cached(top);
    let clock = Instant::now();
    let xs: Vec<u64> = [top - 4, top - 2, top].iter().map(|&k| 1u64 << k).collect();
    let qs = squarefree(qmax);
    let t = parallel_terms(&qs, &xs);
    let fr: u64 = qs.iter().map(|&q| totient(q)).sum();
    let cuts: Vec<u64> = LADDER.iter().copied().filter(|&c| c <= qmax).collect();
    println!("# LADDER");
    println!();
    println!(
        "M_Q = sum_(q <= Q, q squarefree) T_q; {} denominators, {} fractions, {:.1} s; rms and max of r_Q = (Psi - M_Q)/sigma over the 54 distinct non-folded designs and Rudin-Shapiro, at x = 2^k; floor: the iid-sign sd sqrt(1 - theta H_Q/sum Lambda^2), H_Q = sum_(q<=Q) mu^2(q)/phi(q); ratio: rms/floor; excess: rms^2 - floor^2, in sigma^2",
        qs.len(),
        fr,
        clock.elapsed().as_secs_f64()
    );
    println!();
    print!("{:<6} {:<4}", "x", "");
    for c in &cuts {
        print!(" {:>6}", format!("Q{}", c));
    }
    println!();
    let keep: Vec<usize> = (0..COIN)
        .filter(|&d| !["4:12", "5:13", "3:1"].contains(&label(d).as_str()))
        .collect();
    for (xi, &x) in xs.iter().enumerate() {
        let a = at(&rows, x);
        let sg = a.sigma();
        let r: Vec<Vec<f64>> = cuts
            .iter()
            .map(|&c| {
                keep.iter()
                    .map(|&d| {
                        let m: f64 = qs
                            .iter()
                            .zip(&t)
                            .filter(|(q, _)| **q <= c)
                            .map(|(_, v)| v[xi][d])
                            .sum();
                        (a.psi[d] - m) / sg
                    })
                    .collect()
            })
            .collect();
        let floors: Vec<f64> = cuts
            .iter()
            .map(|&c| {
                let hq: f64 = qs
                    .iter()
                    .filter(|&&q| q <= c)
                    .map(|&q| 1.0 / totient(q) as f64)
                    .sum();
                (1.0 - a.theta * hq / a.lam2).sqrt()
            })
            .collect();
        for (stat, f) in [
            ("rms", 0),
            ("max", 1),
            ("floor", 2),
            ("ratio", 3),
            ("excess", 4),
        ] {
            print!("{:<6} {:<5}", format!("2^{}", 63 - x.leading_zeros()), stat);
            for (ci, v) in r.iter().enumerate() {
                let rms = (v.iter().map(|y| y * y).sum::<f64>() / v.len() as f64).sqrt();
                let val = match f {
                    0 => rms,
                    1 => v.iter().fold(0f64, |m, y| m.max(y.abs())),
                    2 => floors[ci],
                    3 => rms / floors[ci],
                    _ => rms * rms - floors[ci] * floors[ci],
                };
                print!(" {:>6.3}", val);
            }
            println!();
        }
    }
    let x = *xs.last().unwrap();
    let sg = at(&rows, x).sigma();
    println!();
    println!("largest T_q at 2^{}, as q:T_q/sigma:g, g the growth exponent of the best orbit through a/q, |S(a/q; b^k)| ~ b^(k g), against 1/2 for the noise:", top);
    for name in ["2:1", "4:13", "6:14", "6:135", "6:23", "6:2", "6:1", "6:5"] {
        let d = find(name);
        let mut v: Vec<(u64, f64)> = qs
            .iter()
            .zip(&t)
            .map(|(q, tt)| (*q, tt[xs.len() - 1][d] / sg))
            .collect();
        v.sort_by(|p, q| q.1.abs().total_cmp(&p.1.abs()));
        print!("{:<6}", name);
        for (q, val) in v.iter().take(10) {
            print!(" {}:{:.2}:{:.3}", q, val, growth(Design::all()[d], *q));
        }
        println!();
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let clock = Instant::now();
    let num = |i: usize, d: u32| args.get(i).and_then(|s| s.parse().ok()).unwrap_or(d);
    match args.get(1).map(String::as_str) {
        Some("census") => summary(num(2, 32)),
        Some("thue") => thue(num(2, 32)),
        Some("law") => law(num(2, 32)),
        Some("ladder") => ladder(num(2, 32), num(3, 4096) as u64),
        Some("exact") => exact(),
        _ => {
            eprintln!("verbs: census [k], thue [k], law [k], ladder [k q], exact");
            std::process::exit(2);
        }
    }
    eprintln!("{:.2} s", clock.elapsed().as_secs_f64());
}

#[cfg(test)]
mod tests {
    use super::*;

    fn twisted(d: Design, a: u64, q: u64, x: u64) -> Cx {
        let tab: Vec<Cx> = (0..q).map(|r| Cx::turn(r, q)).collect();
        let mut out = [Cx::ZERO];
        sums(
            &[d.tile()],
            d.b as usize,
            a,
            q,
            &tab,
            &digits(x, d.b),
            &mut out,
        );
        out[0]
    }

    fn value(d: usize, n: u64) -> i64 {
        match d {
            RS => 1 - 2 * rudin(n) as i64,
            COIN => 1 - 2 * coin(n) as i64,
            _ => {
                let des = Design::all()[d];
                let (mut v, mut p) = (n, 0);
                while v > 0 {
                    p ^= des.has(v % des.b) as i64;
                    v /= des.b;
                }
                1 - 2 * p
            }
        }
    }

    #[test]
    fn reader_matches_the_digits() {
        let reader = Reader::new();
        for v in (0..3000u64).chain((0..3000u64).map(|i| mix(i) >> 32)) {
            let w = reader.word(v);
            for d in 0..STREAMS {
                assert_eq!(sign(w, d), value(d, v), "{} at {}", label(d), v);
            }
        }
    }

    #[test]
    fn census_matches_trial_division() {
        let top = 16;
        let rows = census(top);
        let n = 1usize << top;
        let mut lpf = vec![0u64; n];
        for i in 2..n {
            if lpf[i] == 0 {
                for j in (i..n).step_by(i) {
                    if lpf[j] == 0 {
                        lpf[j] = i as u64;
                    }
                }
            }
        }
        let lam = |m: usize| {
            if m < 2 {
                return 0.0;
            }
            let p = lpf[m];
            let mut r = m as u64;
            while r % p == 0 {
                r /= p;
            }
            if r == 1 {
                (p as f64).ln()
            } else {
                0.0
            }
        };
        let small = |m: u64| (2..=1021u64).find(|p| lpf[*p as usize] == *p && m % p == 0);
        for (x, a) in rows.iter().step_by(7) {
            for d in [0, 5, 20, 40, 56, RS, COIN] {
                let (mut pi, mut psi, mut m, mut u) = (0i64, 0f64, [0i64; BINS], [0f64; BINS]);
                for k in 1..*x {
                    let s = value(d, k);
                    psi += s as f64 * lam(k as usize);
                    if k >= 2 && lpf[k as usize] == k {
                        pi += s;
                    }
                    match small(k) {
                        Some(p) if p <= 7 => {}
                        found => {
                            let bin = found.map_or(BINS, |p| level(p) as usize) - 1;
                            m[bin] += s;
                            if k > 1 {
                                u[bin] += s as f64 / (k as f64).ln();
                            }
                        }
                    }
                }
                assert_eq!(a.pi[d], pi, "{} at {}", label(d), x);
                assert!((a.psi[d] - psi).abs() < 1e-6, "{} at {}", label(d), x);
                for b in 0..BINS {
                    assert_eq!(a.m[b][d], m[b], "{} at {} bin {}", label(d), x, b);
                    assert!((a.u[b][d] - u[b]).abs() < 1e-6);
                }
            }
        }
    }

    #[test]
    fn twisted_sums_match_the_definition() {
        let all = Design::all();
        for (d, a, q, x) in [
            (0usize, 1u64, 3u64, 1000u64),
            (0, 5, 17, 777),
            (20, 2, 7, 3001),
            (40, 4, 37, 5000),
            (56, 3, 10, 1296),
            (RS, 3, 11, 2049),
            (RS, 1, 3, 4096),
        ] {
            let mut want = Cx::ZERO;
            for n in 0..x {
                want = want.add(Cx::turn(a * n % q, q).scale(value(d, n) as f64));
            }
            let got = if d == RS {
                let tab: Vec<Cx> = (0..q).map(|r| Cx::turn(r, q)).collect();
                rudin_sum(a, q, &tab, &digits(x, 2))
            } else {
                twisted(all[d], a, q, x)
            };
            let diff = got.add(want.scale(-1.0));
            assert!(
                diff.re.hypot(diff.im) < 1e-8,
                "{} {}/{} {}",
                label(d),
                a,
                q,
                x
            );
        }
    }

    #[test]
    fn thue_morse_closed_forms() {
        let tm = Design { b: 2, mask: 1 };
        for x in (0..3000u64).map(|i| mix(i) % 100_000) {
            let (s, s4) = (twisted(tm, 1, 3, x), twisted(tm, 1, 3, 4 * x));
            assert!((s4.re - 3.0 * s.re).abs() < 1e-7 && (s4.im - 3.0 * s.im).abs() < 1e-7);
        }
        for k in 2..=20u32 {
            let sum: i64 = (0..1u64 << k)
                .filter(|n| n % 2 == 1 && n % 3 != 0)
                .map(|n| value(0, n))
                .sum();
            assert_eq!(3 * sum, -(3i64.pow(k / 2)));
        }
    }

    #[test]
    fn base_six_fifths_are_fibonacci_and_lucas() {
        for (name, k) in [
            ("6:23", 1u32),
            ("6:23", 4),
            ("6:23", 6),
            ("6:14", 2),
            ("6:14", 5),
        ] {
            let d = find(name);
            let x = 6u64.pow(k);
            let sum: i64 = (1..x)
                .filter(|n| n % 2 == 1 && n % 3 != 0)
                .map(|n| value(d, n) * if n % 5 == 0 { 4 } else { -1 })
                .sum();
            let (f, l) = (lucas(3 * k - 2).0, lucas(3 * k - 1).1);
            let want = if name == "6:23" { l } else { 5.0 * f };
            assert_eq!(sum as f64, want, "{} {}", name, k);
        }
    }
}
