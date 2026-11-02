use std::collections::HashMap;
use std::env;
use std::time::Instant;

// DESIGNS

#[derive(Clone, Copy, PartialEq, Debug)]
enum Design {
    Thue,
    Rudin,
    Ternary,
}

const DESIGNS: [Design; 3] = [Design::Thue, Design::Rudin, Design::Ternary];
const CHUNK: u64 = 59049;
const CHUNK30: u128 = 205_891_132_094_649;

impl Design {
    fn name(self) -> &'static str {
        match self {
            Design::Thue => "thue-morse",
            Design::Rudin => "rudin-shapiro",
            Design::Ternary => "ternary 2s",
        }
    }
    fn base(self) -> i128 {
        match self {
            Design::Ternary => 3,
            _ => 2,
        }
    }
}

struct Kit {
    table: Vec<u8>,
}

impl Kit {
    fn new() -> Kit {
        let mut table = vec![0u8; CHUNK as usize];
        for m in 1..CHUNK as usize {
            table[m] = table[m / 3] ^ u8::from(m % 3 == 2);
        }
        Kit { table }
    }
    fn ternary(&self, mut m: u64) -> u8 {
        let mut p = 0;
        while m > 0 {
            p ^= self.table[(m % CHUNK) as usize];
            m /= CHUNK;
        }
        p
    }
    fn bit(&self, d: Design, m: u128) -> u8 {
        match d {
            Design::Thue => (m.count_ones() & 1) as u8,
            Design::Rudin => ((m & (m >> 1)).count_ones() & 1) as u8,
            Design::Ternary => {
                let (mut m, mut p) = (m, 0);
                while m >> 64 != 0 {
                    p ^= self.ternary((m % CHUNK30) as u64);
                    m /= CHUNK30;
                }
                p ^ self.ternary(m as u64)
            }
        }
    }
    fn sign(&self, d: Design, m: i128) -> i64 {
        assert!(m >= 0, "a design is read at a negative integer");
        1 - 2 * self.bit(d, m as u128) as i64
    }
}

// READINGS

#[derive(Clone)]
enum Reading {
    Poly(&'static str, Vec<i128>),
    Primes,
}

impl Reading {
    fn name(&self) -> &'static str {
        match self {
            Reading::Poly(name, _) => name,
            Reading::Primes => "p_n",
        }
    }
}

fn listed() -> Vec<Reading> {
    vec![
        Reading::Poly("n", vec![0, 1]),
        Reading::Poly("n^2", vec![0, 0, 1]),
        Reading::Poly("n^2 + 1", vec![1, 0, 1]),
        Reading::Poly("n^2 + n", vec![0, 1, 1]),
        Reading::Poly("(2n + 1)^2", vec![1, 4, 4]),
        Reading::Poly("2n^2 + 1", vec![1, 0, 2]),
        Reading::Poly("n^3", vec![0, 0, 0, 1]),
        Reading::Primes,
    ]
}

fn probes() -> Vec<Reading> {
    vec![
        Reading::Poly("n^2 + 2", vec![2, 0, 1]),
        Reading::Poly("n^2 + 3", vec![3, 0, 1]),
        Reading::Poly("n^2 + 4", vec![4, 0, 1]),
        Reading::Poly("n^2 + 5", vec![5, 0, 1]),
        Reading::Poly("(n + 1)^2 + 1", vec![2, 2, 1]),
        Reading::Poly("n^3 + n", vec![0, 1, 0, 1]),
    ]
}

// POLYNOMIALS

fn eval(p: &[i128], n: i128) -> i128 {
    p.iter().rev().fold(0, |acc, &c| acc * n + c)
}

fn compose(p: &[i128], slope: i128, shift: i128) -> Vec<i128> {
    let mut out = vec![0i128; p.len()];
    let mut power = vec![1i128];
    for (i, &c) in p.iter().enumerate() {
        if i > 0 {
            let mut next = vec![0i128; power.len() + 1];
            for (j, &e) in power.iter().enumerate() {
                next[j] += e * shift;
                next[j + 1] += e * slope;
            }
            power = next;
        }
        for (j, &e) in power.iter().enumerate() {
            out[j] += c * e;
        }
    }
    out
}

struct Head {
    level: u32,
    shift: i128,
    tail: i128,
}

fn near_head(p: &[i128], b: i128, phi: (i128, i128), psi: (i128, i128)) -> Option<Head> {
    let f = compose(p, phi.0, phi.1);
    let g = compose(p, psi.0, psi.1);
    let k = f.len() - 1;
    if f[k] == 0 || g[k] % f[k] != 0 {
        return None;
    }
    let (mut ratio, mut level) = (g[k] / f[k], 0u32);
    while ratio % b == 0 {
        ratio /= b;
        level += 1;
    }
    if ratio != 1 || level == 0 {
        return None;
    }
    let scale = b.pow(level);
    if (1..=k).any(|i| g[i] != scale * f[i]) {
        return None;
    }
    let s = g[0] - scale * f[0];
    Some(Head {
        level,
        shift: s.div_euclid(scale),
        tail: s.rem_euclid(scale),
    })
}

// PRIMES

fn sieve(count: usize, mut each: impl FnMut(usize, u64)) {
    if count == 0 {
        return;
    }
    each(0, 2);
    let nf = count.max(6) as f64;
    let limit = nf * (nf.ln() + nf.ln().ln()) + 100.0;
    let root = limit.sqrt() as usize + 2;
    let mut flags = vec![true; root + 1];
    let mut base = Vec::new();
    for i in 2..=root {
        if flags[i] {
            if i > 2 {
                base.push(i as u64);
            }
            let mut j = i * i;
            while j <= root {
                flags[j] = false;
                j += i;
            }
        }
    }
    const SEG: usize = 1 << 18;
    let mut next: Vec<u64> = base.iter().map(|&p| p * p).collect();
    let mut seg = vec![false; SEG];
    let (mut lo, mut found) = (3u64, 1usize);
    loop {
        let hi = lo + 2 * SEG as u64;
        seg.iter_mut().for_each(|x| *x = false);
        for (i, &p) in base.iter().enumerate() {
            let mut m = next[i];
            while m < hi {
                seg[((m - lo) / 2) as usize] = true;
                m += 2 * p;
            }
            next[i] = m;
        }
        for (i, &c) in seg.iter().enumerate() {
            if !c {
                each(found, lo + 2 * i as u64);
                found += 1;
                if found == count {
                    return;
                }
            }
        }
        lo = hi;
    }
}

// STREAMS

struct Stream {
    bits: Vec<u64>,
}

impl Stream {
    fn new(len: usize) -> Stream {
        Stream {
            bits: vec![0; len / 64 + 3],
        }
    }
    fn set(&mut self, i: usize) {
        self.bits[i >> 6] |= 1 << (i & 63);
    }
    fn get(&self, i: usize) -> u64 {
        (self.bits[i >> 6] >> (i & 63)) & 1
    }
    fn window(&self, pos: usize, len: usize) -> u64 {
        let (w, o) = (pos >> 6, pos & 63);
        let mut v = self.bits[w] >> o;
        if o + len > 64 {
            v |= self.bits[w + 1] << (64 - o);
        }
        if len == 64 {
            v
        } else {
            v & ((1u64 << len) - 1)
        }
    }
    fn ones(&self, from: usize, to: usize) -> u64 {
        (from..to).map(|i| self.get(i)).sum()
    }
}

fn build(kit: &Kit, d: Design, p: &[i128], len: usize) -> Stream {
    let mut s = Stream::new(len);
    let small = p.iter().all(|&c| c >= 0) && eval(p, len as i128) < 1i128 << 62;
    let q: Vec<u64> = p.iter().map(|&c| c as u64).collect();
    for n in 0..len {
        let m = if small {
            q.iter().rev().fold(0u64, |acc, &c| acc * n as u64 + c) as u128
        } else {
            eval(p, n as i128) as u128
        };
        if kit.bit(d, m) == 1 {
            s.set(n);
        }
    }
    s
}

fn build_primes(kit: &Kit, len: usize) -> Vec<Stream> {
    let mut out: Vec<Stream> = DESIGNS.iter().map(|_| Stream::new(len)).collect();
    sieve(len, |i, p| {
        for (k, &d) in DESIGNS.iter().enumerate() {
            if kit.bit(d, p as u128) == 1 {
                out[k].set(i);
            }
        }
    });
    out
}

fn build_control(len: usize) -> Stream {
    let mut s = Stream::new(len);
    let mut state = 1u64;
    for n in 0..len {
        state = state.wrapping_add(0x9E3779B97F4A7C15);
        let mut z = state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
        if (z ^ (z >> 31)) >> 63 == 1 {
            s.set(n);
        }
    }
    s
}

fn stream(kit: &Kit, d: Design, r: &Reading, len: usize, primes: &[Stream]) -> Stream {
    match r {
        Reading::Poly(_, p) => build(kit, d, p, len),
        Reading::Primes => {
            let k = DESIGNS.iter().position(|&e| e == d).unwrap();
            Stream {
                bits: primes[k].bits[..len / 64 + 3].to_vec(),
            }
        }
    }
}

// TESTS

fn pearson(n: f64, su: f64, sv: f64, suv: f64) -> Option<f64> {
    let (mu, mv) = (su / n, sv / n);
    let (vu, vv) = (1.0 - mu * mu, 1.0 - mv * mv);
    if vu < 1e-12 || vv < 1e-12 {
        return None;
    }
    Some((suv / n - mu * mv) / (vu * vv).sqrt())
}

fn bit_pearson(n: usize, ones_u: u64, ones_v: u64, apart: u64) -> Option<f64> {
    let nf = n as f64;
    pearson(
        nf,
        nf - 2.0 * ones_u as f64,
        nf - 2.0 * ones_v as f64,
        nf - 2.0 * apart as f64,
    )
}

struct Verdict {
    mean: f64,
    shift: f64,
    shift_at: usize,
    dilation: f64,
    dilation_at: (usize, usize),
    kernel: f64,
    kernel_at: usize,
    square: f64,
    table: Vec<(usize, usize, f64)>,
}

fn verdict(s: &Stream, n: usize) -> Verdict {
    assert!(n % 64 == 0);
    let words = n / 64;
    let ones: u64 = s.bits[..words].iter().map(|w| w.count_ones() as u64).sum();
    let mean = 1.0 - 2.0 * ones as f64 / n as f64;
    let (mut shift, mut shift_at) = (0f64, 0usize);
    for h in 1..=64usize {
        let apart: u64 = (0..words)
            .map(|i| (s.bits[i] ^ s.window(64 * i + h, 64)).count_ones() as u64)
            .sum();
        let ones_v = ones - s.ones(0, h) + s.ones(n, n + h);
        if let Some(p) = bit_pearson(n, ones, ones_v, apart) {
            if p.abs() > shift {
                shift = p.abs();
                shift_at = h;
            }
        }
    }
    let mut ones_v = [[0u64; 9]; 10];
    let mut apart = [[0u64; 9]; 10];
    for i in 0..n {
        let u = s.get(i);
        for big in 2..=9usize {
            let w = s.window(big * i, big);
            for r in 0..big {
                let v = (w >> r) & 1;
                ones_v[big][r] += v;
                apart[big][r] += u ^ v;
            }
        }
    }
    let (mut dilation, mut dilation_at) = (0f64, (0usize, 0usize));
    let mut table = Vec::new();
    for big in 2..=9usize {
        for r in 0..big {
            if let Some(p) = bit_pearson(n, ones, ones_v[big][r], apart[big][r]) {
                table.push((big, r, p));
                if p.abs() > dilation + 1e-12 {
                    dilation = p.abs();
                    dilation_at = (big, r);
                }
            }
        }
    }
    let (mut kernel, mut kernel_at) = (0f64, 0usize);
    for b in [2usize, 3] {
        let (mut ou, mut ov, mut ap) = (0u64, 0u64, 0u64);
        for i in 0..n {
            let (u, v) = (s.get(b * i), s.get(b * b * i));
            ou += u;
            ov += v;
            ap += u ^ v;
        }
        if let Some(p) = bit_pearson(n, ou, ov, ap) {
            if p.abs() > kernel {
                kernel = p.abs();
                kernel_at = b;
            }
        }
    }
    let odd: u64 = (0..n)
        .map(|i| s.get(i) ^ s.get(2 * i) ^ s.get(3 * i) ^ s.get(6 * i))
        .sum();
    Verdict {
        mean,
        shift,
        shift_at,
        dilation,
        dilation_at,
        kernel,
        kernel_at,
        square: 1.0 - 2.0 * odd as f64 / n as f64,
        table,
    }
}

fn lz78(s: &Stream, len: usize) -> f64 {
    let mut dict: HashMap<(u32, u8), u32> = HashMap::new();
    let (mut node, mut next, mut phrases) = (0u32, 1u32, 0usize);
    for i in 0..len {
        let b = s.get(i) as u8;
        match dict.get(&(node, b)) {
            Some(&c) => node = c,
            None => {
                dict.insert((node, b), next);
                next += 1;
                phrases += 1;
                node = 0;
            }
        }
    }
    let c = (phrases + usize::from(node != 0)) as f64;
    c * (c.log2() + 1.0) / len as f64
}

// HEADS

const PERIODS: [usize; 18] = [
    1, 2, 3, 4, 6, 8, 9, 12, 16, 18, 24, 27, 32, 36, 48, 54, 64, 81,
];
const LAMBDAS: [usize; 5] = [2, 3, 4, 8, 9];

struct Identity {
    lambda: usize,
    mu: usize,
    depth: u32,
    period: usize,
    mean: f64,
}

fn heads(s: &Stream, n: usize) -> Vec<Identity> {
    let mut out = Vec::new();
    for &lambda in &LAMBDAS {
        for mu in 0..lambda {
            for depth in 0..3u32 {
                let lj = lambda.pow(depth);
                let off = mu * (lj - 1) / (lambda - 1);
                let w: Vec<u8> = (0..n)
                    .map(|i| {
                        let a = lj * i + off;
                        (s.get(a) ^ s.get(lambda * a + mu)) as u8
                    })
                    .collect();
                if let Some(&period) = PERIODS.iter().find(|&&p| (p..n).all(|i| w[i] == w[i % p])) {
                    let mean = w.iter().map(|&x| 1.0 - 2.0 * x as f64).sum::<f64>() / n as f64;
                    out.push(Identity {
                        lambda,
                        mu,
                        depth,
                        period,
                        mean,
                    });
                }
            }
        }
    }
    out
}

// LOCAL CONSTANTS

fn up(x: f64) -> String {
    if x == 0.0 {
        return "0".into();
    }
    let e = x.log10().floor() as i32;
    let m = (x / 10f64.powi(e - 2) - 1e-9).ceil();
    format!("{:.2}e{}", m / 100.0, e)
}

fn local_constant(kit: &Kit, d: Design, f: &[i128], head: &Head, depth: u32) -> (f64, f64) {
    let b = d.base();
    let window = b.pow(depth);
    let classes = if d == Design::Rudin {
        2 * window
    } else {
        window
    };
    let tail_sign = kit.sign(d, head.tail);
    let tail_top = (head.tail >> (head.level - 1)) & 1;
    let fast = f.iter().all(|&c| c >= 0) && eval(f, classes) < 1i128 << 62;
    let g: Vec<u64> = f.iter().map(|&c| c as u64).collect();
    let (mut total, mut loose) = (0i64, 0i64);
    for n0 in 0..classes {
        let m = if fast {
            (g.iter().rev().fold(0u64, |acc, &c| acc * n0 as u64 + c) % classes as u64) as i128
        } else {
            eval(f, n0).rem_euclid(classes)
        };
        let a = m % window;
        let c = a + head.shift;
        if c < 0 || c >= window {
            loose += 1;
            continue;
        }
        let mut v = tail_sign * kit.sign(d, c) * kit.sign(d, a);
        if d == Design::Rudin {
            let high = m / window;
            let flip = ((c >> (depth - 1)) ^ (a >> (depth - 1))) & 1;
            if high & flip == 1 {
                v = -v;
            }
            if c & 1 & tail_top == 1 {
                v = -v;
            }
        }
        total += v;
    }
    (total as f64 / classes as f64, loose as f64 / classes as f64)
}

fn pair_means(
    kit: &Kit,
    d: Design,
    p: &[i128],
    phi: (i128, i128),
    psi: (i128, i128),
    n: usize,
) -> (f64, Option<f64>) {
    let (mut su, mut sv, mut suv) = (0i64, 0i64, 0i64);
    for i in 0..n as i128 {
        let u = kit.sign(d, eval(p, phi.0 * i + phi.1));
        let v = kit.sign(d, eval(p, psi.0 * i + psi.1));
        su += u;
        sv += v;
        suv += u * v;
    }
    let nf = n as f64;
    (
        suv as f64 / nf,
        pearson(nf, su as f64, sv as f64, suv as f64),
    )
}

fn battery_pairs() -> Vec<(String, (i128, i128), (i128, i128))> {
    let mut out = Vec::new();
    for big in 2..=9i128 {
        for r in 0..big {
            out.push((format!("dilation ({big}, {r})"), (1, 0), (big, r)));
        }
    }
    for b in [2i128, 3] {
        out.push((format!("kernel ({b}n, {}n)", b * b), (b, 0), (b * b, 0)));
    }
    out
}

// BATTERY

fn row(v: &Verdict) -> String {
    format!(
        "mean {:+.4} | shift max|r| {:.4} at h={} | dilation max|r| {:.4} at (B, r)=({}, {}) | kernel max|r| {:.4} at b={} | square mean {:+.4}",
        v.mean, v.shift, v.shift_at, v.dilation, v.dilation_at.0, v.dilation_at.1, v.kernel, v.kernel_at, v.square
    )
}

fn listing(ids: &[&Identity], shown: usize) -> String {
    let mut out: Vec<String> = ids
        .iter()
        .take(shown)
        .map(|i| {
            if i.period == 1 {
                format!("{}n+{} j={} {:+.0}", i.lambda, i.mu, i.depth, i.mean)
            } else {
                format!(
                    "{}n+{} j={} p={} {:+.4}",
                    i.lambda, i.mu, i.depth, i.period, i.mean
                )
            }
        })
        .collect();
    if ids.len() > shown {
        out.push("...".into());
    }
    format!("{} [{}]", ids.len(), out.join("; "))
}

fn battery() {
    let kit = Kit::new();
    let n = 1usize << 15;
    let len = 729 * n + 1024;
    let primes = build_primes(&kit, len);
    println!(
        "battery: N = {n} pairs per test, the design ladder's protocol; null band 4/sqrt(N) = {:.4}",
        4.0 / (n as f64).sqrt()
    );
    println!("identities: y(phi^j n) y(phi^(j+1) n) over n < N for phi(n) = lambda n + mu, lambda in {{2, 3, 4, 8, 9}}, mu < lambda, j = 0, 1, 2, when it is periodic in n with period at most 81: plain at period 1, twisted above, with its mean");
    let mut readings = listed();
    readings.extend(probes());
    for &d in &DESIGNS {
        for r in &readings {
            let s = stream(&kit, d, r, len, &primes);
            let v = verdict(&s, n);
            let lz: Vec<String> = [14u32, 16, 18]
                .iter()
                .map(|&e| format!("{:.3}", lz78(&s, 1 << e)))
                .collect();
            let ids = heads(&s, n);
            let plain: Vec<&Identity> = ids.iter().filter(|i| i.period == 1).collect();
            let twisted: Vec<&Identity> = ids.iter().filter(|i| i.period > 1).collect();
            println!(
                "  {:<13} {:<14} {} | LZ78 at 2^14, 2^16, 2^18: {} | plain {} | twisted {}",
                d.name(),
                r.name(),
                row(&v),
                lz.join(", "),
                listing(&plain, 4),
                listing(&twisted, 4)
            );
        }
    }
    let c = build_control(len);
    let v = verdict(&c, n);
    let ids = heads(&c, n);
    println!(
        "  {:<28} {} | LZ78 at 2^18: {:.3} | identities {}",
        "splitmix64, seed 1",
        row(&v),
        lz78(&c, 1 << 18),
        ids.len()
    );
}

// LOCAL

fn local() {
    let kit = Kit::new();
    println!("near-heads: a battery pair (phi, psi) of a polynomial reading P is a near-head at the design's base b when P(psi(n)) = b^L (P(phi(n)) + q) + s with 0 <= s < b^L; the raw mean of y(psi(n)) y(phi(n)) tends to x(s) kappa, and [kappa_J - eps_J, kappa_J + eps_J] encloses it at J digits; the raw mean at N, a multiple of the class period, must lie in it");
    let mut readings = listed();
    readings.extend(probes());
    let pairs = battery_pairs();
    let (mut heads_seen, mut exact, mut between, mut zero) = (0, 0, 0, 0);
    for &d in &DESIGNS {
        let (depth, n) = if d == Design::Ternary {
            (13u32, 3usize.pow(14))
        } else {
            (20u32, 1usize << 22)
        };
        for r in &readings {
            let p = match r {
                Reading::Poly(_, p) => p,
                Reading::Primes => continue,
            };
            let mut found = 0;
            for (label, phi, psi) in &pairs {
                let Some(h) = near_head(p, d.base(), *phi, *psi) else {
                    continue;
                };
                found += 1;
                let f = compose(p, phi.0, phi.1);
                let (kappa, eps) = local_constant(&kit, d, &f, &h, depth);
                let (raw, pe) = pair_means(&kit, d, p, *phi, *psi, n);
                let inside = (raw - kappa).abs() <= eps + 1e-12;
                heads_seen += 1;
                if eps == 0.0 && kappa.abs() == 1.0 {
                    exact += 1;
                } else if kappa.abs() + eps < 1.0 && kappa.abs() > eps {
                    between += 1;
                } else if kappa.abs() <= eps {
                    zero += 1;
                }
                println!(
                    "  {:<13} {:<14} {:<17} L={} q={:+} s={} | kappa_J {:+.6} eps_J <= {} at J={} | raw mean {:+.6} Pearson {} at N={} | inside {}",
                    d.name(),
                    r.name(),
                    label,
                    h.level,
                    h.shift,
                    h.tail,
                    kappa,
                    up(eps),
                    depth,
                    raw,
                    pe.map_or("const".into(), |x| format!("{x:+.6}")),
                    n,
                    inside
                );
            }
            if found == 0 {
                println!(
                    "  {:<13} {:<14} no near-head among the 46 battery pairs",
                    d.name(),
                    r.name()
                );
            }
        }
    }
    println!("near-head pairs {heads_seen}: exact (kappa = +-1, eps = 0) {exact}; strictly between 0 and 1 {between}; enclosing 0 {zero}");
    println!("the pair law along the squares: the limit of the mean of x(n^2 + c) x(n^2), enclosed at J = 26 digits (16 at base 3), beside the mean at N, the nearest rational p/q with q <= 64 and the distance from the limit to it, at most abs(kappa_J - p/q) + eps_J");
    let square = [0i128, 0, 1];
    for &d in &DESIGNS {
        let (depth, n) = if d == Design::Ternary {
            (16u32, 3usize.pow(14))
        } else {
            (26u32, 1usize << 22)
        };
        let cells: Vec<String> = (1..=8i128)
            .map(|c| {
                let h = Head {
                    level: 1,
                    shift: c,
                    tail: 0,
                };
                let (kappa, eps) = local_constant(&kit, d, &square, &h, depth);
                let raw = (0..n as i128)
                    .map(|i| kit.sign(d, i * i + c) * kit.sign(d, i * i))
                    .sum::<i64>() as f64
                    / n as f64;
                let (p, q) = (1..=64i64)
                    .map(|q| ((kappa * q as f64).round() as i64, q))
                    .min_by(|a, b| {
                        (kappa - a.0 as f64 / a.1 as f64)
                            .abs()
                            .total_cmp(&(kappa - b.0 as f64 / b.1 as f64).abs())
                    })
                    .unwrap();
                let far = (kappa - p as f64 / q as f64).abs() + eps;
                format!(
                    "c={c}: {kappa:+.9} eps <= {} mean {raw:+.9} nearest {p}/{q} within {}",
                    up(eps),
                    up(far)
                )
            })
            .collect();
        println!("  {:<13} {}", d.name(), cells.join(" | "));
    }
}

// FALSIFY

fn free(d: Design, r: &Reading) -> bool {
    match r {
        Reading::Primes => true,
        Reading::Poly(_, p) => battery_pairs()
            .iter()
            .all(|(_, phi, psi)| near_head(p, d.base(), *phi, *psi).is_none()),
    }
}

fn falsify(top: u32) {
    let kit = Kit::new();
    let sizes: Vec<usize> = (20..=top).step_by(2).map(|e| 1usize << e).collect();
    let big = *sizes.last().unwrap();
    let len = 9 * big + 128;
    println!("falsification: the squares, the readings with no near-head among the battery pairs, the primes and the control, at N = {}; z = max|r| sqrt(N), the largest of 64 shift, 44 dilation and 2 kernel correlations, against the control's",
        sizes.iter().map(|n| format!("2^{}", n.trailing_zeros())).collect::<Vec<_>>().join(", "));
    let mut readings = listed();
    readings.extend(probes());
    let squares = Reading::Poly("n^2", vec![0, 0, 1]);
    let report = |name: &str, s: &Stream| {
        let mut path: Vec<(f64, f64, f64, f64)> = Vec::new();
        for &n in &sizes {
            let clock = Instant::now();
            let v = verdict(s, n);
            path.push(((n as f64).log2(), v.shift, v.dilation, v.kernel));
            let z = (n as f64).sqrt();
            let mut top = v.table.clone();
            top.sort_by(|a, b| b.2.abs().total_cmp(&a.2.abs()));
            let lead: Vec<String> = top
                .iter()
                .take(3)
                .map(|&(big, r, p)| format!("({big}, {r}) {p:+.6}"))
                .collect();
            println!(
                "  {:<28} 2^{} {} | z shift {:.2} dilation {:.2} kernel {:.2} square {:.2} | top dilations {}",
                name,
                n.trailing_zeros(),
                row(&v),
                v.shift * z,
                v.dilation * z,
                v.kernel * z,
                v.square.abs() * z,
                lead.join(", ")
            );
            eprintln!(
                "    {name} at 2^{}: {:.1} s",
                n.trailing_zeros(),
                clock.elapsed().as_secs_f64()
            );
        }
        let (first, last) = (path[0], path[path.len() - 1]);
        let slope = |a: f64, b: f64| (b / a).log2() / (last.0 - first.0);
        let steps: Vec<String> = path
            .windows(2)
            .map(|w| format!("{:+.3}", (w[1].2 / w[0].2).log2() / (w[1].0 - w[0].0)))
            .collect();
        println!(
            "  {:<28} slope of log2 max|r| against log2 N from 2^{} to 2^{}: shift {:+.3} dilation {:+.3} kernel {:+.3} (noise -0.500); dilation by step {}",
            name,
            first.0,
            last.0,
            slope(first.1, last.1),
            slope(first.2, last.2),
            slope(first.3, last.3),
            steps.join(", ")
        );
    };
    for &d in &DESIGNS {
        let chosen = std::iter::once(&squares).chain(
            readings
                .iter()
                .filter(|r| matches!(r, Reading::Poly(..)) && free(d, r)),
        );
        for r in chosen {
            let clock = Instant::now();
            let s = stream(&kit, d, r, len, &[]);
            eprintln!(
                "    built {} {} in {:.1} s",
                d.name(),
                r.name(),
                clock.elapsed().as_secs_f64()
            );
            report(&format!("{} {}", d.name(), r.name()), &s);
        }
    }
    let clock = Instant::now();
    let primes = build_primes(&kit, len);
    eprintln!(
        "    sieved {len} primes in {:.1} s",
        clock.elapsed().as_secs_f64()
    );
    for (k, &d) in DESIGNS.iter().enumerate() {
        report(&format!("{} p_n", d.name()), &primes[k]);
    }
    drop(primes);
    report("splitmix64, seed 1", &build_control(len));
}

// LEAK

fn thue_correlations(len: usize) -> Vec<f64> {
    let mut c = vec![0f64; len + 2];
    c[0] = 1.0;
    c[1] = -1.0 / 3.0;
    for h in 2..len + 2 {
        let q = h / 2;
        c[h] = if h % 2 == 0 {
            c[q]
        } else {
            -(c[q] + c[q + 1]) / 2.0
        };
    }
    c
}

fn leak() {
    let m = [[1.0, -0.5, -0.5], [-0.5, -0.5, 1.0], [-0.5, 1.0, -0.5]];
    let mut u = [1.0f64, 0.0, 0.0];
    let check = 24u32;
    let c = thue_correlations(1 << check);
    let mut worst = 0f64;
    let mut sums = vec![u];
    for k in 0..62u32 {
        let r = (1u64 << k) % 3;
        let f = -1.0 / 3.0;
        let mut next = [0f64; 3];
        for i in 0..3 {
            next[i] = (0..3).map(|j| m[i][j] * u[j]).sum();
        }
        if r == 2 {
            next[0] -= f / 2.0;
        } else {
            next[1] -= f / 2.0;
        }
        next[2] += 0.5;
        u = next;
        sums.push(u);
        if k + 1 <= check {
            let h = 1usize << (k + 1);
            for res in 0..3 {
                let direct: f64 = (res..h).step_by(3).map(|i| c[i]).sum();
                worst = worst.max((direct - u[res]).abs());
            }
        }
    }
    println!("thue-morse correlations by residue mod 3: U_r(2H) = sum over r' of M[r][r'] U_r'(H) + b, M = [[1, -1/2, -1/2], [-1/2, -1/2, 1], [-1/2, 1, -1/2]], eigenvalues 0, 3/2, -3/2; the recursion against direct sums of C_t for H = 2^k, k <= {check}: max deviation {worst:.2e}");
    for k in [10u32, 20, 30, 40, 50, 60, 61] {
        let s = sums[k as usize];
        let count = ((1u64 << k) + 2) / 3;
        println!(
            "  k = {k}: U_0(2^k) (2/3)^k = {:+.6}, mean of C_t(h) over 3 | h < 2^k = {:+.4e}, times 2^(k log2(4/3)) = {:+.6}",
            s[0] * (2.0f64 / 3.0).powi(k as i32),
            s[0] / count as f64,
            s[0] / count as f64 * (4.0f64 / 3.0).powi(k as i32)
        );
    }
    let m2 = [[2i128, -1, -1], [-1, -1, 2], [-1, 2, -1]];
    let (mut v, mut holds) = ([6i128, 0, 0], 0u32);
    for k in 0..70u32 {
        let mut next = [0i128; 3];
        for i in 0..3 {
            next[i] = (0..3).map(|j| m2[i][j] * v[j]).sum();
        }
        let p = 1i128 << (k + 1);
        if (1u64 << k) % 3 == 2 {
            next[0] += p;
        } else {
            next[1] += p;
        }
        next[2] += 3 * p;
        v = next;
        let kk = k + 1;
        let delta = if kk % 2 == 0 { 144i128 } else { 126 };
        if 45 * v[0] == 96 * 3i128.pow(kk) + delta * (1i128 << kk) {
            holds += 1;
        }
    }
    println!("  exact: with V_k = 6 . 2^k U(2^k), integer, V_(k+1) = 2M V_k + 2^(k+1) 6b; U_0(2^k) = (16/45)(3/2)^k + 8/15 at even k and + 7/15 at odd k, that is 45 V_k[0] = 96 . 3^k + 144 . 2^k or 126 . 2^k: holds at {holds} of the 70 values k = 1..70");
    let mut big = 0f64;
    for h in (3..1usize << check).step_by(3) {
        big = big.max(c[h].abs());
    }
    println!("  max abs C_t(h) over 3 | h < 2^{check}, h > 0: {big:.6}");
}

// IDENTITIES

fn identities() {
    let kit = Kit::new();
    let top = 1u64 << 24;
    let t = |m: u64| kit.sign(Design::Thue, m as i128);
    let r = |m: u64| kit.sign(Design::Rudin, m as i128);
    let e = |m: u64| kit.sign(Design::Ternary, m as i128);
    let alt = |n: u64| if n % 2 == 0 { 1 } else { -1 };
    let checks: Vec<(&str, Box<dyn Fn(u64) -> bool + '_>)> = vec![
        ("t(4n^2) = t(n^2)", Box::new(|n| t(4 * n * n) == t(n * n))),
        (
            "t(n^2 + 1) = -(-1)^n t(n^2)",
            Box::new(|n| t(n * n + 1) == -alt(n) * t(n * n)),
        ),
        (
            "t(n^2 + 2) = -t(n^2)",
            Box::new(|n| t(n * n + 2) == -t(n * n)),
        ),
        (
            "t(n^2 + 3) = t(n^2)",
            Box::new(|n| t(n * n + 3) == t(n * n)),
        ),
        (
            "t(2n^2 + 1) = -t(n^2)",
            Box::new(|n| t(2 * n * n + 1) == -t(n * n)),
        ),
        (
            "t((2n + 1)^2) = -t(n^2 + n)",
            Box::new(|n| t((2 * n + 1) * (2 * n + 1)) == -t(n * n + n)),
        ),
        (
            "t(n^2 + 4) t(4n^2 + 4) = +1 at n = 0 mod 4, -1 otherwise",
            Box::new(|n| t(n * n + 4) * t(4 * n * n + 4) == if n % 4 == 0 { 1 } else { -1 }),
        ),
        (
            "(n + 1)^2 + 1 under t: y(2n + 1) = y(4n + 3) = -t((n + 1)^2)",
            Box::new(|n| {
                let y = |m: u64| t((m + 1) * (m + 1) + 1);
                y(2 * n + 1) == y(4 * n + 3) && y(2 * n + 1) == -t((n + 1) * (n + 1))
            }),
        ),
        ("r(4n^2) = r(n^2)", Box::new(|n| r(4 * n * n) == r(n * n))),
        (
            "r(n^2 + 1) = r(n^2)",
            Box::new(|n| r(n * n + 1) == r(n * n)),
        ),
        (
            "r(2n^2 + 1) = (-1)^n r(n^2)",
            Box::new(|n| r(2 * n * n + 1) == alt(n) * r(n * n)),
        ),
        (
            "r((2n + 1)^2) = r(n^2 + n)",
            Box::new(|n| r((2 * n + 1) * (2 * n + 1)) == r(n * n + n)),
        ),
        ("t3(9n^2) = t3(n^2)", Box::new(|n| e(9 * n * n) == e(n * n))),
        (
            "t3(n^2 + 1) = -t3(n^2) when 3 does not divide n, t3(n^2) when it does",
            Box::new(|n| e(n * n + 1) == if n % 3 == 0 { e(n * n) } else { -e(n * n) }),
        ),
        (
            "t3 along n^2 + n: y(3n + 1) = -y(n)",
            Box::new(|n| e((3 * n + 1) * (3 * n + 2)) == -e(n * n + n)),
        ),
        (
            "t3 along (2n + 1)^2: y(3n + 1) = y(n)",
            Box::new(|n| e((6 * n + 3) * (6 * n + 3)) == e((2 * n + 1) * (2 * n + 1))),
        ),
        (
            "t3 along 2n^2 + 1: y(3n) = y(9n) = t3(2n^2)",
            Box::new(|n| {
                e(18 * n * n + 1) == e(162 * n * n + 1) && e(18 * n * n + 1) == e(2 * n * n)
            }),
        ),
    ];
    println!("identities between readings, checked at every n < 2^24:");
    for (name, f) in &checks {
        let bad = (0..top).filter(|&n| !f(n)).count();
        println!("  {name}: {bad} violations");
    }
    let tb = |m: u128| kit.bit(Design::Thue, m);
    let eb = |m: u128| kit.bit(Design::Ternary, m);
    let mut link = 0;
    for k in 1..=20u32 {
        for n in 0..1024u128 {
            let x = (n << k) + 1;
            if tb(x * x) != 1 ^ tb((n * n << (k - 1)) + n) {
                link += 1;
            }
            let c = 3u128.pow(k);
            let w = c * n + 1;
            if eb(w * w) != eb(c * n * n + 2 * n) {
                link += 1;
            }
        }
    }
    println!("kernel links: t((2^k n + 1)^2) = 1 + t(2^(k-1) n^2 + n) and t3((3^k n + 1)^2) = t3(3^k n^2 + 2n), k <= 20, n < 1024: {link} violations");
    let (mut pairs, mut split) = (0, 0);
    for k in 3..40u32 {
        let n = (1u128 << (k - 1)) + 1;
        for kk in k + 1..=40u32 {
            pairs += 1;
            if tb((n * n << (k - 1)) + n) != tb((n * n << (kk - 1)) + n) {
                split += 1;
            }
        }
    }
    println!("thue-morse along n^2: g_k(n) = t(2^(k-1) n^2 + n) at n = 2^(k-1) + 1 differs from g_k' for {split} of {pairs} pairs 3 <= k < k' <= 40");
    let (mut pairs, mut split) = (0, 0);
    for k in 1..24u32 {
        let n = (3u128.pow(k) + 1) / 2;
        for kk in k + 1..=24u32 {
            pairs += 1;
            if eb(3u128.pow(k) * n * n + 2 * n) != eb(3u128.pow(kk) * n * n + 2 * n) {
                split += 1;
            }
        }
    }
    println!("ternary 2s along n^2: e_k(n) = t3(3^k n^2 + 2n) at n = (3^k + 1)/2 differs from e_k' for {split} of {pairs} pairs 1 <= k < k' <= 24");
}

fn main() {
    let args: Vec<String> = env::args().collect();
    match args.get(1).map(String::as_str) {
        Some("identities") => identities(),
        Some("battery") => battery(),
        Some("local") => local(),
        Some("falsify") => falsify(args.get(2).and_then(|s| s.parse().ok()).unwrap_or(24)),
        Some("leak") => leak(),
        _ => eprintln!(
            "verbs: identities | battery | local | falsify [even top exponent, default 24] | leak"
        ),
    }
}
