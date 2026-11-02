use crate::gen::{counter_words, transpose, Gen, Words, CTR, LANES, OUT};

// BARS

pub const SAC_LOG: usize = 14;
pub const MAIN_LOG: usize = 20;
pub const BM_WORDS: usize = 1024;
pub const Z: f64 = 4.5;
pub const RUN: usize = 5;
pub const WINDOW: usize = 8;
pub const BEFORE: usize = 8;
pub const SAC_BAR: f64 = 4.0 / 128.0;
pub const WALSH_BAR: f64 = 4.5 / 1024.0;
pub const BM_BAR: i64 = 10;
pub const RANK_BAR: f64 = 27.0;
pub const SCREEN_BAR: f64 = 2048.0 + 4.5 * 64.0;
pub const TESTS: [&str; 6] = ["monobit", "sac", "walsh", "cube", "complexity", "rank"];

// SCREEN

pub fn screen_stat(base: &Words, flip: &[Words]) -> f64 {
    let mut s = 0.0;
    for (i, f) in flip.iter().enumerate() {
        let var = if i < 6 { 128.0 } else { 64.0 };
        for j in 0..OUT {
            let d = 2.0 * (base[j] ^ f[j]).count_ones() as f64 - 64.0;
            s += d * d / var;
        }
    }
    s
}

pub fn screen_null<G: Gen>(g: &G, batches: u64) -> (f64, f64, f64, usize) {
    let mut vals = Vec::new();
    for m in 0..batches {
        let outs: Vec<Words> = (0..=CTR)
            .map(|i| {
                let mut s = g.start(&counter_words(64 * m, i.checked_sub(1)));
                g.step(&mut s, 1);
                let mut o = [0u64; OUT];
                g.out(&s, &mut o);
                o
            })
            .collect();
        vals.push(screen_stat(&outs[0], &outs[1..]));
    }
    let n = vals.len() as f64;
    let mean = vals.iter().sum::<f64>() / n;
    let sd = (vals.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
    let max = vals.iter().copied().fold(0.0, f64::max);
    let over = vals.iter().filter(|&&v| v >= SCREEN_BAR).count();
    (mean, sd, max, over)
}

pub fn screen<G: Gen>(g: &G, cap: usize) -> Option<usize> {
    let mut states: Vec<G::State> = (0..=CTR)
        .map(|i| g.start(&counter_words(0, i.checked_sub(1))))
        .collect();
    let mut outs = vec![[0u64; OUT]; CTR + 1];
    let mut streak = 0;
    for r in 1..=cap + RUN - 1 {
        for (s, o) in states.iter_mut().zip(outs.iter_mut()) {
            g.step(s, r);
            g.out(s, o);
        }
        if screen_stat(&outs[0], &outs[1..]) < SCREEN_BAR {
            streak += 1;
            if streak == RUN {
                return Some(r + 1 - RUN);
            }
        } else {
            streak = 0;
            if r + 1 > cap {
                return None;
            }
        }
    }
    None
}

pub fn trace<G: Gen>(g: &G, rounds: usize) {
    let mut states: Vec<G::State> = (0..=CTR)
        .map(|i| g.start(&counter_words(0, i.checked_sub(1))))
        .collect();
    let mut outs = vec![[0u64; OUT]; CTR + 1];
    for r in 1..=rounds {
        for (s, o) in states.iter_mut().zip(outs.iter_mut()) {
            g.step(s, r);
            g.out(s, o);
        }
        let zero = (1..=CTR)
            .flat_map(|i| (0..OUT).map(move |j| (i, j)))
            .filter(|&(i, j)| outs[0][j] == outs[i][j])
            .count();
        let full = (1..=CTR)
            .flat_map(|i| (0..OUT).map(move |j| (i, j)))
            .filter(|&(i, j)| outs[0][j] ^ outs[i][j] == !0)
            .count();
        println!(
            "round {r}: screen {:.0}, pairs never flipped {zero}, always flipped {full}",
            screen_stat(&outs[0], &outs[1..])
        );
    }
}

pub const LITE_LOG: usize = 14;
pub const LITE_WALSH_BAR: f64 = 4.5 / 128.0;

pub fn lite<G: Gen>(g: &G, from: usize, cap: usize) -> Option<usize> {
    let hi = cap + RUN - 1;
    let mut ones = vec![0u64; hi + 1];
    let mut walsh = vec![[[0i64; LITE_LOG + 1]; OUT]; hi + 1];
    let mut ranks = vec![[0u64; 4]; hi + 1];
    let mut o = [0u64; OUT];
    for m in 0..1usize << (LITE_LOG - 6) {
        let mut s = g.start(&counter_words(64 * m as u64, None));
        for r in 1..=hi {
            g.step(&mut s, r);
            if r < from {
                continue;
            }
            g.out(&s, &mut o);
            for j in 0..OUT {
                let p = o[j].count_ones();
                ones[r] += p as u64;
                let s0 = 64 - 2 * p as i64;
                let w = &mut walsh[r][j];
                w[0] += s0;
                for b in 0..6 {
                    w[1 + b] += 64 - 2 * (o[j] ^ LANES[b]).count_ones() as i64;
                }
                for b in 6..LITE_LOG {
                    w[1 + b] += if m >> (b - 6) & 1 == 1 { -s0 } else { s0 };
                }
            }
            ranks[r][match rank(&o) {
                64 => 0,
                63 => 1,
                62 => 2,
                _ => 3,
            }] += 1;
        }
    }
    let n = (1u64 << LITE_LOG) as f64;
    let bits = n * 64.0;
    let probs = rank_probabilities();
    let k = n / 64.0;
    let ok = |r: usize| {
        let mono = ((2.0 * ones[r] as f64 - bits) / bits.sqrt()).abs();
        let w = walsh[r]
            .iter()
            .flat_map(|w| w.iter())
            .map(|&v| (v as f64 / n).abs())
            .fold(0.0, f64::max);
        let chi: f64 = (0..4)
            .map(|i| (ranks[r][i] as f64 - k * probs[i]).powi(2) / (k * probs[i]))
            .sum();
        mono < Z && w < LITE_WALSH_BAR && chi < RANK_BAR
    };
    (from..=cap).find(|&r| (r..r + RUN).all(ok))
}

// SAC

pub struct Sac {
    pub flips: Vec<Vec<u32>>,
}

impl Sac {
    pub fn run<G: Gen>(g: &G, horizon: usize) -> Sac {
        let mut flips = vec![vec![0u32; CTR * OUT]; horizon + 1];
        let mut base = [0u64; OUT];
        let mut o = [0u64; OUT];
        for m in 0..1u64 << (SAC_LOG - 6) {
            let mut states: Vec<G::State> = (0..=CTR)
                .map(|i| g.start(&counter_words(64 * m, i.checked_sub(1))))
                .collect();
            for r in 1..=horizon {
                for (i, s) in states.iter_mut().enumerate() {
                    g.step(s, r);
                    if i == 0 {
                        g.out(s, &mut base);
                    } else {
                        g.out(s, &mut o);
                        let row = &mut flips[r][(i - 1) * OUT..i * OUT];
                        for j in 0..OUT {
                            row[j] += (base[j] ^ o[j]).count_ones();
                        }
                    }
                }
            }
        }
        Sac { flips }
    }
    pub fn deviation(&self, r: usize) -> f64 {
        let n = (1u64 << SAC_LOG) as f64;
        self.flips[r]
            .iter()
            .map(|&f| (f as f64 / n - 0.5).abs())
            .fold(0.0, f64::max)
    }
    pub fn reached(&self, r: usize) -> bool {
        self.flips[r].iter().all(|&f| f > 0)
    }
}

// MAIN RUN

#[derive(Clone)]
pub struct Acc {
    pub ones: u64,
    pub walsh: Vec<[i64; 21]>,
    pub cubes: [u64; 21],
    chain: [u64; 21],
    pub ranks: [u64; 4],
    pub stream: Vec<u64>,
}

impl Acc {
    fn new() -> Acc {
        Acc {
            ones: 0,
            walsh: vec![[0; 21]; OUT],
            cubes: [0; 21],
            chain: [0; 21],
            ranks: [0; 4],
            stream: Vec::new(),
        }
    }
    fn add(&mut self, m: usize, o: &Words) {
        let mut w6 = 0u64;
        for j in 0..OUT {
            let x = o[j];
            let p = x.count_ones();
            self.ones += p as u64;
            let s0 = 64 - 2 * p as i64;
            let w = &mut self.walsh[j];
            w[0] += s0;
            for b in 0..6 {
                w[1 + b] += 64 - 2 * (x ^ LANES[b]).count_ones() as i64;
            }
            for b in 6..20 {
                w[1 + b] += if m >> (b - 6) & 1 == 1 { -s0 } else { s0 };
            }
            let mut f = x ^ (x >> 1);
            f ^= f >> 2;
            self.cubes[2] += (f & 0x1111_1111_1111_1111).count_ones() as u64;
            f ^= f >> 4;
            self.cubes[3] += (f & 0x0101_0101_0101_0101).count_ones() as u64;
            f ^= f >> 8;
            self.cubes[4] += (f & 0x0001_0001_0001_0001).count_ones() as u64;
            f ^= f >> 16;
            self.cubes[5] += (f & 0x0000_0001_0000_0001).count_ones() as u64;
            f ^= f >> 32;
            self.cubes[6] += f & 1;
            w6 |= (f & 1) << j;
        }
        for d in 7..=20 {
            self.chain[d] ^= w6;
            if (m + 1) % (1 << (d - 6)) == 0 {
                self.cubes[d] += self.chain[d].count_ones() as u64;
                self.chain[d] = 0;
            }
        }
        let r = rank(o);
        self.ranks[match r {
            64 => 0,
            63 => 1,
            62 => 2,
            _ => 3,
        }] += 1;
        if self.stream.len() < BM_WORDS {
            self.stream.extend_from_slice(&transpose(o));
        }
    }
}

pub fn rank(o: &Words) -> usize {
    let mut rows = *o;
    let mut r = 0;
    for bit in 0..64 {
        let Some(p) = (r..64).find(|&i| rows[i] >> bit & 1 == 1) else {
            continue;
        };
        rows.swap(r, p);
        let pivot = rows[r];
        for row in rows.iter_mut().skip(r + 1) {
            if *row >> bit & 1 == 1 {
                *row ^= pivot;
            }
        }
        r += 1;
    }
    r
}

pub fn rank_probability(n: i32, r: i32) -> f64 {
    let mut lg = (r * (2 * n - r) - n * n) as f64 * 2f64.ln();
    for i in 0..r {
        lg += 2.0 * (1.0 - 2f64.powi(i - n)).ln() - (1.0 - 2f64.powi(i - r)).ln();
    }
    lg.exp()
}

pub fn rank_probabilities() -> [f64; 4] {
    let (a, b, c) = (
        rank_probability(64, 64),
        rank_probability(64, 63),
        rank_probability(64, 62),
    );
    [a, b, c, 1.0 - a - b - c]
}

pub fn complexity(words: &[u64]) -> usize {
    let n = words.len() * 64;
    let rev: Vec<u64> = (0..words.len())
        .map(|w| words[words.len() - 1 - w].reverse_bits())
        .chain(std::iter::once(0))
        .collect();
    let cap = n / 64 + 2;
    let mut c = vec![0u64; cap];
    let mut b = vec![0u64; cap];
    c[0] = 1;
    b[0] = 1;
    let (mut l, mut shift, mut lb) = (0usize, 1usize, 0usize);
    let window = |off: usize, w: usize| -> u64 {
        let bit = off + 64 * w;
        let (q, s) = (bit / 64, bit % 64);
        let lo = rev.get(q).copied().unwrap_or(0) >> s;
        if s == 0 {
            lo
        } else {
            lo | rev.get(q + 1).copied().unwrap_or(0) << (64 - s)
        }
    };
    for big in 0..n {
        let off = n - 1 - big;
        let mut d = 0u32;
        for w in 0..=l / 64 {
            d ^= (c[w] & window(off, w)).count_ones();
        }
        if d & 1 == 0 {
            shift += 1;
            continue;
        }
        let old = (2 * l <= big).then(|| c.clone());
        let (q, s) = (shift / 64, shift % 64);
        for w in (0..=lb / 64).rev() {
            if w + q < cap {
                c[w + q] ^= b[w] << s;
            }
            if s != 0 && w + q + 1 < cap {
                c[w + q + 1] ^= b[w] >> (64 - s);
            }
        }
        match old {
            Some(t) => {
                lb = l;
                l = big + 1 - l;
                b = t;
                shift = 1;
            }
            None => shift += 1,
        }
    }
    l
}

// VERDICT

pub struct Verdict {
    pub rstar: Option<usize>,
    pub cone: Option<usize>,
    pub sac_from: Option<usize>,
    pub first: [Option<usize>; 6],
    pub stats: [f64; 6],
    pub horizon: usize,
}

pub struct Rounds {
    pub accs: Vec<Acc>,
}

pub fn main_run<G: Gen>(g: &G, lo: usize, hi: usize) -> Rounds {
    let mut accs = vec![Acc::new(); hi + 1 - lo];
    let mut o = [0u64; OUT];
    for m in 0..1usize << (MAIN_LOG - 6) {
        let mut s = g.start(&counter_words(64 * m as u64, None));
        for r in 1..=hi {
            g.step(&mut s, r);
            if r >= lo {
                g.out(&s, &mut o);
                accs[r - lo].add(m, &o);
            }
        }
    }
    Rounds { accs }
}

pub fn stats(acc: &Acc, sac: &Sac, r: usize) -> [f64; 5] {
    let bits = (1u64 << MAIN_LOG) as f64 * 64.0;
    let mono = (2.0 * acc.ones as f64 - bits) / bits.sqrt();
    let n = (1u64 << MAIN_LOG) as f64;
    let walsh = acc
        .walsh
        .iter()
        .flat_map(|w| w.iter())
        .map(|&v| (v as f64 / n).abs())
        .fold(0.0, f64::max);
    let cube = (2..=20)
        .map(|d| {
            let m = 64.0 * (1u64 << (20 - d)) as f64;
            ((2.0 * acc.cubes[d] as f64 - m) / m.sqrt()).abs()
        })
        .fold(0.0, f64::max);
    let probs = rank_probabilities();
    let k = acc.ranks.iter().sum::<u64>() as f64;
    let chi = (0..4)
        .map(|i| (acc.ranks[i] as f64 - k * probs[i]).powi(2) / (k * probs[i]))
        .sum();
    [mono.abs(), sac.deviation(r), walsh, cube, chi]
}

pub fn passes(st: &[f64; 5]) -> [bool; 5] {
    [
        st[0] < Z,
        st[1] < SAC_BAR,
        st[2] < WALSH_BAR,
        st[3] < Z,
        st[4] < RANK_BAR,
    ]
}

pub fn battery<G: Gen>(g: &G, guess: usize, cap: usize, floor: usize) -> Verdict {
    if floor > cap {
        return Verdict {
            rstar: None,
            cone: None,
            sac_from: None,
            first: [None; 6],
            stats: [f64::NAN; 6],
            horizon: cap,
        };
    }
    let mut horizon = (guess.max(floor) + WINDOW + RUN).min(cap + RUN);
    loop {
        let sac = Sac::run(g, horizon);
        let ok = |r: usize| sac.deviation(r) < SAC_BAR;
        let cone = (1..=horizon).find(|&r| sac.reached(r));
        let sac_from = (1..=horizon + 1 - RUN).find(|&r| (r..r + RUN).all(ok));
        let Some(from) = sac_from else {
            if horizon >= cap + RUN - 1 {
                return Verdict {
                    rstar: None,
                    cone,
                    sac_from,
                    first: [None; 6],
                    stats: [f64::NAN; 6],
                    horizon,
                };
            }
            horizon = (2 * horizon).min(cap + RUN);
            continue;
        };
        let lo = from.saturating_sub(BEFORE).max(1);
        let rounds = main_run(g, lo, horizon);
        let st: Vec<[f64; 5]> = (lo..=horizon)
            .map(|r| stats(&rounds.accs[r - lo], &sac, r))
            .collect();
        let pass: Vec<[bool; 5]> = st.iter().map(passes).collect();
        let mut bm: Vec<Option<i64>> = vec![None; horizon + 1 - lo];
        let mut lc = |r: usize| -> i64 {
            *bm[r - lo].get_or_insert_with(|| {
                complexity(&rounds.accs[r - lo].stream) as i64 - (BM_WORDS * 32) as i64
            })
        };
        let all = |r: usize, lc: &mut dyn FnMut(usize) -> i64| {
            pass[r - lo].iter().all(|&p| p) && lc(r).abs() <= BM_BAR
        };
        let mut rstar = None;
        for r in from.max(floor)..=horizon + 1 - RUN {
            if (r..r + RUN).all(|t| pass[t - lo].iter().all(|&p| p))
                && (r..r + RUN).all(|t| all(t, &mut lc))
            {
                rstar = Some(r);
                break;
            }
        }
        let Some(rs) = rstar else {
            if horizon >= cap + RUN - 1 {
                let s = st[from - lo];
                return Verdict {
                    rstar: None,
                    cone,
                    sac_from,
                    first: [None; 6],
                    stats: [s[0], s[1], s[2], s[3], lc(from) as f64, s[4]],
                    horizon,
                };
            }
            horizon = (2 * horizon).min(cap + RUN);
            continue;
        };
        let end = rs + RUN - 1;
        let mut first = [None; 6];
        for (t, slot) in first.iter_mut().enumerate() {
            let mut f = end + 1;
            while f > lo {
                let r = f - 1;
                let ok = if t == 4 {
                    if r + RUN < rs {
                        break;
                    }
                    lc(r).abs() <= BM_BAR
                } else {
                    pass[r - lo][if t < 4 { t } else { 4 }]
                };
                if !ok {
                    break;
                }
                f = r;
            }
            *slot = Some(f);
        }
        let s = st[rs - lo];
        return Verdict {
            rstar: Some(rs),
            cone,
            sac_from,
            first,
            stats: [s[0], s[1], s[2], s[3], lc(rs) as f64, s[4]],
            horizon,
        };
    }
}

pub fn last(v: &Verdict) -> String {
    let top = v.first.iter().flatten().max().copied().unwrap_or(0);
    TESTS
        .iter()
        .zip(v.first.iter())
        .filter(|(_, f)| **f == Some(top))
        .map(|(n, _)| *n)
        .collect::<Vec<_>>()
        .join("+")
}

// CORRECTIONS

pub const SIDES: [usize; 2] = [8, 16];
pub const BIG_BAR: f64 = 24.0;
pub const BIG_FLOOR: usize = 10;
pub const EXTRA: [&str; 3] = ["full walsh", "rank 512", "rank 1024"];

pub fn full_bar() -> f64 {
    (2.0 * (2.0 * (1u64 << 26) as f64).ln()).sqrt() + 1.0
}

pub fn window<G: Gen>(g: &G, lo: usize, hi: usize) -> Vec<Vec<Words>> {
    let batches = 1usize << (MAIN_LOG - 6);
    let mut outs = vec![vec![[0u64; OUT]; batches]; hi + 1 - lo];
    for m in 0..batches {
        let mut s = g.start(&counter_words(64 * m as u64, None));
        for r in 1..=hi {
            g.step(&mut s, r);
            if r >= lo {
                g.out(&s, &mut outs[r - lo][m]);
            }
        }
    }
    outs
}

pub fn full_spectrum(outs: &[Words]) -> f64 {
    let batches = outs.len();
    let pat: [u64; 64] = std::array::from_fn(|a| {
        (0..6)
            .filter(|b| a >> b & 1 == 1)
            .fold(0, |p, b| p ^ LANES[b])
    });
    let mut f = vec![[0i32; 64]; batches];
    let mut best = 0i32;
    for j in 0..OUT {
        for (row, o) in f.iter_mut().zip(outs) {
            for (v, p) in row.iter_mut().zip(pat) {
                *v = 64 - 2 * (o[j] ^ p).count_ones() as i32;
            }
        }
        let mut h = 1;
        while h < batches {
            for i in (0..batches).step_by(2 * h) {
                let (lo, hi) = f[i..i + 2 * h].split_at_mut(h);
                for (u, v) in lo.iter_mut().zip(hi.iter_mut()) {
                    for a in 0..64 {
                        let (x, y) = (u[a], v[a]);
                        u[a] = x + y;
                        v[a] = x - y;
                    }
                }
            }
            h *= 2;
        }
        for row in &f {
            for &v in row {
                best = best.max(v.abs());
            }
        }
    }
    best as f64 / (64.0 * batches as f64).sqrt()
}

pub fn rank_wide(rows: &mut [u64], w: usize) -> usize {
    let n = rows.len() / w;
    let mut r = 0;
    for col in 0..64 * w {
        let (cw, cb) = (col / 64, col % 64);
        let Some(p) = (r..n).find(|&i| rows[i * w + cw] >> cb & 1 == 1) else {
            continue;
        };
        if p != r {
            for k in 0..w {
                rows.swap(r * w + k, p * w + k);
            }
        }
        for i in r + 1..n {
            if rows[i * w + cw] >> cb & 1 == 1 {
                for k in cw..w {
                    rows[i * w + k] ^= rows[r * w + k];
                }
            }
        }
        r += 1;
    }
    r
}

pub fn big_rank(outs: &[Words], w: usize) -> ([u64; 3], usize) {
    let stream: Vec<u64> = outs.iter().flat_map(transpose).collect();
    let side = 64 * w;
    let mut classes = [0u64; 3];
    let mut low = side;
    for chunk in stream.chunks_exact(side * w) {
        let mut rows = chunk.to_vec();
        let r = rank_wide(&mut rows, w);
        classes[(side - r).min(2)] += 1;
        low = low.min(r);
    }
    (classes, low)
}

pub fn big_chi(classes: &[u64; 3], side: usize) -> f64 {
    let n = side as i32;
    let (a, b) = (rank_probability(n, n), rank_probability(n, n - 1));
    let probs = [a, b, 1.0 - a - b];
    let k = classes.iter().sum::<u64>() as f64;
    (0..3)
        .map(|i| (classes[i] as f64 - k * probs[i]).powi(2) / (k * probs[i]))
        .sum()
}

#[derive(Clone)]
pub struct Extra {
    pub full: f64,
    pub big: [([u64; 3], f64, usize); 2],
}

pub fn extra(outs: &[Words]) -> Extra {
    let big = SIDES.map(|w| {
        let (c, low) = big_rank(outs, w);
        (c, big_chi(&c, 64 * w), low)
    });
    Extra {
        full: full_spectrum(outs),
        big,
    }
}

pub fn extra_pass(e: &Extra) -> [bool; 3] {
    let ok = |i: usize| e.big[i].1 < BIG_BAR && e.big[i].2 + BIG_FLOOR >= 64 * SIDES[i];
    [e.full < full_bar(), ok(0), ok(1)]
}

pub struct Corrected {
    pub v: Verdict,
    pub six: Option<usize>,
    pub extra: Option<Extra>,
    pub blocked: Vec<(usize, &'static str)>,
}

pub fn corrected<G: Gen>(g: &G, guess: usize, cap: usize) -> Corrected {
    let mut floor = 1;
    let mut blocked = Vec::new();
    let mut six = None;
    loop {
        let v = battery(g, guess, cap, floor);
        if floor == 1 {
            six = v.rstar;
        }
        let Some(rs) = v.rstar else {
            return Corrected {
                v,
                six,
                extra: None,
                blocked,
            };
        };
        let outs = window(g, rs, rs + RUN - 1);
        let mut at_rs = None;
        let mut fail = None;
        for (i, o) in outs.iter().enumerate() {
            let e = extra(o);
            let p = extra_pass(&e);
            if i == 0 {
                at_rs = Some(e);
            }
            if let Some(t) = p.iter().position(|&x| !x) {
                fail = Some((rs + i, EXTRA[t]));
                break;
            }
        }
        match fail {
            None => {
                return Corrected {
                    v,
                    six,
                    extra: at_rs,
                    blocked,
                }
            }
            Some((t, name)) => {
                blocked.push((t, name));
                floor = t + 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complexity_of_an_m_sequence_is_its_degree() {
        let mut s: u32 = 1;
        let mut words = vec![0u64; 1024];
        for t in 0..1024 * 64 {
            words[t / 64] |= ((s & 1) as u64) << (t % 64);
            let bit = (s ^ (s >> 1)) & 1;
            s = (s >> 1) | (bit << 6);
        }
        assert_eq!(complexity(&words), 7);
    }

    #[test]
    fn rank_law_matches_every_three_by_three_matrix() {
        let mut counts = [0u64; 4];
        for m in 0u64..512 {
            let mut rows = [0u64; 64];
            for (i, row) in rows.iter_mut().take(3).enumerate() {
                *row = m >> (3 * i) & 7;
            }
            counts[rank(&rows)] += 1;
        }
        for (r, &c) in counts.iter().enumerate() {
            assert!((512.0 * rank_probability(3, r as i32) - c as f64).abs() < 1e-6);
        }
    }
}
