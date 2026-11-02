use crate::design::{anf, offsets, Mess, Spread};

pub const OUT: usize = 64;
pub const CTR: usize = 32;
pub const SEED: u64 = 0x6d61_6368_696e_65;

pub type Words = [u64; OUT];

pub fn splitmix(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

pub fn fill_words(seed: Option<u64>, cells: usize) -> Vec<u64> {
    let Some(seed) = seed else {
        return vec![0; cells];
    };
    let mut state = seed;
    let words: Vec<u64> = (0..cells.div_ceil(64))
        .map(|_| splitmix(&mut state))
        .collect();
    (0..cells)
        .map(|p| {
            if words[p / 64] >> (p % 64) & 1 == 1 {
                !0
            } else {
                0
            }
        })
        .collect()
}

pub const LANES: [u64; 6] = [
    0xaaaa_aaaa_aaaa_aaaa,
    0xcccc_cccc_cccc_cccc,
    0xf0f0_f0f0_f0f0_f0f0,
    0xff00_ff00_ff00_ff00,
    0xffff_0000_ffff_0000,
    0xffff_ffff_0000_0000,
];

pub fn counter_words(base: u64, flip: Option<usize>) -> [u64; CTR] {
    let mut w = [0u64; CTR];
    for (b, word) in w.iter_mut().enumerate() {
        *word = if b < 6 {
            LANES[b]
        } else if base >> b & 1 == 1 {
            !0
        } else {
            0
        };
        if flip == Some(b) {
            *word = !*word;
        }
    }
    w
}

pub fn lane_counters(ctr: &[u64; CTR]) -> [u64; 64] {
    let mut c = [0u64; 64];
    for (b, word) in ctr.iter().enumerate() {
        for (l, v) in c.iter_mut().enumerate() {
            *v |= (word >> l & 1) << b;
        }
    }
    c
}

pub fn transpose(rows: &[u64; 64]) -> [u64; 64] {
    let mut t = [0u64; 64];
    for (i, r) in rows.iter().enumerate() {
        for (j, v) in t.iter_mut().enumerate() {
            *v |= (r >> j & 1) << i;
        }
    }
    t
}

pub trait Gen {
    type State: Clone;
    fn start(&self, ctr: &[u64; CTR]) -> Self::State;
    fn step(&self, s: &mut Self::State, round: usize);
    fn out(&self, s: &Self::State, o: &mut Words);
}

// MACHINE

#[derive(Clone, Copy)]
pub enum Load {
    Fill,
    Zero,
}

pub struct Machine {
    pub n: usize,
    pub k: usize,
    pub fill: Vec<bool>,
    writes: Vec<(u32, u8, [u32; 8])>,
    birth: [u64; 9],
    survive: Option<[u64; 9]>,
    pub src: Vec<[u32; 3]>,
    pub terms: usize,
    stamp: Vec<(usize, usize)>,
    pub out: [u32; OUT],
    pub fold: Option<[u32; OUT]>,
    pub ctr: [u32; CTR],
    seed: Vec<u64>,
}

#[derive(Clone)]
pub struct Grid {
    a: Vec<u64>,
    b: Vec<u64>,
}

impl Grid {
    pub fn cells(&self) -> &[u64] {
        &self.a
    }
}

pub fn block(n: usize, top: usize, left: usize, rows: usize, cols: usize) -> Vec<u32> {
    (0..rows * cols)
        .map(|b| (((top + b / cols) % n) * n + (left + b % cols) % n) as u32)
        .collect()
}

impl Machine {
    pub fn new(fill: Vec<bool>, hood: u32, mess: &Mess, spread: &Spread, load: Load) -> Machine {
        let cells = fill.len();
        let n = (cells as f64).sqrt() as usize;
        let offs = offsets(hood);
        let m = n as i64;
        let writes = (0..cells)
            .filter(|&a| fill[a])
            .map(|a| {
                let (i, j) = ((a / n) as i64, (a % n) as i64);
                let mut r = [0u32; 8];
                let mut c = 0u8;
                for &(di, dj) in &offs {
                    let q = ((i + di).rem_euclid(m) * m + (j + dj).rem_euclid(m)) as usize;
                    if !fill[q] {
                        r[c as usize] = q as u32;
                        c += 1;
                    }
                }
                (a as u32, c, r)
            })
            .collect();
        let (birth, survive) = match mess {
            Mess::Feistel(g) => (anf(g), None),
            Mess::Life { birth, survive } => (anf(birth), Some(anf(survive))),
        };
        let (src, terms) = spread.sources(n);
        let stamp = (0..cells)
            .filter(|&p| fill[p])
            .map(|p| (p / n, p % n))
            .collect();
        let out: [u32; OUT] = block(n, 0, 0, 8, 8).try_into().unwrap();
        let seed = fill_words(
            match load {
                Load::Fill => Some(SEED),
                Load::Zero => None,
            },
            cells,
        );
        let mut machine = Machine {
            n,
            k: offs.len(),
            fill,
            writes,
            birth,
            survive,
            src,
            terms,
            stamp,
            out,
            fold: None,
            ctr: [0; CTR],
            seed,
        };
        machine.place(machine.antipode());
        machine
    }
    pub fn antipode(&self) -> (usize, usize) {
        (self.n / 2 + 2, self.n / 2)
    }
    pub fn place(&mut self, at: (usize, usize)) {
        self.ctr = block(self.n, at.0, at.1, 4, 8).try_into().unwrap();
    }
    pub fn overlap(&self, at: (usize, usize)) -> usize {
        block(self.n, at.0, at.1, 4, 8)
            .iter()
            .filter(|&&p| self.fill[p as usize])
            .count()
    }
    pub fn loads(&self) -> [((usize, usize), usize); 3] {
        let (r0, c0) = self.antipode();
        let n = self.n as i64;
        let mut spots: Vec<(i64, (usize, usize), usize)> = Vec::new();
        for di in -2i64..=2 {
            for dj in -2i64..=2 {
                let at = (
                    (r0 as i64 + di).rem_euclid(n) as usize,
                    (c0 as i64 + dj).rem_euclid(n) as usize,
                );
                spots.push((di * di + dj * dj, at, self.overlap(at)));
            }
        }
        let pick = |key: &dyn Fn(usize) -> i64| {
            let s = spots.iter().min_by_key(|(d, _, o)| (key(*o), *d)).unwrap();
            (s.1, s.2)
        };
        [
            pick(&|o| (o as i64 - 16).abs()),
            pick(&|o| -(o as i64)),
            pick(&|o| o as i64),
        ]
    }
    pub fn structural(&self, cap: usize) -> Option<usize> {
        let cells = self.n * self.n;
        let mut x = vec![0u32; cells];
        let mut y = vec![0u32; cells];
        for (b, &p) in self.ctr.iter().enumerate() {
            x[p as usize] |= 1 << b;
        }
        for r in 1..=cap {
            for &(a, m, ref rd) in &self.writes {
                let mut v = 0u32;
                for &q in &rd[..m as usize] {
                    v |= x[q as usize];
                }
                x[a as usize] |= v;
            }
            for (q, s) in self.src.iter().enumerate() {
                y[q] = s[..self.terms].iter().fold(0, |v, &p| v | x[p as usize]);
            }
            std::mem::swap(&mut x, &mut y);
            let hit =
                |j: usize| x[self.out[j] as usize] | self.fold.map_or(0, |f| x[f[j] as usize]);
            if (0..OUT).all(|j| hit(j) == u32::MAX) {
                return Some(r);
            }
        }
        None
    }
    pub fn spread(&self, s: &mut Grid) {
        if self.terms == 1 {
            for (q, p) in self.src.iter().enumerate() {
                s.b[q] = s.a[p[0] as usize];
            }
        } else {
            for (q, p) in self.src.iter().enumerate() {
                s.b[q] = s.a[p[0] as usize] ^ s.a[p[1] as usize] ^ s.a[p[2] as usize];
            }
        }
        std::mem::swap(&mut s.a, &mut s.b);
    }
    pub fn feistel(&self, x: &mut [u64]) {
        for &(a, m, ref r) in &self.writes {
            let (mut c0, mut c1, mut c2, mut c3) = (0u64, 0u64, 0u64, 0u64);
            for &q in &r[..m as usize] {
                let v = x[q as usize];
                let t0 = c0 & v;
                c0 ^= v;
                let t1 = c1 & t0;
                c1 ^= t0;
                let t2 = c2 & t1;
                c2 ^= t1;
                c3 |= t2;
            }
            let born = eval(&self.birth, c0, c1, c2, c3);
            match &self.survive {
                None => x[a as usize] ^= born,
                Some(sv) => {
                    let kept = eval(sv, c0, c1, c2, c3);
                    let old = x[a as usize];
                    x[a as usize] = born ^ (old & (kept ^ born));
                }
            }
        }
    }
    pub fn stamp(&self, x: &mut [u64], round: usize) {
        let n = self.n;
        let (di, dj) = (round % n, (round / n) % n);
        for &(i, j) in &self.stamp {
            let p = ((i + di) % n) * n + (j + dj) % n;
            x[p] = !x[p];
        }
    }
    pub fn grid(&self, a: Vec<u64>) -> Grid {
        Grid {
            b: vec![0; a.len()],
            a,
        }
    }
    pub fn undo(&self, s: &mut Grid, round: usize) {
        self.stamp(&mut s.a, round);
        for _ in 1..self.n {
            self.spread(s);
        }
        self.feistel(&mut s.a);
    }
    pub fn written(&self) -> Vec<(usize, usize)> {
        self.writes
            .iter()
            .map(|w| (w.0 as usize, w.1 as usize))
            .collect()
    }
    pub fn write_count(&self) -> usize {
        self.writes.len()
    }
    pub fn read_count(&self) -> usize {
        self.writes.iter().map(|w| w.1 as usize).sum()
    }
}

#[inline(always)]
fn eval(sel: &[u64; 9], c0: u64, c1: u64, c2: u64, c3: u64) -> u64 {
    let c01 = c0 & c1;
    sel[0]
        ^ (sel[1] & c0)
        ^ (sel[2] & c1)
        ^ (sel[3] & c01)
        ^ (sel[4] & c2)
        ^ (sel[5] & c0 & c2)
        ^ (sel[6] & c1 & c2)
        ^ (sel[7] & c01 & c2)
        ^ (sel[8] & c3)
}

impl Gen for Machine {
    type State = Grid;
    fn start(&self, ctr: &[u64; CTR]) -> Grid {
        let mut a = self.seed.clone();
        for (b, &p) in self.ctr.iter().enumerate() {
            a[p as usize] ^= ctr[b];
        }
        self.grid(a)
    }
    fn step(&self, s: &mut Grid, round: usize) {
        self.feistel(&mut s.a);
        self.spread(s);
        self.stamp(&mut s.a, round);
    }
    fn out(&self, s: &Grid, o: &mut Words) {
        for (j, &p) in self.out.iter().enumerate() {
            o[j] = s.a[p as usize];
        }
        if let Some(f) = &self.fold {
            for (j, &p) in f.iter().enumerate() {
                o[j] ^= s.a[p as usize];
            }
        }
    }
}

// RULE 30

pub struct Rule30 {
    pub cells: usize,
    pub ctr: [u32; CTR],
    pub out: [u32; OUT],
    nb: Option<Vec<[u32; 3]>>,
    seed: Vec<u64>,
}

impl Rule30 {
    pub fn new(cells: usize, near: bool) -> Rule30 {
        let at = if near { cells - CTR } else { cells / 2 + 16 };
        Rule30 {
            cells,
            ctr: std::array::from_fn(|b| (at + b) as u32),
            out: std::array::from_fn(|j| j as u32),
            nb: None,
            seed: fill_words(Some(SEED), cells),
        }
    }
    pub fn slide(n: usize, at: (usize, usize)) -> Rule30 {
        let cells = n * n;
        let (src, _) = Spread::Map("slide").sources(n);
        let nb = src
            .iter()
            .map(|s| {
                let p = s[0] as usize;
                [
                    ((p + cells - 1) % cells) as u32,
                    p as u32,
                    ((p + 1) % cells) as u32,
                ]
            })
            .collect();
        Rule30 {
            cells,
            ctr: block(n, at.0, at.1, 4, 8).try_into().unwrap(),
            out: block(n, 0, 0, 8, 8).try_into().unwrap(),
            nb: Some(nb),
            seed: fill_words(Some(SEED), cells),
        }
    }
    pub fn grid(&self, a: Vec<u64>) -> Grid {
        Grid {
            b: vec![0; a.len()],
            a,
        }
    }
}

impl Gen for Rule30 {
    type State = Grid;
    fn start(&self, ctr: &[u64; CTR]) -> Grid {
        let mut a = self.seed.clone();
        for (b, &p) in self.ctr.iter().enumerate() {
            a[p as usize] ^= ctr[b];
        }
        self.grid(a)
    }
    fn step(&self, s: &mut Grid, _round: usize) {
        let n = self.cells;
        let (x, y) = (&s.a, &mut s.b);
        match &self.nb {
            Some(nb) => {
                for (q, t) in nb.iter().enumerate() {
                    y[q] = x[t[0] as usize] ^ (x[t[1] as usize] | x[t[2] as usize]);
                }
            }
            None => {
                y[0] = x[n - 1] ^ (x[0] | x[1]);
                for p in 1..n - 1 {
                    y[p] = x[p - 1] ^ (x[p] | x[p + 1]);
                }
                y[n - 1] = x[n - 2] ^ (x[n - 1] | x[0]);
            }
        }
        std::mem::swap(&mut s.a, &mut s.b);
    }
    fn out(&self, s: &Grid, o: &mut Words) {
        for (j, &p) in self.out.iter().enumerate() {
            o[j] = s.a[p as usize];
        }
    }
}

// CHACHA

const SIGMA: [u32; 4] = [0x6170_7865, 0x3320_646e, 0x7962_2d32, 0x6b20_6574];

fn quarter(s: &mut [u32; 16], a: usize, b: usize, c: usize, d: usize) {
    s[a] = s[a].wrapping_add(s[b]);
    s[d] = (s[d] ^ s[a]).rotate_left(16);
    s[c] = s[c].wrapping_add(s[d]);
    s[b] = (s[b] ^ s[c]).rotate_left(12);
    s[a] = s[a].wrapping_add(s[b]);
    s[d] = (s[d] ^ s[a]).rotate_left(8);
    s[c] = s[c].wrapping_add(s[d]);
    s[b] = (s[b] ^ s[c]).rotate_left(7);
}

pub fn chacha_block(key: &[u32; 8], counter: u32, nonce: &[u32; 3], rounds: usize) -> [u32; 16] {
    let mut init = [0u32; 16];
    init[..4].copy_from_slice(&SIGMA);
    init[4..12].copy_from_slice(key);
    init[12] = counter;
    init[13..].copy_from_slice(nonce);
    let mut s = init;
    for _ in 0..rounds / 2 {
        quarter(&mut s, 0, 4, 8, 12);
        quarter(&mut s, 1, 5, 9, 13);
        quarter(&mut s, 2, 6, 10, 14);
        quarter(&mut s, 3, 7, 11, 15);
        quarter(&mut s, 0, 5, 10, 15);
        quarter(&mut s, 1, 6, 11, 12);
        quarter(&mut s, 2, 7, 8, 13);
        quarter(&mut s, 3, 4, 9, 14);
    }
    for (v, i) in s.iter_mut().zip(init) {
        *v = v.wrapping_add(i);
    }
    s
}

pub struct ChaCha {
    pub rounds: usize,
    pub key: [u32; 8],
    pub nonce: [u32; 3],
}

impl ChaCha {
    pub fn new(rounds: usize) -> ChaCha {
        let mut st = SEED;
        let words: Vec<u64> = (0..6).map(|_| splitmix(&mut st)).collect();
        let half = |i: usize| (words[i / 2] >> (32 * (i % 2))) as u32;
        ChaCha {
            rounds,
            key: std::array::from_fn(half),
            nonce: std::array::from_fn(|i| half(8 + i)),
        }
    }
    #[inline(always)]
    pub fn value(&self, c: u32) -> u64 {
        let b = chacha_block(&self.key, c, &self.nonce, self.rounds);
        b[0] as u64 | (b[1] as u64) << 32
    }
}

impl Gen for ChaCha {
    type State = Words;
    fn start(&self, ctr: &[u64; CTR]) -> Words {
        let c = lane_counters(ctr);
        let rows: [u64; 64] = std::array::from_fn(|l| self.value(c[l] as u32));
        transpose(&rows)
    }
    fn step(&self, _s: &mut Words, _round: usize) {}
    fn out(&self, s: &Words, o: &mut Words) {
        *o = *s;
    }
}

// BASELINES

pub struct Counter {
    pub name: &'static str,
    table: Vec<[u64; 4]>,
    jumps: Vec<Vec<[u64; 4]>>,
}

fn xo_step(s: &mut [u64; 4]) -> u64 {
    let out = s[0].wrapping_add(s[3]);
    let t = s[1] << 17;
    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];
    s[2] ^= t;
    s[3] = s[3].rotate_left(45);
    out
}

fn apply(m: &[[u64; 4]], v: &[u64; 4]) -> [u64; 4] {
    let mut out = [0u64; 4];
    for (i, row) in m.iter().enumerate() {
        let bit = (row[0] & v[0]) ^ (row[1] & v[1]) ^ (row[2] & v[2]) ^ (row[3] & v[3]);
        out[i / 64] |= ((bit.count_ones() & 1) as u64) << (i % 64);
    }
    out
}

fn columns_to_rows(cols: &[[u64; 4]]) -> Vec<[u64; 4]> {
    let mut rows = vec![[0u64; 4]; 256];
    for (c, col) in cols.iter().enumerate() {
        for (r, row) in rows.iter_mut().enumerate() {
            row[c / 64] |= (col[r / 64] >> (r % 64) & 1) << (c % 64);
        }
    }
    rows
}

fn compose(a: &[[u64; 4]], b: &[[u64; 4]]) -> Vec<[u64; 4]> {
    let cols: Vec<[u64; 4]> = (0..256)
        .map(|c| {
            let mut e = [0u64; 4];
            e[c / 64] = 1 << (c % 64);
            apply(a, &apply(b, &e))
        })
        .collect();
    columns_to_rows(&cols)
}

impl Counter {
    pub fn xoshiro(log: usize) -> Counter {
        let mut st = SEED;
        let mut s = [0u64; 4];
        for w in s.iter_mut() {
            *w = splitmix(&mut st);
        }
        let mut table = Vec::with_capacity(1 << log);
        for _ in 0..1usize << log {
            table.push(s);
            for _ in 0..64 {
                xo_step(&mut s);
            }
        }
        let cols: Vec<[u64; 4]> = (0..256)
            .map(|c| {
                let mut e = [0u64; 4];
                e[c / 64] = 1 << (c % 64);
                for _ in 0..64 {
                    xo_step(&mut e);
                }
                e
            })
            .collect();
        let mut m = columns_to_rows(&cols);
        let mut jumps = Vec::new();
        for _ in 0..CTR {
            jumps.push(m.clone());
            m = compose(&m, &m);
        }
        Counter {
            name: "xoshiro256+ low bit",
            table,
            jumps,
        }
    }
    pub fn splitmix() -> Counter {
        Counter {
            name: "splitmix64",
            table: Vec::new(),
            jumps: Vec::new(),
        }
    }
    fn value(&self, c: u64) -> u64 {
        if self.table.is_empty() {
            let mut st = SEED ^ c.wrapping_mul(0xd1b5_4a32_d192_ed03);
            return splitmix(&mut st);
        }
        let log = self.table.len().trailing_zeros() as usize;
        let mut s = self.table[(c & ((1 << log) - 1)) as usize];
        for b in log..CTR {
            if c >> b & 1 == 1 {
                s = apply(&self.jumps[b], &s);
            }
        }
        let mut v = 0u64;
        for j in 0..64 {
            v |= (xo_step(&mut s) & 1) << j;
        }
        v
    }
}

impl Gen for Counter {
    type State = Words;
    fn start(&self, ctr: &[u64; CTR]) -> Words {
        let c = lane_counters(ctr);
        let rows: [u64; 64] = std::array::from_fn(|l| self.value(c[l]));
        transpose(&rows)
    }
    fn step(&self, _s: &mut Words, _round: usize) {}
    fn out(&self, s: &Words, o: &mut Words) {
        *o = *s;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quarter_round_matches_rfc_8439_section_2_1_1() {
        let mut s = [0u32; 16];
        s[..4].copy_from_slice(&[0x1111_1111, 0x0102_0304, 0x9b8d_6f43, 0x0123_4567]);
        quarter(&mut s, 0, 1, 2, 3);
        assert_eq!(s[..4], [0xea2a_92f4, 0xcb1c_f8ce, 0x4581_472e, 0x5881_c4bb]);
    }

    #[test]
    fn chacha20_block_matches_rfc_8439_section_2_3_2() {
        let key: [u32; 8] =
            std::array::from_fn(|i| u32::from_le_bytes(std::array::from_fn(|b| (4 * i + b) as u8)));
        let nonce = [0x0900_0000, 0x4a00_0000, 0];
        let want = [
            0xe4e7_f110,
            0x1559_3bd1,
            0x1fdd_0f50,
            0xc471_20a3,
            0xc7f4_d1c7,
            0x0368_c033,
            0x9aaa_2204,
            0x4e6c_d4c3,
            0x4664_82d2,
            0x09aa_9f07,
            0x05d7_c214,
            0xa202_8bd9,
            0xd19c_12b5,
            0xb94e_16de,
            0xe883_d0cb,
            0x4e3c_50a2,
        ];
        assert_eq!(chacha_block(&key, 1, &nonce, 20), want);
    }
}
