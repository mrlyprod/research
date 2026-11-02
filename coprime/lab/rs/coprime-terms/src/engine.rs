use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Instant;

use crate::design::Design;

// CHUNKS

const CHUNK: u32 = 9;
const CHUNK_SPAN: u64 = 19683;

fn chunk_table(invert: bool) -> Vec<u32> {
    let mut out = vec![0u32; CHUNK_SPAN as usize];
    for value in 0..CHUNK_SPAN {
        let mut mask = 0u32;
        let mut rest = value;
        for place in 0..CHUNK {
            let digit = rest % 3;
            rest /= 3;
            if (digit == 1) != invert {
                mask |= 1 << place;
            }
        }
        out[value as usize] = mask;
    }
    out
}

#[inline(always)]
fn mask_of(mut x: u64, table: &[u32], chunks: u32, full: u32) -> u32 {
    let mut mask = 0u32;
    let mut shift = 0u32;
    for _ in 0..chunks {
        mask |= table[(x % CHUNK_SPAN) as usize] << shift;
        x /= CHUNK_SPAN;
        shift += CHUNK;
    }
    mask & full
}

// SIEVE

pub fn mobius(limit: usize) -> Vec<i8> {
    let mut mu = vec![0i8; limit];
    if limit > 1 {
        mu[1] = 1;
    }
    let mut composite = vec![false; limit];
    let mut primes: Vec<usize> = Vec::new();
    for value in 2..limit {
        if !composite[value] {
            primes.push(value);
            mu[value] = -1;
        }
        for &prime in primes.iter() {
            let product = value * prime;
            if product >= limit {
                break;
            }
            composite[product] = true;
            if value % prime == 0 {
                mu[product] = 0;
                break;
            }
            mu[product] = -mu[value];
        }
    }
    mu
}

pub fn primes_to(limit: u64) -> Vec<u64> {
    let size = limit as usize + 1;
    let mut composite = vec![false; size];
    let mut out = Vec::new();
    for value in 2..size {
        if !composite[value] {
            out.push(value as u64);
            let mut multiple = value * value;
            while multiple < size {
                composite[multiple] = true;
                multiple += value;
            }
        }
    }
    out
}

pub fn mobius_range(lo: u64, hi: u64, primes: &[u64], mu: &mut Vec<i8>, rest: &mut Vec<u64>) {
    let width = (hi - lo) as usize;
    mu.clear();
    mu.resize(width, 1);
    rest.clear();
    rest.extend(lo..hi);
    for &prime in primes.iter() {
        if prime * prime >= hi {
            break;
        }
        let mut multiple = lo.div_ceil(prime) * prime;
        if multiple == 0 {
            multiple = prime;
        }
        while multiple < hi {
            let slot = (multiple - lo) as usize;
            mu[slot] = -mu[slot];
            rest[slot] /= prime;
            multiple += prime;
        }
        let square = prime * prime;
        let mut multiple = lo.div_ceil(square) * square;
        if multiple == 0 {
            multiple = square;
        }
        while multiple < hi {
            mu[(multiple - lo) as usize] = 0;
            multiple += square;
        }
    }
    for slot in 0..width {
        if mu[slot] != 0 && rest[slot] > 1 {
            mu[slot] = -mu[slot];
        }
    }
    if lo == 0 {
        mu[0] = 0;
    }
}

// TRANSFORMS

fn pass<T: Copy + std::ops::AddAssign, const BIT: usize>(buf: &mut [T]) {
    for chunk in buf.chunks_exact_mut(2 * BIT) {
        let (low, high) = chunk.split_at_mut(BIT);
        for index in 0..BIT {
            high[index] += low[index];
        }
    }
}

fn zeta<T: Copy + std::ops::AddAssign>(buf: &mut [T], level: u32) {
    for place in 0..level {
        match place {
            0 => pass::<T, 1>(buf),
            1 => pass::<T, 2>(buf),
            2 => pass::<T, 4>(buf),
            3 => pass::<T, 8>(buf),
            4 => pass::<T, 16>(buf),
            _ => {
                let bit = 1usize << place;
                for chunk in buf.chunks_exact_mut(2 * bit) {
                    let (low, high) = chunk.split_at_mut(bit);
                    for index in 0..bit {
                        high[index] += low[index];
                    }
                }
            }
        }
    }
}

// BITS

#[cfg(target_arch = "aarch64")]
#[inline(always)]
fn block_bits(masks: &[u32], probe: u32) -> u64 {
    use std::arch::aarch64::*;
    let quads = masks.len() / 4;
    let mut out = 0u64;
    unsafe {
        let wide = vdupq_n_u32(probe);
        let lanes = vld1q_u32([1u32, 2, 4, 8].as_ptr());
        for quad in 0..quads {
            let block = vld1q_u32(masks.as_ptr().add(4 * quad));
            let hit = vceqzq_u32(vandq_u32(block, wide));
            out |= (vaddvq_u32(vandq_u32(hit, lanes)) as u64) << (4 * quad);
        }
    }
    for (place, &mask) in masks.iter().enumerate().skip(4 * quads) {
        out |= (((mask & probe) == 0) as u64) << place;
    }
    out
}

#[cfg(not(target_arch = "aarch64"))]
#[inline(always)]
fn block_bits(masks: &[u32], probe: u32) -> u64 {
    let mut out = 0u64;
    for (half, group) in masks.chunks(32).enumerate() {
        let mut bits = 0u32;
        for (place, &mask) in group.iter().enumerate() {
            bits |= (((mask & probe) == 0) as u32) << place;
        }
        out |= (bits as u64) << (32 * half);
    }
    out
}

#[inline(always)]
fn above(position: usize) -> u64 {
    (!0u64 << (position & 63)) << 1
}

trait Ring: Copy {
    fn lift(value: u64) -> Self;
    fn zero() -> Self;
    fn plus(self, other: Self) -> Self;
    fn minus(self, other: Self) -> Self;
    fn times(self, other: Self) -> Self;
    fn widen(self) -> u128;
}

impl Ring for u64 {
    fn lift(value: u64) -> Self {
        value
    }
    fn zero() -> Self {
        0
    }
    fn plus(self, other: Self) -> Self {
        self.wrapping_add(other)
    }
    fn minus(self, other: Self) -> Self {
        self.wrapping_sub(other)
    }
    fn times(self, other: Self) -> Self {
        self.wrapping_mul(other)
    }
    fn widen(self) -> u128 {
        self as u128
    }
}

impl Ring for u128 {
    fn lift(value: u64) -> Self {
        value as u128
    }
    fn zero() -> Self {
        0
    }
    fn plus(self, other: Self) -> Self {
        self.wrapping_add(other)
    }
    fn minus(self, other: Self) -> Self {
        self.wrapping_sub(other)
    }
    fn times(self, other: Self) -> Self {
        self.wrapping_mul(other)
    }
    fn widen(self) -> u128 {
        self
    }
}

// THRESHOLDS

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Mode {
    Auto,
    Direct,
    Zeta,
    Convolve,
    Residue,
    Bitset,
    Rows,
    Cube,
}

#[derive(Clone, Copy, Debug)]
pub struct Caps {
    pub residue: u64,
    pub tail: bool,
    pub bitset: usize,
    pub rows: usize,
    pub legacy: Option<Mode>,
}

fn integer_root(value: usize, power: u32) -> usize {
    let mut root = (value as f64).powf(1.0 / power as f64) as usize + 2;
    while root > 1 && root.pow(power) > value {
        root -= 1;
    }
    root
}

pub fn caps(level: u32, dimension: usize) -> Caps {
    let size = 1usize << level;
    let sweep = 2 * level as usize * size;
    if dimension == 2 {
        Caps {
            residue: 64,
            tail: false,
            bitset: 0,
            rows: integer_root(sweep, 2),
            legacy: None,
        }
    } else {
        Caps {
            residue: 16,
            tail: true,
            bitset: integer_root(256 * sweep, 3),
            rows: integer_root(18 * sweep, 2),
            legacy: None,
        }
    }
}

pub fn caps_for(level: u32, dimension: usize, mode: Mode) -> Caps {
    let base = caps(level, dimension);
    match mode {
        Mode::Auto => base,
        Mode::Residue => Caps {
            residue: u64::MAX,
            ..base
        },
        Mode::Direct | Mode::Zeta | Mode::Convolve => {
            let mode = if dimension == 2 && mode == Mode::Convolve {
                Mode::Zeta
            } else {
                mode
            };
            Caps {
                legacy: Some(mode),
                ..base
            }
        }
        Mode::Bitset => Caps {
            tail: false,
            bitset: usize::MAX,
            ..base
        },
        Mode::Rows => Caps {
            tail: false,
            bitset: 0,
            rows: usize::MAX,
            ..base
        },
        Mode::Cube => Caps {
            tail: false,
            bitset: 0,
            rows: 0,
            ..base
        },
    }
}

// RESIDUE AUTOMATON

pub fn residue_count(design: &Design, level: u32, modulus: u64) -> u128 {
    let corners = design.corners();
    let width = modulus as usize;
    let mut step = vec![0usize; width * 3];
    for residue in 0..width {
        for digit in 0..3usize {
            step[residue * 3 + digit] = (3 * residue + digit) % width;
        }
    }
    let states = width.pow(design.dimension as u32);
    let mut cur = vec![0u128; states];
    let mut next = vec![0u128; states];
    cur[0] = 1;
    let coded: Vec<Vec<usize>> = corners
        .iter()
        .map(|v| v.iter().map(|d| *d as usize).collect())
        .collect();
    for _ in 0..level {
        next.iter_mut().for_each(|slot| *slot = 0);
        for (state, &weight) in cur.iter().enumerate() {
            if weight == 0 {
                continue;
            }
            let mut residues = [0usize; 4];
            let mut rest = state;
            for axis in (0..design.dimension).rev() {
                residues[axis] = rest % width;
                rest /= width;
            }
            for corner in coded.iter() {
                let mut target = 0usize;
                for axis in 0..design.dimension {
                    target = target * width + step[residues[axis] * 3 + corner[axis]];
                }
                next[target] += weight;
            }
        }
        std::mem::swap(&mut cur, &mut next);
    }
    cur[0]
}

// CONTEXT

const SHELF: u32 = 6;

struct Ctx {
    level: u32,
    size: usize,
    full: u32,
    chunks: u32,
    dimension: usize,
    count: Vec<u32>,
    masks: Vec<u32>,
    items: Vec<(u32, u32)>,
    imask: Vec<u32>,
    iweight: Vec<u32>,
    sub: Vec<(u32, u32)>,
    zbuf: Vec<u32>,
    zbuf16: Vec<u16>,
    rows: Vec<u64>,
    ranked: Vec<u64>,
    ranked32: Vec<u32>,
    extras: Vec<(u32, u32)>,
    binomial: Vec<Vec<u64>>,
    bucket: Vec<usize>,
}

impl Ctx {
    fn new(level: u32, dimension: usize) -> Ctx {
        let size = 1usize << level;
        let mut binomial = vec![vec![0u64; level as usize + 1]; level as usize + 1];
        for row in 0..=level as usize {
            binomial[row][0] = 1;
            for column in 1..=row {
                binomial[row][column] = binomial[row - 1][column - 1] + binomial[row - 1][column];
            }
        }
        Ctx {
            level,
            size,
            full: if level == 32 {
                u32::MAX
            } else {
                (1u32 << level) - 1
            },
            chunks: level.div_ceil(CHUNK).max(1),
            dimension,
            count: vec![0u32; size],
            masks: Vec::new(),
            items: Vec::new(),
            imask: Vec::new(),
            iweight: Vec::new(),
            sub: Vec::new(),
            zbuf: Vec::new(),
            zbuf16: Vec::new(),
            rows: Vec::new(),
            ranked: Vec::new(),
            ranked32: Vec::new(),
            extras: Vec::new(),
            binomial,
            bucket: Vec::new(),
        }
    }

    fn shelf(&self) -> u32 {
        self.level.min(SHELF)
    }

    fn gather(&mut self, modulus: u64, span: u64, table: &[u32]) {
        self.masks.clear();
        let mut x = 0u64;
        while x < span {
            self.masks.push(mask_of(x, table, self.chunks, self.full));
            x += modulus;
        }
    }

    fn dedup(&mut self) {
        self.items.clear();
        for &mask in self.masks.iter() {
            if self.count[mask as usize] == 0 {
                self.items.push((mask, 0));
            }
            self.count[mask as usize] += 1;
        }
        for slot in self.items.iter_mut() {
            slot.1 = self.count[slot.0 as usize];
        }
        for slot in self.items.iter() {
            self.count[slot.0 as usize] = 0;
        }
        let shelf = self.shelf();
        let shift = self.level - shelf;
        let buckets = 1usize << shelf;
        self.bucket.clear();
        self.bucket.resize(buckets + 1, 0);
        for &(mask, _) in self.items.iter() {
            self.bucket[(mask >> shift) as usize + 1] += 1;
        }
        for slot in 0..buckets {
            self.bucket[slot + 1] += self.bucket[slot];
        }
        let width = self.items.len();
        self.imask.clear();
        self.imask.resize(width, 0);
        self.iweight.clear();
        self.iweight.resize(width, 0);
        let mut cursor = self.bucket.clone();
        for &(mask, weight) in self.items.iter() {
            let slot = &mut cursor[(mask >> shift) as usize];
            self.imask[*slot] = mask;
            self.iweight[*slot] = weight;
            *slot += 1;
        }
    }

    fn direct2(&self) -> u128 {
        let mut total = 0u128;
        for &(a, ca) in self.items.iter() {
            let mut inner = 0u64;
            for &(b, cb) in self.items.iter() {
                inner += cb as u64 * ((a & b == 0) as u64);
            }
            total += ca as u128 * inner as u128;
        }
        total
    }

    fn direct3(&mut self) -> u128 {
        let mut total = 0u128;
        let items = std::mem::take(&mut self.items);
        for &(a, ca) in items.iter() {
            self.sub.clear();
            for &(b, cb) in items.iter() {
                if a & b == 0 {
                    self.sub.push((b, cb));
                }
            }
            let mut inner = 0u128;
            for index in 0..self.sub.len() {
                let (b, cb) = self.sub[index];
                let mut reach = 0u64;
                for &(c, cc) in self.sub.iter() {
                    reach += cc as u64 * ((b & c == 0) as u64);
                }
                inner += cb as u128 * reach as u128;
            }
            total += ca as u128 * inner;
        }
        self.items = items;
        total
    }

    fn spread(&mut self) {
        if self.zbuf.len() < self.size {
            self.zbuf = vec![0u32; self.size];
        }
        for &(mask, weight) in self.items.iter() {
            self.zbuf[mask as usize] = weight;
        }
        zeta(&mut self.zbuf[..self.size], self.level);
    }

    fn clear_spread(&mut self) {
        self.zbuf[..self.size].fill(0);
    }

    fn zeta2(&mut self) -> u128 {
        self.spread();
        let mut total = 0u128;
        for &(a, ca) in self.items.iter() {
            let free = !a & self.full;
            total += ca as u128 * self.zbuf[free as usize] as u128;
        }
        self.clear_spread();
        total
    }

    fn zeta3(&mut self, wide: bool) -> u128 {
        self.spread();
        let mut total = 0u128;
        for &(a, ca) in self.items.iter() {
            let mut inner = 0u128;
            if wide {
                for &(b, cb) in self.items.iter() {
                    if a & b != 0 {
                        continue;
                    }
                    let free = !(a | b) & self.full;
                    inner += cb as u128 * self.zbuf[free as usize] as u128;
                }
            } else {
                let mut narrow = 0u64;
                for &(b, cb) in self.items.iter() {
                    if a & b != 0 {
                        continue;
                    }
                    let free = !(a | b) & self.full;
                    narrow += cb as u64 * self.zbuf[free as usize] as u64;
                }
                inner = narrow as u128;
            }
            total += ca as u128 * inner;
        }
        self.clear_spread();
        total
    }

    fn convolve3(&mut self, wide: bool) -> u128 {
        let ranks = self.level as usize + 1;
        if self.ranked.len() < ranks * self.size {
            self.ranked = vec![0u64; ranks * self.size];
        }
        self.ranked[..ranks * self.size].fill(0);
        for &(mask, weight) in self.items.iter() {
            let rank = mask.count_ones() as usize;
            self.ranked[rank * self.size + mask as usize] += weight as u64;
        }
        for rank in 0..ranks {
            let slice = &mut self.ranked[rank * self.size..(rank + 1) * self.size];
            zeta(slice, self.level);
        }
        let top = self.level as usize;
        let mut poly = vec![0u64; ranks];
        let mut square = vec![0u64; 2 * ranks];
        let mut total = 0i128;
        for set in 0..self.size {
            let rank = (set as u32).count_ones() as usize;
            for (index, slot) in poly.iter_mut().enumerate().take(rank + 1) {
                *slot = self.ranked[index * self.size + set];
            }
            for slot in square[..=(2 * rank).min(top)].iter_mut() {
                *slot = 0;
            }
            for left in 0..=rank {
                let value = poly[left];
                if value == 0 {
                    continue;
                }
                let bound = rank.min(top - left);
                for right in 0..=bound {
                    square[left + right] += value * poly[right];
                }
            }
            let mut acc = 0i128;
            for degree in rank..=top {
                let lower = degree.saturating_sub(rank);
                let upper = degree.min(2 * rank);
                let cube: i128 = if wide {
                    let mut sum = 0u128;
                    for left in lower..=upper {
                        sum += square[left] as u128 * poly[degree - left] as u128;
                    }
                    sum as i128
                } else {
                    let mut sum = 0u64;
                    for left in lower..=upper {
                        sum += square[left] * poly[degree - left];
                    }
                    sum as i128
                };
                let weight = self.binomial[top - rank][degree - rank] as i128;
                if (degree - rank).is_multiple_of(2) {
                    acc += weight * cube;
                } else {
                    acc -= weight * cube;
                }
            }
            total += acc;
        }
        assert!(total >= 0);
        total as u128
    }

    fn tail3(&self) -> u128 {
        let reach = self.masks.len() as u128;
        let zero = self.masks.iter().filter(|m| **m == 0).count() as u128;
        let pair = if reach == 3 {
            (self.masks[1] & self.masks[2] == 0) as u128
        } else {
            0
        };
        6 * pair + 3 * zero * (reach - 1) + zero
    }

    fn bitset3(&mut self) -> u128 {
        let reach = self.masks.len();
        assert!(reach <= 1 << 14);
        let words = reach.div_ceil(64);
        self.rows.clear();
        self.rows.resize(reach * words, 0);
        for i in 0..reach {
            let probe = self.masks[i];
            for (word, block) in self.masks.chunks(64).enumerate() {
                self.rows[i * words + word] = block_bits(block, probe);
            }
        }
        let zero = self.masks.iter().filter(|m| **m == 0).count() as u128;
        let mut triples = 0u128;
        for i in 0..reach {
            let row = &self.rows[i * words..(i + 1) * words];
            let mut acc = 0u64;
            for word in i / 64..words {
                let mut bits = row[word];
                if word == i / 64 {
                    bits &= above(i);
                }
                while bits != 0 {
                    let j = word * 64 + bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    let other = &self.rows[j * words..(j + 1) * words];
                    let head = j / 64;
                    let mut found = (row[head] & other[head] & above(j)).count_ones();
                    for slot in head + 1..words {
                        found += (row[slot] & other[slot]).count_ones();
                    }
                    acc += found as u64;
                }
            }
            triples += acc as u128;
        }
        6 * triples + 3 * zero * (reach as u128 - 1) + zero
    }

    fn rows_pass<G: Copy + Into<u64>>(&self, table: &[G]) -> u128 {
        let width = self.imask.len();
        let full = self.full;
        let shelf = self.shelf();
        let shift = self.level - shelf;
        let top = (1u32 << shelf) - 1;
        let mut total = 0u128;
        let mut diagonal = 0u128;
        for i in 0..width {
            let probe = self.imask[i];
            let wi = self.iweight[i] as u128;
            if probe == 0 {
                diagonal += wi * wi * table[full as usize].into() as u128;
            }
            let own = probe >> shift;
            let room = !own & top;
            let mut acc = 0u64;
            let mut sub = 0u32;
            loop {
                if sub >= own {
                    let lo = self.bucket[sub as usize].max(i + 1);
                    let hi = self.bucket[sub as usize + 1];
                    let mut base = lo;
                    while base < hi {
                        let stop = (base + 64).min(hi);
                        let mut bits = block_bits(&self.imask[base..stop], probe);
                        while bits != 0 {
                            let j = base + bits.trailing_zeros() as usize;
                            bits &= bits - 1;
                            let free = !(probe | self.imask[j]) & full;
                            acc += self.iweight[j] as u64 * table[free as usize].into();
                        }
                        base = stop;
                    }
                }
                if sub == room {
                    break;
                }
                sub = (sub.wrapping_sub(room)) & room;
            }
            total += wi * acc as u128;
        }
        2 * total + diagonal
    }

    fn rows3(&mut self) -> u128 {
        let reach = self.masks.len();
        if reach < 65536 {
            if self.zbuf16.len() < self.size {
                self.zbuf16 = vec![0u16; self.size];
            }
            for &(mask, weight) in self.items.iter() {
                self.zbuf16[mask as usize] = weight as u16;
            }
            zeta(&mut self.zbuf16[..self.size], self.level);
            let total = self.rows_pass(&self.zbuf16[..self.size]);
            self.zbuf16[..self.size].fill(0);
            total
        } else {
            self.spread();
            let total = self.rows_pass(&self.zbuf[..self.size]);
            self.clear_spread();
            total
        }
    }

    fn pick_rank(&self) -> usize {
        let top = self.level as usize;
        let mut histogram = vec![0u64; top + 1];
        for &mask in self.imask.iter() {
            histogram[mask.count_ones() as usize] += 1;
        }
        let slice = (self.size as u64) * 5;
        let mut best = top;
        let mut best_cost = u64::MAX;
        for cap in top / 2..=top {
            let mut cost = (cap as u64 + 1) * slice;
            for (rank, &count) in histogram.iter().enumerate().take(top + 1).skip(cap + 1) {
                cost += count * (1u64 << (top - rank)) * 12;
            }
            if cost < best_cost {
                best_cost = cost;
                best = cap;
            }
        }
        best
    }

    fn cube_pass<R: Ring>(&self, cap: usize) -> u128 {
        let top = self.level as usize;
        let size = self.size;
        let mut left = [0u64; 33];
        let mut right = [0u64; 33];
        let mut table = [R::zero(); 65];
        let mut total = R::zero();
        for set in 0..size / 2 {
            let other = (size - 1) ^ set;
            let rank = (set as u32).count_ones() as usize;
            let degree = rank.min(cap);
            let mut sum = 0u64;
            for (index, slot) in left.iter_mut().enumerate().take(degree + 1) {
                *slot = self.ranked32[index * size + set] as u64;
                sum += *slot;
            }
            if sum == 0 {
                continue;
            }
            let degree_other = (top - rank).min(cap);
            let mut sum = 0u64;
            for (index, slot) in right.iter_mut().enumerate().take(degree_other + 1) {
                *slot = self.ranked32[index * size + other] as u64;
                sum += *slot;
            }
            if sum == 0 {
                continue;
            }
            total = total.plus(pair_term(
                &left[..=degree],
                &right[..=degree_other],
                rank,
                top,
                &mut table,
            ));
            total = total.plus(pair_term(
                &right[..=degree_other],
                &left[..=degree],
                top - rank,
                top,
                &mut table,
            ));
        }
        total.widen()
    }

    fn cube3(&mut self) -> u128 {
        let cap = self.pick_rank();
        let slices = cap + 1;
        let size = self.size;
        if self.ranked32.len() < slices * size {
            self.ranked32 = vec![0u32; slices * size];
        }
        self.ranked32[..slices * size].fill(0);
        self.extras.clear();
        for &(mask, weight) in self.items.iter() {
            let rank = mask.count_ones() as usize;
            if rank <= cap {
                self.ranked32[rank * size + mask as usize] = weight;
            } else {
                self.extras.push((mask, weight));
            }
        }
        for rank in 0..slices {
            zeta(
                &mut self.ranked32[rank * size..(rank + 1) * size],
                self.level,
            );
        }
        let reach = self.masks.len() as u128;
        let total = if reach * reach * reach < 1u128 << 64 {
            self.cube_pass::<u64>(cap)
        } else {
            self.cube_pass::<u128>(cap)
        };
        let mut extra = 0u128;
        for index in 0..self.extras.len() {
            let (mask, weight) = self.extras[index];
            let room = !mask & self.full;
            let mut pairs = 0u128;
            let mut sub = room;
            loop {
                let rank = sub.count_ones() as usize;
                let weight_sub = self.ranked32[rank * size + sub as usize] as u128;
                if weight_sub != 0 {
                    let rest = room & !sub;
                    let bound = (rest.count_ones() as usize).min(cap);
                    let mut inside = 0u64;
                    for slot in 0..=bound {
                        inside += self.ranked32[slot * size + rest as usize] as u64;
                    }
                    pairs += weight_sub * inside as u128;
                }
                if sub == 0 {
                    break;
                }
                sub = (sub - 1) & room;
            }
            extra += weight as u128 * pairs;
        }
        total + 3 * extra
    }

    fn measure(&mut self, modulus: u64, span: u64, table: &[u32], caps: &Caps) -> (u128, usize) {
        self.gather(modulus, span, table);
        let reach = self.masks.len();
        if let Some(mode) = caps.legacy {
            self.dedup();
            let wide_pair = (reach as u128).pow(2) > 1u128 << 62;
            let wide_cube = (reach as u128).pow(3) > 1u128 << 62;
            return match (self.dimension, mode) {
                (2, Mode::Direct) => (self.direct2(), 1),
                (2, _) => (self.zeta2(), 2),
                (_, Mode::Direct) => (self.direct3(), 1),
                (_, Mode::Zeta) => (self.zeta3(wide_pair), 2),
                _ => (self.convolve3(wide_cube), 3),
            };
        }
        if self.dimension == 2 {
            self.dedup();
            return if self.items.len() <= caps.rows {
                (self.direct2(), 1)
            } else {
                (self.zeta2(), 2)
            };
        }
        if caps.tail && reach <= 3 {
            return (self.tail3(), 4);
        }
        if reach <= caps.bitset {
            return (self.bitset3(), 1);
        }
        self.dedup();
        if self.items.len() <= caps.rows {
            (self.rows3(), 2)
        } else {
            (self.cube3(), 3)
        }
    }
}

fn pair_term<R: Ring>(
    poly: &[u64],
    other: &[u64],
    rank: usize,
    top: usize,
    table: &mut [R; 65],
) -> R {
    let degree = poly.len() - 1;
    if 2 * degree < rank {
        return R::zero();
    }
    let length = 2 * degree - rank + 1;
    for (slot, cell) in table.iter_mut().enumerate().take(length) {
        let sum = rank + slot;
        let mut i = sum.saturating_sub(degree);
        let mut k = sum - i;
        let mut half = R::zero();
        while i < k {
            half = half.plus(R::lift(poly[i]).times(R::lift(poly[k])));
            i += 1;
            k -= 1;
        }
        let mut value = half.plus(half);
        if i == k {
            value = value.plus(R::lift(poly[i]).times(R::lift(poly[i])));
        }
        *cell = value;
    }
    let spread = top - rank;
    let lowest = spread + 1 - other.len();
    let mut acc = R::zero();
    for order in 0..=spread {
        if order >= lowest {
            acc = acc.plus(R::lift(other[spread - order]).times(table[0]));
        }
        if order == spread {
            break;
        }
        for slot in 0..length - 1 {
            table[slot] = table[slot].minus(table[slot + 1]);
        }
    }
    acc
}

// PROFILE

pub const METHODS: [&str; 5] = ["residue", "bitset", "rows", "cube", "tail"];

#[derive(Clone, Copy, Debug, Default)]
pub struct Cell {
    pub moduli: u64,
    pub nanos: u64,
}

pub struct Profile {
    pub cells: std::sync::Mutex<Vec<Cell>>,
}

impl Default for Profile {
    fn default() -> Profile {
        Profile {
            cells: std::sync::Mutex::new(vec![Cell::default(); 40 * METHODS.len()]),
        }
    }
}

impl Profile {
    pub fn print(&self, level: u32) {
        let cells = self.cells.lock().unwrap();
        println!(
            "level {} band(log2 Y) method moduli seconds ns/modulus",
            level
        );
        let mut total = 0u64;
        for band in 0..40usize {
            for (index, name) in METHODS.iter().enumerate() {
                let cell = cells[band * METHODS.len() + index];
                if cell.moduli == 0 {
                    continue;
                }
                total += cell.nanos;
                println!(
                    "{} {} {} {:.3} {}",
                    band,
                    name,
                    cell.moduli,
                    cell.nanos as f64 * 1e-9,
                    cell.nanos / cell.moduli
                );
            }
        }
        println!("cpu seconds {:.3}", total as f64 * 1e-9);
    }
}

// LEVELS

#[derive(Clone, Copy, Debug)]
pub struct Level {
    pub level: u32,
    pub value: i128,
    pub seconds: f64,
}

struct Sieve {
    fine_mu: Vec<i8>,
    primes: Vec<u64>,
    table: Vec<u32>,
}

impl Sieve {
    fn new(design: &Design, top: u32) -> Sieve {
        let span = 3u64.pow(top);
        Sieve {
            fine_mu: mobius(std::cmp::min(span, 1 << 16) as usize + 1),
            primes: primes_to((span as f64).sqrt() as u64 + 2),
            table: chunk_table(design.invert),
        }
    }
}

fn weight(
    design: &Design,
    level: u32,
    sieve: &Sieve,
    threads: usize,
    mode: Mode,
    profile: Option<&Profile>,
) -> i128 {
    let (fine_mu, primes, table) = (&sieve.fine_mu[..], &sieve.primes[..], &sieve.table[..]);
    let span = 3u64.pow(level);
    let caps = caps_for(level, design.dimension, mode);
    let origin: i128 = if design.zero_filled() { 1 } else { 0 };
    let fine = std::cmp::min(span, fine_mu.len() as u64 - 1);
    let block = 8192u64;
    let cursor = AtomicUsize::new(0);
    let total: i128 = std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..threads {
            let cursor = &cursor;
            let caps = &caps;
            handles.push(scope.spawn(move || {
                let mut ctx = Ctx::new(level, design.dimension);
                let mut block_mu: Vec<i8> = Vec::new();
                let mut block_rest: Vec<u64> = Vec::new();
                let mut acc = 0i128;
                let mut local = vec![Cell::default(); 40 * METHODS.len()];
                loop {
                    let task = cursor.fetch_add(1, Ordering::Relaxed) as u64;
                    let (lo, hi) = if task < fine {
                        (task + 1, task + 2)
                    } else {
                        let base = fine + 1 + (task - fine) * block;
                        (base, base + block)
                    };
                    if lo >= span {
                        break;
                    }
                    let hi = std::cmp::min(hi, span);
                    if task >= fine {
                        mobius_range(lo, hi, primes, &mut block_mu, &mut block_rest);
                    }
                    for modulus in lo..hi {
                        let sign = if task < fine {
                            fine_mu[modulus as usize]
                        } else {
                            block_mu[(modulus - lo) as usize]
                        };
                        if sign == 0 || modulus % 3 == 0 {
                            continue;
                        }
                        let clock = profile.map(|_| Instant::now());
                        let (found, method) = if modulus <= caps.residue {
                            (residue_count(design, level, modulus), 0)
                        } else {
                            ctx.measure(modulus, span, table, caps)
                        };
                        if let Some(clock) = clock {
                            let band = (64 - (span / modulus).leading_zeros()) as usize;
                            let cell = &mut local[band * METHODS.len() + method];
                            cell.moduli += 1;
                            cell.nanos += clock.elapsed().as_nanos() as u64;
                        }
                        acc += sign as i128 * (found as i128 - origin);
                    }
                }
                if let Some(profile) = profile {
                    let mut cells = profile.cells.lock().unwrap();
                    for (slot, cell) in cells.iter_mut().zip(local.iter()) {
                        slot.moduli += cell.moduli;
                        slot.nanos += cell.nanos;
                    }
                }
                acc
            }));
        }
        handles.into_iter().map(|h| h.join().unwrap()).sum()
    });
    total
}

pub fn terms(design: &Design, top: u32, threads: usize) -> Vec<Level> {
    terms_with(design, top, threads, Mode::Auto)
}

pub fn terms_with(design: &Design, top: u32, threads: usize, mode: Mode) -> Vec<Level> {
    terms_each(design, top, threads, mode, &mut |_| {})
}

pub fn terms_each(
    design: &Design,
    top: u32,
    threads: usize,
    mode: Mode,
    sink: &mut dyn FnMut(&Level),
) -> Vec<Level> {
    let sieve = Sieve::new(design, top);
    let mut previous = 0i128;
    let mut out = Vec::new();
    for level in 1..=top {
        let clock = Instant::now();
        let current = weight(design, level, &sieve, threads, mode, None);
        let value = if design.zero_filled() {
            current - previous
        } else {
            current
        };
        previous = current;
        let entry = Level {
            level,
            value,
            seconds: clock.elapsed().as_secs_f64(),
        };
        sink(&entry);
        out.push(entry);
    }
    out
}

pub fn profile(design: &Design, level: u32, threads: usize) -> i128 {
    let sieve = Sieve::new(design, level);
    let profile = Profile::default();
    let clock = Instant::now();
    let value = weight(design, level, &sieve, threads, Mode::Auto, Some(&profile));
    let seconds = clock.elapsed().as_secs_f64();
    profile.print(level);
    println!("W({}) {} wall {:.3}", level, value, seconds);
    value
}

pub fn count_one(design: &Design, level: u32, modulus: u64, mode: Mode) -> u128 {
    let span = 3u64.pow(level);
    if mode == Mode::Residue {
        return residue_count(design, level, modulus);
    }
    let table = chunk_table(design.invert);
    let mut ctx = Ctx::new(level, design.dimension);
    let caps = caps_for(level, design.dimension, mode);
    let caps = Caps { residue: 0, ..caps };
    ctx.measure(modulus, span, &table, &caps).0
}

pub fn methods(design: &Design, level: u32, modulus: u64) -> Vec<u128> {
    let span = 3u64.pow(level);
    let table = chunk_table(design.invert);
    let mut ctx = Ctx::new(level, design.dimension);
    ctx.gather(modulus, span, &table);
    ctx.dedup();
    let mut out = Vec::new();
    if design.dimension == 2 {
        out.push(ctx.direct2());
        out.push(ctx.zeta2());
    } else {
        out.push(ctx.direct3());
        out.push(ctx.zeta3(true));
        out.push(ctx.zeta3(false));
        out.push(ctx.convolve3(true));
        out.push(ctx.convolve3(false));
        out.push(ctx.bitset3());
        out.push(ctx.rows3());
        out.push(ctx.cube3());
        if ctx.masks.len() <= 3 {
            out.push(ctx.tail3());
        }
    }
    if modulus.pow(design.dimension as u32) <= 4_000_000 {
        out.push(residue_count(design, level, modulus));
    }
    out
}
