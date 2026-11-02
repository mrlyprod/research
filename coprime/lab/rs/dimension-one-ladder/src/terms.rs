use crate::design::Design;
use mrly_num::factor::{gcd, mobius_sieve};
use mrly_num::prime::primes;
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

const CHUNK: u32 = 9;
const CHUNK_SPAN: u64 = 19683;
const FINE: u64 = 1 << 16;
const BLOCK: u64 = 1 << 13;

fn digit_table(invert: bool) -> Vec<u32> {
    (0..CHUNK_SPAN)
        .map(|value| {
            let mut mask = 0u32;
            let mut rest = value;
            for place in 0..CHUNK {
                if (rest % 3 == 1) != invert {
                    mask |= 1 << place;
                }
                rest /= 3;
            }
            mask
        })
        .collect()
}

fn mask_of(mut x: u64, table: &[u32], chunks: u32, full: u32) -> u32 {
    let mut mask = 0u32;
    for chunk in 0..chunks {
        mask |= table[(x % CHUNK_SPAN) as usize] << (chunk * CHUNK);
        x /= CHUNK_SPAN;
    }
    mask & full
}

fn mobius_block(lo: u64, hi: u64, base: &[usize], sign: &mut Vec<i8>, rest: &mut Vec<u64>) {
    let width = (hi - lo) as usize;
    sign.clear();
    sign.resize(width, 1);
    rest.clear();
    rest.extend(lo..hi);
    for &p in base {
        let p = p as u64;
        if p * p >= hi {
            break;
        }
        let mut m = lo.div_ceil(p) * p;
        while m < hi {
            let slot = (m - lo) as usize;
            sign[slot] = -sign[slot];
            rest[slot] /= p;
            m += p;
        }
        let mut m = lo.div_ceil(p * p) * p * p;
        while m < hi {
            sign[(m - lo) as usize] = 0;
            m += p * p;
        }
    }
    for slot in 0..width {
        if sign[slot] != 0 && rest[slot] > 1 {
            sign[slot] = -sign[slot];
        }
    }
}

fn subset_sums(buf: &mut [u64], level: u32) {
    for place in 0..level {
        let bit = 1usize << place;
        for base in (0..buf.len()).step_by(bit << 1) {
            let (low, high) = buf[base..base + (bit << 1)].split_at_mut(bit);
            for index in 0..bit {
                high[index] += low[index];
            }
        }
    }
}

struct Lattice {
    level: u32,
    size: usize,
    full: u32,
    chunks: u32,
    dimension: usize,
    count: Vec<u32>,
    items: Vec<(u32, u32)>,
    sub: Vec<(u32, u32)>,
    sums: Vec<u64>,
    ranked: Vec<u64>,
    choose: Vec<Vec<u64>>,
}

impl Lattice {
    fn new(level: u32, dimension: usize) -> Lattice {
        let size = 1usize << level;
        let ranks = level as usize + 1;
        let mut choose = vec![vec![0u64; ranks]; ranks];
        for row in 0..ranks {
            choose[row][0] = 1;
            for column in 1..=row {
                choose[row][column] = choose[row - 1][column - 1] + choose[row - 1][column];
            }
        }
        Lattice {
            level,
            size,
            full: (1u32 << level) - 1,
            chunks: level.div_ceil(CHUNK).max(1),
            dimension,
            count: vec![0; size],
            items: Vec::new(),
            sub: Vec::new(),
            sums: vec![0; size],
            ranked: if dimension == 3 {
                vec![0; ranks * size]
            } else {
                Vec::new()
            },
            choose,
        }
    }

    fn gather(&mut self, modulus: u64, span: u64, table: &[u32]) {
        self.items.clear();
        for x in (0..span).step_by(modulus as usize) {
            let mask = mask_of(x, table, self.chunks, self.full) as usize;
            if self.count[mask] == 0 {
                self.items.push((mask as u32, 0));
            }
            self.count[mask] += 1;
        }
        for item in self.items.iter_mut() {
            item.1 = self.count[item.0 as usize];
        }
    }

    fn release(&mut self) {
        for &(mask, _) in self.items.iter() {
            self.count[mask as usize] = 0;
        }
    }

    fn spread(&mut self) {
        for &(mask, weight) in self.items.iter() {
            self.sums[mask as usize] = weight as u64;
        }
        subset_sums(&mut self.sums, self.level);
    }

    fn direct2(&self) -> u128 {
        let mut total = 0u128;
        for &(a, ca) in self.items.iter() {
            let inner: u64 = self
                .items
                .iter()
                .filter(|(b, _)| a & b == 0)
                .map(|(_, cb)| *cb as u64)
                .sum();
            total += ca as u128 * inner as u128;
        }
        total
    }

    fn zeta2(&mut self) -> u128 {
        self.spread();
        let mut total = 0u128;
        for &(a, ca) in self.items.iter() {
            total += ca as u128 * self.sums[(!a & self.full) as usize] as u128;
        }
        self.sums.fill(0);
        total
    }

    fn direct3(&mut self) -> u128 {
        let mut total = 0u128;
        for &(a, ca) in self.items.iter() {
            self.sub.clear();
            self.sub
                .extend(self.items.iter().filter(|(b, _)| a & b == 0));
            let mut inner = 0u128;
            for &(b, cb) in self.sub.iter() {
                let reach: u64 = self
                    .sub
                    .iter()
                    .filter(|(c, _)| b & c == 0)
                    .map(|(_, cc)| *cc as u64)
                    .sum();
                inner += cb as u128 * reach as u128;
            }
            total += ca as u128 * inner;
        }
        total
    }

    fn zeta3(&mut self) -> u128 {
        self.spread();
        let width = self.items.len() as u64;
        let mut total = 0u128;
        for &(a, ca) in self.items.iter() {
            let free = !a & self.full;
            let mut inner = 0u64;
            if 1u64 << free.count_ones() <= width {
                let mut b = free;
                loop {
                    inner += self.count[b as usize] as u64 * self.sums[(free ^ b) as usize];
                    if b == 0 {
                        break;
                    }
                    b = (b - 1) & free;
                }
            } else {
                for &(b, cb) in self.items.iter() {
                    if a & b == 0 {
                        inner += cb as u64 * self.sums[(free & !b) as usize];
                    }
                }
            }
            total += ca as u128 * inner as u128;
        }
        self.sums.fill(0);
        total
    }

    fn convolve3(&mut self) -> u128 {
        let top = self.level as usize;
        let size = self.size;
        self.ranked.fill(0);
        for &(mask, weight) in self.items.iter() {
            let rank = mask.count_ones() as usize;
            self.ranked[rank * size + mask as usize] = weight as u64;
        }
        for rank in 0..=top {
            subset_sums(&mut self.ranked[rank * size..(rank + 1) * size], self.level);
        }
        let mut poly = vec![0u128; top + 1];
        let mut square = vec![0u128; top + 1];
        let mut total = 0i128;
        for set in 0..size {
            let rank = (set as u32).count_ones() as usize;
            for (degree, slot) in poly.iter_mut().enumerate().take(rank + 1) {
                *slot = self.ranked[degree * size + set] as u128;
            }
            let reach = (2 * rank).min(top);
            square[..=reach].fill(0);
            for left in 0..=rank.min(reach) {
                for right in 0..=rank.min(reach - left) {
                    square[left + right] += poly[left] * poly[right];
                }
            }
            let mut acc = 0i128;
            for degree in rank..=(3 * rank).min(top) {
                let mut cube = 0u128;
                for left in degree.saturating_sub(rank)..=reach.min(degree) {
                    cube += square[left] * poly[degree - left];
                }
                let weight = self.choose[top - rank][degree - rank] as i128 * cube as i128;
                if (degree - rank).is_multiple_of(2) {
                    acc += weight;
                } else {
                    acc -= weight;
                }
            }
            total += acc;
        }
        assert!(total >= 0);
        total as u128
    }

    fn count(&mut self, modulus: u64, span: u64, table: &[u32]) -> u128 {
        self.gather(modulus, span, table);
        let width = self.items.len();
        let sweep = self.level as usize * self.size;
        let found = if self.dimension == 2 {
            if width * width <= sweep {
                self.direct2()
            } else {
                self.zeta2()
            }
        } else {
            let direct = width * width * width / 40;
            let pair = sweep + width * width;
            let ranked = 2 * self.level as usize * sweep;
            if direct <= pair.min(ranked) {
                self.direct3()
            } else if pair <= ranked {
                self.zeta3()
            } else {
                self.convolve3()
            }
        };
        self.release();
        found
    }
}

fn weight(
    design: &Design,
    level: u32,
    threads: usize,
    table: &[u32],
    small: &[i8],
    base: &[usize],
) -> i128 {
    let span = 3u64.pow(level);
    let fine = span.min(FINE);
    let origin: i128 = if design.origin_filled() { 1 } else { 0 };
    let cursor = AtomicU64::new(0);
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                let cursor = &cursor;
                scope.spawn(move || {
                    let mut lattice = Lattice::new(level, design.dimension);
                    let mut sign = Vec::new();
                    let mut rest = Vec::new();
                    let mut acc = 0i128;
                    loop {
                        let task = cursor.fetch_add(1, Ordering::Relaxed);
                        let (lo, hi) = if task < fine {
                            (task + 1, task + 2)
                        } else {
                            let lo = fine + 1 + (task - fine) * BLOCK;
                            (lo, lo + BLOCK)
                        };
                        if lo >= span {
                            break;
                        }
                        let hi = hi.min(span);
                        if task >= fine {
                            mobius_block(lo, hi, base, &mut sign, &mut rest);
                        }
                        for modulus in lo..hi {
                            let mu = if task < fine {
                                small[modulus as usize]
                            } else {
                                sign[(modulus - lo) as usize]
                            };
                            if mu == 0 || modulus % 3 == 0 {
                                continue;
                            }
                            let found = lattice.count(modulus, span, table) as i128;
                            acc += mu as i128 * (found - origin);
                        }
                    }
                    acc
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).sum()
    })
}

pub struct Term {
    pub level: u32,
    pub value: i128,
    pub seconds: f64,
}

pub fn terms(design: &Design, top: u32, threads: usize) -> Vec<Term> {
    let table = digit_table(design.invert);
    let small = mobius_sieve(FINE as usize);
    let base = primes((3f64.powi(top as i32)).sqrt() as usize + 2);
    let mut previous = 0i128;
    let mut out = Vec::new();
    for level in 1..=top {
        let clock = Instant::now();
        let current = weight(design, level, threads, &table, &small, &base);
        let value = if design.origin_filled() {
            current - previous
        } else {
            current
        };
        previous = current;
        out.push(Term {
            level,
            value,
            seconds: clock.elapsed().as_secs_f64(),
        });
    }
    out
}

fn walk(depth: u32, coords: &mut [u64], corners: &[Vec<u64>], found: &mut u64) {
    if depth == 0 {
        let common = coords.iter().fold(0u128, |g, &c| gcd(g, u128::from(c)));
        if common == 1 {
            *found += 1;
        }
        return;
    }
    for corner in corners {
        for (slot, digit) in coords.iter_mut().zip(corner) {
            *slot = *slot * 3 + digit;
        }
        walk(depth - 1, coords, corners, found);
        for (slot, digit) in coords.iter_mut().zip(corner) {
            *slot = (*slot - digit) / 3;
        }
    }
}

pub fn brute(design: &Design, level: u32) -> u64 {
    let corners = design.corners();
    let mut coords = vec![0u64; design.dimension];
    let mut found = 0;
    walk(level, &mut coords, &corners, &mut found);
    found
}

pub fn stored(design: &Design) -> Vec<(u32, i128)> {
    let path = format!(
        "{}/terms/{}_bfile.txt",
        env!("CARGO_MANIFEST_DIR"),
        design.name
    );
    std::fs::read_to_string(path)
        .unwrap_or_default()
        .lines()
        .filter_map(|line| {
            let mut parts = line.split_whitespace();
            Some((parts.next()?.parse().ok()?, parts.next()?.parse().ok()?))
        })
        .collect()
}
