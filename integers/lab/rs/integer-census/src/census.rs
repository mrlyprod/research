use crate::rows::{Row, Sheet, CEILING};

pub const WINDOWS: [usize; 3] = [1_000, 10_000, 100_000];
pub const CHAMPIONS: usize = 20;
pub const MISSES: usize = 30;

pub struct Census {
    pub counts: Vec<u32>,
    pub incidences: u64,
    pub repeats: u64,
    pub low: u64,
}

pub struct Band {
    pub first: usize,
    pub last: usize,
    pub missed: usize,
}

impl Band {
    pub fn width(&self) -> usize {
        self.last - self.first + 1
    }
    pub fn density(&self) -> f64 {
        self.missed as f64 / self.width() as f64
    }
}

pub fn build(sheet: &Sheet) -> Census {
    let mut counts = vec![0u32; CEILING as usize + 1];
    let mut incidences = 0u64;
    let mut repeats = 0u64;
    let mut low = 0u64;
    for row in &sheet.rows {
        for &term in &row.written {
            counts[term as usize] += 1;
            incidences += 1;
        }
        repeats += row.repeats as u64;
        low += row.low as u64;
    }
    Census {
        counts,
        incidences,
        repeats,
        low,
    }
}

pub fn split(census: &Census, window: usize) -> (usize, usize, usize) {
    let mut never = 0;
    let mut once = 0;
    let mut many = 0;
    for &count in &census.counts[1..=window] {
        match count {
            0 => never += 1,
            1 => once += 1,
            _ => many += 1,
        }
    }
    (never, once, many)
}

pub fn champions(census: &Census, take: usize) -> Vec<(usize, u32)> {
    let mut all: Vec<(usize, u32)> = census
        .counts
        .iter()
        .enumerate()
        .skip(1)
        .map(|(value, &count)| (value, count))
        .collect();
    all.sort_by(|left, right| right.1.cmp(&left.1).then(left.0.cmp(&right.0)));
    all.truncate(take);
    all
}

pub fn misses(census: &Census, take: usize) -> Vec<usize> {
    census
        .counts
        .iter()
        .enumerate()
        .skip(1)
        .filter(|(_, &count)| count == 0)
        .map(|(value, _)| value)
        .take(take)
        .collect()
}

pub fn bands(census: &Census) -> Vec<Band> {
    let mut out = Vec::new();
    let mut first = 1usize;
    while first <= CEILING as usize {
        let last = (first * 10 - 1).min(CEILING as usize);
        let missed = census.counts[first..=last]
            .iter()
            .filter(|&&count| count == 0)
            .count();
        out.push(Band {
            first,
            last,
            missed,
        });
        first *= 10;
    }
    out
}

pub fn writers(sheet: &Sheet, value: i128) -> Vec<&Row> {
    sheet
        .rows
        .iter()
        .filter(|row| row.written.binary_search(&value).is_ok())
        .collect()
}
