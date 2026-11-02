use mrly_num::factor::coprime;
use std::thread;

pub const THREADS: u64 = 8;

fn row_slice(n: u32, rem: u64, step: u64) -> u64 {
    let mut count = 0u64;
    let mut m = 2 * rem + 1;
    while m < 1u64 << n {
        let rest = m - 1;
        let mut sub = rest;
        loop {
            if coprime(sub as usize + 1, m as usize) {
                count += 1;
            }
            if sub == 0 {
                break;
            }
            sub = (sub - 1) & rest;
        }
        m += 2 * step;
    }
    count
}

fn complement_slice(n: u32, rem: u64, step: u64) -> u64 {
    let size = 1u64 << n;
    let mut count = 0u64;
    let mut i = rem;
    while i < size {
        let free = !i & (size - 1);
        let mut j = free;
        while j > i {
            if coprime(i as usize, j as usize) {
                count += 1;
            }
            j = (j - 1) & free;
        }
        i += step;
    }
    count
}

fn spread(n: u32, slice: fn(u32, u64, u64) -> u64) -> u64 {
    let handles: Vec<thread::JoinHandle<u64>> = (0..THREADS)
        .map(|rem| thread::spawn(move || slice(n, rem, THREADS)))
        .collect();
    2 * handles.into_iter().map(|h| h.join().unwrap()).sum::<u64>()
}

pub fn by_rows(n: u32) -> u64 {
    spread(n, row_slice)
}

pub fn by_complement(n: u32) -> u64 {
    spread(n, complement_slice)
}

pub fn by_pascal(n: u32) -> u64 {
    let mut row: Vec<u8> = vec![1];
    let mut total = 0u64;
    for m in 0..1u64 << n {
        for (k, entry) in row.iter().enumerate() {
            if *entry == 1 && coprime(k, m as usize - k) {
                total += 1;
            }
        }
        let mut next: Vec<u8> = Vec::with_capacity(row.len() + 1);
        next.push(1);
        for k in 1..row.len() {
            next.push((row[k - 1] + row[k]) & 1);
        }
        next.push(1);
        row = next;
    }
    total
}

pub fn density(term: u64, n: u32) -> f64 {
    term as f64 / 3f64.powi(n as i32)
}

pub fn limit() -> f64 {
    16.0 / (3.0 * std::f64::consts::PI * std::f64::consts::PI)
}
