use std::sync::atomic::{AtomicUsize, Ordering};

use crate::design::Design;

// GCD

fn gcd(mut a: u64, mut b: u64) -> u64 {
    while b != 0 {
        let rest = a % b;
        a = b;
        b = rest;
    }
    a
}

fn walk(depth: u32, coords: &mut Vec<u64>, corners: &[Vec<u64>], found: &mut u64) {
    if depth == 0 {
        let mut common = coords[0];
        for value in coords.iter().skip(1) {
            common = gcd(common, *value);
        }
        if common == 1 {
            *found += 1;
        }
        return;
    }
    for corner in corners.iter() {
        for (slot, digit) in coords.iter_mut().zip(corner.iter()) {
            *slot = *slot * 3 + digit;
        }
        walk(depth - 1, coords, corners, found);
        for (slot, digit) in coords.iter_mut().zip(corner.iter()) {
            *slot = (*slot - digit) / 3;
        }
    }
}

pub fn count(design: &Design, level: u32, threads: usize) -> u64 {
    let corners = design.corners();
    if level < 2 {
        let mut coords = vec![0u64; design.dimension];
        let mut found = 0u64;
        walk(level, &mut coords, &corners, &mut found);
        return found;
    }
    let fill = corners.len();
    let cursor = AtomicUsize::new(0);
    std::thread::scope(|scope| {
        let mut handles = Vec::new();
        for _ in 0..threads {
            let cursor = &cursor;
            let corners = &corners;
            handles.push(scope.spawn(move || {
                let mut found = 0u64;
                loop {
                    let task = cursor.fetch_add(1, Ordering::Relaxed);
                    if task >= fill * fill {
                        break;
                    }
                    let first = &corners[task / fill];
                    let second = &corners[task % fill];
                    let mut coords: Vec<u64> = first
                        .iter()
                        .zip(second.iter())
                        .map(|(a, b)| a * 3 + b)
                        .collect();
                    walk(level - 2, &mut coords, corners, &mut found);
                }
                found
            }));
        }
        handles.into_iter().map(|h| h.join().unwrap()).sum()
    })
}
