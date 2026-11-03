use mrly_math::bang::baseq::{fill_from_corners, group, orbit, representatives};
use mrly_math::bang::code_to_corners;
use mrly_math::bang::factory::residue_corners;
use mrly_math::bang::universe::permutations;
use mrly_math::bang::Code;
use std::collections::HashSet;
use std::thread;

const SIDE: usize = 3;
const BRUTE: usize = 4;
const BURNSIDE: usize = 6;
const PUBLISHED: [u128; 6] = [4, 12, 80, 3984, 37333248, 25626412338274304];
const SIDES: [usize; 6] = [2, 3, 4, 5, 6, 7];

fn corner_index(cell: &[u8]) -> usize {
    cell.iter()
        .fold(0, |acc, &digit| acc * 2 + (digit & 1) as usize)
}

fn cell_index(cell: &[u8]) -> usize {
    cell.iter()
        .fold(0, |acc, &digit| acc * SIDE + digit as usize)
}

fn corner_of_cell(dimension: usize) -> Vec<usize> {
    residue_corners(dimension, SIDE)
        .iter()
        .map(|cell| corner_index(cell))
        .collect()
}

fn pattern(code: u128, corners: &[usize]) -> u128 {
    corners
        .iter()
        .enumerate()
        .filter(|(_, &corner)| code >> corner & 1 == 1)
        .fold(0, |acc, (cell, _)| acc | 1u128 << cell)
}

fn motions(dimension: usize) -> Vec<Vec<u8>> {
    let cells = residue_corners(dimension, SIDE);
    let last = (SIDE - 1) as u8;
    let mut out = Vec::new();
    for perm in permutations(dimension) {
        for flips in 0..1usize << dimension {
            out.push(
                cells
                    .iter()
                    .map(|cell| {
                        let moved: Vec<u8> = (0..dimension)
                            .map(|i| {
                                let digit = cell[perm[i]];
                                if flips >> i & 1 == 1 {
                                    last - digit
                                } else {
                                    digit
                                }
                            })
                            .collect();
                        cell_index(&moved) as u8
                    })
                    .collect(),
            );
        }
    }
    out
}

fn image(pattern: u128, motion: &[u8]) -> u128 {
    let mut rest = pattern;
    let mut out = 0u128;
    while rest != 0 {
        out |= 1u128 << motion[rest.trailing_zeros() as usize];
        rest &= rest - 1;
    }
    out
}

fn shapes(dimension: usize) -> (usize, usize) {
    let corners = corner_of_cell(dimension);
    let motions = motions(dimension);
    let total = 1u128 << (1usize << dimension);
    let chunk = total.div_ceil(8);
    let found: Vec<(HashSet<u128>, HashSet<u128>)> = thread::scope(|scope| {
        let handles: Vec<_> = (0..8u128)
            .map(|part| {
                let corners = &corners;
                let motions = &motions;
                scope.spawn(move || {
                    let mut raw = HashSet::new();
                    let mut least = HashSet::new();
                    for code in part * chunk..((part + 1) * chunk).min(total) {
                        let drawn = pattern(code, corners);
                        raw.insert(drawn);
                        least.insert(
                            motions
                                .iter()
                                .map(|motion| image(drawn, motion))
                                .min()
                                .expect("the motion group is not empty"),
                        );
                    }
                    (raw, least)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("a part ran"))
            .collect()
    });
    let raw: HashSet<u128> = found
        .iter()
        .flat_map(|(raw, _)| raw.iter().copied())
        .collect();
    let least: HashSet<u128> = found
        .iter()
        .flat_map(|(_, least)| least.iter().copied())
        .collect();
    (raw.len(), least.len())
}

fn cycles(perm: &[usize], dimension: usize) -> u32 {
    let image =
        |corner: usize| (0..dimension).fold(0usize, |acc, i| acc | (corner >> perm[i] & 1) << i);
    let mut seen = vec![false; 1 << dimension];
    let mut count = 0;
    for start in 0..1usize << dimension {
        if !seen[start] {
            count += 1;
            let mut at = start;
            while !seen[at] {
                seen[at] = true;
                at = image(at);
            }
        }
    }
    count
}

fn burnside(dimension: usize) -> u128 {
    let perms = permutations(dimension);
    let total: u128 = perms
        .iter()
        .map(|perm| 1u128 << cycles(perm, dimension))
        .sum();
    total / perms.len() as u128
}

fn fill(code: Code, side: usize, dimension: usize) -> u128 {
    let corners = code_to_corners(code, dimension, 2).expect("the code fits its corners");
    fill_from_corners(&corners, side, dimension)
}

fn agree(dimension: usize) -> usize {
    let corners = corner_of_cell(dimension);
    (0..1u128 << (1usize << dimension))
        .filter(|&code| {
            pattern(code, &corners).count_ones() as u128 != fill(Code::from(code), SIDE, dimension)
        })
        .count()
}

fn split(dimension: usize, side: usize) -> (usize, usize) {
    let cube = group(2, dimension);
    let classes = representatives(2, dimension).expect("the walk fits");
    let mixed = classes
        .iter()
        .filter(|(rep, _)| {
            orbit(&cube, *rep)
                .into_iter()
                .map(|member| fill(member, side, dimension))
                .collect::<HashSet<u128>>()
                .len()
                > 1
        })
        .count();
    (classes.len(), mixed)
}

fn report() {
    println!("SIDE {SIDE}: designs, distinct drawn patterns, classes up to rigid motion, A003180");
    for dimension in 1..=BRUTE {
        let (raw, least) = shapes(dimension);
        println!(
            "dim {dimension}: {} {raw} {least} {} fill mismatches {}",
            1u128 << (1usize << dimension),
            PUBLISHED[dimension - 1],
            agree(dimension)
        );
    }
    println!("BURNSIDE over axis permutations");
    for dimension in 1..=BURNSIDE {
        println!("dim {dimension}: {}", burnside(dimension));
    }
    println!("CUBE CLASSES with a non-constant fill: dim classes mixed");
    for side in SIDES {
        for dimension in 1..=BRUTE {
            let (classes, mixed) = split(dimension, side);
            println!("side {side} dim {dimension}: {classes} {mixed}");
        }
    }
    println!("CLASS of dim 2 code 7: code cells dimension");
    let cube = group(2, 2);
    for member in orbit(&cube, Code::from(7u128)) {
        let cells = fill(member, SIDE, 2);
        println!(
            "{member} {cells} {:.6}",
            (cells as f64).ln() / (SIDE as f64).ln()
        );
    }
}

fn main() {
    report();
}
