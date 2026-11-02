use crate::fills::cell_index;
use crate::tables::write_csv;
use mrly_math::bang::baseq::{distinct_designs, fill_from_corners};
use mrly_math::bang::factory::code_to_corners;
use mrly_math::bang::Code;
use mrly_math::rules::render;
use num_bigint::BigUint;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use std::thread;

const SIDES: [usize; 3] = [3, 5, 7];
const LEVELS: [usize; 3] = [1, 2, 3];
const SEQUENCE: usize = 9;
const TABLE: usize = 8;
const BURNSIDE_LIMIT: usize = 6;
const PUBLISHED: [&str; 16] = [
    "2",
    "4",
    "12",
    "64",
    "700",
    "17424",
    "1053696",
    "160579584",
    "62856336636",
    "63812936890000",
    "168895157342195152",
    "1169048914836855865344",
    "21209591746609937928524800",
    "1010490883477487017627972550656",
    "126641164340871500483202065902080000",
    "41817338589698457759723104703370865147904",
];

fn brute(design: Code, number: usize, dimension: usize, level: usize) -> u128 {
    let filled: Vec<Vec<u8>> = mrly_math::bang::factory::residue_corners(dimension, 2)
        .into_iter()
        .filter(|r| design.get() >> cell_index(r, 2) & 1 == 1)
        .collect();
    let tile = render(&filled, number, dimension, 2).expect("the tile renders");
    u128::from(tile.fractal(level).sum())
}

fn closed(design: Code, number: usize, dimension: usize, level: usize) -> u128 {
    let corners = code_to_corners(design, dimension, 2).expect("the design fits its corners");
    fill_from_corners(&corners, number, dimension).pow(level as u32)
}

fn profile(design: Code, dimension: usize) -> Vec<usize> {
    let mut out = vec![0usize; dimension + 1];
    for corner in 0..1usize << dimension {
        if design.get() >> corner & 1 == 1 {
            out[corner.count_ones() as usize] += 1;
        }
    }
    out
}

fn binomial(n: usize, k: usize) -> u128 {
    (0..k).fold(1u128, |acc, i| acc * (n - i) as u128 / (i as u128 + 1))
}

pub fn a129824(dimension: usize) -> BigUint {
    (0..=dimension)
        .map(|k| BigUint::from(1 + binomial(dimension, k)))
        .product()
}

fn weight_at_most_one(dimension: usize) -> Code {
    Code::from(
        (0..1usize << dimension)
            .filter(|corner| corner.count_ones() <= 1)
            .map(|corner| 1u128 << corner)
            .sum::<u128>(),
    )
}

fn compare(dimension: usize) -> (usize, usize) {
    let designs = 1u128 << (1usize << dimension);
    let workers = thread::available_parallelism().map_or(1, |n| n.get());
    let chunk = (designs as usize).div_ceil(workers) as u128;
    let totals: Vec<(usize, usize)> = thread::scope(|scope| {
        let handles: Vec<_> = (0..workers as u128)
            .map(|w| {
                scope.spawn(move || {
                    let mut checks = 0usize;
                    let mut bad = 0usize;
                    for design in w * chunk..((w + 1) * chunk).min(designs) {
                        for &n in &SIDES {
                            for &level in &LEVELS {
                                checks += 1;
                                if brute(Code::from(design), n, dimension, level)
                                    != closed(Code::from(design), n, dimension, level)
                                {
                                    bad += 1;
                                }
                            }
                        }
                    }
                    (checks, bad)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|handle| handle.join().expect("the worker finishes"))
            .collect()
    });
    totals.iter().fold((0, 0), |(c, b), &(x, y)| (c + x, b + y))
}

fn distinct(dimension: usize) -> (usize, usize) {
    let mut seen: BTreeMap<Vec<u128>, BTreeSet<Vec<usize>>> = BTreeMap::new();
    for design in 0..1u128 << (1usize << dimension) {
        let sequence: Vec<u128> = (1..=SEQUENCE)
            .map(|n| closed(Code::from(design), n, dimension, 1))
            .collect();
        seen.entry(sequence)
            .or_default()
            .insert(profile(Code::from(design), dimension));
    }
    let collisions = seen.values().filter(|profiles| profiles.len() > 1).count();
    (seen.len(), collisions)
}

pub fn report(path: &Path) {
    println!("fill classes at base 2: two generators, the closed form and the published terms");
    let sponge = weight_at_most_one(3);
    let carpet = weight_at_most_one(2);
    println!(
        "Menger sponge n 3: level 1 {} level 2 {}; carpet n 3: {} n 5: {}",
        brute(sponge, 3, 3, 1),
        brute(sponge, 3, 3, 2),
        brute(carpet, 3, 2, 1),
        brute(carpet, 5, 2, 1)
    );
    for dimension in [2usize, 3] {
        let (checks, bad) = compare(dimension);
        println!("D {dimension}: {checks} generator checks, {bad} mismatches");
    }
    for dimension in 1..=4 {
        let (count, collisions) = distinct(dimension);
        let closed = a129824(dimension);
        println!(
            "D {dimension}: distinct fill sequences {count}, closed form {closed}, match {}, profile collisions {collisions}",
            count.to_string() == closed.to_string()
        );
    }
    let agreed = PUBLISHED
        .iter()
        .enumerate()
        .filter(|(d, term)| a129824(*d).to_string() == **term)
        .count();
    let terms: Vec<String> = (0..PUBLISHED.len())
        .map(|d| a129824(d).to_string())
        .collect();
    println!(
        "A129824 closed form D 0..{}: {}",
        PUBLISHED.len() - 1,
        terms.join(", ")
    );
    println!("published terms matched: {agreed} of {}", PUBLISHED.len());
    println!(
        "{:<3}{:>12}{:>18}{:>16}{:>8}",
        "D", "designs", "A000616", "A129824", "ratios"
    );
    let header: Vec<String> = [
        "D",
        "total_designs",
        "cube_classes_A000616",
        "fractal_dim_classes_A129824",
        "limiting_ratio_classes",
    ]
    .iter()
    .map(|name| (*name).to_string())
    .collect();
    let mut records = Vec::new();
    for dimension in 1..=TABLE {
        let designs = BigUint::from(2u32).pow(1 << dimension);
        let cubes = (dimension <= BURNSIDE_LIMIT)
            .then(|| distinct_designs(2, dimension).expect("the Burnside average is an integer"))
            .map(|value| value.to_string())
            .unwrap_or_default();
        let classes = a129824(dimension);
        let ratios = (1u32 << dimension) + 1;
        let shown = if designs.to_string().len() > 12 {
            format!("2^{}", 1 << dimension)
        } else {
            designs.to_string()
        };
        println!("{dimension:<3}{shown:>12}{cubes:>18}{classes:>16}{ratios:>8}");
        records.push(vec![
            dimension.to_string(),
            designs.to_string(),
            cubes,
            classes.to_string(),
            ratios.to_string(),
        ]);
    }
    println!("written counts.csv, A000616 by Burnside to D = {BURNSIDE_LIMIT}");
    write_csv(path, &header, &records);
}
