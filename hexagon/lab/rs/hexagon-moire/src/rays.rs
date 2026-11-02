use crate::sums::{mean, odds};
use mrly_num::series::chi4;

fn triangle(value: f64) -> f64 {
    let rest = value.rem_euclid(2.0);
    1.0 - 2.0 * rest.min(2.0 - rest)
}

fn average(point: f64, cap: usize, twisted: bool) -> f64 {
    let terms: Vec<f64> = odds(cap)
        .map(|n| {
            let sign = if twisted { f64::from(chi4(n)) } else { 1.0 };
            sign * triangle(n as f64 * point)
        })
        .collect();
    mean(&terms)
}

pub fn twisted() {
    println!("(1/M) sum over odd n <= N of chi4(n) T(n x), the untwisted average beside it");
    let probes = [
        0.0,
        1.0 / 3.0,
        2.0 / 3.0,
        0.2,
        1.0 / 7.0,
        0.5,
        0.25,
        1.0 / 9.0,
        2f64.sqrt() - 1.0,
    ];
    for point in probes {
        println!(
            "  x = {point:.6}: N=55 {:+.5}  N=2001 {:+.6}  N=40001 {:+.7}  untwisted N=40001 {:+.6}",
            average(point, 55, true),
            average(point, 2001, true),
            average(point, 40001, true),
            average(point, 40001, false)
        );
    }
}

pub const RESOLUTION: usize = 200003;

fn square(value: f64) -> f64 {
    if (value.floor() as i64) % 2 == 0 {
        1.0
    } else {
        -1.0
    }
}

fn field(cap: usize, first: &[f64], second: &[f64]) -> f64 {
    let mut running = vec![0.0f64; first.len()];
    let mut layers = 0;
    for n in odds(cap) {
        let twist = -f64::from(chi4(n));
        let scale = n as f64;
        for (slot, cell) in running.iter_mut().enumerate() {
            let left = square(scale * first[slot]);
            let right = square(scale * second[slot]);
            *cell += 0.5 + twist / 8.0 + (left + right) / 4.0 - twist * left * right / 8.0;
        }
        layers += 1;
    }
    mean(&running) / layers as f64
}

pub fn crosshairs() {
    let points: Vec<f64> = (0..RESOLUTION)
        .map(|slot| (slot as f64 + 0.5) / RESOLUTION as f64)
        .collect();
    let generic: Vec<f64> = points
        .iter()
        .map(|point| (point * 2f64.sqrt()).rem_euclid(1.0))
        .collect();
    let shifted: Vec<f64> = points
        .iter()
        .map(|point| (point + 1.0 / 3.0).rem_euclid(1.0))
        .collect();
    println!("coarse cut model, {RESOLUTION} samples along a line, excess over the background");
    for cap in [55usize, 555, 5555] {
        let background = field(cap, &points, &generic);
        let diagonal = field(cap, &points, &points) - background;
        let offset = field(cap, &points, &shifted) - background;
        println!("  N = {cap:4}: background {background:.6}  A=C line {diagonal:+.6}  A-C=1/3 line {offset:+.6}");
    }
    let control: Vec<f64> = points
        .iter()
        .map(|point| (point * 3f64.sqrt()).rem_euclid(1.0))
        .collect();
    println!("coarse crosshair at A = a/q, predicted (-1)^a/(4q) for odd q and 0 for even q");
    for cap in [55usize, 5555] {
        let background = field(cap, &control, &points);
        for (a, q) in [(1, 3), (2, 3), (1, 5), (2, 5), (1, 7), (1, 2), (1, 4)] {
            let line = vec![a as f64 / q as f64 + 1e-12; RESOLUTION];
            let excess = field(cap, &line, &points) - background;
            let predicted = if q % 2 == 1 {
                (if a % 2 == 0 { 1.0 } else { -1.0 }) / (4 * q) as f64
            } else {
                0.0
            };
            println!("  N = {cap:4}  A = {a}/{q}: excess {excess:+.6}  predicted {predicted:+.6}");
        }
    }
}
