mod carpet;
mod census;
mod count;
mod design;
mod fit;
mod two;

use count::{arcs, glue, mirrors, series, void_block};
use design::Design;
use std::time::Instant;

// CHECK

fn check() {
    let (mut levels, mut images, mut voids, mut routes, mut bad) = (0, 0, 0, 0, 0);
    for (base, top) in [(2usize, 7usize), (3, 4)] {
        for code in 0..1u128 << (base * base) {
            let d = Design::full(code, base);
            let blocks = series(base, &d.tile, top);
            for (level, &want) in blocks.iter().enumerate() {
                let (side, on) = d.cells(level);
                let (loops, strands) = arcs(side, &on);
                let euler = mirrors(side, &on);
                levels += 1;
                if loops as u128 != want || euler != loops as i64 || strands != 2 * side as u64 {
                    bad += 1;
                    println!("mismatch {} level {level}: arcs {loops}/{strands}, mirrors {euler}, blocks {want}", d.name());
                }
            }
            for g in 1..4 {
                images += 1;
                if series(base, &d.image(g).tile, top) != blocks {
                    bad += 1;
                    println!("symmetry {g} breaks {}", d.name());
                }
            }
        }
    }
    for n in 2..6 {
        let mut s = 1;
        while s * n <= 3000 {
            let glued = glue(n, &vec![true; n * n], &void_block(s), true);
            voids += 1;
            if glued.partner != void_block(s * n).partner || glued.loops != 0 {
                bad += 1;
                println!("void block differs at base {n} side {s}");
            }
            s *= n;
        }
    }
    for (code, n) in [
        (7u128, 3usize),
        (14, 3),
        (11, 3),
        (13, 3),
        (6, 3),
        (9, 3),
        (7, 4),
        (7, 5),
    ] {
        let parity = Design::new(code, n, 2);
        let full = Design::full(parity.full_code(), n);
        for level in 0..=3 {
            routes += 1;
            if parity.cells(level) != full.cells(level) {
                bad += 1;
                println!("parity route differs on {}", parity.name());
            }
        }
    }
    println!("CHECK every code at base 2 to level 7 and at base 3 to level 4");
    println!("  union-find over midpoints = mirror cycle rank = block recursion, strands = 2 side: {levels} levels");
    println!("  half turn, transpose, anti-transpose give the same loop sequence: {images} images");
    println!("  glued void block = formula, no loop: {voids} sides at bases 2 to 5");
    println!("  parity name and full code build the same cells: {routes} levels");
    println!("  mismatches: {bad}");
}

fn main() {
    let verb = std::env::args().nth(1).unwrap_or_else(|| "all".to_string());
    let clock = Instant::now();
    if verb == "check" || verb == "all" {
        check();
    }
    if verb == "two" || verb == "all" {
        two::two();
    }
    if verb == "carpet" || verb == "all" {
        carpet::carpet();
    }
    if verb == "bases" || verb == "all" {
        census::bases();
    }
    if verb == "lengths" || verb == "all" {
        census::lengths(3, 12, &[495, 43, 181, 175, 171]);
    }
    if verb == "gains" || verb == "all" {
        census::gains();
    }
    if verb == "census" || verb == "all" {
        census::census(2, 24, 24, &(0..16).collect::<Vec<u128>>());
        census::census(3, 14, 17, &(0..512).collect::<Vec<u128>>());
    }
    println!("time {:.1} s", clock.elapsed().as_secs_f64());
}
