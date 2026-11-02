mod coprime;
mod design;

use num_bigint::BigUint;
use std::time::Instant;

const COMPLEMENT_N: u32 = 14;
const PASCAL_N: u32 = 12;
const DESIGN_N: usize = 14;
const AFFINE_N: usize = 8;
const ORBIT_N: usize = 2;

fn verdict(agree: bool) -> &'static str {
    if agree {
        "PASS"
    } else {
        "FAIL"
    }
}

fn a396934(top: u32) {
    let clock = Instant::now();
    let rows: Vec<u64> = (0..=top).map(coprime::by_rows).collect();
    println!(
        "A396934 row walk n = 0..{} in {:.1}s",
        top,
        clock.elapsed().as_secs_f64()
    );
    let (low, tiny) = (COMPLEMENT_N.min(top), PASCAL_N.min(top));
    let complements: Vec<u64> = (0..=low).map(coprime::by_complement).collect();
    let pascal: Vec<u64> = (0..=tiny).map(coprime::by_pascal).collect();
    println!(
        "A396934 row walk vs complement walk n = 0..{low}: {}",
        verdict(rows[..=low as usize] == complements[..])
    );
    println!(
        "A396934 row walk vs Pascal mod 2 n = 0..{tiny}: {}",
        verdict(rows[..=tiny as usize] == pascal[..])
    );
    for (n, value) in rows.iter().enumerate() {
        println!("A396934 {n} {value}");
    }
    let density = coprime::density(rows[top as usize], top);
    println!("A396934 a({top})/3^{top} = {density:.10} = {density:.7}");
    println!(
        "A396934 16/(3*Pi^2) = {:.10} = {:.7}",
        coprime::limit(),
        coprime::limit()
    );
    println!("A396934 gap = {:.1e}", density - coprime::limit());
}

fn a398348() {
    let clock = Instant::now();
    let terms: Vec<BigUint> = (1..=DESIGN_N).map(design::by_cycles).collect();
    println!(
        "A398348 cycle walk n = 1..{} in {:.1}s",
        DESIGN_N,
        clock.elapsed().as_secs_f64()
    );
    let clock = Instant::now();
    let affine: Vec<BigUint> = (1..=AFFINE_N).map(design::by_affine).collect();
    println!(
        "A398348 affine powers n = 1..{} in {:.1}s",
        AFFINE_N,
        clock.elapsed().as_secs_f64()
    );
    println!(
        "A398348 cycle walk vs affine powers n = 1..{}: {}",
        AFFINE_N,
        verdict(terms[..AFFINE_N] == affine[..])
    );
    let orbits: Vec<BigUint> = (1..=ORBIT_N)
        .map(|n| BigUint::from(design::by_orbits(n)))
        .collect();
    println!(
        "A398348 cycle walk vs orbit enumeration n = 1..{}: {}",
        ORBIT_N,
        verdict(terms[..ORBIT_N] == orbits[..])
    );
    for (i, value) in terms.iter().enumerate() {
        println!("A398348 {} {}", i + 1, value);
    }
    for n in [7usize, 8] {
        println!(
            "A398348 digits of a({}) = {}",
            n,
            terms[n - 1].to_string().len()
        );
    }
}

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(20);
    println!(
        "domain A396934 n = 0..{top} on {} threads, complement walk to {COMPLEMENT_N}, Pascal to {PASCAL_N}",
        coprime::THREADS
    );
    println!("domain A398348 n = 1..{DESIGN_N}, affine to {AFFINE_N}, orbits to {ORBIT_N}");
    a396934(top);
    a398348();
}
