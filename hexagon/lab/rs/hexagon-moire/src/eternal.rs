use num_rational::Ratio;
use std::collections::BTreeSet;

const DENOMINATORS: i64 = 48;
const CUT: i64 = 55;
const REACH: i64 = 2001;

type Q = Ratio<i64>;

fn weight(n: i64, x: Q, z: Q) -> Option<i64> {
    let cx = ((4 * n * x.numer()).div_euclid(*x.denom())).min(4 * n - 1);
    let cz = (2 * ((2 * n * z.numer()).div_euclid(*z.denom()))).min(4 * n - 2);
    let cy = 6 * n - 2 - cx - cz;
    (0..4 * n).contains(&cy).then(|| {
        [cx, cy, cz]
            .iter()
            .map(|c| c.div_euclid(4) & 1)
            .sum::<i64>()
            & 1
    })
}

fn constant(x: Q, z: Q, reach: i64) -> Option<i64> {
    let mut seen = None;
    for n in (1..=reach).step_by(2) {
        let w = weight(n, x, z)?;
        if seen.is_some_and(|s| s != w) {
            return None;
        }
        seen = Some(w);
    }
    seen
}

fn gap(a: Q, b: Q) -> Q {
    if a > b {
        a - b
    } else {
        b - a
    }
}

fn name(value: Q) -> String {
    format!("{value}")
}

pub fn run() {
    let mut grid = BTreeSet::new();
    for q in 1..=DENOMINATORS {
        for a in 0..=q {
            grid.insert(Ratio::new(a, q));
        }
    }
    let half = Ratio::new(3, 2);
    let (mut early, mut late) = (Vec::new(), Vec::new());
    for x in &grid {
        for z in &grid {
            if let Some(w) = constant(*x, *z, CUT) {
                early.push((*x, *z, w));
                if constant(*x, *z, REACH).is_some() {
                    late.push((*x, *z, w));
                }
            }
        }
    }
    println!("xor pair, code 105 inks the even-weight parity of (floor(x/4), floor(y/4), floor(z/4)) and code 150 the odd");
    println!(
        "  ideal-frame points (X, Z) = a/q, q <= {DENOMINATORS}, {} points a side, whose parity is one value at every odd n: {} up to n = {CUT}, {} up to n = {REACH}",
        grid.len(),
        early.len(),
        late.len()
    );
    for (x, z, w) in &late {
        let y = half - x - z;
        println!(
            "  (X, Y, Z) = ({}, {}, {}): parity {w}, {} at every odd n <= {REACH}",
            name(*x),
            name(y),
            name(*z),
            if *w == 0 {
                "code 105 ink and code 150 paper"
            } else {
                "code 150 ink and code 105 paper"
            }
        );
    }
    let quarter = |v: Q| {
        let scaled = v * 4;
        let near = scaled.round();
        gap(scaled, near) / 4
    };
    let worst = early
        .iter()
        .filter(|(x, z, _)| !late.iter().any(|(a, b, _)| a == x && b == z))
        .map(|(x, z, _)| quarter(*x).max(quarter(*z)))
        .max();
    if let Some(worst) = worst {
        println!(
            "  the {} that fail between n = {CUT} and n = {REACH} sit within {:.4} of the quarter lattice in both coordinates",
            early.len() - late.len(),
            *worst.numer() as f64 / *worst.denom() as f64
        );
    }
    for (label, x, z) in [
        ("(3/4, 3/4, 0)", Ratio::new(3, 4), Ratio::new(0, 1)),
        ("(3/4, 0, 3/4)", Ratio::new(3, 4), Ratio::new(3, 4)),
        ("(0, 3/4, 3/4)", Ratio::new(0, 1), Ratio::new(3, 4)),
    ] {
        let reads: Vec<String> = (1..=15)
            .step_by(2)
            .map(|n| weight(n, x, z).map_or("-".into(), |w| w.to_string()))
            .collect();
        println!(
            "  parity of {label} at n = 1, 3, ..., 15: {}",
            reads.join(" ")
        );
    }
}
