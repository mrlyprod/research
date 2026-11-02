use crate::lattice::{Family, FAMILIES};
use crate::star::{law, slice};
use crate::sums::odds;
use mrly_num::series::{beta, chi4, dirichlet, zeta, CATALAN};
use num_bigint::BigInt;
use num_rational::{BigRational, Ratio};
use std::f64::consts::PI;

const TERMS: usize = 10_000_000;
const LAYERS: usize = 55;

type Q = Ratio<i64>;

struct Shape {
    a: Q,
    b: Q,
    c: Q,
    d: Q,
    e: Q,
    f: Q,
}

fn shape(family: Family) -> Shape {
    let r = |p: i64, q: i64| Ratio::new(p, q);
    match family {
        Family::Carpet => Shape {
            a: r(1, 2),
            b: r(1, 8),
            c: r(1, 2),
            d: r(0, 1),
            e: r(0, 1),
            f: r(-1, 8),
        },
        Family::Net => Shape {
            a: r(1, 2),
            b: r(-1, 8),
            c: r(-1, 2),
            d: r(0, 1),
            e: r(0, 1),
            f: r(1, 8),
        },
        Family::Tree => Shape {
            a: r(1, 4),
            b: r(0, 1),
            c: r(1, 3),
            d: r(-1, 12),
            e: r(1, 6),
            f: r(-1, 6),
        },
        Family::Void => Shape {
            a: r(1, 4),
            b: r(0, 1),
            c: r(0, 1),
            d: r(-1, 4),
            e: r(1, 2),
            f: r(0, 1),
        },
    }
}

fn chi(number: i64) -> i64 {
    -i64::from(chi4(number as usize))
}

fn quasi(shape: &Shape, n: i64) -> Q {
    let (t, square) = (Ratio::from_integer(n), Ratio::from_integer(n * n));
    let s = Ratio::from_integer(chi(n));
    shape.a + shape.b * s + (shape.c + shape.d * s) / t + (shape.e + shape.f * s) / square
}

fn text(value: Q) -> String {
    if *value.denom() == 1 {
        return value.numer().to_string();
    }
    format!("{}/{}", value.numer(), value.denom())
}

fn big(value: Q) -> BigRational {
    BigRational::new(BigInt::from(*value.numer()), BigInt::from(*value.denom()))
}

fn real(value: Q) -> f64 {
    *value.numer() as f64 / *value.denom() as f64
}

fn sums(limit: usize) -> [f64; 4] {
    let mut out = [0.0; 4];
    for n in odds(limit) {
        let x = n as f64;
        let sign = f64::from(chi4(n));
        out[0] += sign;
        out[1] += sign / x;
        out[2] += sign / (x * x);
        out[3] += 1.0 / (x * x);
    }
    out
}

fn shapes_match() {
    for family in FAMILIES {
        let form = shape(family);
        let mut matched = 0;
        let mut layers = 0;
        for n in odds(LAYERS) {
            matched += usize::from(law(family, n) == quasi(&form, n as i64));
            layers += 1;
        }
        println!(
            "  {}: A = {} B = {} c = {} d = {} e = {} f = {}, rebuild the closed form {matched}/{layers}",
            family.name(),
            text(form.a),
            text(form.b),
            text(form.c),
            text(form.d),
            text(form.e),
            text(form.f)
        );
    }
}

fn summed_identity() {
    println!(
        "  the summed identity in exact rationals, counted hexagons against the character sums"
    );
    for family in FAMILIES {
        let form = shape(family);
        let zero = BigRational::from_integer(BigInt::from(0));
        let mut left = zero.clone();
        let (mut s0, mut s1, mut s2, mut s3) = (zero.clone(), zero.clone(), zero.clone(), zero);
        let mut classes = [0usize; 4];
        let mut counts = [0usize; 4];
        for n in odds(LAYERS) {
            let (side, sign) = (n as i64, Ratio::from_integer(i64::from(chi4(n))));
            let ink = slice(family, n).ink();
            left += big(ink - form.a - form.c / Ratio::from_integer(side));
            s0 += big(sign);
            s1 += big(sign / Ratio::from_integer(side));
            s2 += big(sign / Ratio::from_integer(side * side));
            s3 += big(Ratio::new(1, side * side));
            let right =
                -big(form.b) * &s0 - big(form.d) * &s1 + big(form.e) * &s3 - big(form.f) * &s2;
            let slot = (n % 8) / 2;
            counts[slot] += 1;
            classes[slot] += usize::from(left == right);
        }
        println!(
            "  {}: n = 1,3,5,7 mod 8 read {}/{} {}/{} {}/{} {}/{}",
            family.name(),
            classes[0],
            counts[0],
            classes[1],
            counts[1],
            classes[2],
            counts[2],
            classes[3],
            counts[3]
        );
    }
}

fn limits(catalan: f64) -> Vec<(Family, f64, f64)> {
    FAMILIES
        .into_iter()
        .map(|family| {
            let form = shape(family);
            let even =
                -real(form.d) * PI / 4.0 + real(form.e) * PI * PI / 8.0 - real(form.f) * catalan;
            (family, even, even - real(form.b))
        })
        .collect()
}

fn ladder(catalan: f64) {
    let stops = [400usize, 1600, 3200];
    let deepest = stops[stops.len() - 1] + 3;
    for family in FAMILIES {
        let form = shape(family);
        let (b, d, e, f) = (real(form.b), real(form.d), real(form.e), real(form.f));
        let flat = d == 0.0 && e == 0.0;
        let base = -d * PI / 4.0 + e * PI * PI / 8.0 - f * catalan;
        let mut running = 0.0;
        let mut readings: Vec<(usize, f64)> = Vec::new();
        for count in 1..=deepest {
            let n = 2 * count as i64 - 1;
            running += real(law(family, n as usize) - form.a - form.c / Ratio::from_integer(n));
            if !stops.iter().any(|stop| count >= *stop && count <= stop + 3) {
                continue;
            }
            let target = base - if count % 2 == 0 { 0.0 } else { b };
            let scale = count as f64;
            let gap = (running - target) * scale;
            readings.push((count, if flat { gap * scale } else { gap }));
        }
        let power = if flat { "M^2" } else { "M" };
        let even = if flat { f / 8.0 } else { (d - e) / 4.0 };
        let odd = if flat { -f / 8.0 } else { (-d - e) / 4.0 };
        println!(
            "  {}: gap * {power} -> {even:+.8} at even M and {odd:+.8} at odd M",
            family.name()
        );
        for chunk in readings.chunks(4) {
            let cells: Vec<String> = chunk
                .iter()
                .map(|(count, read)| format!("M = {count} {read:+.8}"))
                .collect();
            println!("    {}", cells.join("  "));
        }
    }
}

fn split() {
    let mut cells: Vec<String> = Vec::new();
    for count in [400usize, 1600, 3200] {
        let [_, s1, _, s3] = sums(2 * count - 1);
        let gap = s1 + s3 / 4.0 - (PI / 4.0 + PI * PI / 32.0);
        cells.push(format!("M = {count} {:+.8}", gap * count as f64));
    }
    println!(
        "  carpet split at even M only: gap * M -> {:+.8}, reading {}",
        -5.0 / 16.0,
        cells.join("  ")
    );
}

pub fn run() {
    let catalan = beta(2.0, TERMS);
    let eisenstein = dirichlet(2.0, &[0, 1, -1], TERMS);
    let zeta3 = zeta(3.0, 2_000_000).expect("three is above one");
    println!("constants from their own series");
    println!("  G = {catalan:.10}  L(2, chi_-3) = {eisenstein:.10}  zeta(3) = {zeta3:.10}");
    println!(
        "  carpet split  pi/4 + pi^2/32     = {:.10}",
        PI / 4.0 + PI * PI / 32.0
    );
    println!(
        "  carpet mean   G/8                = {:.10}   G/8 - 1/8 = {:.10}",
        catalan / 8.0,
        catalan / 8.0 - 0.125
    );
    println!(
        "  void mean     (pi + pi^2)/16     = {:.10}",
        (PI + PI * PI) / 16.0
    );
    println!(
        "  flat stack    pi^2 ln2/(7 zeta3) = {:.10}",
        PI * PI * 2f64.ln() / (7.0 * zeta3)
    );
    println!(
        "  tree mean     pi/48+pi^2/48+G/6  = {:.10}",
        PI / 48.0 + PI * PI / 48.0 + catalan / 6.0
    );
    println!("partial character sums over odd n <= N, M layers, from the exact ink laws");
    println!("  the carpet split is an identity at even M only, so it prints at even M only");
    for limit in [53usize, 55] {
        let [s0, s1, s2, s3] = sums(limit);
        let carpet = if limit.div_ceil(2) % 2 == 0 {
            format!("M(I1 - I3 + 1/4) = {:.10}  ", s1 + s3 / 4.0)
        } else {
            String::new()
        };
        println!(
            "  N = {limit}: {carpet}M(mean ink - 1/2 - eps/2) = {:.10}  void M(mean - 1/4) = {:.10}",
            -s0 / 8.0 + s2 / 8.0,
            s1 / 4.0 + s3 / 2.0
        );
    }
    println!("the one-layer law M(mean I - A - c eps) = -B S - d s1 + e s3 - f s2");
    shapes_match();
    summed_identity();
    println!("  the limit -B[M odd] - d pi/4 + e pi^2/8 - f G, both classes of the layer count");
    for (family, even, odd) in limits(CATALAN) {
        println!(
            "  {}: {even:.10} at even M, N = 3 mod 4 and {odd:.10} at odd M, N = 1 mod 4",
            family.name()
        );
    }
    println!("  the approach, every class of M mod 4 so every class of N mod 8");
    ladder(CATALAN);
    split();
}
