use mrly_math::bang::Code;
use mrly_math::counts::fill;
use num_bigint::{BigInt, Sign};
use num_rational::BigRational;

pub type Poly = Vec<BigRational>;

const EXTRA: usize = 5;

fn rational(value: i128) -> BigRational {
    BigRational::from_integer(BigInt::from(value))
}

fn is_zero(value: &BigRational) -> bool {
    value.numer().sign() == Sign::NoSign
}

pub fn trim(mut poly: Poly) -> Poly {
    while poly.last().is_some_and(is_zero) {
        poly.pop();
    }
    poly
}

fn add(a: &[BigRational], b: &[BigRational]) -> Poly {
    let zero = rational(0);
    let width = a.len().max(b.len());
    (0..width)
        .map(|i| a.get(i).unwrap_or(&zero) + b.get(i).unwrap_or(&zero))
        .collect()
}

fn multiply(a: &[BigRational], b: &[BigRational]) -> Poly {
    let mut out = vec![rational(0); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        for (j, y) in b.iter().enumerate() {
            out[i + j] += x * y;
        }
    }
    out
}

pub fn evaluate(poly: &[BigRational], x: i128) -> BigRational {
    let point = rational(x);
    poly.iter()
        .rev()
        .fold(rational(0), |acc, coefficient| acc * &point + coefficient)
}

pub fn interpolate(points: &[(i128, u128)]) -> Poly {
    let mut out = Poly::new();
    for (i, &(xi, yi)) in points.iter().enumerate() {
        let mut term = vec![rational(yi as i128)];
        for (j, &(xj, _)) in points.iter().enumerate() {
            if j != i {
                let scale = rational(xi - xj);
                term = multiply(&term, &[-rational(xj) / &scale, rational(1) / &scale]);
            }
        }
        out = add(&out, &term);
    }
    trim(out)
}

pub fn degree(poly: &[BigRational]) -> Option<usize> {
    trim(poly.to_vec()).len().checked_sub(1)
}

pub fn leading(poly: &[BigRational]) -> BigRational {
    trim(poly.to_vec())
        .last()
        .cloned()
        .unwrap_or_else(|| rational(0))
}

pub fn fit(code: Code, dimension: usize, base: usize) -> (Vec<Poly>, bool) {
    let polys: Vec<Poly> = (0..base)
        .map(|class| {
            let points: Vec<(i128, u128)> = (0..)
                .map(|j| class + j * base)
                .filter(|&n| n >= 1)
                .take(dimension + 1 + EXTRA)
                .map(|n| {
                    let value = fill(code, n, dimension, 1, base).expect("the fill counts");
                    (n as i128, value)
                })
                .collect();
            let poly = interpolate(&points[..=dimension]);
            for &(x, y) in &points {
                assert!(
                    evaluate(&poly, x) == rational(y as i128),
                    "the class polynomial misses the fill"
                );
            }
            poly
        })
        .collect();
    let collapses = polys.iter().all(|poly| *poly == polys[0]);
    (polys, collapses)
}

pub fn fraction(value: &BigRational) -> String {
    if value.denom() == &BigInt::from(1) {
        value.numer().to_string()
    } else {
        format!("{}/{}", value.numer(), value.denom())
    }
}

fn term(magnitude: &BigRational, degree: usize) -> String {
    if degree == 0 {
        return fraction(magnitude);
    }
    let monomial = if degree == 1 {
        String::from("n")
    } else {
        format!("n**{degree}")
    };
    let one = BigInt::from(1);
    match (magnitude.numer() == &one, magnitude.denom() == &one) {
        (true, true) => monomial,
        (true, false) => format!("{monomial}/{}", magnitude.denom()),
        (false, true) => format!("{}*{monomial}", magnitude.numer()),
        (false, false) => format!("{}*{monomial}/{}", magnitude.numer(), magnitude.denom()),
    }
}

pub fn text(poly: &[BigRational]) -> String {
    let trimmed = trim(poly.to_vec());
    if trimmed.is_empty() {
        return String::from("0");
    }
    let mut out = String::new();
    for degree in (0..trimmed.len()).rev() {
        let coefficient = &trimmed[degree];
        if is_zero(coefficient) {
            continue;
        }
        let negative = coefficient.numer().sign() == Sign::Minus;
        let magnitude = if negative {
            -coefficient.clone()
        } else {
            coefficient.clone()
        };
        if out.is_empty() {
            if negative {
                out.push('-');
            }
        } else {
            out.push_str(if negative { " - " } else { " + " });
        }
        out.push_str(&term(&magnitude, degree));
    }
    out
}
