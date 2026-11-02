use std::f64::consts::PI;

const CUT: usize = (1 << 23) - 1;
const FIRST: u32 = 6;

struct Field {
    name: &'static str,
    units: f64,
    disc: f64,
    chi: fn(usize) -> i64,
}

fn chi_gauss(n: usize) -> i64 {
    match n % 4 {
        1 => 1,
        3 => -1,
        _ => 0,
    }
}

fn chi_eisenstein(n: usize) -> i64 {
    match n % 3 {
        1 => 1,
        2 => -1,
        _ => 0,
    }
}

const FIELDS: [Field; 2] = [
    Field {
        name: "Gaussian",
        units: 4.0,
        disc: 4.0,
        chi: chi_gauss,
    },
    Field {
        name: "Eisenstein",
        units: 6.0,
        disc: 3.0,
        chi: chi_eisenstein,
    },
];

impl Field {
    fn rho(&self) -> f64 {
        2.0 * PI / (self.units * self.disc.sqrt())
    }

    fn l_two(&self, terms: usize) -> f64 {
        (1..=terms)
            .rev()
            .map(|n| (self.chi)(n) as f64 / (n as f64 * n as f64))
            .sum()
    }

    fn constant(&self, terms: usize) -> f64 {
        self.rho() / (2.0 * (PI * PI / 6.0) * self.l_two(terms))
    }
}

struct Tables {
    count: Vec<i64>,
    ideals: Vec<i64>,
    norms: Vec<i64>,
    layers: Vec<i64>,
}

fn mobius(limit: usize) -> Vec<i64> {
    let mut mu = vec![1i64; limit + 1];
    let mut composite = vec![false; limit + 1];
    let mut primes = Vec::new();
    mu[0] = 0;
    for n in 2..=limit {
        if !composite[n] {
            primes.push(n);
            mu[n] = -1;
        }
        for &p in &primes {
            if n * p > limit {
                break;
            }
            composite[n * p] = true;
            if n % p == 0 {
                mu[n * p] = 0;
                break;
            }
            mu[n * p] = -mu[n];
        }
    }
    mu
}

fn tables(field: &Field, limit: usize) -> Tables {
    let mu = mobius(limit);
    let chi: Vec<i64> = (0..=limit).map(|n| (field.chi)(n)).collect();
    let mut ideals = vec![0i64; limit + 1];
    let mut mobius_k = vec![0i64; limit + 1];
    for d in 1..=limit {
        for m in (d..=limit).step_by(d) {
            ideals[m] += chi[d];
            mobius_k[m] += mu[d] * mu[m / d] * chi[m / d];
        }
    }
    let weighted: Vec<i64> = ideals
        .iter()
        .enumerate()
        .map(|(n, &a)| n as i64 * a)
        .collect();
    let mut totient = vec![0i64; limit + 1];
    for i in 1..=limit {
        if weighted[i] == 0 {
            continue;
        }
        for k in 1..=limit / i {
            totient[i * k] += weighted[i] * mobius_k[k];
        }
    }
    let prefix = |v: &Vec<i64>| -> Vec<i64> {
        let mut acc = 0;
        v.iter()
            .map(|&x| {
                acc += x;
                acc
            })
            .collect()
    };
    Tables {
        count: prefix(&totient),
        ideals: prefix(&ideals),
        norms: prefix(&weighted),
        layers: mobius_k,
    }
}

#[derive(Default, Clone)]
struct Window {
    n: usize,
    e54: f64,
    lo54: f64,
    e54_max: f64,
    e43_max: f64,
    e43_at: usize,
    e_nln: f64,
    e_nln_max: f64,
    b1: f64,
    np: f64,
    int_p: f64,
    rest: f64,
    cross: f64,
    mean: f64,
}

fn window_of(n: usize) -> usize {
    (usize::BITS - 1 - n.leading_zeros()) as usize
}

fn rms(sum: f64, n: usize) -> f64 {
    (sum / n as f64).sqrt()
}

fn run(field: &Field) {
    let t = tables(field, CUT);
    let rho = field.rho();
    let c = field.constant(CUT);
    println!(
        "{}: rho {:.12} c {:.12} (L(2,chi) summed to {})",
        field.name, rho, c, CUT
    );
    let top = window_of(CUT);
    let mut wins = vec![Window::default(); top + 1];
    let mut sliver = Window::default();
    let sliver_from = (CUT as f64 * 0.98) as usize;
    let mut record = (0.0f64, 0usize);
    let mut last = 0.0;
    for n in 1usize..=CUT {
        let x = n as f64;
        let e = t.count[n] as f64 - c * x * x;
        let b1 = t.norms[n] as f64 - 0.5 * rho * x * x;
        let np = x * (t.ideals[n] as f64 - rho * x);
        let p54 = x.powf(1.25);
        let p43 = x.powf(4.0 / 3.0);
        let r = e / p54;
        let lo = (1u64 << window_of(n)) as f64;
        if r.abs() > record.0 {
            record = (r.abs(), n);
        }
        last = r;
        let acc = |w: &mut Window| {
            w.n += 1;
            w.e54 += r * r;
            w.lo54 += (e / lo.powf(1.25)).powi(2);
            w.mean += r;
            w.e54_max = w.e54_max.max(r.abs());
            if e.abs() / p43 > w.e43_max {
                w.e43_max = e.abs() / p43;
                w.e43_at = n;
            }
            let q = e / (x * x.ln());
            w.e_nln += q * q;
            w.e_nln_max = w.e_nln_max.max(q.abs());
            w.b1 += (b1 / p54).powi(2);
            w.np += (np / p54).powi(2);
            w.int_p += ((np - b1) / p54).powi(2);
            w.rest += ((e - b1) / p54).powi(2);
            w.cross += r * b1 / p54;
        };
        acc(&mut wins[window_of(n)]);
        if n >= sliver_from {
            acc(&mut sliver);
        }
    }
    println!("window  count  rms(E)/N^5/4  rms(E)/lo^5/4  mean(E)/N^5/4  max|E|/N^5/4  max|E|/N^4/3  at  rms(E)/(N ln N)  max|E|/(N ln N)");
    for (k, w) in wins.iter().enumerate().skip(FIRST as usize) {
        println!(
            "2^{:<2} {:>8} {:>10.4} {:>12.4} {:>12.4} {:>12.4} {:>12.4} {:>9} {:>10.4} {:>14.4}",
            k,
            w.n,
            rms(w.e54, w.n),
            rms(w.lo54, w.n),
            w.mean / w.n as f64,
            w.e54_max,
            w.e43_max,
            w.e43_at,
            rms(w.e_nln, w.n),
            w.e_nln_max
        );
    }
    println!(
        "window  rms(E1)/N^5/4  rms(N P)/N^5/4  rms(int P)/N^5/4  rms(E-E1)/N^5/4  corr(E,E1)"
    );
    for (k, w) in wins.iter().enumerate().skip(FIRST as usize) {
        let e = w.e54;
        let corr = w.cross / (e * w.b1).sqrt();
        println!(
            "2^{:<2} {:>12.4} {:>14.4} {:>16.4} {:>16.4} {:>11.4}",
            k,
            rms(w.b1, w.n),
            rms(w.np, w.n),
            rms(w.int_p, w.n),
            rms(w.rest, w.n),
            corr
        );
    }
    let fit = |from: usize| -> f64 {
        let pts: Vec<(f64, f64)> = (from..=top)
            .map(|k| {
                let w = &wins[k];
                (
                    (1.5 * (1u64 << k) as f64).ln(),
                    rms(w.e54, w.n).ln() + 1.25 * (1.5 * (1u64 << k) as f64).ln(),
                )
            })
            .collect();
        let m = pts.len() as f64;
        let sx: f64 = pts.iter().map(|p| p.0).sum();
        let sy: f64 = pts.iter().map(|p| p.1).sum();
        let sxx: f64 = pts.iter().map(|p| p.0 * p.0).sum();
        let sxy: f64 = pts.iter().map(|p| p.0 * p.1).sum();
        (m * sxy - sx * sy) / (m * sxx - sx * sx)
    };
    println!(
        "fit of ln rms(E) on ln N: exponent {:.4} over windows 2^{}..2^{}, {:.4} over 2^{}..2^{}, {:.4} over 2^{}..2^{}",
        fit(FIRST as usize),
        FIRST,
        top,
        fit(14),
        14,
        top,
        fit(top - 3),
        top - 3,
        top
    );
    println!(
        "cut N = {} (every integer N from 1 to the cut counted, top window 2^{} complete)",
        CUT, top
    );
    println!(
        "within 0.02 of the cut (N from {}): count {} rms(E)/N^5/4 {:.4} max|E|/N^4/3 {:.4} at {}, E/N^5/4 at the cut {:+.4}",
        sliver_from,
        sliver.n,
        rms(sliver.e54, sliver.n),
        sliver.e43_max,
        sliver.e43_at,
        last
    );
    println!(
        "record |E|/N^5/4 {:.4} at N = {} ({:.4} of the cut)",
        record.0,
        record.1,
        record.1 as f64 / CUT as f64
    );
}

fn main() {
    for field in &FIELDS {
        run(field);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const BOUNDS: [usize; 9] = [50, 200, 800, 1200, 2000, 3200, 12800, 51200, 102400];
    const GAUSS: [i64; 9] = [
        672, 10608, 168088, 372872, 1045088, 2663864, 42663808, 683283504, 2732257376,
    ];
    const EISENSTEIN: [i64; 9] = [
        630, 9606, 151020, 337026, 945486, 2419950, 38560362, 616981926, 2466558234,
    ];

    #[test]
    fn totient_sums_match_the_class_route() {
        for (field, pinned) in FIELDS.iter().zip([GAUSS, EISENSTEIN]) {
            let t = tables(field, 102400);
            for (b, want) in BOUNDS.iter().zip(pinned) {
                assert_eq!(t.count[*b], want, "{} at {}", field.name, b);
            }
        }
    }

    #[test]
    fn layers_sum_to_the_error() {
        let limit = 400_000;
        for field in &FIELDS {
            let t = tables(field, limit);
            let rho = field.rho();
            let c = field.constant(limit);
            for n in [100usize, 250, 400] {
                let e = t.count[n] as f64 - c * (n * n) as f64;
                let layered: f64 = (1..=limit)
                    .map(|m| {
                        let y = n as f64 / m as f64;
                        let e1 = t.norms[y as usize] as f64 - 0.5 * rho * y * y;
                        t.layers[m] as f64 * e1
                    })
                    .sum();
                assert!((e - layered).abs() < 0.01, "{} at {}", field.name, n);
            }
        }
    }

    #[test]
    fn ideal_counts_match_lattice_points() {
        let limit = 2000i64;
        for field in &FIELDS {
            let t = tables(field, limit as usize);
            let mut points = 0i64;
            for a in -80..=80i64 {
                for b in -80..=80i64 {
                    let norm = if field.units == 4.0 {
                        a * a + b * b
                    } else {
                        a * a - a * b + b * b
                    };
                    if norm <= limit {
                        points += 1;
                    }
                }
            }
            assert_eq!(
                t.ideals[limit as usize],
                (points - 1) / field.units as i64,
                "{}",
                field.name
            );
        }
    }
}
