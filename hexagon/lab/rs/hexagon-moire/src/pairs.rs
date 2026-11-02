use crate::lattice::{cell, Family, Rule};
use crate::sums::pairwise;
use mrly_num::factor::gcd;
use mrly_num::series::chi4;

fn marks(first: i64, second: i64, span: i64) -> Vec<i64> {
    let mut out = Vec::new();
    let (mut a, mut b) = (0i64, 0i64);
    while a <= span || b <= span {
        let next = a.min(b);
        out.push(next);
        if a == next {
            a += first;
        }
        if b == next {
            b += second;
        }
    }
    out
}

pub fn pearson(m: usize, n: usize, rule: &Rule) -> (f64, f64) {
    let g = gcd(m as u128, n as u128) as i64;
    let (unit_m, unit_n) = (n as i64 / g, m as i64 / g);
    let (m, n) = (m as i64, n as i64);
    let (wide, tall) = (4 * m * n / g, 2 * m * n / g);
    let across = marks(unit_m, unit_n, wide);
    let down = marks(unit_m, unit_n, tall);
    let mut strips: Vec<[f64; 4]> = Vec::new();
    for pair in across.windows(2) {
        let width = (pair[1] - pair[0]) as f64 / wide as f64;
        let (xm, xn) = (pair[0] / unit_m, pair[0] / unit_n);
        let mut total = [0.0; 4];
        for strip in down.windows(2) {
            let height = (strip[1] - strip[0]) as f64 / tall as f64;
            let (zm, zn) = (2 * (strip[0] / unit_m), 2 * (strip[0] / unit_n));
            if let (Some(a), Some(b)) = (cell(rule, m, xm, zm), cell(rule, n, xn, zn)) {
                let weight = width * height;
                total[0] += weight;
                total[1] += weight * f64::from(a);
                total[2] += weight * f64::from(b);
                total[3] += weight * f64::from(a && b);
            }
        }
        strips.push(total);
    }
    let sum = |k: usize| pairwise(&strips.iter().map(|row| row[k]).collect::<Vec<f64>>());
    let area = sum(0);
    let (ea, eb, eab) = (sum(1) / area, sum(2) / area, sum(3) / area);
    let covariance = eab - ea * eb;
    (
        covariance / (ea * (1.0 - ea) * eb * (1.0 - eb)).sqrt(),
        covariance,
    )
}

fn det(a: [[f64; 3]; 3]) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

fn richardson(triple: &[(usize, f64)]) -> f64 {
    let row = |m: usize| [1.0, 1.0 / m as f64, 1.0 / (m * m) as f64];
    let a = [row(triple[0].0), row(triple[1].0), row(triple[2].0)];
    let mut b = a;
    for (slot, (_, value)) in triple.iter().enumerate() {
        b[slot][0] = *value;
    }
    det(b) / det(a)
}

fn branch(first: (usize, f64), second: (usize, f64)) -> f64 {
    (second.0 as f64 * second.1 - first.0 as f64 * first.1) / (second.0 - first.0) as f64
}

const SMALL: [(usize, usize); 7] = [(5, 9), (5, 7), (9, 17), (13, 25), (17, 51), (5, 15), (3, 9)];

const DOUBLING: [(usize, usize); 18] = [
    (3, 5),
    (7, 13),
    (7, 15),
    (11, 21),
    (11, 23),
    (19, 39),
    (23, 45),
    (27, 53),
    (103, 205),
    (149, 297),
    (151, 301),
    (157, 313),
    (201, 401),
    (301, 601),
    (401, 801),
    (403, 805),
    (501, 1001),
    (601, 1201),
];

const FITTED: [usize; 6] = [157, 201, 301, 401, 501, 601];

type Row = (&'static str, &'static [(usize, usize)], &'static [usize]);

pub fn small_pairs(rule: &Rule) {
    println!("full-hexagon exact Pearson, carpet, the seven blind pairs");
    for (m, n) in SMALL {
        let (r, covariance) = pearson(m, n, rule);
        println!("  ({m},{n}) r = {r:+.8}  cov = {covariance:+.8}");
    }
}

pub fn doubling(rule: &Rule) {
    println!("doubling pairs (m, 2m+-1), carpet: sign law sign r = -chi4(m) chi4(n)");
    let mut hits = 0;
    let mut fitted: Vec<(usize, f64)> = Vec::new();
    let mut fitted_cov: Vec<(usize, f64)> = Vec::new();
    let mut branches: Vec<(usize, f64)> = Vec::new();
    for (m, n) in DOUBLING {
        let (r, covariance) = pearson(m, n, rule);
        let predicted = -f64::from(chi4(m) * chi4(n));
        let agreed = r.signum() == predicted;
        hits += usize::from(agreed);
        println!(
            "  ({m},{n}) r = {r:+.8}  cov = {covariance:+.8}  predicted sign {predicted:+}  {}",
            if agreed { "ok" } else { "MISS" }
        );
        if FITTED.contains(&m) {
            fitted.push((m, r));
            fitted_cov.push((m, covariance));
        }
        if [103, 403, 501, 601].contains(&m) {
            branches.push((m, r.abs()));
        }
    }
    println!("  sign law: {hits}/{} pairs", DOUBLING.len());
    println!("Richardson r = r_inf + a/m + b/m^2 on sliding triples of m = {FITTED:?}");
    for start in 0..FITTED.len() - 2 {
        let r = richardson(&fitted[start..start + 3]);
        let c = richardson(&fitted_cov[start..start + 3]);
        println!(
            "  m = {:?}: r_inf = {r:+.8}  cov_inf = {c:+.8}",
            &FITTED[start..start + 3]
        );
    }
    let low = branch(branches[2], branches[3]);
    let high = branch(branches[0], branches[1]);
    println!("  branch 1/m extrapolation |r|: (501,601) -> {low:.7}  (103,403) -> {high:.7}");
    println!(
        "  candidates: 253/2160 = {:.8}  19/162 = {:.8}",
        253.0 / 2160.0,
        19.0 / 162.0
    );
}

pub fn persistence() {
    println!("other pairs and families");
    let rules: Vec<(Family, Rule)> = [Family::Carpet, Family::Tree, Family::Void]
        .into_iter()
        .map(|family| (family, Rule::new(family)))
        .collect();
    let table: [Row; 4] = [
        ("doubling", &[(201, 401)], &[1, 2]),
        ("adjacent", &[(199, 201), (249, 251)], &[0, 1, 2]),
        ("gcd echo (m,3m)", &[(67, 201), (99, 297)], &[0, 1, 2]),
        ("coprime control", &[(101, 173), (97, 251)], &[0, 1, 2]),
    ];
    for (label, pairs, chosen) in table {
        for (m, n) in pairs.iter().copied() {
            let mut line = format!("  {label} ({m},{n}):");
            for slot in chosen {
                let (family, rule) = &rules[*slot];
                let (r, _) = pearson(m, n, rule);
                line.push_str(&format!("  {} {r:+.8}", family.name()));
            }
            println!("{line}");
        }
    }
    println!(
        "  candidates: adjacent -11/135 = {:+.7}  echo 29/135 = {:+.7}",
        -11.0 / 135.0,
        29.0 / 135.0
    );
}
