use crate::lattice::{Family, Rule, FAMILIES};
use num_rational::Ratio;

type Q = Ratio<i128>;

fn whole(value: Q) -> i128 {
    value.floor().to_integer()
}

fn wrap(value: Q) -> Q {
    value - (value / 2).floor() * 2
}

fn side(step: i128, drift: i128, point: Q, grain: i128) -> Vec<([i128; 4], Q)> {
    let mut cuts: Vec<Q> = (0..=2 * grain).map(|j| Ratio::new(j, grain)).collect();
    if step != 0 {
        let span = 4 * grain * (step.abs() + drift.abs() + 2);
        for k in -span..=span {
            let v = (Ratio::new(k, grain) - point * drift) / step;
            if v >= Q::from_integer(0) && v <= Q::from_integer(2) {
                cuts.push(v);
            }
        }
    }
    cuts.sort();
    cuts.dedup();
    let mut out: Vec<([i128; 4], Q)> = Vec::new();
    for pair in cuts.windows(2) {
        if pair[0] == pair[1] {
            continue;
        }
        let mid = (pair[0] + pair[1]) / 2;
        let far = wrap(mid * step + point * drift);
        let key = [
            whole(mid).rem_euclid(2),
            whole(mid * grain) - whole(mid) * grain,
            whole(far).rem_euclid(2),
            whole(far * grain) - whole(far) * grain,
        ];
        let weight = (pair[1] - pair[0]) / 2;
        match out.iter_mut().find(|(seen, _)| *seen == key) {
            Some((_, total)) => *total += weight,
            None => out.push((key, weight)),
        }
    }
    out
}

fn overtone(residue: i128, p: i128, q: i128) -> i128 {
    let seed = i128::from((p, q) == (0, 0) || (p, q) == (3, 1));
    if residue == 1 {
        seed
    } else {
        1 - seed
    }
}

fn ink(rule: &Rule, sigma: i128, tau: i128, tone: i128) -> Q {
    let middle = (sigma + tau + tone).rem_euclid(2);
    Q::from_integer(i128::from(rule.filled(
        4 * sigma as i64,
        4 * middle as i64,
        4 * tau as i64,
    )))
}

fn field(rule: &Rule, x: Q, z: Q, near: i128, far: i128, step: i128, drift: i128) -> [Q; 3] {
    let across = side(step, drift, x, 4);
    let down = side(step, drift, z, 2);
    let mut out = [Q::from_integer(0); 3];
    for (a, wide) in &across {
        for (c, tall) in &down {
            let weight = wide * tall;
            let one = ink(rule, a[0], c[0], overtone(near, a[1], c[1]));
            let two = ink(rule, a[2], c[2], overtone(far, a[3], c[3]));
            out[0] += weight * one;
            out[1] += weight * two;
            out[2] += weight * one * two;
        }
    }
    out
}

fn clip(poly: &[(Q, Q)], keep: impl Fn(Q, Q) -> Q) -> Vec<(Q, Q)> {
    let zero = Q::from_integer(0);
    let mut out = Vec::new();
    for index in 0..poly.len() {
        let (p, q) = (poly[index], poly[(index + 1) % poly.len()]);
        let (fp, fq) = (keep(p.0, p.1), keep(q.0, q.1));
        if fp >= zero {
            out.push(p);
        }
        if (fp > zero && fq < zero) || (fp < zero && fq > zero) {
            let t = fp / (fp - fq);
            out.push((p.0 + (q.0 - p.0) * t, p.1 + (q.1 - p.1) * t));
        }
    }
    out
}

fn moments(poly: &[(Q, Q)]) -> [Q; 4] {
    let mut out = [Q::from_integer(0); 4];
    for index in 0..poly.len() {
        let (x0, y0) = poly[index];
        let (x1, y1) = poly[(index + 1) % poly.len()];
        let cross = x0 * y1 - x1 * y0;
        out[0] += cross;
        out[1] += (x0 + x1) * cross;
        out[2] += (y0 + y1) * cross;
        out[3] += (x0 * y1 + x0 * y0 * 2 + x1 * y1 * 2 + x1 * y0) * cross;
    }
    [out[0] / 2, out[1] / 6, out[2] / 6, out[3] / 24]
}

fn knots(drift: i128, grain: i128) -> Vec<Q> {
    if drift == 0 {
        return vec![Q::from_integer(0), Ratio::new(1, 2), Q::from_integer(1)];
    }
    let cells = grain * drift.abs();
    (0..=cells).map(|j| Ratio::new(j, cells)).collect()
}

fn limit(rule: &Rule, near: i128, far: i128, step: i128, drift: i128) -> [Q; 3] {
    let across = knots(drift, 4);
    let down = knots(drift, 2);
    let mut totals = [Q::from_integer(0); 3];
    for xi in 0..across.len() - 1 {
        for zi in 0..down.len() - 1 {
            let (x0, x1) = (across[xi], across[xi + 1]);
            let (z0, z1) = (down[zi], down[zi + 1]);
            let cell = vec![(x0, z0), (x1, z0), (x1, z1), (x0, z1)];
            let cell = clip(&cell, |x, z| x + z - Ratio::new(1, 2));
            if cell.len() < 3 {
                continue;
            }
            let cell = clip(&cell, |x, z| Ratio::new(3, 2) - x - z);
            if cell.len() < 3 {
                continue;
            }
            let [area, mx, mz, mxz] = moments(&cell);
            let corner = |x: Q, z: Q| field(rule, x, z, near, far, step, drift);
            let (v00, v01) = (corner(x0, z0), corner(x0, z1));
            let (v10, v11) = (corner(x1, z0), corner(x1, z1));
            let (dx, dz) = (x1 - x0, z1 - z0);
            for slot in 0..3 {
                let c11 = (v11[slot] - v10[slot] - v01[slot] + v00[slot]) / (dx * dz);
                let c10 = (v10[slot] - v00[slot]) / dx - c11 * z0;
                let c01 = (v01[slot] - v00[slot]) / dz - c11 * x0;
                let c00 = v00[slot] - c10 * x0 - c01 * z0 - c11 * x0 * z0;
                totals[slot] += c00 * area + c10 * mx + c01 * mz + c11 * mxz;
            }
        }
    }
    let hexagon = Ratio::new(3, 4);
    [
        totals[0] / hexagon,
        totals[1] / hexagon,
        totals[2] / hexagon,
    ]
}

const KINDS: [(&str, i128, i128); 4] = [
    ("doubling n = 2m+1", 2, 1),
    ("doubling n = 2m-1", 2, -1),
    ("adjacent n = m+2", 1, 2),
    ("gcd echo n = 3m", 3, 0),
];

fn far_residue(near: i128, step: i128, drift: i128) -> i128 {
    (near * step + drift).rem_euclid(4)
}

pub fn run() {
    println!("the pair limit as an exact rational: layer m of residue r mod 4 against layer n = step*m + drift, the phase (mX mod 2, mZ mod 2) equidistributed on the hexagon");
    println!("  the cut cell law s_y = s_x + s_z + w with w = 1 exactly at (p, q) = (0,0) and (3,1) when n = 1 mod 4 and its complement when n = 3 mod 4, integrated over 1/2 <= X + Z <= 3/2 in exact rationals");
    for family in FAMILIES {
        let rule = Rule::new(family);
        for (label, step, drift) in KINDS {
            let mut line = format!("  {} {label}:", family.name());
            for near in [1i128, 3] {
                let far = far_residue(near, step, drift);
                let [em, en, both] = limit(&rule, near, far, step, drift);
                let cov = both - em * en;
                let spread = em * (Q::from_integer(1) - em);
                let same = spread == en * (Q::from_integer(1) - en);
                let r = cov / spread;
                line.push_str(&format!(
                    "  m={near} mod 4: E {em} {en} cov {cov} r {r} = {:+.10}{}",
                    *r.numer() as f64 / *r.denom() as f64,
                    if same { "" } else { " SPREADS DIFFER" }
                ));
            }
            println!("{line}");
        }
    }
    println!("  the four doubling branches read one magnitude, 253/2160, with the sign law -chi4(m) chi4(n); the adjacent and echo carpet rows are -11/135 and 29/135; the tree and void doubling covariances are exactly 0");
    println!("  the exact full-hexagon Pearson against the limit, carpet: the gap times m is bounded, so the finite-layer reading is the limit plus O(1/m)");
    let carpet = Rule::new(Family::Carpet);
    for (m, n, step, drift) in [
        (301usize, 601usize, 2i128, -1i128),
        (601, 1201, 2, -1),
        (103, 205, 2, -1),
        (203, 405, 2, -1),
        (249, 251, 1, 2),
        (99, 297, 3, 0),
    ] {
        let (measured, _) = crate::pairs::pearson(m, n, &carpet);
        let near = (m as i128).rem_euclid(4);
        let [em, en, both] = limit(&carpet, near, far_residue(near, step, drift), step, drift);
        let cov = both - em * en;
        let value = cov / (em * (Q::from_integer(1) - em));
        let target = *value.numer() as f64 / *value.denom() as f64;
        println!(
            "  ({m},{n}): measured {measured:+.8}  limit {value} = {target:+.8}  gap {:.2e}  gap * m {:+.5}",
            (measured - target).abs(),
            (measured - target) * m as f64
        );
    }
}
