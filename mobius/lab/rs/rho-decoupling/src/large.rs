use num_bigint::BigInt;

use crate::riesz::{alpha, band, families, moment, theta_band, transform_table, Family};
use crate::{bigpow, ratio_f64};

// GRID

pub struct Grid {
    pub l: usize,
    pub n: usize,
    pub k: usize,
    pub fh: Vec<(f64, f64)>,
    pub sq: Vec<f64>,
}

pub fn grid(fam: &Family, l: usize) -> Grid {
    let q = fam.q as usize;
    let n = q.pow(l as u32);
    let g = transform_table(fam.q, &fam.digits, l);
    let mut fh = vec![(0.0f64, 0.0f64); n];
    let mut sq = vec![0.0f64; n];
    for a in 0..n {
        let mut b = a;
        let (mut re, mut im) = (1.0f64, 0.0f64);
        for _ in 0..l {
            let (gr, gi) = g[b];
            let nr = re * gr - im * gi;
            let ni = re * gi + im * gr;
            re = nr;
            im = ni;
            b = (b * q) % n;
        }
        fh[a] = (re, im);
        sq[a] = re * re + im * im;
    }
    Grid {
        l,
        n,
        k: fam.digits.len(),
        fh,
        sq,
    }
}

impl Grid {
    fn kl(&self) -> f64 {
        (self.k as f64).powi(self.l as i32)
    }

    fn l1(&self) -> f64 {
        self.sq.iter().map(|v| v.sqrt()).sum()
    }

    fn second(&self) -> f64 {
        (1..self.n)
            .map(|a| self.sq[a])
            .fold(0.0f64, f64::max)
            .sqrt()
            / self.kl()
    }

    pub fn large_set(&self, eta: f64) -> Vec<usize> {
        let t2 = self.kl().powi(2) * (self.n as f64).powf(-2.0 * eta) * (1.0 - 1e-9);
        (0..self.n).filter(|&a| self.sq[a] >= t2).collect()
    }
}

// ROUNDING

fn read(v: f64, digits: usize) -> String {
    band(v - 1e-12, v + 1e-12, digits)
}

fn down(v: f64, digits: usize) -> String {
    band(v - 1e-12, v - 1e-12, digits)
}

fn up(v: f64, digits: usize) -> String {
    band(v + 1e-12, v + 1e-12, digits)
}

fn sci(v: f64, upward: bool) -> String {
    if v <= 0.0 {
        return "0".to_string();
    }
    let e = v.log10().floor() as i32;
    let m = v / 10f64.powi(e);
    let scaled = m * 1e5 * if upward { 1.0 + 1e-12 } else { 1.0 - 1e-12 };
    let r = if upward {
        scaled.ceil()
    } else {
        scaled.floor()
    };
    if r >= 1e6 {
        format!("1.00000e{}", e + 1)
    } else {
        format!("{:.5}e{}", r / 1e5, e)
    }
}

// STRUCTURE

pub struct Shape {
    pub min_gap: usize,
    pub runs: usize,
    pub max_run: usize,
}

pub fn shape(set: &[usize], n: usize) -> Shape {
    let m = set.len();
    if m < 2 {
        return Shape {
            min_gap: n,
            runs: m,
            max_run: m,
        };
    }
    let gaps: Vec<usize> = (0..m)
        .map(|i| {
            if i + 1 < m {
                set[i + 1] - set[i]
            } else {
                n - set[i] + set[0]
            }
        })
        .collect();
    let min_gap = *gaps.iter().min().unwrap();
    let breaks = gaps.iter().filter(|&&d| d != 1).count();
    if breaks == 0 {
        return Shape {
            min_gap,
            runs: 1,
            max_run: m,
        };
    }
    let start = (0..m).find(|&i| gaps[(i + m - 1) % m] != 1).unwrap();
    let mut max_run = 0usize;
    let mut run = 0usize;
    for t in 0..m {
        run += 1;
        if gaps[(start + t) % m] != 1 {
            max_run = max_run.max(run);
            run = 0;
        }
    }
    Shape {
        min_gap,
        runs: breaks,
        max_run,
    }
}

// GRAM

pub struct Gram {
    pub lo: f64,
    pub hi: f64,
    pub iterations: usize,
}

pub fn gram_top(g: &Grid, set: &[usize], iterations: usize) -> Gram {
    let m = set.len();
    let entry = |i: usize, j: usize| -> (f64, f64) {
        let d = (set[i] + g.n - set[j]) % g.n;
        g.fh[d]
    };
    let mut gersh = 0.0f64;
    let mut frob = 0.0f64;
    for i in 0..m {
        let mut row = 0.0;
        for j in 0..m {
            let (re, im) = entry(i, j);
            let mag = (re * re + im * im).sqrt();
            row += mag;
            frob += mag * mag;
        }
        gersh = gersh.max(row);
    }
    let q = (g.n as f64).powf(1.0 / g.l as f64).round() as usize;
    let mut level = g.l;
    let mut unit = 1usize;
    while level > 0 && set.iter().all(|&a| a % (unit * q) == 0) {
        unit *= q;
        level -= 1;
    }
    let witness = (g.k as f64).powi((g.l - level) as i32) * (q as f64).powi(level as i32);
    let hi = gersh.min(frob.sqrt()).min(g.n as f64).min(witness);
    let mut v: Vec<(f64, f64)> = (0..m)
        .map(|i| (1.0 + 0.1 * (((i * 7919 + 13) % 97) as f64) / 97.0, 0.0))
        .collect();
    let mut lo = 0.0f64;
    let mut used = 0usize;
    for it in 0..iterations {
        let mut w = vec![(0.0f64, 0.0f64); m];
        for (i, slot) in w.iter_mut().enumerate() {
            let (mut sr, mut si) = (0.0f64, 0.0f64);
            for (j, &(vr, vi)) in v.iter().enumerate() {
                let (gr, gi) = entry(i, j);
                sr += gr * vr - gi * vi;
                si += gr * vi + gi * vr;
            }
            *slot = (sr, si);
        }
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for i in 0..m {
            let (vr, vi) = v[i];
            let (wr, wi) = w[i];
            num += vr * wr + vi * wi;
            den += vr * vr + vi * vi;
        }
        let ray = num / den;
        used = it + 1;
        let done = (ray - lo).abs() <= 1e-13 * ray.abs();
        lo = lo.max(ray);
        let norm = w.iter().map(|(a, b)| a * a + b * b).sum::<f64>().sqrt();
        v = w.iter().map(|(a, b)| (a / norm, b / norm)).collect();
        if done {
            break;
        }
    }
    Gram {
        lo,
        hi: hi.max(lo),
        iterations: used,
    }
}

// CELL

pub struct Cell {
    pub eta: f64,
    pub count: usize,
    pub r9: f64,
    pub parseval: f64,
    pub shape: Shape,
    pub floor: f64,
    pub gram: Option<Gram>,
    pub small: f64,
    pub large: f64,
    pub large_l2: f64,
    pub large9: f64,
    pub c_sup: f64,
    pub c: f64,
}

pub fn cell(fam: &Family, g: &Grid, s4: f64, eta: f64, cap: usize) -> Cell {
    let x = g.n as f64;
    let lx = x.ln();
    let al = alpha(fam);
    let kl = g.kl();
    let set = g.large_set(eta);
    let count = set.len();
    let r9 = s4 / kl.powi(4) * x.powf(4.0 * eta);
    let parseval = x / kl * x.powf(2.0 * eta);
    assert!(count as f64 <= r9 * (1.0 + 1e-6), "R.9 count fails");
    assert!(
        count as f64 <= parseval * (1.0 + 1e-6),
        "Parseval count fails"
    );
    let floor = set.iter().map(|&a| g.sq[a]).sum::<f64>() / kl;
    let gram = if count <= cap {
        Some(gram_top(g, &set, 300))
    } else {
        None
    };
    if let Some(gr) = &gram {
        assert!(gr.lo >= floor * (1.0 - 1e-9), "Rayleigh below the floor");
        assert!(gr.lo <= x * (1.0 + 1e-9), "Gram above the full grid");
    }
    let small = (al + 0.5 - eta).min((1.0 + al) / 2.0);
    let large = al + (count as f64).ln() / (2.0 * lx);
    let large_l2 = (floor * kl).ln() / (2.0 * lx);
    let large9 = al + r9.ln() / (2.0 * lx);
    let cap = (1.0 + al) / 2.0;
    Cell {
        eta,
        count,
        r9,
        parseval,
        shape: shape(&set, g.n),
        floor,
        gram,
        small,
        large,
        large_l2,
        large9,
        c_sup: al - small.max(large),
        c: al - small.max(large.min(large_l2).min(cap)),
    }
}

pub fn set_row(fam: &Family, g: &Grid, c: &Cell) -> String {
    let gram = match &c.gram {
        Some(gr) => format!(
            "{}..{} ({})",
            sci(gr.lo, false),
            sci(gr.hi, true),
            gr.iterations
        ),
        None => "skipped".to_string(),
    };
    format!(
        "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
        fam.q,
        fam.label,
        g.l,
        read(c.eta, 6),
        c.count,
        sci(c.r9, true),
        sci(c.parseval, true),
        c.shape.min_gap,
        c.shape.runs,
        c.shape.max_run,
        sci(c.floor, false),
        gram,
        (g.k as u128).pow(g.l as u32),
        g.n
    )
}

pub fn chain_row(fam: &Family, g: &Grid, c: &Cell, c1: f64) -> String {
    format!(
        "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
        fam.q,
        fam.label,
        g.l,
        read(c.eta, 6),
        read(alpha(fam), 6),
        read(c1, 6),
        up(c.small, 6),
        up(c.large, 6),
        up(c.large_l2, 6),
        up(c.large9, 6),
        down(c.c_sup, 6),
        down(c.c, 6)
    )
}

// WITNESS

fn member(q: u64, digits: &[u64], l: usize, mut v: u64) -> bool {
    for _ in 0..l {
        if !digits.contains(&(v % q)) {
            return false;
        }
        v /= q;
    }
    v == 0
}

pub struct Witness {
    pub m: u64,
    pub r: u64,
    pub rq: u64,
    pub cntq: u64,
    pub bal: f64,
}

pub fn witness(fam: &Family, l: usize) -> Witness {
    let q = fam.q;
    let x = q.pow(l as u32);
    let m = ((x as f64).sqrt() / 2.0).floor() as u64;
    let mut r = 0u64;
    let mut rq = 0u64;
    let mut cntq = 0u64;
    for mm in m..2 * m {
        let div = mm % q == 0;
        if div {
            cntq += 1;
        }
        for nn in m..2 * m {
            if member(q, &fam.digits, l, mm * nn) {
                r += 1;
                if div {
                    rq += 1;
                }
            }
        }
    }
    let kl = (fam.digits.len() as f64).powi(l as i32);
    let bal = rq as f64 - kl / x as f64 * cntq as f64 * m as f64;
    Witness {
        m,
        r,
        rq,
        cntq,
        bal,
    }
}

pub fn witness_row(fam: &Family, l: usize, w: &Witness) -> String {
    let x = (fam.q as f64).powi(l as i32);
    let kl = (fam.digits.len() as f64).powi(l as i32);
    let exp = |v: f64| -> String {
        if v > 0.0 {
            read(v.ln() / x.ln(), 6)
        } else {
            "void".to_string()
        }
    };
    format!(
        "| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
        fam.q,
        fam.label,
        l,
        w.m,
        w.r,
        read(w.r as f64 / kl, 6),
        exp(w.r as f64),
        read(alpha(fam), 6),
        w.rq,
        w.cntq,
        read(w.bal, 1),
        read(w.bal / kl, 6),
        exp(w.bal.abs())
    )
}

// RUN

pub fn cells() -> Vec<(u64, &'static str, usize)> {
    vec![
        (3, "01", 8),
        (3, "01", 12),
        (4, "012", 6),
        (4, "012", 9),
        (5, "0123", 6),
        (5, "0123", 8),
        (10, "ex7", 4),
        (10, "ex7", 6),
        (100, "01", 3),
        (100, "0to49", 2),
        (100, "0to49", 3),
    ]
}

pub fn etas(fam: &Family, l: usize) -> (f64, f64, Vec<f64>) {
    let m = moment(fam, 4, l.max(4));
    let (tlo, thi) = theta_band(fam, &m);
    let al = alpha(fam);
    let eta4 = (1.0 + 3.0 * al - (tlo + thi) / 2.0) / 2.0;
    let s4 = ratio_f64(
        &(bigpow(&BigInt::from(fam.q), l) * &m.emod[l]),
        &BigInt::from(1),
    );
    (s4, eta4, vec![0.0, eta4 / 2.0, eta4, 2.0 * eta4, 0.25, 0.5])
}

pub fn cell_row(fam: &Family, g: &Grid, s4: f64, eta4: f64, c1: f64) -> String {
    let x = g.n as f64;
    let second = g.second();
    format!(
        "| {} | {} | {} | {} | {} | {} | {} | {} | {} |",
        fam.q,
        fam.label,
        g.l,
        read(alpha(fam), 6),
        read((s4 / g.kl().powi(4)).ln() / x.ln(), 6),
        read(eta4, 6),
        read(c1, 6),
        read(second, 6),
        read(-second.ln() / x.ln(), 6)
    )
}

pub fn run() {
    let fams = families();
    let cap = 1500usize;
    let mut sets = Vec::new();
    let mut chains = Vec::new();
    let mut wits = Vec::new();
    let mut heads = Vec::new();
    for (q, label, l) in cells() {
        let fam = fams.iter().find(|f| f.q == q && f.label == label).unwrap();
        let g = grid(fam, l);
        let x = g.n as f64;
        let kl = g.kl();
        let (s4, eta4, etas) = etas(fam, l);
        let total2: f64 = g.sq.iter().sum();
        let total4: f64 = g.sq.iter().map(|v| v * v).sum();
        assert!((total2 - x * kl).abs() <= 1e-10 * x * kl, "Parseval fails");
        assert!((total4 - s4).abs() <= 1e-10 * s4, "fourth moment fails");
        let c1 = g.l1().ln() / x.ln();
        heads.push(cell_row(fam, &g, s4, eta4, c1));
        for eta in etas {
            let c = cell(fam, &g, s4, eta, cap);
            sets.push(set_row(fam, &g, &c));
            chains.push(chain_row(fam, &g, &c, c1));
        }
        let w = witness(fam, l);
        wits.push(witness_row(fam, l, &w));
    }
    println!("riesz large values cells");
    println!("| q | F | L | alpha | nu_4(L) | eta_4 | c_1 | second/k^L | eta_1 |");
    for r in heads {
        println!("{r}");
    }
    println!("riesz large values");
    println!("| q | F | L | eta | #A | R.9 | Parseval | min gap | runs | max run | Delta_L floor | Delta_L exact | k^L | q^L |");
    for r in sets {
        println!("{r}");
    }
    println!("riesz large values chain");
    println!(
        "| q | F | L | eta | alpha | c_1 | small | large sup | large l2 | large R.9 | c sup | c |"
    );
    for r in chains {
        println!("{r}");
    }
    println!("riesz large values witness");
    println!("| q | F | L | M | R | R/k^L | log R/log x | alpha | R_q | multiples of q | bal | bal/k^L | log abs bal/log x |");
    for r in wits {
        println!("{r}");
    }
}

// TESTS

#[cfg(test)]
mod tests {
    use super::*;

    fn brute(q: u64, digits: &[u64], l: usize) -> Vec<(f64, f64)> {
        let n = (q as usize).pow(l as u32);
        let mut strings = vec![0u64];
        for _ in 0..l {
            let mut next = Vec::new();
            for &s in &strings {
                for &f in digits {
                    next.push(s * q + f);
                }
            }
            strings = next;
        }
        (0..n)
            .map(|a| {
                let y = a as f64 / n as f64;
                let mut re = 0.0;
                let mut im = 0.0;
                for &s in &strings {
                    let t = 2.0 * std::f64::consts::PI * (s as f64) * y;
                    re += t.cos();
                    im += t.sin();
                }
                (re, im)
            })
            .collect()
    }

    fn fam(q: u64, label: &str) -> Family {
        families()
            .into_iter()
            .find(|f| f.q == q && f.label == label)
            .unwrap()
    }

    #[test]
    fn grid_matches_brute() {
        for (q, label, l) in [(3, "01", 6), (4, "012", 4), (5, "0123", 3), (10, "ex7", 3)] {
            let f = fam(q, label);
            let g = grid(&f, l);
            let b = brute(q, &f.digits, l);
            for (got, want) in g.fh.iter().zip(&b).take(g.n) {
                assert!((got.0 - want.0).abs() < 1e-7 && (got.1 - want.1).abs() < 1e-7);
            }
        }
    }

    #[test]
    fn large_set_matches_brute() {
        let f = fam(3, "01");
        let g = grid(&f, 6);
        let b = brute(3, &f.digits, 6);
        for eta in [0.0, 0.05, 0.1, 0.2, 0.5] {
            let t = 64.0 * 729f64.powf(-eta);
            let direct: Vec<usize> = (0..729)
                .filter(|&a| (b[a].0 * b[a].0 + b[a].1 * b[a].1).sqrt() >= t * (1.0 - 1e-9))
                .collect();
            assert_eq!(g.large_set(eta), direct);
        }
    }

    #[test]
    fn gram_pins() {
        let f = fam(3, "01");
        let g = grid(&f, 6);
        let one = gram_top(&g, &[0], 10);
        assert!((one.lo - 64.0).abs() < 1e-9 && (one.hi - 64.0).abs() < 1e-9);
        let sub: Vec<usize> = (0..9).map(|r| r * 81).collect();
        let w = gram_top(&g, &sub, 300);
        assert!((w.lo - 144.0).abs() < 1e-6 && w.hi >= 144.0 - 1e-6 && w.hi <= 144.0 + 1e-6);
        let part: Vec<usize> = vec![0, 81, 162];
        let pw = gram_top(&g, &part, 300);
        assert!(pw.hi <= 144.0 + 1e-6 && pw.lo <= pw.hi);
        assert_eq!(sci(0.38110351, false), "3.81103e-1");
        assert_eq!(sci(0.38110351, true), "3.81104e-1");
        assert_eq!(sci(99999.96, true), "1.00000e5");
        let full: Vec<usize> = (0..729).collect();
        let fg = gram_top(&g, &full, 300);
        assert!((fg.lo - 729.0).abs() < 1e-6 && fg.hi <= 729.0 * (1.0 + 1e-9));
    }

    #[test]
    fn shape_pins() {
        let s = shape(&[0, 1, 2, 10, 11, 20], 30);
        assert_eq!((s.min_gap, s.runs, s.max_run), (1, 3, 3));
        let w = shape(&[0, 1, 28, 29], 30);
        assert_eq!((w.min_gap, w.runs, w.max_run), (1, 1, 4));
        let e = shape(&[5], 30);
        assert_eq!((e.min_gap, e.runs, e.max_run), (30, 1, 1));
    }

    #[test]
    fn witness_matches_enumeration() {
        let f = fam(3, "01");
        let l = 8;
        let mut strings = vec![0u64];
        for _ in 0..l {
            let mut next = Vec::new();
            for &s in &strings {
                for &d in &f.digits {
                    next.push(s * 3 + d);
                }
            }
            strings = next;
        }
        let w = witness(&f, l);
        let mut r = 0;
        let mut rq = 0;
        for mm in w.m..2 * w.m {
            for nn in w.m..2 * w.m {
                if strings.contains(&(mm * nn)) {
                    r += 1;
                    if mm % 3 == 0 {
                        rq += 1;
                    }
                }
            }
        }
        assert_eq!((w.r, w.rq), (r, rq));
        assert_eq!(w.cntq, (w.m..2 * w.m).filter(|m| m % 3 == 0).count() as u64);
    }
}
