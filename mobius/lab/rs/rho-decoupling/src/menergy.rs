// RANDOM

pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed ^ 0x9e3779b97f4a7c15)
    }

    fn step(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }

    fn wide(&mut self) -> u128 {
        ((self.step() as u128) << 64) | self.step() as u128
    }

    pub fn sign(&mut self) -> f64 {
        if self.step() & 1 == 0 {
            1.0
        } else {
            -1.0
        }
    }
}

// COLUMNS

pub fn column(q: u64, digits: &[u64], l: usize) -> Vec<u128> {
    let mut out = vec![0u128];
    for _ in 0..l {
        let mut next = Vec::with_capacity(out.len() * digits.len());
        for &v in &out {
            for &f in digits {
                next.push(v * q as u128 + f as u128);
            }
        }
        out = next;
    }
    out.retain(|&v| v != 0);
    out.sort_unstable();
    out
}

pub fn random_column(x: u128, count: usize, rng: &mut Rng) -> Vec<u128> {
    let mut out: Vec<u128> = Vec::with_capacity(count);
    while out.len() < count {
        for _ in 0..count - out.len() {
            out.push(1 + rng.wide() % (x - 1));
        }
        out.sort_unstable();
        out.dedup();
    }
    out
}

// GLOBAL ENERGY

pub trait Prod: Copy + Ord {
    fn times(self, other: Self) -> Self;
    fn scramble(self) -> u64;
}

impl Prod for u64 {
    fn times(self, other: u64) -> u64 {
        self * other
    }

    fn scramble(self) -> u64 {
        (self ^ (self >> 33)).wrapping_mul(0xff51afd7ed558ccd) >> 24
    }
}

impl Prod for u128 {
    fn times(self, other: u128) -> u128 {
        self * other
    }

    fn scramble(self) -> u64 {
        let v = (self as u64) ^ ((self >> 64) as u64);
        (v ^ (v >> 33)).wrapping_mul(0xff51afd7ed558ccd) >> 24
    }
}

pub fn energy<T: Prod>(vals: &[T], cap: usize) -> (u128, u32) {
    let n = vals.len() as u128;
    let total = n * n;
    let mut parts = 1usize;
    while total / parts as u128 > cap as u128 {
        parts <<= 1;
    }
    let mask = (parts - 1) as u64;
    let mut e: u128 = 0;
    let mut top: u32 = 0;
    let mut buf: Vec<T> = Vec::with_capacity((total / parts as u128) as usize * 5 / 4 + 64);
    for p in 0..parts as u64 {
        buf.clear();
        for &a in vals {
            for &b in vals {
                let m = a.times(b);
                if parts == 1 || m.scramble() & mask == p {
                    buf.push(m);
                }
            }
        }
        buf.sort_unstable();
        let mut i = 0;
        while i < buf.len() {
            let mut j = i + 1;
            while j < buf.len() && buf[j] == buf[i] {
                j += 1;
            }
            let r = (j - i) as u128;
            e += r * r;
            top = top.max((j - i) as u32);
            i = j;
        }
    }
    (e, top)
}

pub fn energy_of(vals: &[u128], cap: usize) -> (u128, u32) {
    if *vals.last().unwrap() < 1u128 << 32 {
        let small: Vec<u64> = vals.iter().map(|&v| v as u64).collect();
        energy(&small, cap)
    } else {
        energy(vals, cap)
    }
}

// SHIFT SOLUTIONS

fn live(k: u128, m: i64) -> u128 {
    if m >= 1 {
        let mut p = k - 1;
        for _ in 1..m {
            p *= k;
        }
        p
    } else {
        0
    }
}

pub fn shift_excess(k: u128, l: usize) -> u128 {
    let mut t = 0u128;
    for s in 0..2 * l {
        for i in 0..=s {
            for ip in 0..=s {
                if i == ip {
                    continue;
                }
                let a = live(k, l as i64 - i.max(ip) as i64);
                let b = live(k, l as i64 - (s - i).max(s - ip) as i64);
                if a == 0 || b == 0 {
                    continue;
                }
                t += a * b - a.min(b);
            }
        }
    }
    t
}

// MEMBERSHIP

pub struct Bits {
    w: Vec<u64>,
}

impl Bits {
    fn new(x: usize) -> Bits {
        Bits {
            w: vec![0u64; (x >> 6) + 1],
        }
    }

    fn set(&mut self, n: usize) {
        self.w[n >> 6] |= 1u64 << (n & 63);
    }

    pub fn get(&self, n: usize) -> bool {
        self.w[n >> 6] >> (n & 63) & 1 == 1
    }
}

pub fn bits_of(vals: &[u128], x: u128) -> Bits {
    let mut b = Bits::new(x as usize);
    for &v in vals {
        b.set(v as usize);
    }
    b
}

// BOXES

pub struct BoxRow {
    pub m: u64,
    pub n: u64,
    pub r: u128,
    pub e: u128,
    pub ebal: f64,
    pub diag: f64,
    pub triv: f64,
    pub bound: f64,
    pub check: (f64, f64),
}

pub fn box_row(bits: &Bits, mm: u64, nn: u64, delta: f64, cap: usize) -> Option<BoxRow> {
    let mut cs = vec![0u32; nn as usize];
    let mut r: u128 = 0;
    let mut e: u128 = 0;
    let mut spread = 0.0f64;
    for m in mm..2 * mm {
        let mut c = 0u32;
        for l in nn..2 * nn {
            if bits.get((m * l) as usize) {
                c += 1;
                cs[(l - nn) as usize] += 1;
            }
        }
        r += c as u128;
        e += (c as u128) * (c as u128);
        let d = c as f64 - nn as f64 * delta;
        spread += d * d;
    }
    if e > cap as u128 {
        return None;
    }
    let mut keys: Vec<u64> = Vec::with_capacity(e as usize + 8);
    let mut hits: Vec<u32> = Vec::with_capacity(nn as usize);
    for m in mm..2 * mm {
        hits.clear();
        for l in nn..2 * nn {
            if bits.get((m * l) as usize) {
                hits.push((l - nn) as u32);
            }
        }
        for &a in &hits {
            for &b in &hits {
                keys.push(a as u64 * nn + b as u64);
            }
        }
    }
    keys.sort_unstable();
    let base = mm as f64 * delta * delta;
    let mut sorted = cs.clone();
    sorted.sort_unstable();
    let mut dv: Vec<(f64, f64)> = Vec::new();
    let mut i = 0;
    while i < sorted.len() {
        let mut j = i + 1;
        while j < sorted.len() && sorted[j] == sorted[i] {
            j += 1;
        }
        dv.push((sorted[i] as f64, (j - i) as f64));
        i = j;
    }
    let mut ebal = 0.0f64;
    let mut signed = 0.0f64;
    for &(v1, n1) in &dv {
        for &(v2, n2) in &dv {
            ebal += n1 * n2 * (base - delta * (v1 + v2)).abs();
            signed += n1 * n2 * (base - delta * (v1 + v2));
        }
    }
    signed += e as f64;
    let mut diag = 0.0f64;
    for &c in &cs {
        diag += (c as f64 * (1.0 - 2.0 * delta) + base).abs();
    }
    let mut i = 0;
    while i < keys.len() {
        let mut j = i + 1;
        while j < keys.len() && keys[j] == keys[i] {
            j += 1;
        }
        let t = (j - i) as f64;
        let l1 = (keys[i] / nn) as usize;
        let l2 = (keys[i] % nn) as usize;
        let z = base - delta * (cs[l1] as f64 + cs[l2] as f64);
        ebal += (t + z).abs() - z.abs();
        i = j;
    }
    let cells = mm as f64 * nn as f64;
    Some(BoxRow {
        m: mm,
        n: nn,
        r,
        e,
        ebal,
        diag,
        triv: r as f64 * (1.0 - delta) + (cells - r as f64) * delta,
        bound: (mm as f64 * ebal).sqrt(),
        check: (signed, spread),
    })
}

pub fn boxes(x: u128) -> Vec<(u64, u64)> {
    let mut out = Vec::new();
    let mut mm = 8u64;
    while (mm as u128) * 32 <= x {
        let mut nn = 8u64;
        while (nn as u128) * 4 * mm as u128 <= x {
            if (nn as u128) * 32 * mm as u128 >= x {
                out.push((mm, nn));
            }
            nn <<= 1;
        }
        mm <<= 1;
    }
    out
}

pub struct Sweep {
    pub worst: BoxRow,
    pub best: BoxRow,
    pub l5: Option<BoxRow>,
    pub seen: usize,
    pub live: usize,
}

fn sweep(bits: &Bits, x: u128, delta: f64, cap: usize) -> Option<Sweep> {
    let mut worst: Option<BoxRow> = None;
    let mut best: Option<BoxRow> = None;
    let all = boxes(x);
    let mut live = 0;
    let mut l5: Option<BoxRow> = None;
    let edge = (x as f64).powf(0.4);
    for &(mm, nn) in all.iter() {
        let row = match box_row(bits, mm, nn, delta, cap) {
            Some(r) => r,
            None => continue,
        };
        if row.r == 0 {
            continue;
        }
        live += 1;
        if worst.as_ref().map(|b| row.bound > b.bound).unwrap_or(true) {
            worst = Some(row_copy(&row));
        }
        if mm as f64 >= edge
            && nn as f64 >= edge
            && l5.as_ref().map(|b| row.bound > b.bound).unwrap_or(true)
        {
            l5 = Some(row_copy(&row));
        }
        if best.as_ref().map(|b| row.bound < b.bound).unwrap_or(true) {
            best = Some(row);
        }
    }
    match (worst, best) {
        (Some(w), Some(b)) => Some(Sweep {
            worst: w,
            best: b,
            l5,
            seen: all.len(),
            live,
        }),
        _ => None,
    }
}

fn row_copy(r: &BoxRow) -> BoxRow {
    BoxRow {
        m: r.m,
        n: r.n,
        r: r.r,
        e: r.e,
        ebal: r.ebal,
        diag: r.diag,
        triv: r.triv,
        bound: r.bound,
        check: r.check,
    }
}

// SIGMA

fn sigma_max(bits: &Bits, mm: u64, nn: u64, delta: f64, rounds: usize, rng: &mut Rng) -> f64 {
    let mut worst = 0.0f64;
    let mut asign = vec![0.0f64; mm as usize];
    let mut bsign = vec![0.0f64; nn as usize];
    for _ in 0..rounds {
        for a in asign.iter_mut() {
            *a = rng.sign();
        }
        for b in bsign.iter_mut() {
            *b = rng.sign();
        }
        let mut hit = 0.0f64;
        for m in mm..2 * mm {
            let am = asign[(m - mm) as usize];
            for l in nn..2 * nn {
                if bits.get((m * l) as usize) {
                    hit += am * bsign[(l - nn) as usize];
                }
            }
        }
        let asum: f64 = asign.iter().sum();
        let bsum: f64 = bsign.iter().sum();
        worst = worst.max((hit - delta * asum * bsum).abs());
    }
    worst
}

// STUDY

struct Cell {
    q: u64,
    digits: Vec<u64>,
    label: &'static str,
    lmax: usize,
    lrand: usize,
    typeii: Vec<usize>,
}

fn ex(q: u64, e: u64) -> Vec<u64> {
    (0..q).filter(|&f| f != e).collect()
}

fn cells() -> Vec<Cell> {
    vec![
        Cell {
            q: 3,
            digits: vec![0, 1],
            label: "01",
            lmax: 12,
            lrand: 12,
            typeii: vec![10, 12, 14],
        },
        Cell {
            q: 3,
            digits: vec![0, 2],
            label: "02",
            lmax: 10,
            lrand: 99,
            typeii: vec![],
        },
        Cell {
            q: 3,
            digits: vec![1, 2],
            label: "12",
            lmax: 12,
            lrand: 12,
            typeii: vec![12],
        },
        Cell {
            q: 4,
            digits: vec![0, 1, 2],
            label: "012",
            lmax: 8,
            lrand: 8,
            typeii: vec![9],
        },
        Cell {
            q: 5,
            digits: vec![0, 1, 2, 3],
            label: "0123",
            lmax: 6,
            lrand: 6,
            typeii: vec![8],
        },
        Cell {
            q: 5,
            digits: vec![0, 2, 4],
            label: "024",
            lmax: 7,
            lrand: 7,
            typeii: vec![7],
        },
        Cell {
            q: 10,
            digits: ex(10, 7),
            label: "ex7",
            lmax: 4,
            lrand: 4,
            typeii: vec![6],
        },
        Cell {
            q: 10,
            digits: ex(10, 0),
            label: "ex0",
            lmax: 3,
            lrand: 3,
            typeii: vec![4],
        },
        Cell {
            q: 100,
            digits: (0..50).collect(),
            label: "0to49",
            lmax: 2,
            lrand: 2,
            typeii: vec![3],
        },
        Cell {
            q: 100,
            digits: vec![0, 1],
            label: "01",
            lmax: 9,
            lrand: 9,
            typeii: vec![3],
        },
        Cell {
            q: 100,
            digits: ex(100, 37),
            label: "ex37",
            lmax: 1,
            lrand: 1,
            typeii: vec![2],
        },
    ]
}

pub fn qpow(q: u64, l: usize) -> u128 {
    let mut p = 1u128;
    for _ in 0..l {
        p *= q as u128;
    }
    p
}

pub fn exps(v: f64, x: u128) -> f64 {
    if v <= 0.0 {
        0.0
    } else {
        v.ln() / (x as f64).ln()
    }
}

const CAP: usize = 1 << 24;
const PAIRCAP: usize = 8_000_000;

pub fn run() {
    let cs = cells();
    println!("menergy census");
    println!("| q | F | L | K | E_x | theta_x | 2 alpha | E_x/(2K^2 - K) | max r | shift T | shift share | excess share | E_rand | E_rand/(2K^2 - K) | E_x/E_rand |");
    for c in cs.iter() {
        let alpha = (c.digits.len() as f64).ln() / (c.q as f64).ln();
        for l in 1..=c.lmax {
            let x = qpow(c.q, l);
            assert!(x <= u128::MAX / x, "q^(2L) overflows u128");
            let vals = column(c.q, &c.digits, l);
            let k = vals.len() as u128;
            let (e, top) = energy_of(&vals, CAP);
            assert!(e >= 2 * k * k - k, "diagonal floor fails");
            assert!(e <= k * k * k, "trivial ceiling fails");
            assert!(
                e as f64 <= (top as f64) * (k * k) as f64 + 0.5,
                "divisor ceiling fails"
            );
            let diag = (2 * k * k - k) as f64;
            let shift = if c.digits.contains(&0) {
                shift_excess(c.digits.len() as u128, l)
            } else {
                0
            };
            assert!(e >= 2 * k * k - k + shift, "shift floor fails");
            let (rand, randratio, overrand) = if l >= c.lrand {
                let mut rng = Rng::new(0x5eed + l as u64 + c.q * 977);
                let rv = random_column(x, k as usize, &mut rng);
                let (re, _) = energy_of(&rv, CAP);
                (
                    format!("{re}"),
                    format!("{:.4}", re as f64 / diag),
                    format!("{:.4}", e as f64 / re as f64),
                )
            } else {
                ("-".to_string(), "-".to_string(), "-".to_string())
            };
            let excess = e as f64 - diag;
            let exshare = if excess > 0.0 {
                format!("{:.4}", shift as f64 / excess)
            } else {
                "-".to_string()
            };
            println!(
                "| {} | {} | {} | {} | {} | {:.6} | {:.6} | {:.4} | {} | {} | {:.4} | {} | {} | {} | {} |",
                c.q,
                c.label,
                l,
                k,
                e,
                exps(e as f64, x),
                2.0 * alpha,
                e as f64 / diag,
                top,
                shift,
                (diag + shift as f64) / e as f64,
                exshare,
                rand,
                randratio,
                overrand
            );
        }
    }
    println!("menergy type II");
    println!("| q | F | L | x | alpha | kind | box | boxes | M | N | R | E_x(M,N) | E_bal | diag share | bound | bound/triv | bound exp | alpha - exp |");
    for c in cs.iter() {
        let alpha = (c.digits.len() as f64).ln() / (c.q as f64).ln();
        for &l in c.typeii.iter() {
            let x = qpow(c.q, l);
            let vals = column(c.q, &c.digits, l);
            let delta = vals.len() as f64 / x as f64;
            let mut rng = Rng::new(0xb0a7 + l as u64 + c.q * 131);
            let rv = random_column(x, vals.len(), &mut rng);
            for (kind, set) in [("digit", &vals), ("random", &rv)] {
                let bits = bits_of(set, x);
                let sw = match sweep(&bits, x, delta, PAIRCAP) {
                    Some(s) => s,
                    None => {
                        println!(
                            "| {} | {} | {} | {} | {:.6} | {} | none | 0/{} | - | - | - | - | - | - | - | - | - | - |",
                            c.q, c.label, l, x, alpha, kind, boxes(x).len()
                        );
                        continue;
                    }
                };
                let mut shown: Vec<(&str, &BoxRow)> = vec![("worst", &sw.worst)];
                if let Some(r) = sw.l5.as_ref() {
                    shown.push(("l5", r));
                }
                shown.push(("best", &sw.best));
                for (tag, row) in shown {
                    let expo = exps(row.bound, x);
                    println!(
                        "| {} | {} | {} | {} | {:.6} | {} | {} | {}/{} | {} | {} | {} | {} | {:.4e} | {:.4} | {:.4e} | {:.4} | {:.6} | {:.6} |",
                        c.q,
                        c.label,
                        l,
                        x,
                        alpha,
                        kind,
                        tag,
                        sw.live,
                        sw.seen,
                        row.m,
                        row.n,
                        row.r,
                        row.e,
                        row.ebal,
                        row.diag / row.ebal,
                        row.bound,
                        row.bound / row.triv,
                        expo,
                        alpha - expo
                    );
                }
            }
        }
    }
    println!("menergy sigma");
    println!("| q | F | L | M | N | rounds | max Sigma | bound | ratio | signed check |");
    for (q, digits, label, l, mm, nn) in [
        (3u64, vec![0u64, 1], "01", 6usize, 8u64, 8u64),
        (3, vec![0, 1], "01", 8, 8, 32),
        (3, vec![0, 1], "01", 8, 16, 16),
        (3, vec![1, 2], "12", 8, 16, 16),
        (4, vec![0, 1, 2], "012", 6, 16, 32),
        (5, vec![0, 1, 2, 3], "0123", 5, 16, 32),
    ] {
        let x = qpow(q, l);
        let vals = column(q, &digits, l);
        let delta = vals.len() as f64 / x as f64;
        let bits = bits_of(&vals, x);
        let row = box_row(&bits, mm, nn, delta, PAIRCAP).unwrap();
        let mut rng = Rng::new(0x5169 + q * 7 + l as u64);
        let worst = sigma_max(&bits, mm, nn, delta, 40, &mut rng);
        assert!(
            worst <= row.bound * (1.0 + 1e-9),
            "sigma exceeds the energy bound at q={q} L={l}"
        );
        println!(
            "| {} | {} | {} | {} | {} | {} | {:.4e} | {:.4e} | {:.4} | {:.3e} |",
            q,
            label,
            l,
            mm,
            nn,
            40,
            worst,
            row.bound,
            worst / row.bound,
            (row.check.0 - row.check.1).abs() / row.check.1.max(1.0)
        );
    }
}

// TESTS

#[cfg(test)]
mod tests {
    use super::*;

    fn brute_energy(vals: &[u128]) -> u128 {
        let mut e = 0u128;
        for &a in vals {
            for &b in vals {
                for &c in vals {
                    for &d in vals {
                        if a * b == c * d {
                            e += 1;
                        }
                    }
                }
            }
        }
        e
    }

    #[test]
    fn energy_matches_brute() {
        for (q, digits, l) in [
            (3u64, vec![0u64, 1], 4usize),
            (3, vec![1, 2], 4),
            (4, vec![0, 1, 2], 3),
            (5, vec![0, 1, 2, 3], 2),
            (10, ex(10, 7), 2),
            (100, vec![0, 1], 3),
        ] {
            let vals = column(q, &digits, l);
            let (e, _) = energy_of(&vals, CAP);
            assert_eq!(e, brute_energy(&vals), "energy q={q} L={l}");
        }
    }

    #[test]
    fn partition_invariance() {
        let vals = column(3, &[0, 1], 8);
        let (a, ta) = energy_of(&vals, CAP);
        let (b, tb) = energy_of(&vals, 64);
        let (c, tc) = energy(&vals, 7);
        assert_eq!((a, ta), (b, tb));
        assert_eq!((a, ta), (c, tc));
    }

    #[test]
    fn diagonal_and_ceiling() {
        for (q, digits, l) in [
            (3u64, vec![0u64, 1], 9usize),
            (4, vec![0, 1, 2], 6),
            (10, ex(10, 0), 3),
            (100, (0..50).collect::<Vec<u64>>(), 2),
        ] {
            let vals = column(q, &digits, l);
            let k = vals.len() as u128;
            let (e, top) = energy_of(&vals, CAP);
            assert!(e >= 2 * k * k - k);
            assert!(e <= k * k * k);
            assert!(e <= top as u128 * k * k);
        }
    }

    #[test]
    fn scaling_invariance() {
        for l in 1..=9 {
            let a = energy_of(&column(3, &[0, 1], l), CAP);
            let b = energy_of(&column(3, &[0, 2], l), CAP);
            assert_eq!(a, b, "scaling q=3 L={l}");
        }
        for l in 1..=5 {
            let a = energy_of(&column(5, &[0, 1, 2], l), CAP);
            let b = energy_of(&column(5, &[0, 2, 4], l), CAP);
            assert_eq!(a, b, "scaling q=5 L={l}");
        }
        let a = energy_of(&column(3, &[0, 1], 6), CAP);
        let b = energy_of(&column(3, &[1, 2], 6), CAP);
        assert!(a != b, "translation is not a multiplicative invariance");
    }

    #[test]
    fn restricted_energy_matches_brute() {
        let (q, digits, l, mm, nn) = (3u64, vec![0u64, 1], 8usize, 8u64, 16u64);
        let x = qpow(q, l);
        let vals = column(q, &digits, l);
        let delta = vals.len() as f64 / x as f64;
        let bits = bits_of(&vals, x);
        let row = box_row(&bits, mm, nn, delta, PAIRCAP).unwrap();
        let mut r = 0u128;
        let mut e = 0u128;
        for m in mm..2 * mm {
            for l1 in nn..2 * nn {
                if bits.get((m * l1) as usize) {
                    r += 1;
                }
                for l2 in nn..2 * nn {
                    if bits.get((m * l1) as usize) && bits.get((m * l2) as usize) {
                        e += 1;
                    }
                }
            }
        }
        assert_eq!(row.r, r);
        assert_eq!(row.e, e);
        let mut ebal = 0.0f64;
        for l1 in nn..2 * nn {
            for l2 in nn..2 * nn {
                let mut w = 0.0f64;
                for m in mm..2 * mm {
                    let p1 = if bits.get((m * l1) as usize) {
                        1.0
                    } else {
                        0.0
                    } - delta;
                    let p2 = if bits.get((m * l2) as usize) {
                        1.0
                    } else {
                        0.0
                    } - delta;
                    w += p1 * p2;
                }
                ebal += w.abs();
            }
        }
        assert!(
            (row.ebal - ebal).abs() < 1e-6 * ebal,
            "E_bal {} against {ebal}",
            row.ebal
        );
        assert!((row.check.0 - row.check.1).abs() < 1e-6 * row.check.1);
    }

    #[test]
    fn sigma_stays_under_the_bound() {
        for (q, digits, l, mm, nn) in [
            (3u64, vec![0u64, 1], 8usize, 8u64, 16u64),
            (3, vec![1, 2], 8, 16, 16),
            (4, vec![0, 1, 2], 6, 8, 32),
        ] {
            let x = qpow(q, l);
            let vals = column(q, &digits, l);
            let delta = vals.len() as f64 / x as f64;
            let bits = bits_of(&vals, x);
            let row = box_row(&bits, mm, nn, delta, PAIRCAP).unwrap();
            let mut rng = Rng::new(11 + q + l as u64);
            let worst = sigma_max(&bits, mm, nn, delta, 60, &mut rng);
            assert!(
                worst <= row.bound,
                "sigma {worst} over bound {} at q={q}",
                row.bound
            );
        }
    }

    #[test]
    fn shift_solutions_are_solutions() {
        for (q, digits, l) in [
            (3u64, vec![0u64, 1], 6usize),
            (4, vec![0, 1, 2], 4),
            (5, vec![0, 2, 4], 4),
            (10, ex(10, 7), 3),
        ] {
            let k = digits.len() as u128;
            let vals = column(q, &digits, l);
            let (e, _) = energy_of(&vals, CAP);
            let kk = vals.len() as u128;
            let t = shift_excess(k, l);
            assert!(t > 0);
            assert!(e >= 2 * kk * kk - kk + t, "shift floor q={q} L={l}");
        }
        let mut seen = std::collections::HashSet::new();
        let (q, digits, l) = (3u64, vec![0u64, 1], 5usize);
        let vals = column(q, &digits, l);
        let set: std::collections::HashSet<u128> = vals.iter().copied().collect();
        let mut count = 0u128;
        for s in 0..2 * l {
            for i in 0..=s {
                for ip in 0..=s {
                    if i == ip {
                        continue;
                    }
                    for &u in vals.iter() {
                        for &v in vals.iter() {
                            if u == v || u % q as u128 == 0 || v % q as u128 == 0 {
                                continue;
                            }
                            let p = |e: usize| (q as u128).pow(e as u32);
                            let quad = (p(i) * u, p(s - i) * v, p(ip) * u, p(s - ip) * v);
                            if !set.contains(&quad.0)
                                || !set.contains(&quad.1)
                                || !set.contains(&quad.2)
                                || !set.contains(&quad.3)
                            {
                                continue;
                            }
                            assert_eq!(quad.0 * quad.1, quad.2 * quad.3);
                            if seen.insert(quad) {
                                count += 1;
                            }
                        }
                    }
                }
            }
        }
        assert_eq!(count, shift_excess(2, l), "shift count L={l}");
    }

    #[test]
    fn random_column_is_a_set() {
        let mut rng = Rng::new(4);
        let v = random_column(1000, 400, &mut rng);
        assert_eq!(v.len(), 400);
        assert!(v.windows(2).all(|w| w[0] < w[1]));
        assert!(v[0] >= 1 && *v.last().unwrap() < 1000);
    }
}
