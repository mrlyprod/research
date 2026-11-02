use std::collections::HashSet;

// CATALOG

pub const NEIGHBOURHOODS: [(&str, u32); 3] = [("carpet", 495), ("net", 186), ("void", 341)];
pub const SOURCES: [(&str, &[u32]); 4] = [
    ("primes", &[2, 3, 5, 7]),
    ("binary", &[1, 2, 4, 8]),
    ("fibonacci", &[0, 1, 2, 3, 5, 8]),
    ("odds", &[1, 3, 5, 7]),
];
pub const SPREADS: [&str; 6] = ["id", "shift", "drop", "slide", "turn", "flip"];

pub fn spread_at(name: &str, n: usize, i: usize, j: usize) -> (usize, usize) {
    match name {
        "id" => (i, j),
        "shift" => (i, (j + 1) % n),
        "drop" => ((i + 1) % n, j),
        "slide" => ((i + 1) % n, (j + 1) % n),
        "turn" => (j, n - 1 - i),
        "flip" => (j, i),
        _ => unreachable!(),
    }
}

pub fn offsets(code: u32) -> Vec<(i64, i64)> {
    (0..9)
        .filter(|&b| code >> b & 1 == 1 && b != 4)
        .map(|b| (b as i64 / 3 - 1, b as i64 % 3 - 1))
        .collect()
}

pub fn rule_of(values: &[u32]) -> u32 {
    values
        .iter()
        .filter(|&&v| v < 32)
        .fold(0, |r, &v| r | 1 << v)
}

pub fn word(k1: u32, k2: u32) -> u32 {
    let mut mask = 0;
    for i in 0..4 {
        for j in 0..4 {
            let outer = k1 >> (2 * (i / 2) + j / 2) & 1;
            let inner = k2 >> (2 * (i % 2) + j % 2) & 1;
            mask |= (outer & inner) << (4 * i + j);
        }
    }
    mask
}

pub fn three() -> Vec<u32> {
    (1..511).collect()
}

pub fn four() -> Vec<u32> {
    (1..16)
        .flat_map(|a| (1..16).map(move |b| word(a, b)))
        .filter(|&m| m != 0xffff)
        .collect()
}

// TORUS

pub fn symmetries(n: usize) -> Vec<Vec<usize>> {
    let mats: [[i64; 4]; 8] = [
        [1, 0, 0, 1],
        [0, -1, 1, 0],
        [-1, 0, 0, -1],
        [0, 1, -1, 0],
        [0, 1, 1, 0],
        [0, -1, -1, 0],
        [-1, 0, 0, 1],
        [1, 0, 0, -1],
    ];
    let m = n as i64;
    let mut out = Vec::new();
    for mat in mats {
        for ti in 0..m {
            for tj in 0..m {
                let perm = (0..n * n)
                    .map(|p| {
                        let (i, j) = ((p / n) as i64, (p % n) as i64);
                        let a = (mat[0] * i + mat[1] * j + ti).rem_euclid(m);
                        let b = (mat[2] * i + mat[3] * j + tj).rem_euclid(m);
                        (a * m + b) as usize
                    })
                    .collect();
                out.push(perm);
            }
        }
    }
    out
}

pub fn image(perm: &[usize], mask: u32) -> u32 {
    (0..perm.len())
        .filter(|&p| mask >> p & 1 == 1)
        .fold(0, |m, p| m | 1 << perm[p])
}

fn commutes(a: &[usize], b: &[usize]) -> bool {
    (0..a.len()).all(|p| a[b[p]] == b[a[p]])
}

fn reps(masks: &[u32], group: &[&Vec<usize>]) -> Vec<u32> {
    let mut seen = HashSet::new();
    masks
        .iter()
        .copied()
        .filter(|&m| {
            let canon = group.iter().map(|s| image(s, m)).min().unwrap();
            seen.insert(canon)
        })
        .collect()
}

// ROUND

pub struct Round {
    pub bits: usize,
    pub fill: u32,
    pub reads: Vec<(u32, u32)>,
    pub rule: u32,
    pub perm: Vec<usize>,
    lo: [u32; 256],
    hi: [u32; 256],
}

impl Round {
    pub fn new(n: usize, fill: u32, hood: u32, rule: u32, spread: &str) -> Round {
        let perm: Vec<usize> = (0..n * n)
            .map(|p| {
                let (i, j) = spread_at(spread, n, p / n, p % n);
                i * n + j
            })
            .collect();
        Round::with_perm(n, fill, hood, rule, perm)
    }

    pub fn with_perm(n: usize, fill: u32, hood: u32, rule: u32, perm: Vec<usize>) -> Round {
        let bits = n * n;
        let void = ((1u64 << bits) - 1) as u32 & !fill;
        let offs = offsets(hood);
        let m = n as i64;
        let reads = (0..bits)
            .filter(|&a| fill >> a & 1 == 1)
            .map(|a| {
                let (i, j) = ((a / n) as i64, (a % n) as i64);
                let mask = offs.iter().fold(0u32, |mk, &(di, dj)| {
                    let q = (i + di).rem_euclid(m) * m + (j + dj).rem_euclid(m);
                    mk | 1 << q
                }) & void;
                (a as u32, mask)
            })
            .collect();
        let table = |shift: usize| {
            let mut t = [0u32; 256];
            for (byte, slot) in t.iter_mut().enumerate() {
                for b in 0..8 {
                    let p = shift + b;
                    if p < bits && byte >> b & 1 == 1 {
                        *slot |= 1 << perm[p];
                    }
                }
            }
            t
        };
        Round {
            bits,
            fill,
            rule,
            lo: table(0),
            hi: table(8),
            reads,
            perm,
        }
    }

    pub fn void(&self) -> u32 {
        ((1u64 << self.bits) - 1) as u32 & !self.fill
    }

    pub fn mess(&self, x: u32) -> u32 {
        self.reads.iter().fold(0, |m, &(a, mask)| {
            m | ((self.rule >> (x & mask).count_ones()) & 1) << a
        })
    }

    pub fn spread(&self, x: u32) -> u32 {
        self.lo[(x & 255) as usize] | self.hi[(x >> 8 & 255) as usize]
    }

    pub fn table(&self) -> Vec<u32> {
        (0..1u32 << self.bits)
            .map(|x| self.spread(x ^ self.mess(x)))
            .collect()
    }

    pub fn deps(&self) -> Vec<(u32, u32)> {
        self.reads
            .iter()
            .map(|&(a, mask)| {
                let k = mask.count_ones();
                let live = (0..k).any(|c| (self.rule >> c & 1) != (self.rule >> (c + 1) & 1));
                (a, if live { mask } else { 0 })
            })
            .collect()
    }

    pub fn closed(&self) -> bool {
        let all = ((1u64 << self.bits) - 1) as u32;
        let orbit = |s: usize| {
            let (mut p, mut o) = (s, 0u32);
            loop {
                o |= 1 << p;
                p = self.perm[p];
                if p == s {
                    return o;
                }
            }
        };
        let deps = self.deps();
        (0..self.bits).any(|s| {
            let mut set = orbit(s);
            loop {
                let grown = deps
                    .iter()
                    .filter(|&&(_, mask)| mask & set != 0)
                    .fold(set, |g, &(a, _)| g | orbit(a as usize));
                if grown == set {
                    break;
                }
                set = grown;
            }
            set != all
        })
    }
}

// CELLS

pub struct Spec {
    pub n: usize,
    pub fill: u32,
    pub hood: usize,
    pub source: usize,
    pub spread: usize,
}

impl Spec {
    pub fn round(&self) -> Round {
        Round::new(
            self.n,
            self.fill,
            NEIGHBOURHOODS[self.hood].1,
            rule_of(SOURCES[self.source].1),
            SPREADS[self.spread],
        )
    }

    pub fn linear(&self) -> bool {
        SOURCES[self.source].0 == "odds"
    }

    pub fn name(&self) -> String {
        format!(
            "n {} fill {:#06x} {} {} {}",
            self.n,
            self.fill,
            NEIGHBOURHOODS[self.hood].0,
            SOURCES[self.source].0,
            SPREADS[self.spread]
        )
    }
}

pub fn cells(n: usize, masks: &[u32]) -> Vec<Spec> {
    let syms = symmetries(n);
    let mut out = Vec::new();
    for spread in 0..SPREADS.len() {
        let probe = Round::new(n, 0, 0, 0, SPREADS[spread]);
        let central: Vec<&Vec<usize>> = syms.iter().filter(|s| commutes(s, &probe.perm)).collect();
        for hood in 0..NEIGHBOURHOODS.len() {
            for source in 0..SOURCES.len() {
                for fill in reps(masks, &central) {
                    out.push(Spec {
                        n,
                        fill,
                        hood,
                        source,
                        spread,
                    });
                }
            }
        }
    }
    out
}

pub fn cycles(perm: &[u32]) -> Vec<usize> {
    let mut seen = vec![false; perm.len()];
    let mut out = Vec::new();
    for s in 0..perm.len() {
        if seen[s] {
            continue;
        }
        let (mut x, mut len) = (s, 0);
        while !seen[x] {
            seen[x] = true;
            x = perm[x] as usize;
            len += 1;
        }
        out.push(len);
    }
    out
}
