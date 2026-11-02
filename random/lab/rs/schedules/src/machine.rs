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

pub fn rule_of(values: &[u32]) -> u64 {
    values
        .iter()
        .filter(|&&v| v < 64)
        .fold(0, |r, &v| r | 1 << v)
}

pub fn word(k1: u32, k2: u32) -> u128 {
    let mut mask = 0u128;
    for i in 0..4 {
        for j in 0..4 {
            let outer = k1 >> (2 * (i / 2) + j / 2) & 1;
            let inner = k2 >> (2 * (i % 2) + j % 2) & 1;
            mask |= ((outer & inner) as u128) << (4 * i + j);
        }
    }
    mask
}

pub fn full(cells: usize) -> u128 {
    if cells == 128 {
        u128::MAX
    } else {
        (1u128 << cells) - 1
    }
}

pub fn ones(mask: u128) -> Vec<usize> {
    (0..128).filter(|&p| mask >> p & 1 == 1).collect()
}

// SYMMETRIC RULES

pub fn values(rule: u64, m: usize) -> Vec<u8> {
    (0..=m).map(|t| (rule >> t & 1) as u8).collect()
}

pub fn delta(g: &[u8], k: usize) -> Vec<u8> {
    let mut out = g.to_vec();
    for _ in 0..k {
        out = out.windows(2).map(|w| w[0] ^ w[1]).collect();
    }
    out
}

pub fn constant(g: &[u8]) -> bool {
    g.iter().all(|&b| b == g[0])
}

pub fn affine(g: &[u8]) -> bool {
    g.windows(3).all(|w| w[0] == w[2])
}

pub fn selfcomplementary(g: &[u8]) -> bool {
    let m = g.len() - 1;
    (0..=m).all(|t| g[t] == g[m - t])
}

pub fn good(g: &[u8]) -> bool {
    g.len() >= 3 && !affine(g) && !selfcomplementary(g)
}

// SPLITS

#[derive(Clone)]
pub struct Split {
    pub n: usize,
    pub rule: u64,
    pub perm: Vec<usize>,
    pub reads: Vec<(usize, u128)>,
}

impl Split {
    pub fn new(n: usize, fill: u128, hood: u32, rule: u64, spread: &str) -> Split {
        let cells = n * n;
        let void = full(cells) & !fill;
        let offs = offsets(hood);
        let m = n as i64;
        let reads = ones(fill)
            .into_iter()
            .map(|a| {
                let (i, j) = ((a / n) as i64, (a % n) as i64);
                let mask = offs.iter().fold(0u128, |mk, &(di, dj)| {
                    let q = (i + di).rem_euclid(m) * m + (j + dj).rem_euclid(m);
                    mk | 1u128 << q
                }) & void;
                (a, mask)
            })
            .collect();
        let perm = (0..cells)
            .map(|p| {
                let (i, j) = spread_at(spread, n, p / n, p % n);
                i * n + j
            })
            .collect();
        Split {
            n,
            rule,
            perm,
            reads,
        }
    }

    pub fn cells(&self) -> usize {
        self.n * self.n
    }

    pub fn window(&self, mask: u128) -> Vec<u8> {
        values(self.rule, mask.count_ones() as usize)
    }

    pub fn live(&self, mask: u128) -> bool {
        !constant(&self.window(mask))
    }

    pub fn live_reads(&self) -> Vec<(usize, u128)> {
        self.reads
            .iter()
            .copied()
            .filter(|&(_, mask)| self.live(mask))
            .collect()
    }
}

// MACHINE

pub struct Machine {
    pub bits: usize,
    reads: Vec<(u32, u32)>,
    rule: u64,
    tables: Vec<[u32; 256]>,
}

impl Machine {
    pub fn new(s: &Split) -> Machine {
        let bits = s.cells();
        assert!(bits <= 32);
        let tables = (0..4)
            .map(|k| {
                let mut t = [0u32; 256];
                for (byte, slot) in t.iter_mut().enumerate() {
                    for b in 0..8 {
                        let p = 8 * k + b;
                        if p < bits && byte >> b & 1 == 1 {
                            *slot |= 1 << s.perm[p];
                        }
                    }
                }
                t
            })
            .collect();
        Machine {
            bits,
            reads: s.reads.iter().map(|&(a, m)| (a as u32, m as u32)).collect(),
            rule: s.rule,
            tables,
        }
    }

    pub fn mess(&self, x: u32) -> u32 {
        self.reads.iter().fold(0, |m, &(a, mask)| {
            m | (((self.rule >> (x & mask).count_ones()) & 1) as u32) << a
        })
    }

    pub fn spread(&self, x: u32) -> u32 {
        self.tables[0][(x & 255) as usize]
            | self.tables[1][(x >> 8 & 255) as usize]
            | self.tables[2][(x >> 16 & 255) as usize]
            | self.tables[3][(x >> 24 & 255) as usize]
    }

    pub fn apply(&self, x: u32) -> u32 {
        self.spread(x ^ self.mess(x))
    }

    pub fn table(&self) -> Vec<u32> {
        (0..1u32 << self.bits).map(|x| self.apply(x)).collect()
    }
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

pub fn image(perm: &[usize], mask: u128) -> u128 {
    ones(mask).into_iter().fold(0, |m, p| m | 1u128 << perm[p])
}

pub fn commutes(a: &[usize], b: &[usize]) -> bool {
    (0..a.len()).all(|p| a[b[p]] == b[a[p]])
}

pub fn reps(masks: &[u128], group: &[&Vec<usize>], all: Option<u128>) -> Vec<u128> {
    let mut seen = HashSet::new();
    masks
        .iter()
        .copied()
        .filter(|&m| {
            let canon = group
                .iter()
                .flat_map(|s| [image(s, m), all.map_or(u128::MAX, |f| image(s, f & !m))])
                .min()
                .unwrap();
            seen.insert(canon)
        })
        .collect()
}

pub fn orbits(perm: &[usize]) -> Vec<u128> {
    let mut seen = 0u128;
    let mut out = Vec::new();
    for s in 0..perm.len() {
        if seen >> s & 1 == 1 {
            continue;
        }
        let (mut p, mut orbit) = (s, 0u128);
        while orbit >> p & 1 == 0 {
            orbit |= 1u128 << p;
            p = perm[p];
        }
        seen |= orbit;
        out.push(orbit);
    }
    out
}
