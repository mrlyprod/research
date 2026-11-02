// MACHINE

pub fn rot(i: usize, level: usize) -> usize {
    (i >> 1) | ((i & 1) << (level - 1))
}

pub fn perm(level: usize, name: &str, t: usize) -> Vec<usize> {
    let n = 1usize << level;
    (0..n * n)
        .map(|p| {
            let (i, j) = (p / n, p % n);
            let (a, b) = match name {
                "slide" => (i + 1, j + 1),
                "shift" => (i, j + 1),
                "scales" => (i + (1 << (t % level)), j + (1 << (t % level))),
                "rotate" => (rot(i, level) + 1, rot(j, level) + 1),
                "hyper" => (2 * i + j + 1, 3 * i + 2 * j),
                "shear" => (i + j + 1, j + 1),
                "row" => (rot(i, level) + 1, j + 1),
                "cat" => (2 * i + j + 1, i + j + 1),
                _ => unreachable!(),
            };
            (a % n) * n + b % n
        })
        .collect()
}

pub struct Grid {
    pub rule: u32,
    pub level: usize,
    pub n: usize,
    pub cells: usize,
    pub fill: Vec<bool>,
    pub reads: Vec<(usize, Vec<usize>)>,
}

impl Grid {
    pub fn new(level: usize) -> Grid {
        Grid::split(level, "carpet", 0b1010_1100)
    }

    pub fn split(level: usize, design: &str, rule: u32) -> Grid {
        let n = 1usize << level;
        let fill: Vec<bool> = (0..n * n)
            .map(|p| match design {
                "carpet" => (p / n) & (p % n) == 0,
                "board" => (p / n + p % n) % 2 == 0,
                _ => unreachable!(),
            })
            .collect();
        let reads = (0..n * n)
            .filter(|&a| fill[a])
            .map(|a| {
                let (i, j) = (a / n, a % n);
                let mut r = Vec::new();
                for di in [n - 1, 0, 1] {
                    for dj in [n - 1, 0, 1] {
                        let q = ((i + di) % n) * n + (j + dj) % n;
                        if q != a && !fill[q] {
                            r.push(q);
                        }
                    }
                }
                (a, r)
            })
            .collect();
        Grid {
            rule,
            level,
            n,
            cells: n * n,
            fill,
            reads,
        }
    }

    pub fn mask(&self) -> u32 {
        (0..self.cells)
            .filter(|&p| self.fill[p])
            .fold(0, |m, p| m | 1 << p)
    }

    pub fn hidden(&self, schedule: &str) -> usize {
        let p = perm(self.level, schedule, 0);
        let mut seen = vec![false; self.cells];
        let mut count = 0;
        for s in 0..self.cells {
            if seen[s] {
                continue;
            }
            let (mut q, mut inside) = (s, true);
            while !seen[q] {
                seen[q] = true;
                inside &= self.fill[q];
                q = p[q];
            }
            count += inside as usize;
        }
        count
    }

    fn perms(&self, schedule: &str) -> Vec<Vec<usize>> {
        let period = if schedule == "scales" { self.level } else { 1 };
        (0..period).map(|t| perm(self.level, schedule, t)).collect()
    }

    pub fn stuck(&self, schedule: &str) -> usize {
        let perms = self.perms(schedule);
        let mut readers = vec![Vec::new(); self.cells];
        for (a, r) in &self.reads {
            for &q in r {
                readers[q].push(*a);
            }
        }
        (0..self.cells)
            .filter(|&start| {
                let mut set = vec![false; self.cells];
                let mut stack = vec![start];
                set[start] = true;
                let mut size = 1;
                while let Some(q) = stack.pop() {
                    for x in perms.iter().map(|p| p[q]).chain(readers[q].iter().copied()) {
                        if !set[x] {
                            set[x] = true;
                            size += 1;
                            stack.push(x);
                        }
                    }
                }
                size < self.cells
            })
            .count()
    }

    pub fn cover(&self, schedule: &str, most: usize) -> Vec<Option<usize>> {
        let perms = self.perms(schedule);
        let mut out = vec![None; self.cells];
        for base in (0..self.cells).step_by(64) {
            let lanes = (self.cells - base).min(64);
            let all = if lanes == 64 {
                u64::MAX
            } else {
                (1u64 << lanes) - 1
            };
            let mut set = vec![0u64; self.cells];
            for k in 0..lanes {
                set[base + k] |= 1 << k;
            }
            let mut done = 0u64;
            for t in 0..most {
                let full = set.iter().fold(all, |m, &x| m & x);
                for k in 0..lanes {
                    if (full & !done) >> k & 1 == 1 {
                        out[base + k] = Some(t);
                    }
                }
                done |= full;
                if done == all {
                    break;
                }
                let mut grown = set.clone();
                for (a, r) in &self.reads {
                    grown[*a] |= r.iter().fold(0, |m, &q| m | set[q]);
                }
                let p = &perms[t % perms.len()];
                for q in 0..self.cells {
                    set[p[q]] = grown[q];
                }
            }
        }
        out
    }
}

// SAMPLER

type Bits = [u64; 4];

struct Fast {
    rule: u32,
    reads: Vec<(usize, Bits)>,
    perms: Vec<Vec<usize>>,
}

fn bit(x: &Bits, p: usize) -> bool {
    x[p >> 6] >> (p & 63) & 1 == 1
}

fn mess(f: &Fast, y: &Bits) -> Bits {
    let mut out = [0u64; 4];
    for (a, m) in &f.reads {
        let c = (0..4).map(|w| (y[w] & m[w]).count_ones()).sum::<u32>();
        if f.rule >> c & 1 == 1 {
            out[a >> 6] |= 1 << (a & 63);
        }
    }
    out
}

pub struct Rng(pub u64);

impl Rng {
    pub fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

pub struct Sampled {
    pub bias: Vec<f64>,
    pub collisions: Vec<u64>,
}

pub fn sample(
    grid: &Grid,
    schedule: &str,
    start: usize,
    runs: usize,
    rounds: usize,
    seed: u64,
) -> Sampled {
    let cells = grid.cells;
    let period = if schedule == "scales" { grid.level } else { 1 };
    let f = Fast {
        rule: grid.rule,
        reads: grid
            .reads
            .iter()
            .map(|(a, r)| {
                let mut m = [0u64; 4];
                for &q in r {
                    m[q >> 6] |= 1 << (q & 63);
                }
                (*a, m)
            })
            .collect(),
        perms: (0..period).map(|t| perm(grid.level, schedule, t)).collect(),
    };
    let words = cells.div_ceil(64);
    let top = if cells % 64 == 0 {
        u64::MAX
    } else {
        (1u64 << (cells % 64)) - 1
    };
    let mut rng = Rng(seed);
    let mut state: Vec<Bits> = vec![
        {
            let mut d = [0u64; 4];
            d[start >> 6] |= 1 << (start & 63);
            d
        };
        runs
    ];
    let full = (1u128 << cells.min(127)) as f64;
    let half = if cells >= 64 {
        0.5
    } else {
        (full / 2.0) / (full - 1.0)
    };
    let mut bias = vec![1.0 - half];
    let mut collisions = vec![(runs * (runs - 1) / 2) as u64];
    for t in 0..rounds - 1 {
        let p = &f.perms[t % period];
        let mut ones = vec![0u32; cells];
        for d in state.iter_mut() {
            let mut y = [0u64; 4];
            for (w, slot) in y.iter_mut().enumerate().take(words) {
                *slot = rng.next() & if w + 1 == words { top } else { u64::MAX };
            }
            let mut z = y;
            for w in 0..4 {
                z[w] ^= d[w];
            }
            let (m1, m2) = (mess(&f, &y), mess(&f, &z));
            let mut pre = *d;
            for w in 0..4 {
                pre[w] ^= m1[w] ^ m2[w];
            }
            let mut next = [0u64; 4];
            for q in 0..cells {
                if bit(&pre, q) {
                    next[p[q] >> 6] |= 1 << (p[q] & 63);
                }
            }
            for (q, o) in ones.iter_mut().enumerate() {
                *o += bit(&next, q) as u32;
            }
            *d = next;
        }
        bias.push(
            ones.iter()
                .map(|&o| (o as f64 / runs as f64 - half).abs())
                .fold(0.0, f64::max),
        );
        let mut sorted = state.clone();
        sorted.sort();
        let mut coll = 0u64;
        let mut run = 1u64;
        for i in 1..sorted.len() {
            if sorted[i] == sorted[i - 1] {
                run += 1;
            } else {
                coll += run * (run - 1) / 2;
                run = 1;
            }
        }
        coll += run * (run - 1) / 2;
        collisions.push(coll);
    }
    Sampled { bias, collisions }
}
