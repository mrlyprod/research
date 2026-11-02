use std::collections::{BTreeMap, HashSet, VecDeque};

#[derive(Clone, Copy, PartialEq)]
pub enum Seed {
    Cell,
    Tm,
}

impl Seed {
    pub fn name(self) -> &'static str {
        match self {
            Seed::Cell => "one cell",
            Seed::Tm => "thue-morse",
        }
    }

    fn bit(self, n: i64) -> bool {
        match self {
            Seed::Cell => n == 0,
            Seed::Tm => n >= 0 && (n as u64).count_ones() % 2 == 1,
        }
    }
}

// RULES

fn aff_rule(i: u8) -> u8 {
    let (a, b, d, e) = (i >> 3 & 1, i >> 2 & 1, i >> 1 & 1, i & 1);
    (0..8u8).fold(0, |acc, p| {
        let v = (a & (p >> 2)) ^ (b & (p >> 1 & 1)) ^ (d & (p & 1)) ^ e;
        acc | (v & 1) << p
    })
}

fn anf(rule: u8) -> [u64; 8] {
    let mut f = [0u8; 8];
    for (p, slot) in f.iter_mut().enumerate() {
        *slot = rule >> p & 1;
    }
    for bit in [1usize, 2, 4] {
        for s in 0..8 {
            if s & bit != 0 {
                f[s] ^= f[s ^ bit];
            }
        }
    }
    f.map(|x| if x == 1 { !0u64 } else { 0 })
}

fn apply(k: &[u64; 8], l: u64, c: u64, r: u64) -> u64 {
    k[0] ^ (k[1] & r)
        ^ (k[2] & c)
        ^ (k[3] & c & r)
        ^ (k[4] & l)
        ^ (k[5] & l & r)
        ^ (k[6] & l & c)
        ^ (k[7] & l & c & r)
}

fn window_name(mask: u8) -> String {
    let names: Vec<String> = (0..8u8)
        .filter(|p| mask >> p & 1 == 1)
        .map(|p| format!("{}{}{}", p >> 2 & 1, p >> 1 & 1, p & 1))
        .collect();
    format!("{{{}}}", names.join(","))
}

fn fits(rule: u8, g: u8, mask: u8) -> bool {
    (rule ^ g) & mask == 0
}

// ORBIT

pub struct Orbit {
    pub masks: Vec<u8>,
    pub first: [Option<(usize, i64)>; 8],
    pub hash: (u64, u64),
    pub counts: [u64; 8],
}

#[derive(Clone, Copy)]
pub enum Region {
    All,
    Bulk(i64),
    Interior,
}

fn range_mask(i: usize, lo: i64, hi: i64) -> u64 {
    let base = 64 * i as i64;
    let a = (lo - base).clamp(0, 64);
    let b = (hi - base + 1).clamp(0, 64);
    if b <= a {
        return 0;
    }
    let width = (b - a) as u32;
    let ones = if width == 64 {
        !0u64
    } else {
        (1u64 << width) - 1
    };
    ones << a
}

pub fn orbit(rule: u8, seed: Seed, t_max: usize, region: Region, tally_from: usize) -> Orbit {
    let off = t_max + 8;
    let right = match seed {
        Seed::Cell => t_max + 8,
        Seed::Tm => 4 * t_max + 64,
    };
    let nw = (off + right + 1).div_ceil(64);
    let pend = nw * 64 - 1;
    let mut cur = vec![0u64; nw];
    let mut next = vec![0u64; nw];
    for p in 0..=pend {
        if seed.bit(p as i64 - off as i64) {
            cur[p / 64] |= 1 << (p % 64);
        }
    }
    let k = anf(rule);
    let mut bg = 0u64;
    let mut masks = Vec::with_capacity(t_max + 1);
    let mut first: [Option<(usize, i64)>; 8] = [None; 8];
    let mut counts = [0u64; 8];
    let (mut h1, mut h2) = (0xcbf29ce484222325u64, 0x9e3779b97f4a7c15u64);
    for t in 0..=t_max {
        let imax = match seed {
            Seed::Cell => nw - 1,
            Seed::Tm => (pend - t) / 64 - 1,
        };
        let (lo_p, hi_p) = match region {
            Region::All => (0i64, pend as i64),
            Region::Bulk(n_min) => (off as i64 + n_min, pend as i64),
            Region::Interior => (off as i64 - t as i64 + 2, off as i64 + t as i64 - 2),
        };
        let mut found = 0u8;
        for i in 0..=imax {
            let c = cur[i];
            let prev = if i == 0 { bg } else { cur[i - 1] };
            let nxt = match (i + 1 == nw, seed) {
                (false, _) => cur[i + 1],
                (true, Seed::Cell) => bg,
                (true, Seed::Tm) => 0,
            };
            let l = (c << 1) | (prev >> 63);
            let r = (c >> 1) | (nxt << 63);
            let keep = range_mask(i, lo_p, hi_p);
            if t >= tally_from && keep != 0 {
                for pat in 0..8u8 {
                    let m = (if pat & 4 != 0 { l } else { !l })
                        & (if pat & 2 != 0 { c } else { !c })
                        & (if pat & 1 != 0 { r } else { !r })
                        & keep;
                    counts[pat as usize] += u64::from(m.count_ones());
                    found |= u8::from(m != 0) << pat;
                }
            } else if found != 0xff && keep != 0 {
                for pat in 0..8u8 {
                    if found >> pat & 1 == 1 {
                        continue;
                    }
                    let m = (if pat & 4 != 0 { l } else { !l })
                        & (if pat & 2 != 0 { c } else { !c })
                        & (if pat & 1 != 0 { r } else { !r })
                        & keep;
                    if m != 0 {
                        found |= 1 << pat;
                        if first[pat as usize].is_none() {
                            let n = 64 * i as i64 + m.trailing_zeros() as i64 - off as i64;
                            first[pat as usize] = Some((t, n));
                        }
                    }
                }
            }
            h1 = (h1 ^ c).wrapping_mul(0x100000001b3);
            h2 = h2.rotate_left(7) ^ c.wrapping_mul(0x9e3779b97f4a7c15);
        }
        masks.push(found);
        if t == t_max {
            break;
        }
        for i in 0..nw {
            let c = cur[i];
            let prev = if i == 0 { bg } else { cur[i - 1] };
            let nxt = match (i + 1 == nw, seed) {
                (false, _) => cur[i + 1],
                (true, Seed::Cell) => bg,
                (true, Seed::Tm) => 0,
            };
            let l = (c << 1) | (prev >> 63);
            let r = (c >> 1) | (nxt << 63);
            next[i] = apply(&k, l, c, r);
        }
        std::mem::swap(&mut cur, &mut next);
        bg = if rule >> (if bg == 0 { 0 } else { 7 }) & 1 == 1 {
            !0
        } else {
            0
        };
    }
    Orbit {
        masks,
        first,
        hash: (h1, h2),
        counts,
    }
}

fn union(o: &Orbit) -> u8 {
    o.masks.iter().fold(0, |a, &m| a | m)
}

// EXACT

type Vec15 = u16;
type Comp = (Vec15, i32);
type State = (Vec<Comp>, u8);

fn slot(fam: usize, i: i32) -> u32 {
    (fam as i32 * 5 + i + 2) as u32
}

fn eval(v: Vec15, n: i64) -> u8 {
    let mut s = 0u8;
    for fam in 0..3 {
        for i in -2..=2i32 {
            if v >> slot(fam, i) & 1 == 1 {
                let m = n - i as i64;
                s ^= match fam {
                    0 => u8::from(m == 0),
                    1 => u8::from(m >= 0 && (m as u64).count_ones() % 2 == 1),
                    _ => u8::from(m >= 0),
                };
            }
        }
    }
    s
}

fn flip(v: &mut Vec15, fam: usize, i: i32) {
    assert!((-2..=2).contains(&i));
    *v ^= 1 << slot(fam, i);
}

fn phi(v: Vec15, a: u8, b: u8, d: u8) -> Vec15 {
    let mut out = 0;
    for fam in 0..3 {
        for i in -2..=2i32 {
            if v >> slot(fam, i) & 1 == 1 {
                if a == 1 {
                    flip(&mut out, fam, i + 1);
                }
                if b == 1 {
                    flip(&mut out, fam, i);
                }
                if d == 1 {
                    flip(&mut out, fam, i - 1);
                }
            }
        }
    }
    out
}

fn section(v: Vec15, f: i32) -> Vec15 {
    let mut out = 0;
    for fam in 0..3 {
        for i in -2..=2i32 {
            if v >> slot(fam, i) & 1 == 1 {
                let s = f - i;
                let par = s.rem_euclid(2);
                let shift = -s.div_euclid(2);
                match (fam, par) {
                    (0, 1) => {}
                    (1, 1) => {
                        flip(&mut out, 1, shift);
                        flip(&mut out, 2, shift);
                    }
                    (fam, _) => flip(&mut out, fam, shift),
                }
            }
        }
    }
    out
}

fn step(s: &[Comp], g: (u8, u8, u8), e: i32, f: i32) -> Vec<Comp> {
    s.iter()
        .map(|&(v, o)| {
            let w = if e == 1 { phi(v, g.0, g.1, g.2) } else { v };
            let t = f + o;
            (section(w, t.rem_euclid(2)), t.div_euclid(2))
        })
        .collect()
}

fn pattern(s: &[Comp], m: i64, c: u8) -> u16 {
    s.iter().enumerate().fold(0u16, |acc, (j, &(v, o))| {
        acc | u16::from(eval(v, m + o as i64) ^ c) << j
    })
}

fn base(s: &[Comp], c: u8) -> u8 {
    let mut mask = 0u8;
    for m in [0i64, -1] {
        let w = pattern(s, m, c);
        mask |= 1 << (4 * (w & 1) + 2 * (w >> 1 & 1) + (w >> 2 & 1));
    }
    mask
}

pub fn exact_mask(i: u8, seed: Seed) -> u8 {
    let (a, b, d, e) = (i >> 3 & 1, i >> 2 & 1, i >> 1 & 1, i & 1);
    let sigma = a ^ b ^ d;
    let u: Vec15 = match seed {
        Seed::Cell => 1 << slot(0, 0),
        Seed::Tm => 1 << slot(1, 0),
    };
    let init: Vec<Comp> = vec![(u, -1), (u, 0), (u, 1)];
    let mut mask = base(&init, 0);
    let mut seen: HashSet<State> = HashSet::new();
    let mut queue: VecDeque<State> = VecDeque::new();
    for e0 in 0..2 {
        for f0 in 0..2 {
            let q = if e == 1 { e0 as u8 } else { 0 };
            let s = (step(&init, (a, b, d), e0, f0), q);
            if seen.insert(s.clone()) {
                queue.push_back(s);
            }
        }
    }
    while let Some((s, q)) = queue.pop_front() {
        mask |= base(&s, if e == 1 { q } else { 0 });
        for e1 in 0..2 {
            for f1 in 0..2 {
                let q1 = if e == 0 {
                    0
                } else if sigma == 1 {
                    q
                } else {
                    q | e1 as u8
                };
                let n = (step(&s, (a, b, d), e1, f1), q1);
                if seen.insert(n.clone()) {
                    queue.push_back(n);
                }
            }
        }
    }
    mask
}

// SHADOW

pub fn shadow() {
    let horizon = 4096usize;
    println!("shadow: exact window sets of the 16 affine orbits over all t and all n, against the simulation t <= {horizon}");
    let mut exact = [[0u8; 16]; 2];
    let mut bad = 0;
    for (si, seed) in [Seed::Cell, Seed::Tm].into_iter().enumerate() {
        for i in 0..16u8 {
            let g = aff_rule(i);
            let ex = exact_mask(i, seed);
            let sim = union(&orbit(g, seed, horizon, Region::All, usize::MAX));
            exact[si][i as usize] = ex;
            let verdict = if ex == sim {
                "equal"
            } else {
                bad += 1;
                "DIFFERENT"
            };
            println!(
                "  {:<10} g = rule {:<3} exact {} sim {} {}",
                seed.name(),
                g,
                window_name(ex),
                window_name(sim),
                verdict
            );
        }
    }
    println!("  exact against simulation: {bad} differences in 32 pairs");
    for (si, seed) in [Seed::Cell, Seed::Tm].into_iter().enumerate() {
        census(seed, &exact[si], horizon);
    }
}

fn census(seed: Seed, exact: &[u8; 16], horizon: usize) {
    println!("census, seed {}, rules 0..255", seed.name());
    let mut shadowed: BTreeMap<u8, Vec<u8>> = BTreeMap::new();
    for f in 0..=255u8 {
        let gs: Vec<u8> = (0..16u8)
            .filter(|&i| fits(f, aff_rule(i), exact[i as usize]))
            .map(aff_rule)
            .collect();
        if !gs.is_empty() {
            shadowed.insert(f, gs);
        }
    }
    println!("  linear on the reach (exact): {} rules", shadowed.len());
    let mut by_set: BTreeMap<Vec<u8>, Vec<u8>> = BTreeMap::new();
    for (f, gs) in &shadowed {
        by_set.entry(gs.clone()).or_default().push(*f);
    }
    for (gs, fs) in &by_set {
        println!("    shadows {:?}: rules {:?}", gs, fs);
    }
    let mut late = Vec::new();
    let mut sim_linear = 0;
    let mut classes: BTreeMap<(u64, u64), Vec<u8>> = BTreeMap::new();
    let mut last_new = 0usize;
    for f in 0..=255u8 {
        let o = orbit(f, seed, horizon, Region::All, usize::MAX);
        for pat in 0..8 {
            if let Some((t, _)) = o.first[pat] {
                last_new = last_new.max(t);
            }
        }
        let w = union(&o);
        let sim = (0..16u8).any(|i| fits(f, aff_rule(i), w));
        sim_linear += usize::from(sim);
        if sim != shadowed.contains_key(&f) {
            late.push(f);
        }
        classes.entry(o.hash).or_default().push(f);
    }
    println!(
        "  linear on the reach by simulation to t <= {horizon}: {sim_linear}; disagreeing with the exact set: {:?}; the last new window of any rule appears at t = {last_new}",
        late
    );
    for fs in classes.values().filter(|fs| fs.len() > 1) {
        println!("  identical diagrams to t <= {horizon}: rules {:?}", fs);
    }
}

// COLUMNS

pub fn columns(rule: u8, seed: Seed, steps: usize, centre: i64, probes: &[i64]) -> Vec<Vec<u8>> {
    let half = steps + 8;
    let pmax = 2 * half;
    let nw = pmax / 64 + 3;
    let mut cur = vec![0u64; nw];
    let mut next = vec![0u64; nw];
    for p in 0..=pmax {
        if seed.bit(p as i64 - half as i64 + centre) {
            cur[p / 64 + 1] |= 1 << (p % 64);
        }
    }
    let k = anf(rule);
    let mut cols: Vec<Vec<u8>> = probes.iter().map(|_| Vec::with_capacity(steps)).collect();
    for t in 0..steps {
        for (j, &n) in probes.iter().enumerate() {
            let p = (half as i64 + n) as usize;
            cols[j].push((cur[p / 64 + 1] >> (p % 64) & 1) as u8);
        }
        if t + 1 == steps {
            break;
        }
        let lo = 1 + (t + 1) / 64;
        let hi = 1 + (pmax - t - 1) / 64;
        for i in lo..=hi {
            let c = cur[i];
            let l = (c << 1) | (cur[i - 1] >> 63);
            let r = (c >> 1) | (cur[i + 1] << 63);
            next[i] = apply(&k, l, c, r);
        }
        std::mem::swap(&mut cur, &mut next);
    }
    cols
}

fn period(col: &[u8]) -> Option<(usize, usize)> {
    let n = col.len();
    for p in 1..=512usize {
        let last = (0..n - p).rev().find(|&t| col[t] != col[t + p]);
        let pre = last.map_or(0, |t| t + 1);
        if pre <= n / 8 {
            return Some((p, pre));
        }
    }
    None
}

fn last_exercise(f: u8, o: &Orbit) -> i64 {
    (0..16u8)
        .map(|i| {
            let g = aff_rule(i);
            o.masks
                .iter()
                .enumerate()
                .filter(|&(_, &m)| !fits(f, g, m))
                .map(|(t, _)| t as i64)
                .max()
                .unwrap_or(-1)
        })
        .min()
        .unwrap()
}

fn density(f: u8, o: &Orbit) -> f64 {
    let total: u64 = o.counts.iter().sum();
    let best = (0..16u8)
        .map(|i| {
            let g = aff_rule(i);
            (0..8usize)
                .filter(|&p| (f ^ g) >> p & 1 == 1)
                .map(|p| o.counts[p])
                .sum::<u64>()
        })
        .min()
        .unwrap();
    best as f64 / total.max(1) as f64
}

fn rounds_exercised(f: u8, o: &Orbit) -> usize {
    o.masks
        .iter()
        .filter(|&&m| !(0..16u8).any(|i| fits(f, aff_rule(i), m)))
        .count()
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum Kind {
    Periodic,
    Automatic,
    Grows,
}

fn plateau(col: &[u8], ell: usize) -> (usize, bool) {
    let counts = super::kernel_counts(col, 2, ell);
    let last = counts[counts.len() - 1];
    (last, last == counts[counts.len() - 4])
}

fn kind(col: &[u8]) -> (Kind, usize, bool) {
    let (last, stops) = plateau(col, 24);
    let kind = match (stops, period(col)) {
        (true, Some(_)) => Kind::Periodic,
        (true, None) => Kind::Automatic,
        (false, _) => Kind::Grows,
    };
    (kind, last, plateau(col, 48).1)
}

struct Probe {
    label: &'static str,
    seed: Seed,
    steps: usize,
    centre: i64,
    offsets: [i64; 3],
    region: Region,
    inner: Option<Region>,
}

pub fn escape() {
    let horizon = 4096usize;
    println!("escape: three columns per rule; kernel words of length 24 at base 2; a column stops when the last depth adds nothing to three depths before; periodic means period <= 512 after a preperiod of at most an eighth of the column");
    for steps in [1usize << 15, 1 << 14] {
        let mut state = 1u64;
        let ctrl: Vec<u8> = (0..steps)
            .map(|_| (super::splitmix(&mut state) >> 63) as u8)
            .collect();
        println!(
            "  control splitmix64 seed 1, t < {steps}: {:?}",
            super::kernel_counts(&ctrl, 2, 24)
        );
    }
    let probes = [
        Probe {
            label: "one cell, columns -1, 0, 1, t < 2^15",
            seed: Seed::Cell,
            steps: 1 << 15,
            centre: 0,
            offsets: [-1, 0, 1],
            region: Region::All,
            inner: Some(Region::Interior),
        },
        Probe {
            label: "thue-morse, columns -1, 0, 1 at the seed edge, t < 2^15",
            seed: Seed::Tm,
            steps: 1 << 15,
            centre: 0,
            offsets: [-1, 0, 1],
            region: Region::All,
            inner: None,
        },
        Probe {
            label: "thue-morse, columns 2^14, +1, +2 in the bulk, t < 2^14",
            seed: Seed::Tm,
            steps: 1 << 14,
            centre: 1 << 14,
            offsets: [0, 1, 2],
            region: Region::Bulk(horizon as i64 + 72),
            inner: None,
        },
    ];
    for pr in &probes {
        println!("{}", pr.label);
        let mut table: BTreeMap<(&str, Kind), usize> = BTreeMap::new();
        let mut lrs_max = 0usize;
        let mut lrs_n = 0usize;
        let mut bad: Vec<String> = Vec::new();
        let mut aperiodic: Vec<String> = Vec::new();
        let mut reference = Vec::new();
        let mut lists: BTreeMap<(&str, Kind), Vec<u8>> = BTreeMap::new();
        let mut rhos: BTreeMap<(&str, Kind), Vec<f64>> = BTreeMap::new();
        let mut twice: BTreeMap<&str, usize> = BTreeMap::new();
        let mut burn_max = -1i64;
        let mut past_min = i64::MAX;
        let mut tops: BTreeMap<(&str, Kind), Vec<usize>> = BTreeMap::new();
        for f in 0..=255u8 {
            let o = orbit(f, pr.seed, horizon, pr.region, horizon / 2);
            let sim_fit =
                (0..16u8).any(|i| fits(f, aff_rule(i), o.masks.iter().fold(0, |a, &m| a | m)));
            let ls = last_exercise(f, &o);
            let class = if sim_fit {
                "linear on the reach"
            } else if ls <= 512 {
                "burn-in, last exercise t <= 512"
            } else {
                "exercised past t = 512"
            };
            let oi = pr.inner.map(|r| orbit(f, pr.seed, horizon, r, horizon / 2));
            let inner = oi
                .as_ref()
                .map(|oi| (rounds_exercised(f, oi), last_exercise(f, oi)));
            let rho = match (&oi, pr.region) {
                (Some(oi), _) => Some(density(f, oi)),
                (None, Region::Bulk(_)) => Some(density(f, &o)),
                _ => None,
            };
            let cols = columns(f, pr.seed, pr.steps, pr.centre, &pr.offsets);
            let kinds: Vec<(Kind, usize, bool)> = cols.iter().map(|c| kind(c)).collect();
            let worst = kinds.iter().map(|k| k.0).max().unwrap();
            let top = kinds.iter().map(|k| k.1).max().unwrap();
            if class.starts_with("burn-in") {
                burn_max = burn_max.max(ls);
            }
            if class.starts_with("exercised") {
                past_min = past_min.min(ls);
            }
            if class.starts_with("burn-in") && worst == Kind::Grows {
                let late = orbit(f, pr.seed, pr.steps, pr.region, usize::MAX);
                println!(
                    "    burn-in rule {f} with a growing column: last exercise t = {} over t <= {}, kernel counts {:?}, depth profile {:?}",
                    last_exercise(f, &late),
                    pr.steps,
                    kinds.iter().map(|k| k.1).collect::<Vec<_>>(),
                    super::kernel_counts(&cols[0], 2, 24)
                );
            }
            if worst == Kind::Automatic && kinds.iter().all(|k| k.2) {
                *twice.entry(class).or_default() += 1;
            }
            *table.entry((class, worst)).or_default() += 1;
            lists.entry((class, worst)).or_default().push(f);
            tops.entry((class, worst)).or_default().push(top);
            if let Some(r) = rho {
                rhos.entry((class, worst)).or_default().push(r);
            }
            if sim_fit {
                lrs_n += 1;
                lrs_max = lrs_max.max(top);
                if worst == Kind::Grows {
                    bad.push(format!("rule {f} grows {:?}", kinds));
                }
            } else if worst == Kind::Automatic {
                aperiodic.push(format!(
                    "rule {f} last exercise t = {ls}, rounds exercised {}, interior {:?}, density {:?}, kernel counts {:?}",
                    rounds_exercised(f, &o),
                    inner,
                    rho,
                    kinds.iter().map(|k| k.1).collect::<Vec<_>>()
                ));
            }
            if [18, 22, 30, 90, 110, 126, 146, 150, 182].contains(&f) {
                reference.push(format!(
                    "rule {f}: {:?} counts {:?} last exercise {ls}",
                    kinds.iter().map(|k| k.0).collect::<Vec<_>>(),
                    kinds.iter().map(|k| k.1).collect::<Vec<_>>()
                ));
            }
        }
        for ((class, kind), n) in &table {
            println!("  {class:<34} worst column {kind:?}: {n}");
            if *kind != Kind::Periodic && !class.starts_with("linear") {
                println!("    rules {:?}", lists[&(*class, *kind)]);
            }
            if let Some(v) = tops.get(&(*class, *kind)) {
                println!(
                    "    largest last-depth count per rule: {} to {}",
                    v.iter().min().unwrap(),
                    v.iter().max().unwrap()
                );
            }
            if *kind == Kind::Automatic {
                println!(
                    "    plateau also at word length 48: {}",
                    twice.get(class).copied().unwrap_or(0)
                );
            }
            if let Some(v) = rhos.get_mut(&(*class, *kind)) {
                v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                println!(
                    "    exercise density over t > {}: min {:.4} median {:.4} max {:.4}",
                    horizon / 2,
                    v[0],
                    v[v.len() / 2],
                    v[v.len() - 1]
                );
            }
        }
        println!("  last exercise: burn-in at most {burn_max}, exercised at least {past_min}");
        println!(
            "  linear on the reach: {lrs_n} rules, largest last-depth count {lrs_max}, with a growing column {}",
            bad.len()
        );
        for b in &bad {
            println!("    {b}");
        }
        println!(
            "  not linear on the reach, aperiodic, kernel stops: {}",
            aperiodic.len()
        );
        for a in &aperiodic {
            println!("    {a}");
        }
        for r in &reference {
            println!("  {r}");
        }
    }
}

// LEMMAS

fn row_windows(radius: i32) -> HashSet<(u16, u8)> {
    let g = (1u8, 0u8, 1u8);
    let u: Vec15 = 1 << slot(0, 0);
    let init: Vec<Comp> = (-radius..=radius).map(|o| (u, o)).collect();
    let mut out: HashSet<(u16, u8)> = HashSet::new();
    out.insert((pattern(&init, 0, 0), 0));
    out.insert((pattern(&init, -1, 0), 1));
    let mut seen: HashSet<State> = HashSet::new();
    let mut queue: VecDeque<State> = VecDeque::new();
    for e0 in 0..2 {
        for f0 in 0..2 {
            let s = (step(&init, g, e0, f0), (e0 ^ f0) as u8);
            if seen.insert(s.clone()) {
                queue.push_back(s);
            }
        }
    }
    while let Some((s, q)) = queue.pop_front() {
        for m in [0i64, -1] {
            out.insert((pattern(&s, m, 0), q));
        }
        for e1 in 0..2 {
            for f1 in 0..2 {
                let n = (step(&s, g, e1, f1), q);
                if seen.insert(n.clone()) {
                    queue.push_back(n);
                }
            }
        }
    }
    println!(
        "  closure of the width-{} windows of the rule 90 rows: {} states, {} tagged windows",
        2 * radius + 1,
        seen.len(),
        out.len()
    );
    out
}

fn cell(rule: u8, l: u8, c: u8, r: u8) -> u8 {
    rule >> (4 * l + 2 * c + r) & 1
}

fn pascal(s: i64, n: i64) -> u8 {
    if s < 0 || n.abs() > s || (s + n) % 2 != 0 {
        return 0;
    }
    let k = (s + n) / 2;
    u8::from(k & (s - k) == 0)
}

fn psi22(row: &dyn Fn(i64) -> u8, even: bool, n: i64) -> u8 {
    if even {
        row(n - 1) ^ row(n + 1)
    } else {
        row(n) & (1 ^ row(n - 2)) & (1 ^ row(n + 2))
    }
}

pub fn lemma() {
    let horizon = 4096usize;
    println!("lemma: nonlinear rules that draw automatic diagrams from one cell");
    let w = 2 * horizon + 16;
    let off = (horizon + 8) as i64;
    let mut x182 = 0usize;
    let mut x22 = 0usize;
    let mut x126 = 0usize;
    let mut rows: Vec<Vec<u8>> = vec![vec![0u8; w]; 4];
    rows[0][off as usize] = 1;
    rows[1][off as usize] = 1;
    rows[2][off as usize] = 1;
    rows[3][off as usize - 1] = 1;
    rows[3][off as usize] = 1;
    let rules = [182u8, 22, 126, 18];
    for t in 0..=horizon as i64 {
        for n in 0..w as i64 {
            let nn = n - off;
            let got182 = rows[0][n as usize];
            let want182 = u8::from(nn.abs() <= t) & (1 ^ pascal(t + 1, nn));
            x182 += usize::from(got182 != want182);
            if t >= 1 {
                let row = |m: i64| pascal(t - 1, m);
                let want22 = psi22(&row, (t + nn) % 2 == 0, nn);
                x22 += usize::from(rows[1][n as usize] != want22);
            }
            if (n as usize) + 1 < w {
                let d = rows[2][n as usize] ^ rows[2][n as usize + 1];
                x126 += usize::from(d != rows[3][n as usize]);
            }
        }
        for (i, rule) in rules.iter().enumerate() {
            let cur = rows[i].clone();
            for n in 0..w {
                let l = if n > 0 { cur[n - 1] } else { 0 };
                let r = if n + 1 < w { cur[n + 1] } else { 0 };
                rows[i][n] = cell(*rule, l, cur[n], r);
            }
        }
    }
    println!("  rule 182 against [|n| <= t] (1 + P(t+1, n)), t <= {horizon}: {x182} cells differ");
    println!("  rule 22 against the block code of the rule 90 row t-1, t <= {horizon}: {x22} cells differ");
    println!("  rule 126, difference field against rule 18 from the pair 11 at (-1, 0), t <= {horizon}: {x126} cells differ");
    let o126 = orbit(126, Seed::Cell, horizon, Region::Interior, usize::MAX);
    let off90: Vec<usize> = (0..=horizon)
        .filter(|&t| !fits(126, 90, o126.masks[t]))
        .collect();
    println!(
        "  rule 126, rounds t <= {horizon} whose interior windows rule 90 does not fit: {} of them, {:?}",
        off90.len(),
        off90
    );
    let wins = row_windows(3);
    let mut sim: HashSet<(u16, u8)> = HashSet::new();
    for row in 0..=horizon as i64 {
        for n in -row - 4..=row + 4 {
            let pat = (-3..=3i64).fold(0u16, |a, j| a | u16::from(pascal(row, n + j)) << (j + 3));
            sim.insert((pat, (row + n).rem_euclid(2) as u8));
        }
    }
    println!(
        "  simulated windows of the rows s <= {horizon}: {}, inside the closure: {}",
        sim.len(),
        sim.iter().all(|w| wins.contains(w))
    );
    let mut bad = 0usize;
    let mut checked = 0usize;
    for &(pat, q) in &wins {
        let p = |j: i32| (pat >> (j + 3) & 1) as u8;
        let x = |j: i32| {
            let even = (q as i32 + 1 + j).rem_euclid(2) == 0;
            if even {
                p(j - 1) ^ p(j + 1)
            } else {
                p(j) & (1 ^ p(j - 2)) & (1 ^ p(j + 2))
            }
        };
        let lhs = cell(22, x(-1), x(0), x(1));
        let rhs = if q == 0 {
            p(-2) ^ p(2)
        } else {
            let pt = |k: i32| p(k - 1) ^ p(k + 1);
            pt(0) & (1 ^ pt(-2)) & (1 ^ pt(2))
        };
        checked += 1;
        bad += usize::from(lhs != rhs);
    }
    println!("  rule 22 step lemma over every window of every row of rule 90: {checked} tagged windows, {bad} violations");
    let mut first = true;
    for n in -6..=6i64 {
        let row = |m: i64| pascal(0, m);
        let want = psi22(&row, (1 + n) % 2 == 0, n);
        let got = u8::from(n.abs() <= 1);
        if want != got {
            first = false;
        }
    }
    println!("  rule 22 base case x_1 = 111 from the block code of the row 0: {first}");
    pascal_images(horizon);
}

fn rho(t: i64, n: i64) -> i64 {
    if (t + n).rem_euclid(2) == 0 {
        n
    } else if n >= 0 {
        n + 1
    } else {
        n - 1
    }
}

fn step_bg(rule: u8, cur: &[u8], bg: u8) -> Vec<u8> {
    let w = cur.len();
    (0..w)
        .map(|n| {
            let l = if n > 0 { cur[n - 1] } else { bg };
            let r = if n + 1 < w { cur[n + 1] } else { bg };
            cell(rule, l, cur[n], r)
        })
        .collect()
}

fn pascal_images(horizon: usize) {
    let w = 2 * horizon + 16;
    let off = (horizon + 8) as i64;
    let mut seed = vec![0u8; w];
    seed[off as usize] = 1;
    let (mut a126, mut a129, mut a161) = (seed.clone(), seed.clone(), seed);
    let mut bg = 0u8;
    let (mut form, mut comp, mut same, mut lone) = (0usize, 0usize, 0usize, 0usize);
    for t in 0..=horizon as i64 {
        for n in 0..w {
            let nn = n as i64 - off;
            form += usize::from(a126[n] != pascal(t, rho(t, nn)));
            let want = if t == 0 {
                u8::from(nn == 0)
            } else {
                1 ^ a126[n]
            };
            comp += usize::from(a129[n] != want);
            same += usize::from(a161[n] != a129[n]);
            if t >= 1 && n > 0 && n + 1 < w {
                lone += usize::from(a126[n - 1] == 0 && a126[n] == 1 && a126[n + 1] == 0);
            }
        }
        a126 = step_bg(126, &a126, 0);
        a129 = step_bg(129, &a129, bg);
        a161 = step_bg(161, &a161, bg);
        bg = cell(129, bg, bg, bg);
    }
    println!("  rule 126 against P(t, rho_t(n)), rho_t(n) = n if t + n is even and n + sgn(n) if odd, t <= {horizon}: {form} cells differ");
    println!("  rule 129 against 1 + rule 126 for t >= 1 (the seed at t = 0), background tracked, t <= {horizon}: {comp} cells differ");
    println!("  rule 161 against rule 129, t <= {horizon}: {same} cells differ");
    println!("  rule 126 rows t >= 1 holding an isolated one (window 010), t <= {horizon}: {lone}");
    let steps = 1usize << 15;
    let probes: Vec<i64> = (-8..=8).collect();
    let cols = columns(126, Seed::Cell, steps, 0, &probes);
    let mut diff = 0usize;
    let mut worst = 0usize;
    let mut counts = Vec::new();
    for (j, &n0) in probes.iter().enumerate() {
        let sig = if n0 >= 0 { 1 } else { -1 };
        for t in 0..steps {
            let want = pascal(t as i64, n0) ^ pascal(t as i64, n0 + sig);
            diff += usize::from(cols[j][t] != want);
        }
        let last = *super::kernel_counts(&cols[j], 2, 24).last().unwrap();
        worst = worst.max(last);
        counts.push(last);
    }
    println!("  rule 126 columns -8..8 against P(t, n0) + P(t, n0 + sgn(n0)), t < 2^15: {diff} cells differ");
    println!("  rule 126 kernel words at depth 10 over columns -8..8: {counts:?}, largest {worst}");
}
