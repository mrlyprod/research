use std::collections::{BTreeMap, HashSet};
use std::env;
use std::time::Instant;

// CATALOG

const NEIGHBOURHOODS: [(&str, u32); 3] = [("carpet", 495), ("net", 186), ("void", 341)];
const SOURCES: [(&str, &[u32]); 4] = [
    ("primes", &[2, 3, 5, 7]),
    ("binary", &[1, 2, 4, 8]),
    ("fibonacci", &[0, 1, 2, 3, 5, 8]),
    ("odds", &[1, 3, 5, 7]),
];
const SPREADS: [&str; 6] = ["id", "shift", "drop", "slide", "turn", "flip"];
const LAMBDA: f64 = 0.624_329_988_543_550_9;
const TRIES: u32 = 512;

fn spread_at(name: &str, n: usize, i: usize, j: usize) -> (usize, usize) {
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

fn offsets(code: u32) -> Vec<(i64, i64)> {
    (0..9)
        .filter(|&b| code >> b & 1 == 1 && b != 4)
        .map(|b| (b as i64 / 3 - 1, b as i64 % 3 - 1))
        .collect()
}

fn rule_of(values: &[u32]) -> u32 {
    values
        .iter()
        .filter(|&&v| v < 32)
        .fold(0, |r, &v| r | 1 << v)
}

fn word(k1: u32, k2: u32) -> u32 {
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

// TORUS

fn symmetries(n: usize) -> Vec<Vec<usize>> {
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

fn image(perm: &[usize], mask: u32) -> u32 {
    (0..perm.len())
        .filter(|&p| mask >> p & 1 == 1)
        .fold(0, |m, p| m | 1 << perm[p])
}

fn commutes(a: &[usize], b: &[usize]) -> bool {
    (0..a.len()).all(|p| a[b[p]] == b[a[p]])
}

// ROUND

struct Round {
    bits: usize,
    fill: u32,
    reads: Vec<(u32, u32)>,
    rule: u32,
    perm: Vec<usize>,
    lo: [u32; 256],
    hi: [u32; 256],
    blo: [u32; 256],
    bhi: [u32; 256],
}

impl Round {
    fn new(n: usize, fill: u32, hood: u32, rule: u32, spread: &str) -> Round {
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
        let perm: Vec<usize> = (0..bits)
            .map(|p| {
                let (i, j) = spread_at(spread, n, p / n, p % n);
                i * n + j
            })
            .collect();
        let mut back = vec![0; bits];
        for p in 0..bits {
            back[perm[p]] = p;
        }
        let table = |map: &[usize], shift: usize| {
            let mut t = [0u32; 256];
            for (byte, slot) in t.iter_mut().enumerate() {
                for b in 0..8 {
                    let p = shift + b;
                    if p < bits && byte >> b & 1 == 1 {
                        *slot |= 1 << map[p];
                    }
                }
            }
            t
        };
        Round {
            bits,
            fill,
            rule,
            lo: table(&perm, 0),
            hi: table(&perm, 8),
            blo: table(&back, 0),
            bhi: table(&back, 8),
            reads,
            perm,
        }
    }

    fn deps(&self) -> Vec<(u32, u32)> {
        self.reads
            .iter()
            .map(|&(a, mask)| {
                let k = mask.count_ones();
                let live = (0..k).any(|c| (self.rule >> c & 1) != (self.rule >> (c + 1) & 1));
                (a, if live { mask } else { 0 })
            })
            .collect()
    }

    fn complement_free(&self) -> bool {
        self.reads.iter().all(|&(_, mask)| {
            let k = mask.count_ones();
            (0..=k).all(|c| (self.rule >> c & 1) == (self.rule >> (k - c) & 1))
        })
    }

    fn mess(&self, x: u32) -> u32 {
        self.reads.iter().fold(0, |m, &(a, mask)| {
            m | ((self.rule >> (x & mask).count_ones()) & 1) << a
        })
    }

    fn feistel(&self, x: u32) -> u32 {
        x ^ self.mess(x)
    }

    fn spread(&self, x: u32) -> u32 {
        self.lo[(x & 255) as usize] | self.hi[(x >> 8 & 255) as usize]
    }

    fn unspread(&self, x: u32) -> u32 {
        self.blo[(x & 255) as usize] | self.bhi[(x >> 8 & 255) as usize]
    }

    fn table(&self) -> Vec<u32> {
        (0..1u32 << self.bits)
            .map(|x| self.spread(self.feistel(x)))
            .collect()
    }

    fn weight(&self) -> usize {
        (0..1u32 << self.bits)
            .filter(|&x| x & self.fill == 0 && self.mess(x) != 0)
            .count()
    }

    fn weight_law(&self) -> usize {
        let k = self.reads[0].1.count_ones() as usize;
        (1 << (self.bits - 1 - k)) * ones(self.rule, k)
    }
}

fn binomial(k: usize, s: usize) -> usize {
    (0..s).fold(1, |c, i| c * (k - i) / (i + 1))
}

fn ones(rule: u32, k: usize) -> usize {
    (0..=k)
        .filter(|&s| rule >> s & 1 == 1)
        .map(|s| binomial(k, s))
        .sum()
}

fn weights(n: usize) {
    let rest = n * n - 1;
    for (hood, code) in NEIGHBOURHOODS {
        let k = Round::new(n, 1, code, 0, "id").reads[0].1.count_ones() as usize;
        let sums: Vec<String> = SOURCES
            .iter()
            .map(|(name, s)| format!("{name} {}", ones(rule_of(s), k)))
            .collect();
        let odd: Vec<&str> = SOURCES
            .iter()
            .filter(|(_, s)| k == rest && ones(rule_of(s), k) % 2 == 1)
            .map(|(name, _)| *name)
            .collect();
        println!(
            "  one written cell, {hood} reads {k} of {rest}: sum of C({k}, s) over S: {}; weight 2^{} times it, odd for {odd:?}",
            sums.join(", "),
            rest - k
        );
    }
}

// PERMUTATIONS

fn cycles(perm: &[u32]) -> Vec<usize> {
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

fn is_permutation(perm: &[u32]) -> bool {
    let mut seen = vec![false; perm.len()];
    perm.iter()
        .all(|&y| !std::mem::replace(&mut seen[y as usize], true))
}

fn even(lengths: &[usize]) -> bool {
    lengths.iter().map(|l| l - 1).sum::<usize>() % 2 == 0
}

fn is_prime(p: usize) -> bool {
    p >= 2 && (2..).take_while(|d| d * d <= p).all(|d| p % d != 0)
}

fn lone_prime(lengths: &[usize], size: usize) -> Option<usize> {
    let mut count: BTreeMap<usize, usize> = BTreeMap::new();
    for &l in lengths {
        *count.entry(l).or_default() += 1;
    }
    count
        .iter()
        .filter(|&(&l, &c)| c == 1 && is_prime(l) && l + 3 <= size)
        .map(|(&l, _)| l)
        .find(|&p| lengths.iter().filter(|&&l| l % p == 0).count() == 1)
}

// SUBSPACES

#[derive(Clone)]
struct Basis {
    rows: [u32; 32],
    dim: usize,
    gens: Vec<u32>,
}

impl Basis {
    fn new() -> Basis {
        Basis {
            rows: [0; 32],
            dim: 0,
            gens: Vec::new(),
        }
    }

    fn reduce(&self, mut v: u32) -> u32 {
        while v != 0 {
            let top = 31 - v.leading_zeros() as usize;
            if self.rows[top] == 0 {
                return v;
            }
            v ^= self.rows[top];
        }
        0
    }

    fn insert(&mut self, v: u32) -> bool {
        let r = self.reduce(v);
        if r == 0 {
            return false;
        }
        self.rows[31 - r.leading_zeros() as usize] = r;
        self.dim += 1;
        self.gens.push(v);
        true
    }
}

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

fn block(gens: &[&[u32]], v: u32, bits: usize, full: &[bool], rng: &mut Rng) -> Option<Basis> {
    let size = 1usize << bits;
    let mut u = Basis::new();
    u.insert(v);
    for _ in 0..4 * bits {
        let g = u.gens[rng.next() as usize % u.gens.len()];
        let x = rng.next() as usize % size;
        let h = gens[rng.next() as usize % gens.len()];
        let w = h[x ^ g as usize] ^ h[x];
        if full[w as usize] {
            return None;
        }
        u.insert(w);
        if u.dim == bits {
            return None;
        }
    }
    let mut next = 0;
    while next < u.gens.len() {
        let g = u.gens[next] as usize;
        next += 1;
        for h in gens {
            for x in 0..size {
                let w = h[x ^ g] ^ h[x];
                u.insert(w);
                if u.dim == bits || full[w as usize] {
                    return None;
                }
            }
        }
    }
    Some(u)
}

fn first_block(gens: &[&[u32]], bits: usize) -> Option<Basis> {
    let size = 1usize << bits;
    let mut full = vec![false; size];
    let mut rng = Rng(0x9e37_79b9_7f4a_7c15);
    for v in 1..size {
        match block(gens, v as u32, bits, &full, &mut rng) {
            None => full[v] = true,
            Some(u) => return Some(u),
        }
    }
    None
}

fn perp(u: &Basis, bits: usize) -> Vec<u32> {
    let mut out = Basis::new();
    for phi in 1..1u32 << bits {
        if u.gens.iter().all(|&g| (g & phi).count_ones() % 2 == 0) {
            out.insert(phi);
        }
    }
    out.gens
}

fn is_block(gens: &[&[u32]], u: &Basis, bits: usize) -> bool {
    gens.iter().all(|h| {
        u.gens
            .iter()
            .all(|&g| (0..1usize << bits).all(|x| u.reduce(h[x ^ g as usize] ^ h[x]) == 0))
    })
}

// LAWS

fn orbit_of(perm: &[usize], s: usize) -> u32 {
    let (mut p, mut orbit) = (s, 0u32);
    loop {
        orbit |= 1 << p;
        p = perm[p];
        if p == s {
            return orbit;
        }
    }
}

fn closed_set(round: &Round) -> Option<u32> {
    let all = ((1u64 << round.bits) - 1) as u32;
    (0..round.bits).find_map(|s| {
        let mut set = orbit_of(&round.perm, s);
        loop {
            let grown = round
                .deps()
                .iter()
                .filter(|&&(_, mask)| mask & set != 0)
                .fold(set, |g, &(a, _)| g | orbit_of(&round.perm, a as usize));
            if grown == set {
                break;
            }
            set = grown;
        }
        (set != all).then_some(set)
    })
}

fn subgroups(n: usize) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = Vec::new();
    for g in 0..n * n {
        for k in 0..n * n {
            let mut set = vec![0usize];
            let mut i = 0;
            while i < set.len() {
                let p = set[i];
                i += 1;
                for q in [g, k] {
                    let r = ((p / n + q / n) % n) * n + (p % n + q % n) % n;
                    if !set.contains(&r) {
                        set.push(r);
                    }
                }
            }
            set.sort();
            if set.len() > 1 && set.len() < n * n && !out.contains(&set) {
                out.push(set);
            }
        }
    }
    out
}

fn product_law(round: &Round, n: usize) -> Option<usize> {
    let deps = round.deps();
    subgroups(n).into_iter().find_map(|h| {
        let label: Vec<usize> = (0..n * n)
            .map(|p| {
                h.iter()
                    .map(|&q| ((p / n + q / n) % n) * n + (p % n + q % n) % n)
                    .min()
                    .unwrap()
            })
            .collect();
        let permuted = (0..n * n).all(|p| {
            (0..n * n).all(|q| label[p] != label[q] || label[round.perm[p]] == label[round.perm[q]])
        });
        let local = deps.iter().all(|&(a, mask)| {
            (0..n * n).all(|b| mask >> b & 1 == 0 || label[b] == label[a as usize])
        });
        (permuted && local).then_some(h.len())
    })
}

fn partitions(n: usize) -> Vec<(String, Vec<usize>)> {
    let mut out: Vec<(String, Vec<usize>)> = vec![("torus".to_string(), vec![0; n * n])];
    for h in subgroups(n) {
        let label = (0..n * n)
            .map(|p| {
                h.iter()
                    .map(|&q| ((p / n + q / n) % n) * n + (p % n + q % n) % n)
                    .min()
                    .unwrap()
            })
            .collect();
        out.push((format!("cosets of a subgroup of order {}", h.len()), label));
    }
    if n == 4 {
        out.push((
            "2 x 2 tiles".to_string(),
            (0..16).map(|p| (p / 8) * 2 + (p % 4) / 2).collect(),
        ));
    }
    out
}

fn parity_law(h: &[u32], n: usize) -> Option<String> {
    partitions(n).into_iter().find_map(|(name, label)| {
        let mut u = Basis::new();
        for p in 0..n * n {
            if let Some(q) = (p + 1..n * n).find(|&q| label[q] == label[p]) {
                u.insert(1 << p | 1 << q);
            }
        }
        (u.dim > 0 && is_block(&[h], &u, n * n)).then_some(name)
    })
}

fn module(perm: &[usize], seeds: &[u32]) -> Basis {
    let mut u = Basis::new();
    let mut queue: Vec<u32> = seeds.to_vec();
    while let Some(v) = queue.pop() {
        if u.insert(v) {
            queue.push(image(perm, v));
        }
    }
    u
}

fn classes(round: &Round) -> Vec<u32> {
    let mut by: BTreeMap<u32, u32> = BTreeMap::new();
    for (a, mask) in round.deps() {
        if mask != 0 {
            *by.entry(mask).or_default() |= 1 << a;
        }
    }
    by.into_values().collect()
}

fn mess_module(round: &Round) -> Basis {
    module(&round.perm, &classes(round))
}

// CELLS

struct Cell {
    source: usize,
    spread: usize,
    cycles: usize,
    longest: usize,
    fixed: usize,
    even: bool,
    sym: bool,
    group: String,
    law: String,
    prime: usize,
    word: usize,
    laws: [bool; 5],
}

struct Verdict {
    group: &'static str,
    block: Option<Basis>,
    prime: usize,
    word: usize,
}

fn verdict(
    gens: &[&[u32]],
    bits: usize,
    odd: bool,
    words: &dyn Fn(u32) -> Vec<u32>,
    known: &'static str,
) -> Verdict {
    let size = 1usize << bits;
    let alt = if odd { "symmetric" } else { "alternating" };
    if let Some(u) = first_block(gens, bits) {
        assert!(is_block(gens, &u, bits), "a closure that is not a block");
        return Verdict {
            group: "imprimitive",
            block: Some(u),
            prime: 0,
            word: 0,
        };
    }
    for c in 0..TRIES.min(size as u32) {
        if let Some(p) = lone_prime(&cycles(&words(c)), size) {
            return Verdict {
                group: alt,
                block: None,
                prime: p,
                word: 1,
            };
        }
    }
    let mut rng = Rng(0x2545_f491_4f6c_dd1d);
    for t in 0..TRIES {
        let len = 2 + t as usize % 3;
        let mut w: Vec<u32> = (0..size as u32).collect();
        for _ in 0..len {
            let h = gens[rng.next() as usize % gens.len()];
            let c = rng.next() as u32 % size as u32;
            for y in w.iter_mut() {
                *y = h[*y as usize] ^ c;
            }
        }
        if let Some(p) = lone_prime(&cycles(&w), size) {
            return Verdict {
                group: alt,
                block: None,
                prime: p,
                word: len,
            };
        }
    }
    Verdict {
        group: if known.is_empty() { "undecided" } else { known },
        block: None,
        prime: 0,
        word: 0,
    }
}

fn reps(masks: &[u32], group: &[&Vec<usize>], full: Option<u32>) -> Vec<u32> {
    let mut seen = HashSet::new();
    masks
        .iter()
        .copied()
        .filter(|&m| {
            let canon = group
                .iter()
                .flat_map(|s| [image(s, m), full.map_or(u32::MAX, |f| image(s, f & !m))])
                .min()
                .unwrap();
            seen.insert(canon)
        })
        .collect()
}

fn census(n: usize, masks: &[u32], verbose: bool) -> (Vec<Cell>, usize) {
    let bits = n * n;
    let size = 1usize << bits;
    let syms = symmetries(n);
    let mut out = Vec::new();
    let mut faults = 0usize;
    for spread in 0..SPREADS.len() {
        let probe = Round::new(n, 0, 0, 0, SPREADS[spread]);
        let central: Vec<&Vec<usize>> = syms.iter().filter(|s| commutes(s, &probe.perm)).collect();
        for hood in 0..NEIGHBOURHOODS.len() {
            for source in 0..SOURCES.len() {
                for fill in reps(masks, &central, None) {
                    let round = Round::new(
                        n,
                        fill,
                        NEIGHBOURHOODS[hood].1,
                        rule_of(SOURCES[source].1),
                        SPREADS[spread],
                    );
                    let h = round.table();
                    let r: Vec<u32> = h.iter().map(|y| y ^ fill).collect();
                    let undo = (0..size as u32).all(|x| {
                        let f = round.feistel(x);
                        round.feistel(f) == x
                            && round.feistel(round.unspread(r[x as usize] ^ fill)) == x
                    });
                    let lengths = cycles(&r);
                    let weight = if fill.count_ones() == 1 {
                        round.weight()
                    } else {
                        0
                    };
                    let predicted = fill.count_ones() >= 2 || weight % 2 == 0;
                    if !is_permutation(&r)
                        || !undo
                        || even(&lengths) != predicted
                        || fill.count_ones() == 1 && weight != round.weight_law()
                    {
                        faults += 1;
                    }
                    let nonlinear = SOURCES[source].0 != "odds";
                    let product = product_law(&round, n).is_some();
                    let known = if !nonlinear {
                        "affine"
                    } else if product {
                        "product"
                    } else {
                        ""
                    };
                    let words = |c: u32| h.iter().map(|y| y ^ c).collect::<Vec<u32>>();
                    let v = verdict(&[&h], bits, !predicted, &words, known);
                    let closed = closed_set(&round).is_some();
                    let mess = mess_module(&round).dim < bits;
                    let complement = round.complement_free();
                    let parity = v.block.is_none() && parity_law(&h, n).is_some();
                    if (closed || mess || complement || parity) && v.block.is_none()
                        || !known.is_empty() && v.prime > 0
                    {
                        faults += 1;
                    }
                    let name = format!(
                        "n {n} fill {fill:#06x} {} {} {}",
                        NEIGHBOURHOODS[hood].0, SOURCES[source].0, SPREADS[spread]
                    );
                    let law = match (&v.block, nonlinear) {
                        (None, _) => String::new(),
                        (Some(_), false) => "affine".to_string(),
                        _ if closed => "closed".to_string(),
                        _ if mess => "mess".to_string(),
                        _ if complement => "complement".to_string(),
                        _ => parity_law(&h, n)
                            .map_or("unexplained".to_string(), |p| format!("parity of {p}")),
                    };
                    if let Some(u) = &v.block {
                        if verbose || law.starts_with("parity") || law == "unexplained" {
                            println!(
                                "  {law}: {name}: block dim {}, kept parities {:04x?}",
                                u.dim,
                                perp(u, bits)
                            );
                        }
                    }
                    if v.group == "undecided" {
                        println!("  undecided: {name}");
                    }
                    if verbose {
                        println!(
                            "  {} {law}: {name}: cycles {}, longest {}, fixed {}, prime {}",
                            v.group,
                            lengths.len(),
                            lengths.iter().max().unwrap(),
                            lengths.iter().filter(|&&l| l == 1).count(),
                            v.prime
                        );
                    }
                    out.push(Cell {
                        source,
                        spread,
                        cycles: lengths.len(),
                        longest: *lengths.iter().max().unwrap(),
                        fixed: lengths.iter().filter(|&&l| l == 1).count(),
                        even: even(&lengths),
                        sym: central.iter().filter(|s| image(s, fill) == fill).count() > 1,
                        group: v.group.to_string(),
                        law,
                        prime: v.prime,
                        word: v.word,
                        laws: [closed, mess, complement, parity, product],
                    });
                }
            }
        }
    }
    println!("  n {n}: {} cells up to the torus symmetries that commute with the spread, {faults} faults in bijection, involution, undo, sign, weight and law", out.len());
    (out, faults)
}

// EXPECTATIONS

fn harmonic(m: usize) -> (f64, f64) {
    let h: f64 = (1..=m).map(|k| 1.0 / k as f64).sum();
    let h2: f64 = (1..=m).map(|k| 1.0 / (k * k) as f64).sum();
    (h, (h - h2).sqrt())
}

fn longest_exact(m: usize) -> f64 {
    let mut total = 0.0;
    for k in 0..m {
        let mut p = vec![0.0f64; m + 1];
        let mut prefix = vec![0.0f64; m + 2];
        p[0] = 1.0;
        prefix[1] = 1.0;
        for j in 1..=m {
            let lo = j.saturating_sub(k);
            p[j] = (prefix[j] - prefix[lo]) / j as f64;
            prefix[j + 1] = prefix[j] + p[j];
        }
        total += 1.0 - p[m];
    }
    total
}

fn null() {
    for m in 2..=8usize {
        let mut perm: Vec<usize> = (0..m).collect();
        let (mut count, mut cyc, mut fix) = (0u64, 0u64, 0u64);
        loop {
            let table: Vec<u32> = perm.iter().map(|&p| p as u32).collect();
            let lengths = cycles(&table);
            if even(&lengths) {
                count += 1;
                cyc += lengths.len() as u64;
                fix += lengths.iter().filter(|&&l| l == 1).count() as u64;
            }
            let Some(i) = (0..m - 1).rev().find(|&i| perm[i] < perm[i + 1]) else {
                break;
            };
            let j = (i + 1..m).rev().find(|&j| perm[j] > perm[i]).unwrap();
            perm.swap(i, j);
            perm[i + 1..].reverse();
        }
        let (h, _) = harmonic(m);
        let formula = h + if m % 2 == 0 { 1.0 } else { -1.0 } / (m * (m - 1)) as f64;
        println!("  M {m}: {count} even permutations, mean cycles {:.12} against H_M + (-1)^M/(M(M-1)) = {formula:.12}, mean fixed points {:.12}", cyc as f64 / count as f64, fix as f64 / count as f64);
    }
    let (h, _) = harmonic(512);
    println!(
        "  M 512: H_M {h:.10}, even correction {:.3e}",
        1.0 / (512.0 * 511.0)
    );
}

// REPORTS

fn summary(n: usize, cells: &[Cell], expected_longest: f64) {
    let size = 1usize << (n * n);
    let (h, sd) = harmonic(size);
    println!("  uniform permutation of {size}: cycles {h:.4} (sd {sd:.4}), longest/size {:.5}, fixed points 1", expected_longest / size as f64);
    let linear = |c: &Cell| SOURCES[c.source].0 == "odds";
    let rows: [(&str, Box<dyn Fn(&Cell) -> bool>); 5] = [
        ("all", Box::new(|_| true)),
        ("odds", Box::new(move |c| linear(c))),
        ("nonlinear", Box::new(move |c| !linear(c))),
        (
            "nonlinear symmetric",
            Box::new(move |c| !linear(c) && c.sym),
        ),
        (
            "nonlinear asymmetric",
            Box::new(move |c| !linear(c) && !c.sym),
        ),
    ];
    for (name, keep) in rows.iter() {
        let picked: Vec<&Cell> = cells.iter().filter(|c| keep(c)).collect();
        let k = picked.len() as f64;
        let mean = |f: &dyn Fn(&Cell) -> f64| picked.iter().map(|c| f(c)).sum::<f64>() / k;
        let cyc = mean(&|c| c.cycles as f64);
        let mut sorted: Vec<usize> = picked.iter().map(|c| c.cycles).collect();
        sorted.sort();
        println!(
            "  {name:<21} cells {:>5}  cycles mean {cyc:>9.2} median {:>6}  longest/size {:.5}  fixed {:>8.2}  odd {}",
            picked.len(),
            sorted[sorted.len() / 2],
            mean(&|c| c.longest as f64 / size as f64),
            mean(&|c| c.fixed as f64),
            picked.iter().filter(|c| !c.even).count()
        );
    }
}

fn groups(cells: &[Cell]) {
    let mut table: BTreeMap<String, BTreeMap<String, usize>> = BTreeMap::new();
    for c in cells {
        let tag = if SOURCES[c.source].0 == "odds" {
            "odds"
        } else {
            SPREADS[c.spread]
        };
        *table
            .entry(tag.to_string())
            .or_default()
            .entry(c.group.clone())
            .or_default() += 1;
    }
    for (tag, row) in &table {
        println!("  {tag:<6} {row:?}");
    }
    let mut laws: BTreeMap<String, usize> = BTreeMap::new();
    for c in cells.iter().filter(|c| c.group == "imprimitive") {
        *laws
            .entry(c.law.split(" of ").next().unwrap().to_string())
            .or_default() += 1;
    }
    println!("  imprimitive by first law: {laws:?}");
    for (name, linear) in [("nonlinear", false), ("odds", true)] {
        let kind: Vec<&Cell> = cells
            .iter()
            .filter(|c| (SOURCES[c.source].0 == "odds") == linear)
            .collect();
        let primitive: Vec<&&Cell> = kind.iter().filter(|c| c.group != "imprimitive").collect();
        let held: Vec<usize> = (0..5)
            .map(|l| primitive.iter().filter(|c| c.laws[l]).count())
            .collect();
        println!(
            "  {name}: product law in {} cells; primitive {}, laws I-V hold in {} (I {}, II {}, III {}, IV {}, V {})",
            kind.iter().filter(|c| c.laws[4]).count(),
            primitive.len(),
            primitive.iter().filter(|c| c.laws.iter().any(|&b| b)).count(),
            held[0],
            held[1],
            held[2],
            held[3],
            held[4]
        );
    }
    let mut primes: BTreeMap<usize, usize> = BTreeMap::new();
    let mut lens: BTreeMap<usize, usize> = BTreeMap::new();
    for c in cells.iter().filter(|c| c.prime > 0) {
        *primes.entry(c.prime).or_default() += 1;
        *lens.entry(c.word).or_default() += 1;
    }
    println!("  lone primes {primes:?}");
    println!("  word lengths {lens:?}");
}

fn pairs(n: usize, masks: &[u32], verbose: bool) {
    let bits = n * n;
    let full = ((1u64 << bits) - 1) as u32;
    let syms = symmetries(n);
    let all: Vec<&Vec<usize>> = syms.iter().collect();
    let reps = reps(masks, &all, Some(full));
    let mut table: BTreeMap<String, BTreeMap<String, usize>> = BTreeMap::new();
    for hood in 0..NEIGHBOURHOODS.len() {
        for source in 0..SOURCES.len() {
            for &fill in &reps {
                let rule = rule_of(SOURCES[source].1);
                let one = Round::new(n, fill, NEIGHBOURHOODS[hood].1, rule, "id");
                let two = Round::new(n, full & !fill, NEIGHBOURHOODS[hood].1, rule, "id");
                let (h1, h2) = (one.table(), two.table());
                let odd = (fill.count_ones() == 1 && one.weight() % 2 == 1)
                    ^ ((full & !fill).count_ones() == 1 && two.weight() % 2 == 1);
                let words = |c: u32| h1.iter().map(|&y| h2[y as usize] ^ c).collect::<Vec<u32>>();
                let known = if SOURCES[source].0 == "odds" {
                    "affine"
                } else {
                    ""
                };
                let v = verdict(&[&h1, &h2], bits, odd, &words, known);
                if let Some(u) = &v.block {
                    if verbose && known.is_empty() {
                        println!("  pair n {n} fill {fill:#06x} {} {}: block dim {}, kept parities {:04x?}", NEIGHBOURHOODS[hood].0, SOURCES[source].0, u.dim, perp(u, bits));
                    }
                }
                *table
                    .entry(SOURCES[source].0.to_string())
                    .or_default()
                    .entry(v.group.to_string())
                    .or_default() += 1;
            }
        }
    }
    println!(
        "  n {n}: {} splits up to torus symmetry and complement",
        reps.len()
    );
    for (tag, row) in &table {
        println!("  {tag:<10} {row:?}");
    }
}

fn main() {
    let verb = env::args().nth(1).unwrap_or_else(|| "all".to_string());
    let start = Instant::now();
    let three: Vec<u32> = (1..511).collect();
    let four: Vec<u32> = (1..16)
        .flat_map(|a| (1..16).map(move |b| word(a, b)))
        .filter(|&m| m != 0xffff)
        .collect();
    let (mut total, mut faults) = (0, 0);
    if verb == "all" || verb == "three" || verb == "blocks" {
        println!("# 3 x 3: every base-3 level-1 split, {} masks", three.len());
        weights(3);
        let (cells, f) = census(3, &three, verb == "blocks");
        summary(3, &cells, longest_exact(512));
        groups(&cells);
        total += cells.len();
        faults += f;
    }
    if verb == "all" || verb == "four" {
        println!(
            "# 4 x 4: every base-2 two-letter word, {} masks",
            four.len()
        );
        weights(4);
        let (cells, f) = census(4, &four, false);
        summary(4, &cells, LAMBDA * 65536.0);
        groups(&cells);
        total += cells.len();
        faults += f;
    }
    if verb == "all" {
        println!("# both catalogs: {total} cells, {faults} faults");
    }
    if verb == "word" {
        let k: Vec<u32> = env::args().skip(2).map(|a| a.parse().unwrap()).collect();
        census(4, &[word(k[0], k[1])], true);
    }
    if verb == "null" {
        null();
    }
    if verb == "all" || verb == "pairs" {
        println!("# a split and its complement in turn, identity spread");
        pairs(3, &three, false);
        pairs(4, &four, false);
    }
    println!("# {:.1} s", start.elapsed().as_secs_f64());
}
