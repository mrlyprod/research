use mrly_num::gauss::Ring;
use mrly_math::radix::{flowsnake, koch, terdragon, twindragon, Base, Radix};
use std::collections::{HashMap, HashSet};
use std::time::Instant;

type Z = (i64, i64);

const BASES: [(Ring, Z); 7] = [
    (Ring::Gaussian, (2, 0)),
    (Ring::Gaussian, (1, 1)),
    (Ring::Gaussian, (2, 1)),
    (Ring::Eisenstein, (2, 0)),
    (Ring::Eisenstein, (2, 1)),
    (Ring::Eisenstein, (3, 0)),
    (Ring::Eisenstein, (3, 1)),
];

// RING

fn add(z: Z, w: Z) -> Z {
    (z.0 + w.0, z.1 + w.1)
}

fn sub(z: Z, w: Z) -> Z {
    (z.0 - w.0, z.1 - w.1)
}

#[derive(Clone)]
struct Lab {
    ring: Ring,
    base: Base,
    b: Z,
    n: usize,
    units: Vec<Z>,
    residues: Vec<Z>,
}

impl Lab {
    fn new(ring: Ring, value: Z) -> Lab {
        let base = Base::new(ring, value).unwrap();
        Lab {
            ring,
            base,
            b: value,
            n: ring.units(),
            units: ring.associates(1, 0),
            residues: base.residues().unwrap(),
        }
    }
    fn mul(&self, z: Z, w: Z) -> Z {
        self.ring.mul(z, w)
    }
    fn norm(&self, z: Z) -> u64 {
        self.ring.norm(z.0, z.1)
    }
    fn turn(&self, z: Z, k: usize) -> Z {
        self.mul(z, self.units[k % self.n])
    }
    fn back(&self, k: usize) -> usize {
        (self.n - k % self.n) % self.n
    }
    fn conj(&self, z: Z) -> Z {
        self.ring.conjugate(z.0, z.1)
    }
    fn mirror(&self) -> Option<usize> {
        let c = self.conj(self.b);
        (0..self.n).find(|&k| self.turn(self.b, k) == c)
    }
    fn sum(&self, walk: &[usize]) -> Z {
        walk.iter().fold((0, 0), |acc, &k| add(acc, self.units[k]))
    }
    fn power(&self, level: usize) -> Z {
        self.base.power(level)
    }
    fn modulus(&self) -> f64 {
        (self.base.norm() as f64).sqrt()
    }
}

fn mul_wide(ring: Ring, (a, b): (i128, i128), (c, d): (i128, i128)) -> (i128, i128) {
    match ring {
        Ring::Gaussian => (a * c - b * d, a * d + b * c),
        Ring::Eisenstein => (a * c - b * d, a * d + b * c - b * d),
    }
}

fn norm_wide(ring: Ring, (a, b): (i128, i128)) -> i128 {
    match ring {
        Ring::Gaussian => a * a + b * b,
        Ring::Eisenstein => a * a - a * b + b * b,
    }
}

fn wide(z: Z) -> (i128, i128) {
    (z.0 as i128, z.1 as i128)
}

fn spell(ring: Ring, z: Z) -> String {
    let (re, im, unit) = match ring {
        Ring::Gaussian => (z.0, z.1, "i"),
        Ring::Eisenstein => (z.0 - z.1, z.1, "w"),
    };
    match (re, im) {
        (a, 0) => format!("{a}"),
        (0, 1) => unit.to_string(),
        (0, -1) => format!("-{unit}"),
        (0, b) => format!("{b}{unit}"),
        (a, 1) => format!("{a}+{unit}"),
        (a, -1) => format!("{a}-{unit}"),
        (a, b) if b > 0 => format!("{a}+{b}{unit}"),
        (a, b) => format!("{a}{b}{unit}"),
    }
}

fn ring_word(ring: Ring) -> &'static str {
    match ring {
        Ring::Gaussian => "Z[i]",
        Ring::Eisenstein => "Z[w]",
    }
}

fn spell_list(ring: Ring, list: &[Z]) -> String {
    list.iter()
        .map(|&z| spell(ring, z))
        .collect::<Vec<_>>()
        .join(" ")
}

fn spell_turns(turns: &[usize]) -> String {
    turns
        .iter()
        .map(|k| k.to_string())
        .collect::<Vec<_>>()
        .join("")
}

// DESIGN

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
struct Design {
    digits: Vec<Z>,
    twists: Vec<usize>,
}

impl Design {
    fn radix(&self, lab: &Lab) -> Radix {
        let twists = self.twists.iter().map(|&k| lab.units[k]).collect();
        Radix::new(lab.base, self.digits.clone(), twists).unwrap()
    }
    fn size(&self) -> usize {
        self.digits.len()
    }
}

fn canonical(lab: &Lab, code: usize, twists: &[usize]) -> Design {
    let digits: Vec<Z> = (0..lab.residues.len())
        .filter(|i| (code >> i) & 1 == 1)
        .map(|i| lab.residues[i])
        .collect();
    Design {
        digits,
        twists: twists.to_vec(),
    }
}

fn walk_design(lab: &Lab, walk: &[usize]) -> Design {
    let mut digits = Vec::with_capacity(walk.len());
    let mut at = (0, 0);
    for &k in walk {
        digits.push(at);
        at = add(at, lab.units[k]);
    }
    Design {
        digits,
        twists: walk.to_vec(),
    }
}

fn points(lab: &Lab, design: &Design, level: usize) -> Vec<Z> {
    let mut out = vec![(0i64, 0i64)];
    for step in 1..=level {
        let shift = lab.power(step - 1);
        let mut next = Vec::with_capacity(out.len() * design.size());
        for (&d, &t) in design.digits.iter().zip(design.twists.iter()) {
            let head = lab.mul(d, shift);
            for &p in &out {
                next.push(add(head, lab.turn(p, t)));
            }
        }
        out = next;
    }
    out
}

fn incongruent(lab: &Lab, digits: &[Z]) -> bool {
    (0..digits.len()).all(|i| (0..i).all(|j| !lab.base.congruent(digits[i], digits[j])))
}

// STEP

fn junction_key(lab: &Lab, design: &Design, a: usize) -> Z {
    let k = design.size();
    let (d0, dk) = (design.digits[0], design.digits[k - 1]);
    let f0 = sub(lab.b, lab.units[design.twists[0]]);
    let fk = sub(lab.b, lab.units[design.twists[k - 1]]);
    let gap = sub(design.digits[a + 1], design.digits[a]);
    let one = lab.mul(lab.mul(gap, f0), fk);
    let two = lab.mul(lab.mul(lab.units[design.twists[a + 1]], d0), fk);
    let three = lab.mul(lab.mul(lab.units[design.twists[a]], dk), f0);
    sub(add(one, two), three)
}

fn joints(lab: &Lab, design: &Design, top: usize) -> Vec<Vec<(i128, i128)>> {
    let ring = lab.ring;
    let k = design.size();
    let (d0, dk) = (wide(design.digits[0]), wide(design.digits[k - 1]));
    let (u0, uk) = (
        wide(lab.units[design.twists[0]]),
        wide(lab.units[design.twists[k - 1]]),
    );
    let b = wide(lab.b);
    let mut out = vec![Vec::with_capacity(top); k - 1];
    let (mut first, mut last, mut scale) = ((0i128, 0i128), (0i128, 0i128), (1i128, 0i128));
    for _ in 1..=top {
        for (a, row) in out.iter_mut().enumerate() {
            let gap = wide(sub(design.digits[a + 1], design.digits[a]));
            let one = mul_wide(ring, gap, scale);
            let two = mul_wide(ring, wide(lab.units[design.twists[a + 1]]), first);
            let three = mul_wide(ring, wide(lab.units[design.twists[a]]), last);
            row.push((one.0 + two.0 - three.0, one.1 + two.1 - three.1));
        }
        let f = mul_wide(ring, u0, first);
        let l = mul_wide(ring, uk, last);
        let df = mul_wide(ring, d0, scale);
        let dl = mul_wide(ring, dk, scale);
        first = (df.0 + f.0, df.1 + f.1);
        last = (dl.0 + l.0, dl.1 + l.1);
        scale = mul_wide(ring, scale, b);
    }
    out
}

fn steps_by_level(lab: &Lab, design: &Design, top: usize) -> Vec<bool> {
    let table = joints(lab, design, top);
    (0..top)
        .map(|l| table.iter().all(|row| norm_wide(lab.ring, row[l]) == 1))
        .collect()
}

fn lemma_step(lab: &Lab, design: &Design) -> bool {
    let k = design.size();
    k >= 2
        && (0..k - 1).all(|a| junction_key(lab, design, a) == (0, 0))
        && steps_by_level(lab, design, lab.n).iter().all(|&s| s)
}

fn first_break(lab: &Lab, design: &Design) -> Option<usize> {
    if lemma_step(lab, design) {
        return None;
    }
    let top = 36;
    let found = steps_by_level(lab, design, top).iter().position(|&s| !s);
    assert!(
        found.is_some(),
        "a design outside the lemma kept unit steps to level {top}"
    );
    found.map(|l| l + 1)
}

fn brute_steps(lab: &Lab, design: &Design, level: usize) -> bool {
    let p = points(lab, design, level);
    p.windows(2).all(|w| lab.norm(sub(w[1], w[0])) == 1)
}

// MACHINE

struct Machine {
    points: Vec<Z>,
    turn: Vec<u32>,
    near: Vec<bool>,
    zero: u32,
    k: usize,
    n: usize,
    head: Vec<u32>,
    edges: Vec<(u8, u8, u32)>,
    gaps: Vec<Option<u32>>,
    ends: Vec<Option<u32>>,
}

fn machine(lab: &Lab, digits: &[Z]) -> Machine {
    let n = lab.n;
    let k = digits.len();
    let reach = digits.iter().map(|&d| lab.norm(d)).max().unwrap_or(0) as f64;
    let rho = 2.0 * reach.sqrt() / (lab.modulus() - 1.0);
    let cap = (rho * rho + 1e-9).max(1.0);
    let side = (2.0 * cap).sqrt().ceil() as i64 + 1;
    let mut points = Vec::new();
    for a in -side..=side {
        for b in -side..=side {
            let norm = lab.norm((a, b));
            if norm <= 1 || (norm as f64) < cap {
                points.push((a, b));
            }
        }
    }
    let width = 2 * side + 1;
    let mut grid = vec![u32::MAX; (width * width) as usize];
    for (i, &(a, b)) in points.iter().enumerate() {
        grid[((a + side) * width + b + side) as usize] = i as u32;
    }
    let find = |z: Z| -> Option<u32> {
        if z.0.abs() > side || z.1.abs() > side {
            return None;
        }
        let i = grid[((z.0 + side) * width + z.1 + side) as usize];
        (i != u32::MAX).then_some(i)
    };
    let mut turn = Vec::with_capacity(points.len() * n);
    for &z in &points {
        for j in 0..n {
            turn.push(find(lab.turn(z, j)).unwrap());
        }
    }
    let near = points.iter().map(|&z| lab.norm(z) <= 1).collect();
    let zero = find((0, 0)).unwrap();
    let mut head = vec![0u32];
    let mut edges = Vec::new();
    for &delta in &points {
        let scaled = lab.mul(lab.b, delta);
        for r in 0..n {
            for (e, &de) in digits.iter().enumerate() {
                for (f, &df) in digits.iter().enumerate() {
                    let x = sub(add(scaled, de), lab.turn(df, r));
                    if let Some(i) = find(x) {
                        edges.push((e as u8, f as u8, i));
                    }
                }
            }
            head.push(edges.len() as u32);
        }
    }
    let mut gaps = vec![None; k * k];
    for d in 0..k {
        for f in 0..k {
            if d != f {
                gaps[d * k + f] = find(sub(digits[d], digits[f]));
            }
        }
    }
    let ends = digits.iter().map(|&d| find(sub(lab.b, d))).collect();
    Machine {
        points,
        turn,
        near,
        zero,
        k,
        n,
        head,
        edges,
        gaps,
        ends,
    }
}

impl Machine {
    fn states(&self) -> usize {
        self.points.len() * self.n
    }
    fn back(&self, k: usize) -> usize {
        (self.n - k % self.n) % self.n
    }
    fn hop(&self, t: &[usize], s: usize, edge: (u8, u8, u32)) -> usize {
        let (e, f, x) = (edge.0 as usize, edge.1 as usize, edge.2 as usize);
        let r = s % self.n;
        let p = self.turn[x * self.n + self.back(t[e])] as usize;
        p * self.n + (r + t[f] + self.n - t[e]) % self.n
    }
    fn start(&self, t: &[usize], d: usize, f: usize) -> Option<usize> {
        self.gaps[d * self.k + f].map(|g| {
            let p = self.turn[g as usize * self.n + self.back(t[d])] as usize;
            p * self.n + (t[f] + self.n - t[d]) % self.n
        })
    }
    fn end_start(&self, t: &[usize], f: usize) -> Option<usize> {
        self.ends[f].map(|g| g as usize * self.n + t[f] % self.n)
    }
    fn out(&self, s: usize) -> &[(u8, u8, u32)] {
        &self.edges[self.head[s] as usize..self.head[s + 1] as usize]
    }
    fn search(&self, t: &[usize], seeds: &[usize], maps: bool, tail: bool) -> Option<usize> {
        let mut dist = vec![u32::MAX; self.states()];
        let mut queue = std::collections::VecDeque::new();
        for &s in seeds {
            if dist[s] == u32::MAX {
                dist[s] = 0;
                queue.push_back(s);
            }
        }
        while let Some(s) = queue.pop_front() {
            let p = s / self.n;
            if p as u32 == self.zero && (!maps || s % self.n == 0) {
                return Some(dist[s] as usize + 1);
            }
            for &edge in self.out(s) {
                if tail && edge.0 != 0 {
                    continue;
                }
                let next = self.hop(t, s, edge);
                if dist[next] == u32::MAX {
                    dist[next] = dist[s] + 1;
                    queue.push_back(next);
                }
            }
        }
        None
    }
    fn glue(&self, t: &[usize], maps: bool) -> Option<usize> {
        let mut seeds = Vec::new();
        for d in 0..self.k {
            for f in d + 1..self.k {
                if let Some(s) = self.start(t, d, f) {
                    seeds.push(s);
                }
            }
        }
        self.search(t, &seeds, maps, false)
    }
    fn ending(&self, t: &[usize]) -> Option<usize> {
        let seeds: Vec<usize> = (0..self.k).filter_map(|f| self.end_start(t, f)).collect();
        self.search(t, &seeds, false, true)
    }
}

fn linked(k: usize, edges: impl Iterator<Item = (usize, usize)>) -> bool {
    let mut root: Vec<usize> = (0..k).collect();
    fn find(root: &mut [usize], mut x: usize) -> usize {
        while root[x] != x {
            root[x] = root[root[x]];
            x = root[x];
        }
        x
    }
    let mut parts = k;
    for (a, b) in edges {
        let (x, y) = (find(&mut root, a), find(&mut root, b));
        if x != y {
            root[x] = y;
            parts -= 1;
        }
    }
    parts <= 1
}

impl Machine {
    fn piece(&self, t: &[usize]) -> Option<usize> {
        let k = self.k;
        if k <= 1 {
            return None;
        }
        let mut pairs = Vec::new();
        for d in 0..k {
            for f in d + 1..k {
                if let Some(s) = self.start(t, d, f) {
                    pairs.push((d, f, s));
                }
            }
        }
        let mut local = vec![u32::MAX; self.states()];
        let mut order: Vec<usize> = Vec::new();
        for &(_, _, s) in &pairs {
            if local[s] == u32::MAX {
                local[s] = order.len() as u32;
                order.push(s);
            }
        }
        let mut head = vec![0u32];
        let mut succ: Vec<u32> = Vec::new();
        let mut i = 0;
        while i < order.len() {
            let s = order[i];
            for &edge in self.out(s) {
                let next = self.hop(t, s, edge);
                if local[next] == u32::MAX {
                    local[next] = order.len() as u32;
                    order.push(next);
                }
                succ.push(local[next]);
            }
            head.push(succ.len() as u32);
            i += 1;
        }
        let size = order.len();
        let words = size.div_ceil(64);
        let bit = |w: &[u64], j: usize| (w[j / 64] >> (j % 64)) & 1 == 1;
        let mut w = vec![0u64; words];
        for (j, &s) in order.iter().enumerate() {
            if self.near[s / self.n] {
                w[j / 64] |= 1 << (j % 64);
            }
        }
        let mut history: Vec<Vec<u64>> = Vec::new();
        for level in 1.. {
            let edges = pairs
                .iter()
                .filter(|&&(_, _, s)| bit(&w, local[s] as usize))
                .map(|&(d, f, _)| (d, f));
            if !linked(k, edges) {
                return Some(level);
            }
            history.push(w.clone());
            let mut next = vec![0u64; words];
            for j in 0..size {
                let list = &succ[head[j] as usize..head[j + 1] as usize];
                if list.iter().any(|&x| bit(&w, x as usize)) {
                    next[j / 64] |= 1 << (j % 64);
                }
            }
            if history.contains(&next) {
                return None;
            }
            assert!(history.len() < 4096, "the adjacency sets have not cycled");
            w = next;
        }
        None
    }
}

fn brute_piece(lab: &Lab, design: &Design, level: usize) -> bool {
    let set: HashSet<Z> = points(lab, design, level).into_iter().collect();
    let start = match set.iter().next() {
        Some(&z) => z,
        None => return true,
    };
    let mut seen = HashSet::from([start]);
    let mut stack = vec![start];
    while let Some(z) = stack.pop() {
        for &u in &lab.units {
            let w = add(z, u);
            if set.contains(&w) && seen.insert(w) {
                stack.push(w);
            }
        }
    }
    seen.len() == set.len()
}

fn brute_glue(lab: &Lab, design: &Design, level: usize) -> (bool, bool) {
    let mut place = HashSet::new();
    let mut map = HashSet::new();
    let mut word_twist = vec![0usize];
    for _ in 0..level {
        let mut next = Vec::with_capacity(word_twist.len() * design.size());
        for &t in &design.twists {
            for &w in &word_twist {
                next.push((t + w) % lab.n);
            }
        }
        word_twist = next;
    }
    let p = points(lab, design, level);
    for (z, &u) in p.iter().zip(word_twist.iter()) {
        place.insert(*z);
        map.insert((*z, u));
    }
    (place.len() < p.len(), map.len() < p.len())
}

// GROUP

#[derive(Clone)]
struct Symmetry {
    perm: Vec<usize>,
    twist: Vec<usize>,
    mirror: bool,
    unit: usize,
}

fn residue_index(lab: &Lab) -> HashMap<Z, usize> {
    lab.residues
        .iter()
        .enumerate()
        .map(|(i, &z)| (z, i))
        .collect()
}

fn stabiliser(lab: &Lab) -> Vec<Symmetry> {
    let index = residue_index(lab);
    let mut out = Vec::new();
    for unit in 0..lab.n {
        let perm: Option<Vec<usize>> = lab
            .residues
            .iter()
            .map(|&z| index.get(&lab.turn(z, unit)).copied())
            .collect();
        if let Some(perm) = perm {
            out.push(Symmetry {
                perm,
                twist: (0..lab.n).collect(),
                mirror: false,
                unit,
            });
        }
    }
    if let Some(eps) = lab.mirror() {
        for unit in 0..lab.n {
            let perm: Option<Vec<usize>> = lab
                .residues
                .iter()
                .map(|&z| index.get(&lab.turn(lab.conj(z), unit)).copied())
                .collect();
            if let Some(perm) = perm {
                let twist = (0..lab.n)
                    .map(|k| (lab.back(k) + lab.back(eps)) % lab.n)
                    .collect();
                out.push(Symmetry {
                    perm,
                    twist,
                    mirror: true,
                    unit,
                });
            }
        }
    }
    out
}

fn act(sym: &Symmetry, code: usize, twists: &[usize]) -> (usize, Vec<usize>) {
    let q = sym.perm.len();
    let mut slot = vec![usize::MAX; q];
    let mut j = 0;
    let mut image = 0usize;
    for c in 0..q {
        if (code >> c) & 1 == 1 {
            slot[sym.perm[c]] = sym.twist[twists[j]];
            image |= 1 << sym.perm[c];
            j += 1;
        }
    }
    (
        image,
        slot.into_iter().filter(|&t| t != usize::MAX).collect(),
    )
}

fn pack(n: usize, twists: &[usize]) -> u64 {
    twists
        .iter()
        .rev()
        .fold(0u64, |acc, &t| acc * n as u64 + t as u64)
}

fn unpack(n: usize, k: usize, mut index: u64) -> Vec<usize> {
    let mut out = Vec::with_capacity(k);
    for _ in 0..k {
        out.push((index % n as u64) as usize);
        index /= n as u64;
    }
    out
}

fn key(n: usize, code: usize, twists: &[usize]) -> u64 {
    ((code as u64) << 40) | pack(n, twists)
}

fn fixed(lab: &Lab, sym: &Symmetry) -> u64 {
    let q = lab.residues.len();
    let mut total = 0u64;
    for code in 0..1usize << q {
        if code.count_ones() < 2 || act(sym, code, &vec![0; code.count_ones() as usize]).0 != code {
            continue;
        }
        let mut seen = 0usize;
        let mut count = 1u64;
        for c in 0..q {
            if (code >> c) & 1 == 0 || (seen >> c) & 1 == 1 {
                continue;
            }
            let mut len = 0;
            let mut at = c;
            while (seen >> at) & 1 == 0 {
                seen |= 1 << at;
                at = sym.perm[at];
                len += 1;
            }
            let stay = (0..lab.n)
                .filter(|&t| {
                    let mut x = t;
                    for _ in 0..len {
                        x = sym.twist[x];
                    }
                    x == t
                })
                .count();
            count *= stay as u64;
        }
        total += count;
    }
    total
}

fn burnside(lab: &Lab, group: &[Symmetry]) -> u64 {
    let total: u64 = group.iter().map(|sym| fixed(lab, sym)).sum();
    assert_eq!(total % group.len() as u64, 0, "Burnside is not an integer");
    total / group.len() as u64
}

fn same_up_to_turn(lab: &Lab, left: &[Z], right: &[Z]) -> bool {
    let target: HashSet<Z> = right.iter().copied().collect();
    (0..lab.n).any(|v| left.iter().all(|&z| target.contains(&lab.turn(z, v))))
        && left.len() == right.len()
}

// CENSUS

#[derive(Clone, Default)]
struct Cell {
    pairs: u64,
    orbits: u64,
    step: u64,
    curve: u64,
    piece: u64,
    piece_orbits: u64,
    plane: u64,
    plane_orbits: u64,
    distinct: u64,
    distinct_orbits: u64,
}

impl Cell {
    fn merge(&mut self, o: &Cell) {
        self.pairs += o.pairs;
        self.orbits += o.orbits;
        self.step += o.step;
        self.curve += o.curve;
        self.piece += o.piece;
        self.piece_orbits += o.piece_orbits;
        self.plane += o.plane;
        self.plane_orbits += o.plane_orbits;
        self.distinct += o.distinct;
        self.distinct_orbits += o.distinct_orbits;
    }
}

type Mark = (usize, usize, Vec<usize>);

#[derive(Clone, Default)]
struct Survey {
    cells: Vec<Cell>,
    steps: Vec<(usize, Vec<usize>, Option<usize>)>,
    held: Mark,
    split: Mark,
    overlap: Mark,
}

fn raise(mark: &mut Mark, level: usize, code: usize, twists: &[usize]) {
    if level > mark.0 || (level == mark.0 && (code, twists) < (mark.1, &mark.2[..])) {
        *mark = (level, code, twists.to_vec());
    }
}

fn merge_mark(mark: &mut Mark, other: &Mark) {
    if other.0 > 0 {
        raise(mark, other.0, other.1, &other.2);
    }
}

fn unit_path(lab: &Lab, digits: &[Z]) -> bool {
    digits.windows(2).all(|w| lab.norm(sub(w[1], w[0])) == 1)
}

fn unit_linked(lab: &Lab, digits: &[Z]) -> bool {
    let k = digits.len();
    let mut edges = Vec::new();
    for a in 0..k {
        for b in a + 1..k {
            if lab.norm(sub(digits[a], digits[b])) == 1 {
                edges.push((a, b));
            }
        }
    }
    linked(k, edges.into_iter())
}

fn survey_chunk(lab: &Lab, group: &[Symmetry], code: usize, lo: u64, hi: u64) -> Survey {
    let q = lab.residues.len();
    let n = lab.n;
    let mut design = canonical(lab, code, &[]);
    let k = design.size();
    let m = machine(lab, &design.digits);
    let path = unit_path(lab, &design.digits);
    let joined = unit_linked(lab, &design.digits);
    let mut out = Survey {
        cells: vec![Cell::default(); q + 1],
        ..Survey::default()
    };
    let cell = &mut out.cells[k];
    for index in lo..hi {
        let t = unpack(n, k, index);
        cell.pairs += 1;
        design.twists = t.clone();
        if path {
            if lemma_step(lab, &design) {
                assert_eq!(
                    lab.sum(&t),
                    lab.b,
                    "a unit-step pair whose twists do not sum to the base"
                );
                cell.step += 1;
                let glue = m.glue(&t, false);
                if glue.is_none() {
                    cell.curve += 1;
                }
                out.steps.push((code, t.clone(), glue));
            } else {
                let level = first_break(lab, &design).unwrap();
                raise(&mut out.held, level - 1, code, &t);
            }
        }
        let own = key(n, code, &t);
        let mut images = Vec::with_capacity(group.len());
        let mut rep = true;
        for sym in group {
            let (c, u) = act(sym, code, &t);
            let image = key(n, c, &u);
            if image < own {
                rep = false;
                break;
            }
            images.push(image);
        }
        if !rep {
            continue;
        }
        images.sort_unstable();
        images.dedup();
        let size = images.len() as u64;
        cell.orbits += 1;
        if joined {
            match m.piece(&t) {
                None => {
                    cell.piece += size;
                    cell.piece_orbits += 1;
                }
                Some(level) => raise(&mut out.split, level, code, &t),
            }
        }
        if k == q {
            match m.glue(&t, true) {
                None => {
                    cell.plane += size;
                    cell.plane_orbits += 1;
                }
                Some(level) => raise(&mut out.overlap, level, code, &t),
            }
            if m.glue(&t, false).is_none() {
                cell.distinct += size;
                cell.distinct_orbits += 1;
            }
        }
    }
    out
}

fn survey(lab: &Lab, group: &[Symmetry]) -> Survey {
    let q = lab.residues.len();
    let n = lab.n as u64;
    let chunk = 1u64 << 15;
    let mut work = Vec::new();
    for code in 0..1usize << q {
        let k = code.count_ones();
        if k < 2 {
            continue;
        }
        let total = n.pow(k);
        let mut lo = 0;
        while lo < total {
            work.push((code, lo, (lo + chunk).min(total)));
            lo += chunk;
        }
    }
    let next = std::sync::atomic::AtomicUsize::new(0);
    let threads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(4);
    let parts: Vec<Survey> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                scope.spawn(|| {
                    let mut acc = Survey {
                        cells: vec![Cell::default(); q + 1],
                        ..Survey::default()
                    };
                    loop {
                        let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if i >= work.len() {
                            break;
                        }
                        let (code, lo, hi) = work[i];
                        let part = survey_chunk(lab, group, code, lo, hi);
                        for (a, b) in acc.cells.iter_mut().zip(part.cells.iter()) {
                            a.merge(b);
                        }
                        acc.steps.extend(part.steps);
                        merge_mark(&mut acc.held, &part.held);
                        merge_mark(&mut acc.split, &part.split);
                        merge_mark(&mut acc.overlap, &part.overlap);
                    }
                    acc
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut out = Survey {
        cells: vec![Cell::default(); q + 1],
        ..Survey::default()
    };
    for part in parts {
        for (a, b) in out.cells.iter_mut().zip(part.cells.iter()) {
            a.merge(b);
        }
        out.steps.extend(part.steps);
        merge_mark(&mut out.held, &part.held);
        merge_mark(&mut out.split, &part.split);
        merge_mark(&mut out.overlap, &part.overlap);
    }
    out.steps.sort();
    out
}

fn code_action(lab: &Lab) -> Vec<Symmetry> {
    let perms = lab.base.group().unwrap();
    let eps = lab.mirror();
    let mut out = Vec::new();
    let mut i = 0;
    for unit in 0..lab.n {
        for flip in [false, true] {
            if flip && eps.is_none() {
                continue;
            }
            let twist = match (flip, eps) {
                (true, Some(e)) => (0..lab.n)
                    .map(|k| (lab.back(k) + lab.back(e)) % lab.n)
                    .collect(),
                _ => (0..lab.n).collect(),
            };
            out.push(Symmetry {
                perm: perms[i].clone(),
                twist,
                mirror: flip,
                unit,
            });
            i += 1;
        }
    }
    out
}

fn spell_symmetry(lab: &Lab, sym: &Symmetry) -> String {
    let v = spell(lab.ring, lab.units[sym.unit]);
    if sym.mirror {
        format!("x -> {v} conj(x)")
    } else {
        format!("x -> {v} x")
    }
}

type Breach = (usize, Vec<usize>, usize, Vec<usize>, String);

fn breach(lab: &Lab, action: &[Symmetry]) -> Option<Breach> {
    let q = lab.residues.len();
    for code in 0..1usize << q {
        let k = code.count_ones() as usize;
        if k < 2 {
            continue;
        }
        let m = machine(lab, &canonical(lab, code, &[]).digits);
        for index in 0..(lab.n as u64).pow(k as u32) {
            let t = unpack(lab.n, k, index);
            let (piece, glue) = (m.piece(&t), m.glue(&t, false));
            for sym in action {
                let (c, u) = act(sym, code, &t);
                let other = machine(lab, &canonical(lab, c, &[]).digits);
                let (p2, g2) = (other.piece(&u), other.glue(&u, false));
                let by = spell_symmetry(lab, sym);
                if piece.is_none() != p2.is_none() {
                    return Some((
                        code,
                        t,
                        c,
                        u,
                        format!(
                            "by {by}, one piece {} against {}",
                            piece.is_none(),
                            p2.is_none()
                        ),
                    ));
                }
                if glue.is_none() != g2.is_none() {
                    return Some((
                        code,
                        t,
                        c,
                        u,
                        format!(
                            "by {by}, no two words on one point {} against {}",
                            glue.is_none(),
                            g2.is_none()
                        ),
                    ));
                }
            }
        }
    }
    None
}

fn run_census() {
    println!("CENSUS  canonical pairs (code, twist), |F| >= 2, digits in canonical order, every property decided at every level");
    println!("  step: unit steps between consecutive words; curve: step and no two words on one point; piece: connected under unit adjacency; plane: |F| = N(b) and no two words with one place map; distinct: |F| = N(b) and no two words on one point");
    for (ring, value) in BASES {
        let clock = Instant::now();
        let lab = Lab::new(ring, value);
        let q = lab.residues.len();
        let n = lab.n as u64;
        let group = stabiliser(&lab);
        let action = code_action(&lab);
        let s = survey(&lab, &group);
        let mut total = Cell::default();
        for cell in &s.cells {
            total.merge(cell);
        }
        let expect = (n + 1).pow(q as u32) - 1 - q as u64 * n;
        assert_eq!(
            total.pairs, expect,
            "the pair count is not sum_(k>=2) binom(q,k) n^k"
        );
        assert_eq!(
            total.orbits,
            burnside(&lab, &group),
            "orbit walk and Burnside disagree"
        );
        println!(
            "  {} base {}  N(b) {q}  residues {}",
            ring_word(ring),
            spell(ring, value),
            spell_list(ring, &lab.residues)
        );
        println!(
            "    stabiliser of the residue system, order {}: {}",
            group.len(),
            group
                .iter()
                .map(|g| spell_symmetry(&lab, g))
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!(
            "    cut: {} pairs, {} orbits under the stabiliser; the code action of order {} would give {} orbits",
            total.pairs,
            total.orbits,
            action.len(),
            burnside(&lab, &action)
        );
        println!("    |F|  pairs  orbits  step  curve  piece  piece-orbits  plane  plane-orbits  distinct  distinct-orbits");
        for (k, c) in s.cells.iter().enumerate() {
            if c.pairs == 0 {
                continue;
            }
            println!(
                "    {k}  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}",
                c.pairs,
                c.orbits,
                c.step,
                c.curve,
                c.piece,
                c.piece_orbits,
                c.plane,
                c.plane_orbits,
                c.distinct,
                c.distinct_orbits
            );
        }
        println!(
            "    all  {}  {}  {}  {}  {}  {}  {}  {}  {}  {}",
            total.pairs,
            total.orbits,
            total.step,
            total.curve,
            total.piece,
            total.piece_orbits,
            total.plane,
            total.plane_orbits,
            total.distinct,
            total.distinct_orbits
        );
        let show = |m: &Mark| format!("level {} at code {} twists {}", m.0, m.1, spell_turns(&m.2));
        println!(
            "    deepest unit steps outside the lemma: {}",
            show(&s.held)
        );
        println!("    deepest first split: {}", show(&s.split));
        println!("    deepest first shared place map: {}", show(&s.overlap));
        for (code, t, glue) in &s.steps {
            let d = canonical(&lab, *code, t);
            println!(
                "    step pair: code {code} digits {} twists {} first glue {}",
                spell_list(ring, &d.digits),
                spell_turns(t),
                glue.map_or("never".to_string(), |g| format!("level {g}"))
            );
        }
        if q <= 5 {
            match breach(&lab, &action) {
                Some((c, t, c2, u, why)) => println!(
                    "    the code action is not a symmetry: code {c} twists {} goes to code {c2} twists {}, {why}",
                    spell_turns(&t),
                    spell_turns(&u)
                ),
                None => println!("    the code action keeps every class on this base"),
            }
        }
        println!("    {:.2} s", clock.elapsed().as_secs_f64());
    }
}

fn main() {
    let verbs: Vec<String> = std::env::args().skip(1).collect();
    let verbs = if verbs.is_empty() {
        ["census"].iter().map(|s| s.to_string()).collect()
    } else {
        verbs
    };
    for verb in verbs {
        let clock = Instant::now();
        match verb.as_str() {
            "census" => run_census(),
            "curves" => run_curves(),
            "junction" => run_junction(),
            "named" => run_named(),
            other => panic!("unknown verb {other}"),
        }
        println!("  {verb} {:.2} s", clock.elapsed().as_secs_f64());
    }
}

// WALKS

#[derive(Clone)]
struct Walk {
    turns: Vec<usize>,
    radix: bool,
    glue: Option<usize>,
    end: Option<usize>,
    plane: bool,
    canonical: Option<usize>,
}

impl Walk {
    fn crossing(&self) -> Option<usize> {
        match (self.glue, self.end) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (a, b) => a.or(b),
        }
    }
}

fn walks(lab: &Lab, k: usize) -> Vec<Vec<usize>> {
    fn grow(lab: &Lab, k: usize, at: Z, turns: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
        let left = (k - turns.len()) as u64;
        if lab.norm(sub(lab.b, at)) > left * left {
            return;
        }
        if left == 0 {
            out.push(turns.clone());
            return;
        }
        for t in 0..lab.n {
            turns.push(t);
            grow(lab, k, add(at, lab.units[t]), turns, out);
            turns.pop();
        }
    }
    let mut out = Vec::new();
    grow(lab, k, (0, 0), &mut Vec::new(), &mut out);
    out
}

fn images(lab: &Lab, turns: &[usize]) -> Vec<Vec<usize>> {
    let mut out = vec![turns.to_vec()];
    let mut back = turns.to_vec();
    back.reverse();
    out.push(back);
    if let Some(eps) = lab.mirror() {
        for w in out.clone() {
            out.push(
                w.iter()
                    .map(|&t| (lab.back(t) + lab.back(eps)) % lab.n)
                    .collect(),
            );
        }
    }
    out
}

fn class_of(lab: &Lab, turns: &[usize]) -> Vec<usize> {
    images(lab, turns).into_iter().min().unwrap()
}

fn seen_by_census(lab: &Lab, design: &Design) -> Option<usize> {
    let index = residue_index(lab);
    for z in 0..lab.n {
        let classes: Option<Vec<usize>> = design
            .digits
            .iter()
            .map(|&d| index.get(&lab.turn(d, z)).copied())
            .collect();
        if let Some(c) = classes {
            if c.windows(2).all(|w| w[0] < w[1]) {
                return Some(c.iter().fold(0, |acc, &i| acc | 1 << i));
            }
        }
    }
    None
}

fn classify(lab: &Lab, turns: &[usize]) -> Walk {
    let design = walk_design(lab, turns);
    assert!(lemma_step(lab, &design), "a walk outside the step lemma");
    let m = machine(lab, &design.digits);
    let plane = turns.len() == lab.residues.len() && m.glue(turns, true).is_none();
    Walk {
        turns: turns.to_vec(),
        radix: incongruent(lab, &design.digits),
        glue: m.glue(turns, false),
        end: m.ending(turns),
        plane,
        canonical: seen_by_census(lab, &design),
    }
}

fn classify_all(lab: &Lab, list: &[Vec<usize>]) -> Vec<Walk> {
    let next = std::sync::atomic::AtomicUsize::new(0);
    let threads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(4);
    let chunk = 256;
    let mut parts: Vec<Vec<(usize, Walk)>> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                scope.spawn(|| {
                    let mut acc = Vec::new();
                    loop {
                        let i = next.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                        if i * chunk >= list.len() {
                            break;
                        }
                        for (j, item) in list.iter().enumerate().skip(i * chunk).take(chunk) {
                            acc.push((j, classify(lab, item)));
                        }
                    }
                    acc
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let mut all: Vec<(usize, Walk)> = parts.drain(..).flatten().collect();
    all.sort_by_key(|x| x.0);
    all.into_iter().map(|x| x.1).collect()
}

fn run_curves() {
    println!("CURVES  every walk of k unit steps from 0 to the base: digits the partial sums, twists the steps; every such design has unit steps at every level");
    println!("  arc: the k^L + 1 vertices of level L distinct at every level; plane: k = N(b) and no two words with one place map; classes up to reversal and, where conj(b) is an associate, the mirror");
    for (ring, value) in BASES {
        let clock = Instant::now();
        let lab = Lab::new(ring, value);
        let q = lab.residues.len();
        println!(
            "  {} base {}  N(b) {q}",
            ring_word(ring),
            spell(ring, value)
        );
        println!("    k  dim  walks  radix  arcs  arc-classes  plane  plane-arcs  canonical  deepest-crossing");
        let mut listed = Vec::new();
        for k in 2..=q {
            let list = walks(&lab, k);
            let all = classify_all(&lab, &list);
            let radix: Vec<&Walk> = all.iter().filter(|w| w.radix).collect();
            let arcs: Vec<&&Walk> = radix.iter().filter(|w| w.crossing().is_none()).collect();
            let classes: HashSet<Vec<usize>> =
                arcs.iter().map(|w| class_of(&lab, &w.turns)).collect();
            let plane = radix.iter().filter(|w| w.plane).count();
            let plane_arcs = arcs.iter().filter(|w| w.plane).count();
            let canonical = radix.iter().filter(|w| w.canonical.is_some()).count();
            let mut deep = (0usize, Vec::new());
            for w in &radix {
                if let Some(c) = w.crossing() {
                    if c > deep.0 || (c == deep.0 && w.turns < deep.1) {
                        deep = (c, w.turns.clone());
                    }
                }
            }
            let other_arcs = all
                .iter()
                .filter(|w| !w.radix && w.crossing().is_none())
                .count();
            let widest = all.iter().filter_map(|w| w.crossing()).max().unwrap_or(0);
            println!(
                "    {k}  {:.4}  {}  {}  {}  {}  {}  {}  {}  level {} at {}  (over every walk: {other_arcs} more arcs, deepest first revisit level {widest})",
                2.0 * (k as f64).ln() / (q as f64).ln(),
                all.len(),
                radix.len(),
                arcs.len(),
                classes.len(),
                plane,
                plane_arcs,
                canonical,
                deep.0,
                spell_turns(&deep.1)
            );
            let mut reps: Vec<Vec<usize>> = classes.into_iter().collect();
            reps.sort();
            for r in reps {
                let w = classify(&lab, &r);
                listed.push((k, r, w.plane, w.canonical));
            }
        }
        for (k, r, plane, canonical) in &listed {
            println!(
                "    arc class k {k}: turns {}{}{}",
                spell_turns(r),
                if *plane { ", plane-filling" } else { "" },
                canonical.map_or(String::new(), |c| format!(
                    ", canonical code {c} after a turn"
                ))
            );
        }
        println!("    {:.2} s", clock.elapsed().as_secs_f64());
    }
}

// JUNCTION

const BRUTE: usize = 4096;

#[derive(Default)]
struct Tests {
    designs: u64,
    checks: u64,
    misses: u64,
}

fn brute_top(k: usize) -> usize {
    let mut top = 1;
    while k.pow(top as u32 + 1) <= BRUTE {
        top += 1;
    }
    top.max(2)
}

fn brute_arc(lab: &Lab, design: &Design, level: usize) -> bool {
    let mut p = points(lab, design, level);
    p.push(lab.power(level));
    p.iter().collect::<HashSet<_>>().len() == p.len()
}

fn test_design(lab: &Lab, design: &Design, walk: bool, tests: &mut Tests) {
    let t = &design.twists;
    let m = machine(lab, &design.digits);
    let broken = first_break(lab, design);
    let (glue, maps, split) = (m.glue(t, false), m.glue(t, true), m.piece(t));
    let end = if walk { m.ending(t) } else { None };
    tests.designs += 1;
    for level in 1..=brute_top(design.size()) {
        let mut verdicts = vec![
            (
                brute_steps(lab, design, level),
                broken.is_none_or(|b| level < b),
            ),
            (
                brute_piece(lab, design, level),
                split.is_none_or(|s| level < s),
            ),
        ];
        let (glued, mapped) = brute_glue(lab, design, level);
        verdicts.push((glued, glue.is_some_and(|g| g <= level)));
        verdicts.push((mapped, maps.is_some_and(|g| g <= level)));
        if walk {
            let crossing = [glue, end].iter().flatten().min().copied();
            verdicts.push((
                brute_arc(lab, design, level),
                crossing.is_none_or(|c| level < c),
            ));
        }
        for (seen, said) in verdicts {
            tests.checks += 1;
            if seen != said {
                tests.misses += 1;
            }
        }
    }
    if incongruent(lab, &design.digits) {
        tests.checks += 1;
        if design.radix(lab).words(3) != points(lab, design, 3) {
            tests.misses += 1;
        }
    }
}

fn test_symmetry(lab: &Lab, group: &[Symmetry], code: usize, t: &[usize]) -> bool {
    let level = 3;
    let mine = points(lab, &canonical(lab, code, t), level);
    group.iter().all(|sym| {
        let (c, u) = act(sym, code, t);
        let theirs = points(lab, &canonical(lab, c, &u), level);
        let moved: Vec<Z> = if sym.mirror {
            mine.iter().map(|&z| lab.conj(z)).collect()
        } else {
            mine.clone()
        };
        same_up_to_turn(lab, &moved, &theirs)
    })
}

fn run_junction() {
    println!("JUNCTION  the step law, the difference automaton and the stabiliser against brute force over the words of each level");
    println!("  per design and per level up to |F|^level <= {BRUTE}: unit steps, one piece, two words on one point, two words with one place map, and for walks the arc with its end vertex");
    let mut total = Tests::default();
    for (ring, value) in BASES {
        let clock = Instant::now();
        let lab = Lab::new(ring, value);
        let q = lab.residues.len();
        let group = stabiliser(&lab);
        let mut all: Vec<(usize, u64)> = Vec::new();
        for code in 0..1usize << q {
            let k = code.count_ones();
            if k >= 2 {
                for index in 0..(lab.n as u64).pow(k) {
                    all.push((code, index));
                }
            }
        }
        let stride = (all.len() / 3000).max(1);
        let mut tests = Tests::default();
        let mut turned = 0u64;
        for (i, &(code, index)) in all.iter().enumerate() {
            let k = code.count_ones() as usize;
            let t = unpack(lab.n, k, index);
            let design = canonical(&lab, code, &t);
            let path = unit_path(&lab, &design.digits) && lemma_step(&lab, &design);
            if i % stride != 0 && !path {
                continue;
            }
            test_design(&lab, &design, false, &mut tests);
            if i % (stride * 4) == 0 {
                assert!(
                    test_symmetry(&lab, &group, code, &t),
                    "the stabiliser moved a level set"
                );
                turned += 1;
            }
        }
        let mut walked = 0;
        for k in 2..=q {
            for turns in walks(&lab, k) {
                let design = walk_design(&lab, &turns);
                if incongruent(&lab, &design.digits) {
                    test_design(&lab, &design, true, &mut tests);
                    walked += 1;
                }
            }
        }
        println!(
            "  {} base {}: {} canonical pairs of {} (every {stride}th and every step pair) and {walked} radix walks, {} level checks, {} misses; stabiliser images checked on {turned} pairs  {:.2} s",
            ring_word(ring),
            spell(ring, value),
            tests.designs - walked,
            all.len(),
            tests.checks,
            tests.misses,
            clock.elapsed().as_secs_f64()
        );
        assert_eq!(tests.misses, 0, "the lemma and the brute force disagree");
        total.designs += tests.designs;
        total.checks += tests.checks;
        total.misses += tests.misses;
    }
    println!(
        "  all bases: {} designs, {} checks, {} misses",
        total.designs, total.checks, total.misses
    );
    deep_witnesses();
}

fn deep_witnesses() {
    println!("  deep witnesses, brute force level by level");
    let lab = Lab::new(Ring::Gaussian, (1, 1));
    let held = canonical(&lab, 3, &[0, 2]);
    let row: Vec<bool> = (1..=5).map(|l| brute_steps(&lab, &held, l)).collect();
    println!(
        "    Z[i] base 1+i digits {} twists {}: unit steps at levels 1..5 {:?}, lemma first break level {:?}",
        spell_list(lab.ring, &held.digits),
        spell_turns(&held.twists),
        row,
        first_break(&lab, &held)
    );
    let lab = Lab::new(Ring::Gaussian, (2, 1));
    let point = Design {
        digits: vec![(0, 1), (1, 1), (1, 2)],
        twists: vec![0, 1, 2],
    };
    let fixed = point
        .digits
        .iter()
        .zip(point.twists.iter())
        .all(|(&d, &t)| lab.mul(sub(lab.b, lab.units[t]), (1, 1)) == add(d, d));
    println!(
        "    Z[i] base 2+i digits {} twists {}: incongruent {}, step law {}, unit steps at levels 1..8 {}, twists sum to {}, every map fixes (1+i)/2 {fixed}",
        spell_list(lab.ring, &point.digits),
        spell_turns(&point.twists),
        incongruent(&lab, &point.digits),
        lemma_step(&lab, &point),
        (1..=8).all(|l| brute_steps(&lab, &point, l)),
        spell(lab.ring, lab.sum(&point.twists))
    );
    let lab = Lab::new(Ring::Eisenstein, (2, 1));
    let gloss = canonical(&lab, 6, &[0, 1]);
    let level_two: Vec<u64> = points(&lab, &gloss, 2)
        .windows(2)
        .map(|w| lab.norm(sub(w[1], w[0])))
        .collect();
    println!(
        "    Z[w] base 1+w digits {} twists {}: junction identity {}, step norms at level 2 {:?}",
        spell_list(lab.ring, &gloss.digits),
        spell_turns(&gloss.twists),
        junction_key(&lab, &gloss, 0) == (0, 0),
        level_two
    );
    let middle = Design {
        digits: vec![(-1, 0), (0, 0), (1, 0)],
        twists: vec![0, 2, 0],
    };
    let m = machine(&lab, &middle.digits);
    println!(
        "    Z[w] base 1+w digits {} twists {}: incongruent {}, step law {}, twists sum to {}, first revisit {:?}, first shared map {:?}, first split {:?}, unit steps and distinct points by brute force at levels 1..8 {}",
        spell_list(lab.ring, &middle.digits),
        spell_turns(&middle.twists),
        incongruent(&lab, &middle.digits),
        lemma_step(&lab, &middle),
        spell(lab.ring, lab.sum(&middle.twists)),
        m.glue(&middle.twists, false),
        m.glue(&middle.twists, true),
        m.piece(&middle.twists),
        (1..=8).all(|l| brute_steps(&lab, &middle, l) && !brute_glue(&lab, &middle, l).0)
    );
    let walk = walk_design(&lab, &[0, 1]);
    let row: Vec<bool> = (1..=6).map(|l| brute_arc(&lab, &walk, l)).collect();
    let m = machine(&lab, &walk.digits);
    println!(
        "    Z[w] base {} walk 01: arc at levels 1..6 {:?}, automaton first revisit level {:?}",
        spell(lab.ring, lab.b),
        row,
        [m.glue(&walk.twists, false), m.ending(&walk.twists)]
            .iter()
            .flatten()
            .min()
    );
    for (ring, value, code, turns, top) in [
        (Ring::Gaussian, (2, 0), 14usize, vec![1usize, 3, 1], 6usize),
        (Ring::Gaussian, (2, 1), 15, vec![0, 1, 1, 0], 7),
        (Ring::Eisenstein, (2, 0), 11, vec![0, 3, 2], 10),
        (Ring::Eisenstein, (2, 1), 6, vec![2, 4], 8),
        (Ring::Eisenstein, (3, 0), 223, vec![0, 4, 3, 5, 2, 0, 4], 7),
        (Ring::Eisenstein, (3, 1), 63, vec![4, 1, 1, 5, 2, 0], 7),
    ] {
        let lab = Lab::new(ring, value);
        let design = canonical(&lab, code, &turns);
        let m = machine(&lab, &design.digits);
        let row: Vec<bool> = (1..=top).map(|l| brute_piece(&lab, &design, l)).collect();
        println!(
            "    {} base {} code {code} twists {}: one piece at levels 1..{top} {:?}, automaton first split level {:?}",
            ring_word(ring),
            spell(ring, value),
            spell_turns(&turns),
            row,
            m.piece(&turns)
        );
    }
}

// NAMED

fn as_design(lab: &Lab, radix: &Radix) -> Design {
    Design {
        digits: radix.digits().to_vec(),
        twists: radix
            .twists()
            .iter()
            .map(|&u| lab.units.iter().position(|&v| v == u).unwrap())
            .collect(),
    }
}

fn locate(name: &str, lab: &Lab, design: &Design) {
    let t = &design.twists;
    let m = machine(lab, &design.digits);
    let step = lemma_step(lab, design);
    let q = lab.residues.len();
    let plane = design.size() == q && m.glue(t, true).is_none();
    let distinct = m.glue(t, false).is_none();
    let piece = m.piece(t).is_none();
    let canonical = design.digits.iter().all(|d| lab.residues.contains(d));
    let walk = step && design.digits[0] == (0, 0);
    let crossing = if walk {
        [m.glue(t, false), m.ending(t)]
            .iter()
            .flatten()
            .min()
            .copied()
    } else {
        None
    };
    println!(
        "  {name}: {} base {}, digits {}, twists {}, dim {:.6}",
        ring_word(lab.ring),
        spell(lab.ring, lab.b),
        spell_list(lab.ring, &design.digits),
        spell_turns(t),
        2.0 * (design.size() as f64).ln() / (q as f64).ln()
    );
    println!(
        "    canonical digits {canonical}, seen by the canonical census {}, step {step}, piece {piece}, no two words on one point {distinct}, plane {plane}{}",
        canonical && design.digits.windows(2).all(|w| {
            lab.residues.iter().position(|&r| r == w[0]) < lab.residues.iter().position(|&r| r == w[1])
        }),
        if walk {
            format!(
                ", walk {} with first revisit {}",
                spell_turns(t),
                crossing.map_or("never, an arc".to_string(), |c| format!("at level {c}"))
            )
        } else if step {
            String::new()
        } else {
            format!(", unit steps fail at level {}", first_break(lab, design).unwrap())
        }
    );
}

fn hilbert(order: usize) -> Vec<Z> {
    let side = 1i64 << order;
    (0..side * side)
        .map(|d| {
            let (mut x, mut y, mut t) = (0i64, 0i64, d);
            let mut s = 1i64;
            while s < side {
                let rx = 1 & (t / 2);
                let ry = 1 & (t ^ rx);
                if ry == 0 {
                    if rx == 1 {
                        x = s - 1 - x;
                        y = s - 1 - y;
                    }
                    std::mem::swap(&mut x, &mut y);
                }
                x += s * rx;
                y += s * ry;
                t /= 4;
                s *= 2;
            }
            (x, y)
        })
        .collect()
}

fn steps_of(list: &[Z]) -> Vec<Z> {
    list.windows(2).map(|w| sub(w[1], w[0])).collect()
}

fn run_named() {
    println!("NAMED  the classical curves located in the census; w = e^(i pi/3) in every printed element");
    let lab = Lab::new(Ring::Eisenstein, (3, 0));
    let k = as_design(&lab, &koch().unwrap());
    assert_eq!(
        k,
        walk_design(&lab, &[0, 1, 5, 0]),
        "the Koch design is not the walk 0150"
    );
    locate("Koch", &lab, &k);
    let lab = Lab::new(Ring::Eisenstein, (2, 1));
    let t = as_design(&lab, &terdragon().unwrap());
    assert_eq!(
        t,
        canonical(&lab, 7, &[0, 2, 0]),
        "the terdragon is not code 7 twisted 020"
    );
    assert_eq!(
        t,
        walk_design(&lab, &[0, 2, 0]),
        "the terdragon is not the walk 020"
    );
    locate("terdragon", &lab, &t);
    let lab = Lab::new(Ring::Gaussian, (1, 1));
    let d = as_design(&lab, &twindragon().unwrap());
    assert_eq!(
        d,
        canonical(&lab, 3, &[0, 0]),
        "the twindragon is not code 3 untwisted"
    );
    locate("twindragon", &lab, &d);
    locate("walk 01 at 1+i", &lab, &walk_design(&lab, &[0, 1]));
    let lab = Lab::new(Ring::Eisenstein, (3, 1));
    let d = as_design(&lab, &flowsnake().unwrap());
    assert_eq!(
        d,
        canonical(&lab, 127, &[0; 7]),
        "the flowsnake is not code 127 untwisted"
    );
    locate("flowsnake", &lab, &d);
    let gosper = [0usize, 1, 3, 2, 0, 0, 5];
    assert_eq!(
        lab.sum(&gosper),
        lab.b,
        "the mirrored Gosper generator does not reach the base"
    );
    locate(
        "Gosper generator 0132005 with every flag F",
        &lab,
        &walk_design(&lab, &gosper),
    );
    let s: Vec<Z> = gosper.iter().map(|&g| lab.units[g]).collect();
    let flags = [false, true, true, false, false, false, true];
    let copy = |j: usize| -> Vec<Z> {
        let order: Vec<usize> = if flags[j] {
            (0..7).rev().collect()
        } else {
            (0..7).collect()
        };
        order[..6].iter().map(|&i| lab.mul(s[j], s[i])).collect()
    };
    let inner: Vec<Z> = s[..6].to_vec();
    let fits: Vec<usize> = (0..7)
        .filter(|&j| {
            (0..lab.n).any(|v| inner.iter().map(|&z| lab.turn(z, v)).collect::<Vec<_>>() == copy(j))
        })
        .collect();
    println!(
        "  Gosper with flags FRRFFFR: the copies whose inner steps are a turn of the level-1 inner steps are {:?} of 0..6; a radix word order needs all seven",
        fits
    );
    let lab = Lab::new(Ring::Gaussian, (2, 0));
    let one: Vec<Z> = hilbert(1);
    let two: Vec<Z> = hilbert(2);
    let first = steps_of(&one);
    let quarter = steps_of(&two[..4]);
    let turns: Vec<usize> = (0..lab.n)
        .filter(|&v| first.iter().map(|&z| lab.turn(z, v)).collect::<Vec<_>>() == quarter)
        .collect();
    println!(
        "  Hilbert: level 1 {} has steps {}, the first quarter of level 2 has steps {}; units carrying one to the other: {:?}",
        spell_list(lab.ring, &one),
        spell_list(lab.ring, &first),
        spell_list(lab.ring, &quarter),
        turns
    );
    let walk: Vec<usize> = first
        .iter()
        .map(|&z| lab.units.iter().position(|&u| u == z).unwrap())
        .chain(std::iter::once(0))
        .collect();
    assert_eq!(lab.sum(&walk), lab.b);
    locate(
        "Hilbert level 1 as a walk with every flag F",
        &lab,
        &walk_design(&lab, &walk),
    );
}
