use crate::machine::{delta, good, image, ones, selfcomplementary, Split};

// GRAPHS

pub fn strongly_connected(edges: &[Vec<usize>]) -> bool {
    let n = edges.len();
    let mut back = vec![Vec::new(); n];
    for (p, out) in edges.iter().enumerate() {
        for &q in out {
            back[q].push(p);
        }
    }
    let reach = |adj: &[Vec<usize>]| {
        let mut seen = vec![false; n];
        let mut stack = vec![0];
        seen[0] = true;
        while let Some(p) = stack.pop() {
            for &q in &adj[p] {
                if !seen[q] {
                    seen[q] = true;
                    stack.push(q);
                }
            }
        }
        seen.iter().all(|&s| s)
    };
    reach(edges) && reach(&back)
}

pub fn depths(edges: &[Vec<usize>], start: usize) -> Vec<usize> {
    let mut depth = vec![usize::MAX; edges.len()];
    let mut queue = std::collections::VecDeque::from([start]);
    depth[start] = 0;
    while let Some(p) = queue.pop_front() {
        for &q in &edges[p] {
            if depth[q] == usize::MAX {
                depth[q] = depth[p] + 1;
                queue.push_back(q);
            }
        }
    }
    depth
}

// SECOND DIFFERENCES

pub struct Probe {
    pub readers: Vec<usize>,
    pub span: Vec<bool>,
    pub zero: bool,
}

fn reduce(rows: &[u32], mut v: u32) -> u32 {
    for &r in rows {
        if v ^ r < v {
            v ^= r;
        }
    }
    v
}

fn insert(rows: &mut Vec<u32>, v: u32) {
    let r = reduce(rows, v);
    if r != 0 {
        rows.push(r);
        rows.sort_by(|a, b| b.cmp(a));
    }
}

pub fn probe(s: &Split, p: usize) -> Probe {
    let live = s.live_reads();
    let readers: Vec<(usize, u128)> = live
        .iter()
        .copied()
        .filter(|&(_, mask)| mask >> p & 1 == 1)
        .collect();
    let window: Vec<usize> = ones(readers.iter().fold(0u128, |w, &(_, m)| w | m) & !(1u128 << p));
    assert!(window.len() <= 22, "a window of {} cells", window.len());
    let local: Vec<(u32, Vec<u8>)> = readers
        .iter()
        .map(|&(_, mask)| {
            let bits = window
                .iter()
                .enumerate()
                .filter(|&(_, &q)| mask >> q & 1 == 1)
                .fold(0u32, |b, (k, _)| b | 1 << k);
            (bits, delta(&s.window(mask), 1))
        })
        .collect();
    let at = |z: u32| {
        local.iter().enumerate().fold(0u32, |v, (i, (bits, d))| {
            v | (d[(z & bits).count_ones() as usize] as u32) << i
        })
    };
    let base = at(0);
    let mut rows = Vec::new();
    for z in 1..1u32 << window.len() {
        insert(&mut rows, at(z) ^ base);
    }
    Probe {
        readers: readers.iter().map(|&(a, _)| a).collect(),
        span: (0..readers.len())
            .map(|i| reduce(&rows, 1 << i) == 0)
            .collect(),
        zero: reduce(&rows, base) == 0,
    }
}

// CONDITIONS

pub struct Conditions {
    pub live: bool,
    pub reads: bool,
    pub spans: bool,
    pub gamma: bool,
    pub windows: bool,
    pub isolated: bool,
    pub edges: Vec<Vec<usize>>,
}

impl Conditions {
    pub fn theorem(&self) -> bool {
        self.live && self.gamma && self.isolated
    }

    pub fn brief(&self) -> bool {
        self.reads && self.spans && self.windows
    }
}

pub fn isolated(s: &Split, q: usize) -> bool {
    let live = s.live_reads();
    live.iter().any(|&(a, mask)| {
        if mask >> q & 1 == 0 {
            return false;
        }
        let rest = ones(mask & !(1u128 << q));
        let m = mask.count_ones() as usize;
        (1u32..1 << rest.len()).any(|pick| {
            let sub = rest
                .iter()
                .enumerate()
                .filter(|&(k, _)| pick >> k & 1 == 1)
                .fold(0u128, |w, (_, &c)| w | 1u128 << c);
            let k = sub.count_ones() as usize;
            good(&delta(&s.window(mask), k)[..=m - k])
                && live
                    .iter()
                    .filter(|&&(_, other)| other & sub == sub)
                    .all(|&(b, _)| b == a)
        })
    })
}

pub fn conditions(schedule: &[Split]) -> Conditions {
    let cells = schedule[0].cells();
    let mut reads = vec![Vec::new(); cells];
    let mut gamma = vec![Vec::new(); cells];
    let (mut live, mut spans, mut windows, mut iso) = (false, true, true, true);
    for s in schedule {
        let lr = s.live_reads();
        let read: u128 = lr.iter().fold(0, |w, &(_, m)| w | m);
        live |= read != 0;
        for &(_, mask) in &lr {
            windows &= good(&s.window(mask));
        }
        for p in 0..cells {
            reads[p].push(s.perm[p]);
            if read >> p & 1 == 0 {
                gamma[p].push(s.perm[p]);
                continue;
            }
            let pr = probe(s, p);
            for (i, &a) in pr.readers.iter().enumerate() {
                reads[p].push(s.perm[a]);
                if pr.span[i] {
                    gamma[p].push(s.perm[a]);
                } else {
                    spans = false;
                }
            }
            if pr.zero {
                gamma[p].push(s.perm[p]);
            }
            iso &= isolated(s, p);
        }
    }
    Conditions {
        live,
        reads: strongly_connected(&reads),
        spans,
        gamma: strongly_connected(&gamma),
        windows,
        isolated: iso,
        edges: gamma,
    }
}

// LAWS

fn module(schedule: &[Split]) -> usize {
    let mut rows: Vec<u128> = Vec::new();
    let reduce = |rows: &[u128], mut v: u128| {
        for &r in rows {
            if v ^ r < v {
                v ^= r;
            }
        }
        v
    };
    let mut queue: Vec<u128> = Vec::new();
    for s in schedule {
        let mut by: std::collections::BTreeMap<u128, u128> = Default::default();
        for (a, mask) in s.live_reads() {
            *by.entry(mask).or_default() |= 1u128 << a;
        }
        queue.extend(by.into_values());
    }
    while let Some(v) = queue.pop() {
        let r = reduce(&rows, v);
        if r != 0 {
            rows.push(r);
            rows.sort_by(|a, b| b.cmp(a));
            for s in schedule {
                queue.push(image(&s.perm, v));
            }
        }
    }
    rows.len()
}

pub fn law(schedule: &[Split], c: &Conditions, gens: &[&[u32]]) -> String {
    let cells = schedule[0].cells();
    if !c.reads {
        return "closed".to_string();
    }
    if module(schedule) < cells {
        return "mess".to_string();
    }
    if schedule.iter().all(|s| {
        s.reads
            .iter()
            .all(|&(_, m)| selfcomplementary(&s.window(m)))
    }) {
        return "complement".to_string();
    }
    parity(schedule[0].n, gens).map_or("other".to_string(), |p| format!("parity of {p}"))
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

fn parity(n: usize, gens: &[&[u32]]) -> Option<String> {
    let mut parts: Vec<(String, Vec<usize>)> = vec![("torus".to_string(), vec![0; n * n])];
    for h in subgroups(n) {
        let label = (0..n * n)
            .map(|p| {
                h.iter()
                    .map(|&q| ((p / n + q / n) % n) * n + (p % n + q % n) % n)
                    .min()
                    .unwrap()
            })
            .collect();
        parts.push((format!("cosets of order {}", h.len()), label));
    }
    if n == 4 {
        parts.push((
            "2 x 2 tiles".to_string(),
            (0..16).map(|p| (p / 8) * 2 + (p % 4) / 2).collect(),
        ));
    }
    parts.into_iter().find_map(|(name, label)| {
        let mut u = crate::span::Span::new(n * n);
        for p in 0..n * n {
            if let Some(q) = (p + 1..n * n).find(|&q| label[q] == label[p]) {
                u.insert(1 << p | 1 << q);
            }
        }
        (u.dim > 0 && crate::span::is_block(gens, &u)).then_some(name)
    })
}
