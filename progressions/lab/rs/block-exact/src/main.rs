use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use std::time::Instant;

type Mask = u128;

const THREADS: usize = 8;
const KEEP: usize = 4096;

// THE GRID

struct Grid {
    m: usize,
    n: usize,
    k: usize,
    free: bool,
    mid: Vec<u8>,
    sub: Vec<u8>,
    ext: Vec<u8>,
    line_pts: Vec<u8>,
    line_mask: Vec<Mask>,
    mats: Vec<Vec<u8>>,
    class_pts: Vec<u8>,
    maxpeel: Vec<u8>,
}

fn bits(mask: Mask, out: &mut [u8; 128]) -> usize {
    let mut k = 0;
    let mut rest = mask;
    while rest != 0 {
        out[k] = rest.trailing_zeros() as u8;
        rest &= rest - 1;
        k += 1;
    }
    k
}

fn inner_1d(m: usize, k: usize, live: usize) -> usize {
    let inv2 = m.div_ceil(2);
    let mut inner = 0usize;
    for x in 0..m {
        for z in x + 1..m {
            if live >> x & 1 == 0 || live >> z & 1 == 0 {
                continue;
            }
            if k == 3 {
                inner |= 1 << ((x + z) * inv2 % m);
            } else if live >> ((2 * x + m - z) % m) & 1 == 1
                && live >> ((2 * z + m - x) % m) & 1 == 1
            {
                inner |= 1 << x | 1 << z;
            }
        }
    }
    inner & live
}

fn peelable_1d(m: usize, k: usize, free: bool, set: usize) -> bool {
    let mut live = set;
    while live != 0 {
        let inner = inner_1d(m, k, live);
        if inner == live || (free && inner != 0) {
            return false;
        }
        live = inner;
    }
    true
}

impl Grid {
    fn new(m: usize, k: usize, free: bool) -> Grid {
        assert!(
            m % 2 == 1 && m >= 3 && (2..m).all(|d| !m.is_multiple_of(d)),
            "m must be an odd prime"
        );
        assert!(k == 3 || (k == 4 && m >= 5), "k is 3, or 4 at m >= 5");
        let n = m * m;
        let inv2 = m.div_ceil(2);
        let at = |a: usize, b: usize| (a % m) * m + (b % m);
        let coord = |i: usize| (i / m, i % m);
        let mut mid = vec![0u8; n * n];
        let mut sub = vec![0u8; n * n];
        let mut ext = vec![0u8; n * n];
        for x in 0..n {
            for z in 0..n {
                let (xa, xb) = coord(x);
                let (za, zb) = coord(z);
                mid[x * n + z] = at((xa + za) * inv2, (xb + zb) * inv2) as u8;
                sub[x * n + z] = at(xa + m - za, xb + m - zb) as u8;
                ext[x * n + z] = at(2 * xa + m - za, 2 * xb + m - zb) as u8;
            }
        }
        let mut line_pts = vec![0u8; n * n * m];
        let mut line_mask = vec![0 as Mask; n * n];
        for t in 0..n {
            for u in 0..n {
                if t == u {
                    continue;
                }
                let (ta, tb) = coord(t);
                let (da, db) = coord(sub[u * n + t] as usize);
                for c in 0..m {
                    let p = at(ta + c * da, tb + c * db);
                    line_pts[(t * n + u) * m + c] = p as u8;
                    line_mask[t * n + u] |= 1 << p;
                }
            }
        }
        let mut mats = vec![Vec::new(); n];
        for p in 0..m {
            for q in 0..m {
                for r in 0..m {
                    for s in 0..m {
                        if (p * s + m * m - q * r).is_multiple_of(m) {
                            continue;
                        }
                        for (d, row) in mats.iter_mut().enumerate().skip(1) {
                            let (da, db) = coord(d);
                            if (p * da + q * db) % m == 0 && (r * da + s * db) % m == 1 {
                                for v in 0..n {
                                    let (a, b) = coord(v);
                                    row.push(at(p * a + q * b, r * a + s * b) as u8);
                                }
                            }
                        }
                    }
                }
            }
        }
        for row in mats.iter().skip(1) {
            assert_eq!(row.len(), (n - m) * n);
        }
        let mut class_pts = vec![0u8; (m + 1) * m * m];
        let dirs: Vec<(usize, usize)> = std::iter::once((0, 1))
            .chain((0..m).map(|j| (1, j)))
            .collect();
        for (ci, &(da, db)) in dirs.iter().enumerate() {
            for l in 0..m {
                let (ba, bb) = if da == 0 { (l, 0) } else { (0, l) };
                for c in 0..m {
                    class_pts[(ci * m + l) * m + c] = at(ba + c * da, bb + c * db) as u8;
                }
            }
        }
        let full = 1usize << m;
        let mut maxpeel = vec![0u8; full * full];
        for p in 0..full {
            if !peelable_1d(m, k, free, p) {
                continue;
            }
            let k = p.count_ones() as u8;
            let comp = !p & (full - 1);
            let mut f = p;
            loop {
                let mut s = comp;
                loop {
                    let e = &mut maxpeel[f * full + (p | s)];
                    if *e < k {
                        *e = k;
                    }
                    if s == 0 {
                        break;
                    }
                    s = (s - 1) & comp;
                }
                if f == 0 {
                    break;
                }
                f = (f - 1) & p;
            }
        }
        Grid {
            m,
            n,
            k,
            free,
            mid,
            sub,
            ext,
            line_pts,
            line_mask,
            mats,
            class_pts,
            maxpeel,
        }
    }

    fn group(&self) -> u128 {
        (self.n * (self.n - 1) * (self.n - self.m)) as u128
    }

    fn inner(&self, live: Mask) -> Mask {
        let mut pts = [0u8; 128];
        let c = bits(live, &mut pts);
        let mut inner: Mask = 0;
        for i in 0..c {
            let x = pts[i] as usize;
            for &z in &pts[i + 1..c] {
                let z = z as usize;
                if self.k == 3 {
                    inner |= 1 << self.mid[x * self.n + z];
                } else if live >> self.ext[x * self.n + z] & 1 == 1
                    && live >> self.ext[z * self.n + x] & 1 == 1
                {
                    inner |= 1 << x | 1 << z;
                }
            }
        }
        inner & live
    }

    fn peelable(&self, set: Mask) -> bool {
        let mut live = set;
        while live != 0 {
            let inner = self.inner(live);
            if inner == live || (self.free && inner != 0) {
                return false;
            }
            live = inner;
        }
        true
    }

    fn rounds(&self, set: Mask) -> Option<Vec<u8>> {
        let mut pts = [0u8; 128];
        let mut level = vec![255u8; self.n];
        let mut live = set;
        let mut r = 0u8;
        while live != 0 {
            let inner = self.inner(live);
            let rem = live & !inner;
            if rem == 0 || (self.free && inner != 0) {
                return None;
            }
            let kr = bits(rem, &mut pts);
            for i in 0..kr {
                level[pts[i] as usize] = r;
            }
            live = inner;
            r += 1;
        }
        Some(level)
    }

    fn canonical(&self, t: Mask, aut: &mut u64) -> bool {
        let (m, n) = (self.m, self.n);
        let lo: Mask = (1 << m) - 1;
        let t_line = (t & lo) as usize;
        let t_off = t & !lo;
        let mut pts = [0u8; 128];
        let mut offpts = [0u8; 128];
        let k = bits(t, &mut pts);
        if k < 2 {
            *aut = (self.group() / if k == 0 { 1 } else { n as u128 }) as u64;
            return k == 0 || t == 1;
        }
        let per = n - m;
        for i in 0..k {
            let a = pts[i] as usize;
            for (j, &b) in pts[..k].iter().enumerate() {
                if i == j {
                    continue;
                }
                let b = b as usize;
                let pair = a * n + b;
                let lp = &self.line_pts[pair * m..(pair + 1) * m];
                let mut f = 0usize;
                for (c, &p) in lp.iter().enumerate() {
                    f |= (((t >> p) & 1) as usize) << c;
                }
                if f != t_line {
                    let diff = f ^ t_line;
                    if f >> diff.trailing_zeros() & 1 == 1 {
                        return false;
                    }
                    continue;
                }
                let d = self.sub[b * n + a] as usize;
                let ko = bits(t & !self.line_mask[pair], &mut offpts);
                let tab = &self.mats[d];
                for x in 0..per {
                    let perm = &tab[x * n..(x + 1) * n];
                    let mut h: Mask = 0;
                    for v in 0..ko {
                        h |= 1 << perm[self.sub[offpts[v] as usize * n + a] as usize];
                    }
                    if h == t_off {
                        *aut += 1;
                        continue;
                    }
                    let diff = h ^ t_off;
                    if h >> diff.trailing_zeros() & 1 == 1 {
                        return false;
                    }
                }
            }
        }
        true
    }

    fn line_bound(&self, s: Mask, u: Mask) -> usize {
        let m = self.m;
        let mut best = usize::MAX;
        for ci in 0..m + 1 {
            let mut total = 0usize;
            for l in 0..m {
                let base = (ci * m + l) * m;
                let (mut f, mut a) = (0usize, 0usize);
                for c in 0..m {
                    let p = self.class_pts[base + c];
                    f |= (((s >> p) & 1) as usize) << c;
                    a |= (((u >> p) & 1) as usize) << c;
                }
                total += self.maxpeel[(f << m) | a] as usize;
            }
            best = best.min(total);
        }
        best
    }

    fn cert_holds(&self, set: Mask, g: &[u64]) -> bool {
        let mut pts = [0u8; 128];
        let c = bits(set, &mut pts);
        for i in 0..c {
            for j in i + 1..c {
                let (x, z) = (pts[i] as usize, pts[j] as usize);
                let ok = if self.k == 3 {
                    let y = self.mid[x * self.n + z] as usize;
                    set >> y & 1 == 0 || g[x] + g[z] > 2 * g[y]
                } else {
                    let u = self.ext[x * self.n + z] as usize;
                    let w = self.ext[z * self.n + x] as usize;
                    set >> u & 1 == 0 || set >> w & 1 == 0 || g[u] + g[w] > g[x] + g[z]
                };
                if !ok {
                    return false;
                }
            }
        }
        true
    }

    fn bitmap(&self, set: Mask, level: Option<&[u8]>) -> Vec<String> {
        let m = self.m;
        (0..m)
            .rev()
            .map(|b| {
                (0..m)
                    .map(|a| {
                        let v = a * m + b;
                        match (set >> v & 1 == 1, level) {
                            (false, _) => '.',
                            (true, Some(lv)) => char::from_digit(lv[v] as u32, 36).unwrap(),
                            (true, None) => '#',
                        }
                    })
                    .collect()
            })
            .collect()
    }
}

// THE SEARCH

#[derive(Default, Clone)]
struct Stats {
    nodes: u64,
    tests: u64,
    hist: Vec<u64>,
    found: Vec<(Mask, u64)>,
    orbits: u128,
    sets: u128,
    origin: u128,
}

struct Frontier {
    depth: usize,
    items: Vec<(Mask, usize, u64, Vec<u8>)>,
}

struct Run<'a> {
    grid: &'a Grid,
    target: usize,
    census: bool,
    stop: &'a AtomicBool,
    deadline: Option<Instant>,
}

impl<'a> Run<'a> {
    fn visit(&self, s: Mask, size: usize, aut: u64, st: &mut Stats) {
        st.nodes += 1;
        if st.nodes & 4095 == 0 {
            if let Some(d) = self.deadline {
                if Instant::now() >= d {
                    self.stop.store(true, Ordering::Relaxed);
                }
            }
        }
        if size >= self.target {
            st.hist[size] += 1;
            if st.found.len() < KEEP {
                st.found.push((s, aut));
            }
        }
        if self.census {
            let orbit = self.grid.group() / aut as u128;
            st.orbits += 1;
            st.sets += orbit;
            st.origin += orbit * size as u128 / self.grid.n as u128;
        }
    }

    fn children(
        &self,
        s: Mask,
        size: usize,
        cands: &[u8],
        mut f: impl FnMut(Mask, usize, u64, Vec<u8>),
    ) {
        for (i, &v) in cands.iter().enumerate() {
            let rest = cands.len() - 1 - i;
            if size + 1 + rest < self.target {
                break;
            }
            if self.stop.load(Ordering::Relaxed) {
                return;
            }
            let t = s | (1 << v);
            let mut aut = 0u64;
            if !self.grid.canonical(t, &mut aut) {
                continue;
            }
            let mut newc = Vec::with_capacity(rest);
            let mut u = t;
            for &c in &cands[i + 1..] {
                if self.grid.peelable(t | (1 << c)) {
                    newc.push(c);
                    u |= 1 << c;
                }
            }
            if size + 1 + newc.len() < self.target || self.grid.line_bound(t, u) < self.target {
                continue;
            }
            f(t, size + 1, aut, newc);
        }
    }

    fn dfs(&self, s: Mask, size: usize, aut: u64, cands: &[u8], st: &mut Stats) {
        self.visit(s, size, aut, st);
        self.children(s, size, cands, |t, sz, au, newc| {
            st.tests += 1;
            self.dfs(t, sz, au, &newc, st);
        });
    }

    fn frontier(
        &self,
        s: Mask,
        size: usize,
        aut: u64,
        cands: &[u8],
        st: &mut Stats,
        fr: &mut Frontier,
    ) {
        if size == fr.depth {
            fr.items.push((s, size, aut, cands.to_vec()));
            return;
        }
        self.visit(s, size, aut, st);
        self.children(s, size, cands, |t, sz, au, newc| {
            self.frontier(t, sz, au, &newc, st, fr);
        });
    }
}

fn merge(into: &mut Stats, from: Stats) {
    into.nodes += from.nodes;
    into.tests += from.tests;
    for (a, b) in into.hist.iter_mut().zip(from.hist) {
        *a += b;
    }
    for item in from.found {
        if into.found.len() < KEEP {
            into.found.push(item);
        }
    }
    into.orbits += from.orbits;
    into.sets += from.sets;
    into.origin += from.origin;
}

fn search(grid: &Grid, target: usize, census: bool, seconds: Option<f64>) -> (Stats, bool) {
    let stop = AtomicBool::new(false);
    let deadline = seconds.map(|s| Instant::now() + std::time::Duration::from_secs_f64(s));
    let run = Run {
        grid,
        target,
        census,
        stop: &stop,
        deadline,
    };
    let n = grid.n;
    let mut root_aut = 0u64;
    assert!(grid.canonical(3, &mut root_aut));
    let mut st = Stats {
        hist: vec![0; n + 1],
        ..Default::default()
    };
    if census {
        run.visit(0, 0, grid.group() as u64, &mut st);
        run.visit(1, 1, (grid.group() / n as u128) as u64, &mut st);
    }
    let cands: Vec<u8> = (2..n as u8)
        .filter(|&c| grid.peelable(3 | (1 << c)))
        .collect();
    let mut fr = Frontier {
        depth: 3,
        items: Vec::new(),
    };
    loop {
        let mut probe = Stats {
            hist: vec![0; n + 1],
            ..Default::default()
        };
        fr.items.clear();
        run.frontier(3, 2, root_aut, &cands, &mut probe, &mut fr);
        if fr.items.len() >= 64 * THREADS || fr.items.is_empty() || fr.depth >= 8 {
            st = Stats {
                hist: vec![0; n + 1],
                ..Default::default()
            };
            if census {
                run.visit(0, 0, grid.group() as u64, &mut st);
                run.visit(1, 1, (grid.group() / n as u128) as u64, &mut st);
            }
            merge(&mut st, probe);
            break;
        }
        fr.depth += 1;
    }
    let items = fr.items;
    let next = AtomicUsize::new(0);
    let parts: Vec<Stats> = std::thread::scope(|sc| {
        let handles: Vec<_> = (0..THREADS)
            .map(|_| {
                sc.spawn(|| {
                    let mut local = Stats {
                        hist: vec![0; n + 1],
                        ..Default::default()
                    };
                    loop {
                        let i = next.fetch_add(1, Ordering::Relaxed);
                        if i >= items.len() || stop.load(Ordering::Relaxed) {
                            break;
                        }
                        let (s, size, aut, cands) = &items[i];
                        run.dfs(*s, *size, *aut, cands, &mut local);
                    }
                    local
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for p in parts {
        merge(&mut st, p);
    }
    (st, stop.load(Ordering::Relaxed))
}

fn greedy(grid: &Grid, tries: u64) -> usize {
    let n = grid.n;
    let mut best = 0;
    let mut seed = 0x9e3779b97f4a7c15u64;
    for _ in 0..tries {
        let mut order: Vec<u8> = (0..n as u8).collect();
        for i in (1..n).rev() {
            seed ^= seed << 13;
            seed ^= seed >> 7;
            seed ^= seed << 17;
            order.swap(i, (seed % (i as u64 + 1)) as usize);
        }
        let mut set: Mask = 0;
        let mut size = 0;
        for &v in &order {
            if grid.peelable(set | (1 << v)) {
                set |= 1 << v;
                size += 1;
            }
        }
        best = best.max(size);
    }
    best
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

// THE CLIMB

struct Cube {
    m: usize,
    d: usize,
    n: usize,
    k: usize,
    coords: Vec<u8>,
    sum: Vec<u16>,
    ext: Vec<u16>,
    half_off: Vec<u32>,
    half: Vec<u16>,
}

impl Cube {
    fn new(m: usize, d: usize, k: usize) -> Cube {
        let n = m.pow(d as u32);
        assert!(m >= 2 && d >= 1 && n <= 1 << 16, "m^d must fit in 16 bits");
        assert!(k == 3 || (k == 4 && m >= 4), "k is 3, or 4 at m >= 4");
        let mut coords = vec![0u8; n * d];
        for v in 0..n {
            let mut r = v;
            for i in (0..d).rev() {
                coords[v * d + i] = (r % m) as u8;
                r /= m;
            }
        }
        let place: Vec<usize> = (0..d).map(|i| m.pow((d - 1 - i) as u32)).collect();
        let add = |x: usize, z: usize, k: usize| {
            (0..d)
                .map(|i| {
                    (k * coords[x * d + i] as usize + coords[z * d + i] as usize) % m * place[i]
                })
                .sum::<usize>()
        };
        let mut sum = vec![0u16; n * n];
        let mut ext = vec![0u16; n * n];
        let neg: Vec<usize> = (0..n).map(|z| add(z, 0, m - 1)).collect();
        for x in 0..n {
            for z in 0..n {
                sum[x * n + z] = add(x, z, 1) as u16;
                ext[x * n + z] = add(x, neg[z], 2) as u16;
            }
        }
        let mut buckets: Vec<Vec<u16>> = vec![Vec::new(); n];
        for y in 0..n {
            buckets[add(y, y, 1)].push(y as u16);
        }
        let mut half_off = vec![0u32];
        let mut half = Vec::new();
        for b in &buckets {
            half.extend_from_slice(b);
            half_off.push(half.len() as u32);
        }
        Cube {
            m,
            d,
            n,
            k,
            coords,
            sum,
            ext,
            half_off,
            half,
        }
    }

    fn mids(&self, x: usize, z: usize) -> &[u16] {
        let v = self.sum[x * self.n + z] as usize;
        &self.half[self.half_off[v] as usize..self.half_off[v + 1] as usize]
    }

    fn ext(&self, x: usize, z: usize) -> usize {
        self.ext[x * self.n + z] as usize
    }

    fn distinct(t: [usize; 4]) -> bool {
        t[0] != t[1] && t[0] != t[2] && t[0] != t[3] && t[1] != t[2] && t[1] != t[3] && t[2] != t[3]
    }

    fn inner_pair(&self, a: usize, b: usize, has: impl Fn(usize) -> bool) -> bool {
        let (u, w) = (self.ext(a, b), self.ext(b, a));
        has(u) && has(w) && Cube::distinct([u, a, b, w])
    }

    fn rounds(&self, set: &[usize]) -> Result<Vec<u32>, Vec<usize>> {
        let c = set.len();
        let mut pos = vec![usize::MAX; self.n];
        for (i, &v) in set.iter().enumerate() {
            pos[v] = i;
        }
        let mut cnt = vec![0u32; c];
        for a in 0..c {
            for b in a + 1..c {
                if self.k == 3 {
                    for &w in self.mids(set[a], set[b]) {
                        if pos[w as usize] != usize::MAX {
                            cnt[pos[w as usize]] += 1;
                        }
                    }
                } else if self.inner_pair(set[a], set[b], |v| pos[v] != usize::MAX) {
                    cnt[a] += 1;
                    cnt[b] += 1;
                }
            }
        }
        let mut live = vec![true; c];
        let mut round = vec![u32::MAX; c];
        let mut left = c;
        let mut r = 0;
        while left > 0 {
            let out: Vec<usize> = (0..c).filter(|&i| live[i] && cnt[i] == 0).collect();
            if out.is_empty() {
                return Err((0..c).filter(|&i| live[i]).collect());
            }
            for &i in &out {
                round[i] = r;
            }
            for &i in &out {
                live[i] = false;
                left -= 1;
                let alive = |v: usize| pos[v] != usize::MAX && live[pos[v]];
                for j in 0..c {
                    if !live[j] {
                        continue;
                    }
                    let (x, y) = (set[i], set[j]);
                    if self.k == 3 {
                        for &w in self.mids(x, y) {
                            if pos[w as usize] != usize::MAX {
                                cnt[pos[w as usize]] -= 1;
                            }
                        }
                        continue;
                    }
                    if self.inner_pair(x, y, alive) {
                        cnt[j] -= 1;
                    }
                    let b = self.ext(y, x);
                    let e = self.ext(b, y);
                    if alive(b) && alive(e) && Cube::distinct([x, y, b, e]) {
                        cnt[j] -= 1;
                        cnt[pos[b]] -= 1;
                    }
                }
            }
            r += 1;
        }
        Ok(round)
    }

    fn touched(&self, q: usize, set: &[usize], ins: &[bool], mut f: impl FnMut(usize)) {
        for &c in set {
            if self.k == 3 {
                for &w in self.mids(q, c) {
                    f(w as usize);
                }
                continue;
            }
            let w = self.ext(q, c);
            let x = self.ext(w, q);
            if ins[x] && Cube::distinct([x, w, q, c]) {
                f(w);
            }
            for &w in self.mids(q, c) {
                let w = w as usize;
                let x = self.ext(c, w);
                if ins[x] && Cube::distinct([q, w, c, x]) {
                    f(w);
                }
            }
            let w = self.ext(c, q);
            let x = self.ext(w, c);
            if ins[x] && Cube::distinct([x, w, c, q]) {
                f(w);
            }
        }
    }

    fn point(&self, v: usize) -> String {
        let c: Vec<String> = (0..self.d)
            .map(|i| self.coords[v * self.d + i].to_string())
            .collect();
        format!("({})", c.join(","))
    }
}

struct Climber<'a> {
    cube: &'a Cube,
    ins: Vec<bool>,
    cov: Vec<u32>,
    set: Vec<usize>,
    stamp: Vec<u32>,
    tick: u32,
    seed: u64,
}

impl<'a> Climber<'a> {
    fn new(cube: &'a Cube, seed: u64) -> Climber<'a> {
        let n = cube.n;
        Climber {
            cube,
            ins: vec![false; n],
            cov: vec![0; n],
            set: Vec::new(),
            stamp: vec![0; n],
            tick: 0,
            seed: seed.wrapping_mul(0x9e3779b97f4a7c15) | 1,
        }
    }

    fn rand(&mut self) -> u64 {
        self.seed ^= self.seed << 13;
        self.seed ^= self.seed >> 7;
        self.seed ^= self.seed << 17;
        self.seed
    }

    fn reset(&mut self, set: &[usize]) {
        self.ins.iter_mut().for_each(|b| *b = false);
        self.cov.iter_mut().for_each(|c| *c = 0);
        self.set.clear();
        for &p in set {
            self.insert(p);
        }
    }

    fn insert(&mut self, p: usize) {
        let cov = &mut self.cov;
        self.cube.touched(p, &self.set, &self.ins, |w| cov[w] += 1);
        self.ins[p] = true;
        self.set.push(p);
    }

    fn remove(&mut self, i: usize) {
        let p = self.set.swap_remove(i);
        self.ins[p] = false;
        let cov = &mut self.cov;
        self.cube.touched(p, &self.set, &self.ins, |w| cov[w] -= 1);
    }

    fn score(&mut self, p: usize) -> u32 {
        self.tick += 1;
        let (ins, cov, stamp, tick) = (&self.ins, &self.cov, &mut self.stamp, self.tick);
        let mut c = 0;
        self.cube.touched(p, &self.set, ins, |w| {
            if !ins[w] && cov[w] == 0 && w != p && stamp[w] != tick {
                stamp[w] = tick;
                c += 1;
            }
        });
        c
    }

    fn free_extend(&mut self, temp: f64) {
        loop {
            let mut best = None;
            let mut key = f64::INFINITY;
            for p in 0..self.cube.n {
                if !self.ins[p] && self.cov[p] == 0 {
                    let k = self.score(p) as f64
                        + temp * (self.rand() >> 11) as f64 / (1u64 << 53) as f64;
                    if k < key {
                        key = k;
                        best = Some(p);
                    }
                }
            }
            match best {
                Some(p) => self.insert(p),
                None => return,
            }
        }
    }

    fn deep_extend(&mut self, temp: f64, tries: usize) {
        loop {
            let mut cands: Vec<usize> = (0..self.cube.n).filter(|&p| !self.ins[p]).collect();
            for i in (1..cands.len()).rev() {
                let j = (self.rand() % (i as u64 + 1)) as usize;
                cands.swap(i, j);
            }
            let mut added = false;
            for &p in cands.iter().take(tries) {
                let mut trial = self.set.clone();
                trial.push(p);
                if self.cube.rounds(&trial).is_ok() {
                    self.insert(p);
                    self.free_extend(temp);
                    added = true;
                    break;
                }
            }
            if !added {
                return;
            }
        }
    }

    fn kick(&mut self, temp: f64, tries: usize, big: bool) {
        let out = if big {
            4 + (self.rand() % 9) as usize
        } else {
            1 + (self.rand() % 3) as usize
        };
        if !big && self.rand().is_multiple_of(2) {
            let cands: Vec<usize> = (0..self.cube.n).filter(|&p| !self.ins[p]).collect();
            if !cands.is_empty() {
                let p = cands[(self.rand() % cands.len() as u64) as usize];
                self.insert(p);
                loop {
                    match self.cube.rounds(&self.set) {
                        Ok(_) => break,
                        Err(core) => {
                            let core: Vec<usize> =
                                core.into_iter().filter(|&i| self.set[i] != p).collect();
                            let i = if core.is_empty() {
                                self.set.iter().position(|&v| v == p).unwrap()
                            } else {
                                core[(self.rand() % core.len() as u64) as usize]
                            };
                            self.remove(i);
                        }
                    }
                }
            }
        } else {
            for _ in 0..out.min(self.set.len()) {
                let i = (self.rand() % self.set.len() as u64) as usize;
                self.remove(i);
            }
        }
        self.free_extend(temp);
        self.deep_extend(temp, tries);
    }
}

fn climb_thread(
    cube: &Cube,
    seed: u64,
    deadline: Instant,
    stale_cap: usize,
    tries: usize,
) -> (Vec<usize>, u64) {
    let mut c = Climber::new(cube, seed);
    let mut best: Vec<usize> = Vec::new();
    let mut moves = 0u64;
    let mut round = 0u64;
    while Instant::now() < deadline {
        let temp = 0.5 + (c.rand() % 4) as f64;
        if round.is_multiple_of(8) || best.is_empty() {
            c.reset(&[]);
            c.free_extend(temp);
            c.deep_extend(temp, tries);
        } else {
            c.reset(&best);
            c.kick(temp, tries, true);
        }
        round += 1;
        let mut cur = c.set.clone();
        let mut stale = 0;
        while stale < stale_cap && Instant::now() < deadline {
            c.reset(&cur);
            c.kick(temp, tries, false);
            moves += 1;
            if c.set.len() > cur.len() {
                stale = 0;
            } else {
                stale += 1;
            }
            if c.set.len() >= cur.len() {
                cur = c.set.clone();
            }
            if cur.len() > best.len() {
                best = cur.clone();
            }
        }
        if cur.len() > best.len() {
            best = cur.clone();
        }
    }
    (best, moves)
}

fn lift_step(cube: &Cube, blocks: usize, a: usize, b: usize) -> usize {
    let n = cube.n;
    let (mut a, mut b, mut v, mut place) = (a, b, 0, 1);
    for _ in 0..blocks {
        v += cube.ext(b % n, a % n) * place;
        place *= n;
        a /= n;
        b /= n;
    }
    v
}

fn lift_point(cube: &Cube, blocks: usize, v: usize) -> String {
    let mut r = v;
    (0..blocks)
        .map(|_| {
            let p = cube.point(r % cube.n);
            r /= cube.n;
            p
        })
        .collect()
}

fn lift_check(
    cube: &Cube,
    set: &[usize],
    g: &[u128],
    blocks: usize,
) -> (usize, usize, usize, Option<String>) {
    let (c, n, k) = (set.len(), cube.n, cube.k);
    let total = c.pow(blocks as u32);
    let size = n.pow(blocks as u32);
    assert!(size <= 1 << 28, "the lift must fit in memory");
    let tuple = |t: usize| {
        let (mut r, mut s, mut v, mut place) = (t, 0u128, 0usize, 1usize);
        for _ in 0..blocks {
            s += g[r % c];
            v += set[r % c] * place;
            place *= n;
            r /= c;
        }
        (s, v)
    };
    let mut keys: Vec<u128> = (0..total).map(|t| tuple(t).0).collect();
    keys.sort();
    keys.dedup();
    let mut slot = vec![u32::MAX; size];
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); keys.len()];
    for t in 0..total {
        let (s, v) = tuple(t);
        let id = keys.binary_search(&s).unwrap();
        slot[v] = id as u32;
        members[id].push(v);
    }
    let mut found = 0;
    let mut first = None;
    for (id, pts) in members.iter().enumerate() {
        for &x0 in pts {
            for &x1 in pts {
                if x0 == x1 {
                    continue;
                }
                let mut t = [x0, x1, 0, 0];
                let mut ok = true;
                for j in 2..k {
                    t[j] = lift_step(cube, blocks, t[j - 2], t[j - 1]);
                    if slot[t[j]] != id as u32 {
                        ok = false;
                        break;
                    }
                }
                let distinct = (0..k).all(|i| (i + 1..k).all(|j| t[i] != t[j]));
                if ok && distinct && t[0] < t[k - 1] {
                    found += 1;
                    if first.is_none() {
                        let terms: Vec<String> = t[..k]
                            .iter()
                            .map(|&v| lift_point(cube, blocks, v))
                            .collect();
                        first = Some(format!("slice {}, terms {}", keys[id], terms.join(" ")));
                    }
                }
            }
        }
    }
    let largest = members.iter().map(Vec::len).max().unwrap_or(0);
    (keys.len(), largest, found, first)
}

fn report(cube: &Cube, best: &[usize], blocks: usize) -> usize {
    let (m, d, k) = (cube.m, cube.d, best.len());
    let level = cube.rounds(best).expect("the set peels");
    let r = level.iter().max().map_or(0, |&l| l + 1);
    assert!(r <= 80, "the certificate must fit in 128 bits");
    let g: Vec<u128> = level.iter().map(|&l| 3u128.pow(r - 1 - l)).collect();
    let pos: std::collections::HashMap<usize, usize> =
        best.iter().enumerate().map(|(i, &v)| (v, i)).collect();
    let mut progs = 0u64;
    for a in 0..k {
        for b in a + 1..k {
            let (x, z) = (best[a], best[b]);
            if cube.k == 3 {
                for &w in cube.mids(x, z) {
                    if let Some(&c) = pos.get(&(w as usize)) {
                        progs += 1;
                        assert!(g[a] + g[b] > 2 * g[c], "certificate fails");
                    }
                }
            } else if cube.inner_pair(x, z, |v| pos.contains_key(&v)) {
                progs += 1;
                let (u, w) = (pos[&cube.ext(x, z)], pos[&cube.ext(z, x)]);
                assert!(g[u] + g[w] > g[a] + g[b], "certificate fails");
            }
        }
    }
    let what = if cube.k == 3 {
        "triples x + z = 2y with x != z"
    } else {
        "four-term progressions with distinct terms"
    };
    println!("  peel rounds {r}, certificate g = 3^({} - round) checked in exact integers on all {progs} {what} inside the set", r.saturating_sub(1));
    let pts: Vec<String> = best.iter().map(|&v| cube.point(v)).collect();
    println!("  points {}", pts.join(" "));
    let rs: Vec<String> = level.iter().map(|l| l.to_string()).collect();
    println!("  rounds {}", rs.join(" "));
    if k.pow(blocks as u32) > 4_000_000 {
        return 0;
    }
    let (slices, largest, found, first) = lift_check(cube, best, &g, blocks);
    println!(
        "  lift to Z_{m}^{}: {slices} slices of T^{blocks} by the sum of g, the largest {largest} points, {found} {}-term progressions with distinct terms inside the slices by brute force over pairs",
        d * blocks,
        cube.k
    );
    if let Some(w) = first {
        println!("  first: {w}");
    }
    found
}

fn climb(m: usize, d: usize, k: usize, seconds: f64, seed: u64, blocks: usize) {
    let t0 = Instant::now();
    let cube = Cube::new(m, d, k);
    let n = cube.n;
    let deadline = t0 + std::time::Duration::from_secs_f64(seconds);
    let tries = if n <= 400 { n } else { 48 };
    let results: Vec<(Vec<usize>, u64)> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..THREADS as u64)
            .map(|t| {
                let cube = &cube;
                sc.spawn(move || climb_thread(cube, seed * 1000 + t, deadline, 40, tries))
            })
            .collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let sizes: Vec<String> = results.iter().map(|(s, _)| s.len().to_string()).collect();
    let moves: u64 = results.iter().map(|(_, k)| k).sum();
    let mut best = results
        .into_iter()
        .map(|(s, _)| s)
        .max_by_key(|s| s.len())
        .unwrap();
    best.sort();
    let c = best.len();
    let what = if k == 3 {
        format!("peelable set, against the digit base {}", m / 2 + 1)
    } else {
        "4-peelable set".to_string()
    };
    println!(
        "m = {m}, d = {d}: best {what}: {c} of {n} points, constant card^(1/d) = {:.4}; per-thread bests {}, {moves} kicks, {THREADS} threads, {:.1} s, seed {seed}",
        (c as f64).powf(1.0 / d as f64),
        sizes.join(" "),
        t0.elapsed().as_secs_f64()
    );
    let found = report(&cube, &best, blocks);
    if k == 3 || !m.is_multiple_of(3) {
        assert_eq!(found, 0, "a lifted slice holds a progression");
    }
}

fn control() {
    let cube = Cube::new(9, 1, 4);
    println!("m = 9, d = 1, k = 4: the set {{0, 1, 2, 3, 4, 7}}, lifted to Z_9^4");
    let found = report(&cube, &[0, 1, 2, 3, 4, 7], 4);
    assert!(found > 0, "the negative control must break");
}

// THE VERBS

fn optimum(m: usize, k: usize, free: bool, target: Option<usize>, seconds: Option<f64>) {
    let t0 = Instant::now();
    let grid = Grid::new(m, k, free);
    let what = match (k, free) {
        (3, false) => "peelable".to_string(),
        (_, false) => format!("{k}-peelable"),
        (_, true) => format!("{k}-term-progression-free"),
    };
    let n = grid.n;
    let lower = greedy(&grid, 200);
    let mut target = target.unwrap_or(lower).max(lower);
    println!(
        "m = {m}: {what} subsets of Z_{m}^2, n = {n}, group {}, a line holds at most {}, greedy lower bound {lower}, target {target}",
        grid.group(),
        grid.maxpeel[(1 << m) - 1]
    );
    loop {
        let (st, stopped) = search(&grid, target, false, seconds);
        let top = (0..=n).rev().find(|&k| st.hist[k] > 0);
        if stopped {
            println!(
                "  aborted after {:.1} s at {} nodes and {} canonicity tests: no value at m = {m}",
                t0.elapsed().as_secs_f64(),
                st.nodes,
                st.tests
            );
            return;
        }
        let Some(top) = top else {
            println!(
                "  no {what} set of {target} or more points, {} nodes, {} tests; target lowered",
                st.nodes, st.tests
            );
            target -= 1;
            continue;
        };
        let g = gcd(top, n);
        println!("  every {what} set of {target} or more points enumerated up to affine symmetry: {} nodes, {} canonicity tests, {} threads, {:.1} s", st.nodes, st.tests, THREADS, t0.elapsed().as_secs_f64());
        for k in target..=top {
            if st.hist[k] > 0 {
                let kept: Vec<&(Mask, u64)> = st
                    .found
                    .iter()
                    .filter(|(s, _)| s.count_ones() as usize == k)
                    .collect();
                let orbits: u128 = kept.iter().map(|(_, a)| grid.group() / *a as u128).sum();
                if kept.len() as u64 == st.hist[k] {
                    println!(
                        "  {k} points: {} sets up to symmetry, {orbits} sets",
                        st.hist[k]
                    );
                } else {
                    println!(
                        "  {k} points: {} sets up to symmetry, orbit total not computed",
                        st.hist[k]
                    );
                }
            }
        }
        let alpha = top as f64 / n as f64;
        let c = 2.0 * (1.0 / alpha).log2().sqrt();
        if k == 3 && !free {
            println!(
                "  alpha_{m} = {top}/{n} = {}/{} = {alpha:.5}, C(alpha) = {c:.4}, sqrt(card T) = {:.4}",
                top / g,
                n / g,
                (top as f64).sqrt()
            );
        } else {
            println!(
                "  the largest {what} subset of Z_{m}^2 has {top} points, sqrt(card T) = {:.4}",
                (top as f64).sqrt()
            );
        }
        let mut maxes: Vec<(Mask, u64)> = st
            .found
            .iter()
            .copied()
            .filter(|(s, _)| s.count_ones() as usize == top)
            .collect();
        maxes.sort();
        for (s, aut) in &maxes {
            let level = grid.rounds(*s).expect("maximiser peels");
            let rounds = *level.iter().filter(|&&l| l != 255).max().unwrap() as u32 + 1;
            let g: Vec<u64> = level
                .iter()
                .map(|&l| {
                    if l == 255 {
                        0
                    } else {
                        3u64.pow(rounds - 1 - l as u32)
                    }
                })
                .collect();
            assert!(grid.cert_holds(*s, &g));
            let mut pts = [0u8; 128];
            let k = bits(*s, &mut pts);
            let coords: Vec<String> = pts[..k]
                .iter()
                .map(|&v| format!("({},{})", v as usize / m, v as usize % m))
                .collect();
            println!(
                "  maximiser {} rows b = {}..0, automorphisms {aut}, orbit {}",
                grid.bitmap(*s, None).join("/"),
                m - 1,
                grid.group() / *aut as u128
            );
            println!("    points {}", coords.join(" "));
            println!("    peel rounds {} (digit = round removed), certificate g = 3^({} - round) checked on every progression", grid.bitmap(*s, Some(&level)).join("/"), rounds - 1);
        }
        return;
    }
}

fn witness(m: usize, seconds: Option<f64>) {
    let grid = Grid::new(m, 3, false);
    let mut target = greedy(&grid, 200);
    loop {
        let (st, stopped) = search(&grid, target, false, seconds);
        assert!(!stopped, "no value at m = {m}");
        let Some(top) = (0..=grid.n).rev().find(|&k| st.hist[k] > 0) else {
            target -= 1;
            continue;
        };
        let maxes: Vec<Mask> = st
            .found
            .iter()
            .map(|&(s, _)| s)
            .filter(|s| s.count_ones() as usize == top)
            .collect();
        assert_eq!(
            maxes.len() as u64,
            st.hist[top],
            "maximisers beyond the keep"
        );
        let mask = maxes.into_iter().min().unwrap();
        let level = grid.rounds(mask).expect("maximiser peels");
        let rounds = *level.iter().filter(|&&l| l != 255).max().unwrap() as u32 + 1;
        let g: Vec<u64> = level
            .iter()
            .map(|&l| {
                if l == 255 {
                    0
                } else {
                    3u64.pow(rounds - 1 - l as u32)
                }
            })
            .collect();
        assert!(grid.cert_holds(mask, &g));
        let mut pts = [0u8; 128];
        let k = bits(mask, &mut pts);
        let coords: Vec<String> = pts[..k]
            .iter()
            .map(|&v| format!("[{},{}]", v as usize / m, v as usize % m))
            .collect();
        println!(
            "{{\"m\":{m},\"rounds\":{rounds},\"points\":[{}]}}",
            coords.join(",")
        );
        return;
    }
}

fn census(m: usize, seconds: Option<f64>) {
    let t0 = Instant::now();
    let grid = Grid::new(m, 3, false);
    let (st, stopped) = search(&grid, 0, true, seconds);
    if stopped {
        println!(
            "m = {m}: census aborted after {:.1} s at {} nodes, nothing counted",
            t0.elapsed().as_secs_f64(),
            st.nodes
        );
        return;
    }
    let top = (0..=grid.n).rev().find(|&k| st.hist[k] > 0).unwrap();
    let by_size: Vec<String> = (0..=top).map(|k| st.hist[k].to_string()).collect();
    println!("m = {m}: every peelable set enumerated up to affine symmetry: {} orbits, {} sets, {} containing the origin, largest {top}, {} nodes, {:.1} s", st.orbits, st.sets, st.origin, st.nodes, t0.elapsed().as_secs_f64());
    println!("  orbits by size 0..{top}: {}", by_size.join(" "));
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let verb = args.first().map(String::as_str).unwrap_or("optimum");
    let mut ms = Vec::new();
    let mut target = None;
    let mut seconds = None;
    let mut dim = 2;
    let mut seed = 1;
    let mut k = 3;
    let mut free = false;
    let mut blocks = 2;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--k" => {
                k = args[i + 1].parse().unwrap();
                i += 1;
            }
            "--blocks" => {
                blocks = args[i + 1].parse().unwrap();
                i += 1;
            }
            "--free" => free = true,
            "--target" => {
                target = Some(args[i + 1].parse().unwrap());
                i += 1;
            }
            "--seconds" => {
                seconds = Some(args[i + 1].parse().unwrap());
                i += 1;
            }
            "--dim" => {
                dim = args[i + 1].parse().unwrap();
                i += 1;
            }
            "--seed" => {
                seed = args[i + 1].parse().unwrap();
                i += 1;
            }
            v => ms.push(v.parse::<usize>().unwrap()),
        }
        i += 1;
    }
    if verb == "control" {
        control();
        return;
    }
    if ms.is_empty() {
        ms = vec![3, 5, 7];
    }
    for &m in &ms {
        match verb {
            "optimum" => optimum(m, k, free, target, seconds),
            "census" => census(m, seconds),
            "witness" => witness(m, seconds),
            "climb" => climb(m, dim, k, seconds.unwrap_or(30.0), seed, blocks),
            other => panic!("unknown verb {other}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reflect_free(grid: &Grid, set: Mask) -> bool {
        let (m, n) = (grid.m, grid.n);
        let at = |a: usize, b: usize| (a % m) * m + (b % m);
        let mut live = set;
        loop {
            if live == 0 {
                return true;
            }
            let mut rem: Mask = 0;
            for y in 0..n {
                if live >> y & 1 == 0 {
                    continue;
                }
                let (ya, yb) = (y / m, y % m);
                let blocked = (0..n).any(|x| {
                    let z = at(2 * ya + m - x / m, 2 * yb + m - x % m);
                    x != y && live >> x & 1 == 1 && live >> z & 1 == 1
                });
                if !blocked {
                    rem |= 1 << y;
                }
            }
            if rem == 0 {
                return false;
            }
            live &= !rem;
        }
    }

    fn random_sets(n: usize, k: usize, count: usize) -> Vec<Mask> {
        let mut seed = 12345u64;
        (0..count)
            .map(|_| {
                let mut s: Mask = 0;
                while (s.count_ones() as usize) < k {
                    seed ^= seed << 13;
                    seed ^= seed >> 7;
                    seed ^= seed << 17;
                    s |= 1 << (seed % n as u64);
                }
                s
            })
            .collect()
    }

    #[test]
    fn peel_agrees_with_reflection() {
        for m in [3, 5, 7] {
            let grid = Grid::new(m, 3, false);
            for k in [4, 8, 12] {
                for s in random_sets(grid.n, k.min(grid.n - 1), 300) {
                    assert_eq!(
                        grid.peelable(s),
                        reflect_free(&grid, s),
                        "m = {m}, set {s:b}"
                    );
                }
            }
        }
    }

    #[test]
    fn one_dimension_is_half_plus() {
        for m in [3, 5, 7, 11] {
            let grid = Grid::new(m, 3, false);
            let full = (1usize << m) - 1;
            assert_eq!(grid.maxpeel[full] as usize, m.div_ceil(2));
        }
    }

    #[test]
    fn canonical_sets_match_brute_force_at_three() {
        let grid = Grid::new(3, 3, false);
        let mut orbits = 0u128;
        let mut sets = 0u128;
        let mut origin = 0u128;
        for s in 0..(1 as Mask) << 9 {
            let mut aut = 0;
            if grid.peelable(s) && grid.canonical(s, &mut aut) {
                orbits += 1;
                sets += grid.group() / aut as u128;
                origin += grid.group() / aut as u128 * s.count_ones() as u128 / 9;
            }
        }
        let (st, stopped) = search(&grid, 0, true, None);
        assert!(!stopped);
        assert_eq!((st.orbits, st.sets, st.origin), (orbits, sets, origin));
        assert_eq!(origin, 57);
        let all = (0..(1 as Mask) << 9).filter(|&s| grid.peelable(s)).count() as u128;
        assert_eq!(sets, all);
    }

    #[test]
    fn line_bound_is_sound_on_random_sets() {
        let grid = Grid::new(5, 3, false);
        for s in random_sets(grid.n, 10, 200) {
            if !grid.peelable(s) {
                continue;
            }
            let mut pts = [0u8; 128];
            let k = bits(s, &mut pts);
            let half: Mask = pts[..k / 2].iter().fold(0, |a, &v| a | 1 << v);
            assert!(k <= half.count_ones() as usize + grid.line_bound(half, s));
        }
    }

    #[test]
    fn cube_peel_agrees_with_brute_force() {
        for (m, d) in [(4, 2), (5, 2), (6, 2), (7, 2), (3, 3), (4, 3)] {
            let cube = Cube::new(m, d, 3);
            let n = cube.n;
            let at = |v: usize, i: usize| cube.coords[v * d + i] as usize;
            let mid = |x: usize, y: usize, z: usize| {
                (0..d).all(|i| (at(x, i) + at(z, i) + 2 * m - 2 * at(y, i)) % m == 0)
            };
            for k in [4, 7, 10] {
                for s in random_sets(n, k.min(n - 1), 200) {
                    let set: Vec<usize> = (0..n).filter(|&v| s >> v & 1 == 1).collect();
                    let mut live = set.clone();
                    let brute = loop {
                        if live.is_empty() {
                            break true;
                        }
                        let free = live.iter().position(|&y| {
                            !live
                                .iter()
                                .any(|&x| live.iter().any(|&z| x < z && mid(x, y, z)))
                        });
                        match free {
                            Some(i) => {
                                live.remove(i);
                            }
                            None => break false,
                        }
                    };
                    assert_eq!(
                        cube.rounds(&set).is_ok(),
                        brute,
                        "m = {m}, d = {d}, set {set:?}"
                    );
                }
            }
        }
    }

    fn four_inner(m: usize, d: usize, live: &[Vec<usize>]) -> Vec<bool> {
        let mut inner = vec![false; live.len()];
        for x in live {
            for (j, y) in live.iter().enumerate() {
                let t: Vec<Vec<usize>> = (0..4)
                    .map(|c| (0..d).map(|i| (x[i] + c * (y[i] + m - x[i])) % m).collect())
                    .collect();
                let distinct = (0..4).all(|a| (a + 1..4).all(|b| t[a] != t[b]));
                if distinct && live.contains(&t[2]) && live.contains(&t[3]) {
                    inner[j] = true;
                    inner[live.iter().position(|v| *v == t[2]).unwrap()] = true;
                }
            }
        }
        inner
    }

    #[test]
    fn four_peel_agrees_with_brute_force() {
        for (m, d) in [
            (4, 2),
            (5, 1),
            (5, 2),
            (6, 2),
            (7, 1),
            (7, 2),
            (9, 1),
            (9, 2),
            (4, 3),
        ] {
            let cube = Cube::new(m, d, 4);
            let grid = (m % 2 == 1 && m != 9 && d == 2).then(|| Grid::new(m, 4, false));
            let n = cube.n;
            for k in [5, 9, 14] {
                for s in random_sets(n, k.min(n - 1), 200) {
                    let set: Vec<usize> = (0..n).filter(|&v| s >> v & 1 == 1).collect();
                    let mut live: Vec<Vec<usize>> = set
                        .iter()
                        .map(|&v| {
                            cube.coords[v * d..(v + 1) * d]
                                .iter()
                                .map(|&c| c as usize)
                                .collect()
                        })
                        .collect();
                    let brute = loop {
                        if live.is_empty() {
                            break true;
                        }
                        match four_inner(m, d, &live).iter().position(|&b| !b) {
                            Some(i) => {
                                live.remove(i);
                            }
                            None => break false,
                        }
                    };
                    assert_eq!(
                        cube.rounds(&set).is_ok(),
                        brute,
                        "m = {m}, d = {d}, set {set:?}"
                    );
                    if let Some(g) = &grid {
                        assert_eq!(g.peelable(s), brute, "m = {m}, set {set:?}");
                    }
                }
            }
        }
    }

    #[test]
    fn four_free_agrees_with_brute_force() {
        for m in [5, 7] {
            let grid = Grid::new(m, 4, true);
            for k in [4, 8, 12, 16] {
                for s in random_sets(grid.n, k, 300) {
                    let live: Vec<Vec<usize>> = (0..grid.n)
                        .filter(|&v| s >> v & 1 == 1)
                        .map(|v| vec![v / m, v % m])
                        .collect();
                    let brute = !four_inner(m, 2, &live).contains(&true);
                    assert_eq!(grid.peelable(s), brute, "m = {m}, set {s:b}");
                }
            }
        }
    }

    #[test]
    fn four_peel_lift_breaks_at_nine() {
        let cube = Cube::new(9, 1, 4);
        let set = [0, 1, 2, 3, 4, 7];
        let level = cube.rounds(&set).unwrap();
        let r = level.iter().max().unwrap() + 1;
        let g: Vec<u128> = level.iter().map(|&l| 3u128.pow(r - 1 - l)).collect();
        assert!(lift_check(&cube, &set, &g, 4).2 > 0);
    }
}
