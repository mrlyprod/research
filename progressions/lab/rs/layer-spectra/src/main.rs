use std::collections::HashMap;
use std::env;
use std::time::Instant;

const P: i32 = 7;

fn md(x: i32) -> u8 {
    x.rem_euclid(P) as u8
}

// VECTORS OVER F_7

struct Space {
    size: usize,
    add: Vec<u32>,
}

impl Space {
    fn new(dim: usize) -> Space {
        let size = 7usize.pow(dim as u32);
        let digits = |mut i: usize| -> Vec<u8> {
            (0..dim)
                .map(|_| {
                    let d = (i % 7) as u8;
                    i /= 7;
                    d
                })
                .collect()
        };
        let mut add = vec![0u32; size * size];
        for i in 0..size {
            let a = digits(i);
            for j in 0..size {
                let b = digits(j);
                let mut idx = 0usize;
                for k in (0..dim).rev() {
                    idx = idx * 7 + ((a[k] + b[k]) % 7) as usize;
                }
                add[i * size + j] = idx as u32;
            }
        }
        Space { size, add }
    }
    #[inline(always)]
    fn plus(&self, i: usize, j: usize) -> usize {
        self.add[i * self.size + j] as usize
    }
}

// SPECTRA

fn pair_bit(d0: u8, d1: u8) -> u64 {
    1u64 << (d0 as u32 * 7 + d1 as u32)
}

fn spectrum(space: &Space, f: &[u8]) -> u64 {
    let mut mask = 0u64;
    for w in 0..space.size {
        for c in 1..space.size {
            let x1 = space.plus(w, c);
            let x2 = space.plus(x1, c);
            let x3 = space.plus(x2, c);
            let (a, b, cc, d) = (f[w] as i32, f[x1] as i32, f[x2] as i32, f[x3] as i32);
            mask |= pair_bit(md(a - 2 * b + cc), md(b - 2 * cc + d));
        }
    }
    mask
}

fn sumset(a: u64, b: u64) -> u64 {
    let mut out = 0u64;
    for x in 0..7u8 {
        for y in 0..7u8 {
            if a & pair_bit(x, y) == 0 {
                continue;
            }
            for u in 0..7u8 {
                for v in 0..7u8 {
                    if b & pair_bit(u, v) != 0 {
                        out |= pair_bit((x + u) % 7, (y + v) % 7);
                    }
                }
            }
        }
    }
    out
}

fn swapped(mask: u64) -> u64 {
    let mut out = 0u64;
    for x in 0..7u8 {
        for y in 0..7u8 {
            if mask & pair_bit(x, y) != 0 {
                out |= pair_bit(y, x);
            }
        }
    }
    out
}

fn diagonal() -> u64 {
    (0..7u8).map(|t| pair_bit(t, t)).fold(0, |m, b| m | b)
}

fn show(mask: u64) -> String {
    let mut parts = Vec::new();
    for x in 0..7u8 {
        for y in 0..7u8 {
            if mask & pair_bit(x, y) != 0 {
                parts.push(format!("({x},{y})"));
            }
        }
    }
    parts.join(" ")
}

// UNIVARIATE

fn degree(h: &[u8]) -> usize {
    let mut deg = 0;
    for k in 1..7usize {
        let mut s = 0i32;
        for x in 0..7i32 {
            s += h[x as usize] as i32 * x.pow(6 - k as u32);
        }
        if md(-s) != 0 {
            deg = k;
        }
    }
    deg
}

fn univariate() -> Vec<u64> {
    let space = Space::new(1);
    let diag = diagonal();
    let mut good_by_degree = [0u64; 7];
    let mut cubic_difference_by_degree = [0u64; 7];
    let mut size_hist: HashMap<u32, u64> = HashMap::new();
    let mut spectra: HashMap<u64, (u64, usize, [u8; 7])> = HashMap::new();
    let mut h = [0u8; 7];
    for code in 0..823543usize {
        let mut c = code;
        for k in 0..7 {
            h[k] = (c % 7) as u8;
            c /= 7;
        }
        let mask = spectrum(&space, &h);
        assert_eq!(mask, swapped(mask));
        if mask & 1 != 0 {
            continue;
        }
        let deg = degree(&h);
        good_by_degree[deg] += 1;
        if mask & diag == 0 {
            cubic_difference_by_degree[deg] += 1;
        }
        *size_hist.entry(mask.count_ones()).or_insert(0) += 1;
        let e = spectra.entry(mask).or_insert((0, deg, h));
        e.0 += 1;
        if deg < e.1 {
            e.1 = deg;
            e.2 = h;
        }
    }
    println!(
        "good functions F_7 -> F_7: {} by degree {:?}",
        good_by_degree.iter().sum::<u64>(),
        good_by_degree
    );
    println!(
        "with nonzero third difference on every progression: {} by degree {:?}",
        cubic_difference_by_degree.iter().sum::<u64>(),
        cubic_difference_by_degree
    );
    let mut quintics = std::collections::HashSet::new();
    for a in 1..7i32 {
        for t in 0..7i32 {
            for c2 in 0..7i32 {
                for c1 in 0..7i32 {
                    for c0 in 0..7i32 {
                        let mut g = [0u8; 7];
                        for x in 0..7i32 {
                            g[x as usize] = md(a * (x + t).pow(5) + c2 * x * x + c1 * x + c0);
                        }
                        let mask = spectrum(&space, &g);
                        assert!(mask & 1 == 0 && mask & diag == 0 && degree(&g) == 5);
                        quintics.insert(g);
                    }
                }
            }
        }
    }
    println!(
        "the quintics a (x + t)^5 + quadratic: {} distinct, all with nonzero third differences",
        quintics.len()
    );
    assert_eq!(quintics.len() as u64, cubic_difference_by_degree[5]);
    let mut sizes: Vec<_> = size_hist.into_iter().collect();
    sizes.sort();
    let shown: Vec<String> = sizes.iter().map(|(s, n)| format!("{s}:{n}")).collect();
    println!("spectrum size:count {}", shown.join(" "));
    println!("distinct spectra {}", spectra.len());
    let mut small: Vec<_> = spectra
        .iter()
        .filter(|(m, _)| m.count_ones() <= 9)
        .collect();
    small.sort_by_key(|(m, _)| (m.count_ones(), **m));
    for (m, (n, deg, h)) in small {
        println!(
            "spectrum of {} points, {} functions, least degree {}, example {:?}: {}",
            m.count_ones(),
            n,
            deg,
            h,
            show(*m)
        );
    }
    spectra.keys().copied().collect()
}

// LAYERS OF UNIVARIATE FUNCTIONS

fn layers(spectra: &[u64], want: usize) {
    let mut masks: Vec<u64> = spectra.to_vec();
    masks.sort_by_key(|m| (m.count_ones(), *m));
    let n = masks.len();
    let mut best = 0usize;
    let mut best_choice: Vec<usize> = Vec::new();
    let mut count_at = vec![0u64; want + 1];
    fn dfs(
        masks: &[u64],
        start: usize,
        reach: u64,
        chosen: &mut Vec<usize>,
        want: usize,
        best: &mut usize,
        best_choice: &mut Vec<usize>,
        count_at: &mut Vec<u64>,
        budget: &mut u64,
    ) {
        if chosen.len() > *best {
            *best = chosen.len();
            *best_choice = chosen.clone();
        }
        if chosen.len() == want || *budget == 0 {
            return;
        }
        for i in start..masks.len() {
            *budget -= 1;
            if *budget == 0 {
                return;
            }
            let s = masks[i];
            let next = reach | s | sumset(reach, s);
            if next & 1 != 0 {
                continue;
            }
            chosen.push(i);
            count_at[chosen.len()] += 1;
            dfs(
                masks,
                i,
                next,
                chosen,
                want,
                best,
                best_choice,
                count_at,
                budget,
            );
            chosen.pop();
        }
    }
    let mut chosen = Vec::new();
    let mut budget = 200_000_000u64;
    dfs(
        &masks,
        0,
        0,
        &mut chosen,
        want,
        &mut best,
        &mut best_choice,
        &mut count_at,
        &mut budget,
    );
    println!(
        "univariate layers from {n} spectra: deepest zero-sum-free stack {best} of {want} wanted, stacks found by depth {:?}, budget left {budget}",
        &count_at[1..]
    );
    for &i in &best_choice {
        println!(
            "  layer spectrum {} points: {}",
            masks[i].count_ones(),
            show(masks[i])
        );
    }
}

// LIFT

fn lift(m: usize) {
    let space = Space::new(m);
    let size = space.size;
    let digit = |i: usize, k: usize| (i / 7usize.pow(k as u32)) % 7;
    let mut pivots: Vec<(usize, Vec<u8>)> = Vec::new();
    let mut rows = 0usize;
    for c in 1..size {
        let k = (0..m).find(|&k| digit(c, k) != 0).unwrap();
        if digit(c, k) != 1 {
            continue;
        }
        for w in 0..size {
            if digit(w, k) != 0 {
                continue;
            }
            for j in 0..3u32 {
                let mut row = vec![0u8; size];
                let mut x = w;
                for t in 0..7i32 {
                    row[x] = md(row[x] as i32 + t.pow(j));
                    x = space.plus(x, c);
                }
                rows += 1;
                for (col, prow) in &pivots {
                    let f = row[*col];
                    if f != 0 {
                        for i in 0..size {
                            row[i] = md(row[i] as i32 - f as i32 * prow[i] as i32);
                        }
                    }
                }
                if let Some(col) = (0..size).find(|&i| row[i] != 0) {
                    let inv = (1..7).find(|&v| (v * row[col] as i32) % 7 == 1).unwrap();
                    for i in 0..size {
                        row[i] = md(row[i] as i32 * inv);
                    }
                    pivots.push((col, row));
                }
            }
        }
    }
    let dim = size - pivots.len();
    let binom = (m + 1) * (m + 2) * (m + 3) / 6;
    println!(
        "lift m = {m}: {rows} line constraints, rank {}, the functions cubic on every line form a space of dimension {dim}, the cubic polynomials one of dimension {binom}",
        pivots.len()
    );
    assert_eq!(dim, binom);
}

// EXTENSION FIELDS

struct Field {
    q: usize,
    mul: Vec<u16>,
    log: Vec<usize>,
    exp: Vec<usize>,
    trace: Vec<u8>,
}

impl Field {
    fn new(deg: usize) -> Field {
        let q = 7usize.pow(deg as u32);
        let digits = |mut i: usize| -> [i32; 3] {
            let mut d = [0i32; 3];
            for k in 0..deg {
                d[k] = (i % 7) as i32;
                i /= 7;
            }
            d
        };
        let index = |d: [i32; 3]| -> usize {
            let mut i = 0usize;
            for k in (0..deg).rev() {
                i = i * 7 + md(d[k]) as usize;
            }
            i
        };
        let mut mul = vec![0u16; q * q];
        for i in 0..q {
            let x = digits(i);
            for j in 0..q {
                let y = digits(j);
                let prod = if deg == 2 {
                    [x[0] * y[0] + 3 * x[1] * y[1], x[0] * y[1] + x[1] * y[0], 0]
                } else {
                    let (a, b, c) = (x[0], x[1], x[2]);
                    let (d, e, f) = (y[0], y[1], y[2]);
                    [
                        a * d + 2 * (b * f + c * e),
                        a * e + b * d + 2 * c * f,
                        a * f + b * e + c * d,
                    ]
                };
                mul[i * q + j] = index(prod) as u16;
            }
        }
        let mut gen = 0;
        let mut exp = Vec::new();
        for g in 2..q {
            let mut seen = vec![false; q];
            let mut x = 1usize;
            let mut order = 0;
            loop {
                x = mul[x * q + g] as usize;
                order += 1;
                if seen[x] {
                    break;
                }
                seen[x] = true;
                if x == 1 {
                    break;
                }
            }
            if order == q - 1 {
                gen = g;
                break;
            }
        }
        assert!(gen != 0);
        let mut log = vec![0usize; q];
        let mut x = 1usize;
        for i in 0..q - 1 {
            exp.push(x);
            log[x] = i;
            x = mul[x * q + gen] as usize;
        }
        let power = |v: usize, k: usize| -> usize {
            if v == 0 {
                0
            } else {
                exp[(log[v] * k) % (q - 1)]
            }
        };
        let mut trace = vec![0u8; q];
        let space = Space::new(deg);
        for v in 0..q {
            let mut t = 0usize;
            for i in 0..deg {
                t = space.plus(t, power(v, 7usize.pow(i as u32)));
            }
            assert!(t < 7);
            trace[v] = t as u8;
        }
        Field {
            q,
            mul,
            log,
            exp,
            trace,
        }
    }
    fn power(&self, v: usize, k: usize) -> usize {
        if v == 0 {
            0
        } else {
            self.exp[(self.log[v] * k) % (self.q - 1)]
        }
    }
    fn times(&self, a: usize, b: usize) -> usize {
        self.mul[a * self.q + b] as usize
    }
}

fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

struct Monomial {
    k: usize,
    class: usize,
    mask: u64,
}

fn monomials(deg: usize) -> Vec<Monomial> {
    let t0 = Instant::now();
    let field = Field::new(deg);
    let space = Space::new(deg);
    let q = field.q;
    let mut good = Vec::new();
    let mut tried = 0;
    for k in 1..q {
        for class in 0..gcd(k, q - 1) {
            let a = field.exp[class];
            let f: Vec<u8> = (0..q)
                .map(|v| field.trace[field.times(a, field.power(v, k))])
                .collect();
            tried += 1;
            let mask = spectrum(&space, &f);
            if mask & 1 == 0 {
                good.push(Monomial { k, class, mask });
            }
        }
    }
    let mut hist: HashMap<u32, u64> = HashMap::new();
    for g in &good {
        *hist.entry(g.mask.count_ones()).or_insert(0) += 1;
    }
    let mut sizes: Vec<_> = hist.into_iter().collect();
    sizes.sort();
    let shown: Vec<String> = sizes.iter().map(|(s, n)| format!("{s}:{n}")).collect();
    println!(
        "F_{q}: {tried} functions Tr(a v^k) up to v -> lambda v, {} good, spectrum size:count {}, {:.1}s",
        good.len(),
        shown.join(" "),
        t0.elapsed().as_secs_f64()
    );
    let norm = (q - 1) / 6;
    for g in good
        .iter()
        .filter(|g| (g.k == norm && g.class == 0) || g.mask.count_ones() <= 12)
    {
        println!(
            "  k = {} class {} spectrum {} points{}: {}",
            g.k,
            g.class,
            g.mask.count_ones(),
            if g.k == norm { " (the norm)" } else { "" },
            show(g.mask)
        );
    }
    good
}

fn stacks(label: &str, parts: &[&[Monomial]]) {
    let mut found = 0u64;
    let mut example = None;
    fn rec(
        parts: &[&[Monomial]],
        depth: usize,
        reach: u64,
        chosen: &mut Vec<(usize, usize)>,
        found: &mut u64,
        example: &mut Option<Vec<(usize, usize)>>,
    ) {
        if depth == parts.len() {
            *found += 1;
            if example.is_none() {
                *example = Some(chosen.clone());
            }
            return;
        }
        for g in parts[depth] {
            let next = reach | g.mask | sumset(reach, g.mask);
            if next & 1 != 0 {
                continue;
            }
            chosen.push((g.k, g.class));
            rec(parts, depth + 1, next, chosen, found, example);
            chosen.pop();
        }
    }
    rec(parts, 0, 0, &mut Vec::new(), &mut found, &mut example);
    println!(
        "{label}: {found} zero-sum-free stacks{}",
        match example {
            Some(e) => format!(", for instance (k, class) = {e:?}"),
            None => String::new(),
        }
    );
}

// THE PLANE WITH A SEED LINE

struct Csp {
    n: usize,
    windows: Vec<[usize; 4]>,
    at: Vec<Vec<usize>>,
    value: Vec<u8>,
    domain: Vec<u8>,
    trail: Vec<(usize, u8)>,
    nodes: u64,
}

impl Csp {
    fn new(dim: usize) -> Csp {
        let space = Space::new(dim);
        let n = space.size;
        let mut windows = Vec::new();
        for w in 0..n {
            for c in 1..n {
                let x1 = space.plus(w, c);
                let x2 = space.plus(x1, c);
                let x3 = space.plus(x2, c);
                if w < x3 {
                    windows.push([w, x1, x2, x3]);
                }
            }
        }
        assert_eq!(windows.len(), n * (n - 1) / 2);
        let mut at = vec![Vec::new(); n];
        for (i, win) in windows.iter().enumerate() {
            for &cell in win {
                at[cell].push(i);
            }
        }
        Csp {
            n,
            windows,
            at,
            value: vec![255; n],
            domain: vec![127; n],
            trail: Vec::new(),
            nodes: 0,
        }
    }
    fn assign(&mut self, cell: usize, v: u8) -> bool {
        self.value[cell] = v;
        let coef = [1i32, -3, 3, -1];
        for &wi in &self.at[cell].clone() {
            let win = self.windows[wi];
            let mut free = None;
            let mut sum = 0i32;
            let mut count = 0;
            for j in 0..4 {
                let x = self.value[win[j]];
                if x == 255 {
                    free = Some(j);
                } else {
                    count += 1;
                    sum += coef[j] * x as i32;
                }
            }
            if count == 4 {
                if md(sum) == 0 {
                    return false;
                }
                continue;
            }
            if count != 3 {
                continue;
            }
            let j = free.unwrap();
            let inv = if coef[j].rem_euclid(7) == 1 {
                1
            } else if coef[j].rem_euclid(7) == 6 {
                6
            } else if coef[j].rem_euclid(7) == 3 {
                5
            } else {
                2
            };
            let forbidden = md(-sum * inv);
            let other = win[j];
            let bit = 1u8 << forbidden;
            if self.domain[other] & bit != 0 {
                self.trail.push((other, self.domain[other]));
                self.domain[other] &= !bit;
                if self.domain[other] == 0 {
                    return false;
                }
            }
        }
        true
    }
    fn search(&mut self, limit: u64, solutions: &mut Vec<Vec<u8>>, cap: usize) -> bool {
        self.nodes += 1;
        if self.nodes > limit {
            return false;
        }
        let mut pick = None;
        for cell in 0..self.n {
            if self.value[cell] == 255 {
                let size = self.domain[cell].count_ones();
                if pick.map_or(true, |(_, s)| size < s) {
                    pick = Some((cell, size));
                }
            }
        }
        let Some((cell, _)) = pick else {
            solutions.push(self.value.clone());
            return solutions.len() >= cap;
        };
        let dom = self.domain[cell];
        for v in 0..7u8 {
            if dom >> v & 1 == 0 {
                continue;
            }
            let mark = self.trail.len();
            let ok = self.assign(cell, v);
            let done = ok && self.search(limit, solutions, cap);
            while self.trail.len() > mark {
                let (c, d) = self.trail.pop().unwrap();
                self.domain[c] = d;
            }
            self.value[cell] = 255;
            if done || self.nodes > limit {
                return done;
            }
        }
        false
    }
}

fn plane(dim: usize, seed: &[u8; 7], label: &str) {
    let t0 = Instant::now();
    let mut csp = Csp::new(dim);
    let mut alive = true;
    for x in 0..7 {
        alive = alive && csp.assign(x, seed[x]);
    }
    let mut solutions = Vec::new();
    let limit = 150_000_000u64;
    let found = if alive {
        csp.search(limit, &mut solutions, 1)
    } else {
        false
    };
    let status = if found {
        "extends".to_string()
    } else if csp.nodes > limit {
        "undecided at the node cap".to_string()
    } else {
        "extends to no function".to_string()
    };
    println!(
        "{label} seed {:?} on the axis of F_7^{dim}: {status}, {} nodes, {:.1}s",
        seed,
        csp.nodes,
        t0.elapsed().as_secs_f64()
    );
    if let Some(sol) = solutions.first() {
        println!("  completion {:?}", sol);
    }
}

fn seeds() -> Vec<([u8; 7], String)> {
    let space = Space::new(1);
    let mut out: Vec<([u8; 7], String)> = Vec::new();
    let mut seen: Vec<[u8; 7]> = Vec::new();
    let mut h = [0u8; 7];
    for code in 0..823543usize {
        let mut c = code;
        for k in 0..7 {
            h[k] = (c % 7) as u8;
            c /= 7;
        }
        if h[0] != 0 || h[1] != 0 || h[2] != 0 {
            continue;
        }
        let mask = spectrum(&space, &h);
        if mask & 1 != 0 || mask & diagonal() != 0 {
            continue;
        }
        let mut canon = h;
        for a in 1..7i32 {
            for b in 0..7i32 {
                for lam in 1..7i32 {
                    let g: Vec<u8> = (0..7i32)
                        .map(|x| md(lam * h[md(a * x + b) as usize] as i32))
                        .collect();
                    let mut q = [0u8; 7];
                    let (g0, g1, g2) = (g[0] as i32, g[1] as i32, g[2] as i32);
                    let c2 = md((g0 - 2 * g1 + g2) * 4);
                    let c1 = md(g1 - g0 - c2 as i32);
                    for x in 0..7i32 {
                        q[x as usize] =
                            md(g[x as usize] as i32 - c2 as i32 * x * x - c1 as i32 * x - g0);
                    }
                    if q < canon {
                        canon = q;
                    }
                }
            }
        }
        if !seen.contains(&canon) {
            seen.push(canon);
            out.push((canon, format!("degree {}", degree(&canon))));
        }
    }
    out
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let verb = args.first().map(|s| s.as_str()).unwrap_or("all");
    let t0 = Instant::now();
    if verb == "univariate" || verb == "layers" || verb == "all" {
        let spectra = univariate();
        println!("{:.1}s", t0.elapsed().as_secs_f64());
        if verb != "univariate" {
            layers(&spectra, 6);
            println!("{:.1}s", t0.elapsed().as_secs_f64());
        }
    }
    if verb == "lift" || verb == "all" {
        lift(2);
        lift(3);
        println!("{:.1}s", t0.elapsed().as_secs_f64());
    }
    if verb == "plane" || verb == "all" {
        let list = seeds();
        println!(
            "third-difference seeds up to affine maps, scaling and quadratics: {}",
            list.len()
        );
        for (seed, label) in &list {
            plane(2, seed, label);
        }
        println!("{:.1}s", t0.elapsed().as_secs_f64());
    }
    if verb == "monomial" || verb == "all" {
        let two = monomials(2);
        let three = monomials(3);
        stacks("F_49 + F_343", &[&two, &three]);
        stacks("F_343 + F_343", &[&three, &three]);
        stacks("F_49 + F_49 + F_49", &[&two, &two, &two]);
        stacks("F_49 + F_49", &[&two, &two]);
        println!("{:.1}s", t0.elapsed().as_secs_f64());
    }
}
