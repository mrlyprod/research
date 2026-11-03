use mrly_math::bang::factory::create;
use mrly_math::bang::Code;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::env;
use std::hash::{BuildHasherDefault, Hasher};
use std::time::Instant;

#[derive(Default)]
struct Fx(u64);

impl Hasher for Fx {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for &byte in bytes {
            self.write_u64(byte as u64);
        }
    }
    fn write_u64(&mut self, n: u64) {
        self.0 = (self.0.rotate_left(5) ^ n).wrapping_mul(0x517c_c1b7_2722_0a95);
    }
}

type Fast = BuildHasherDefault<Fx>;

#[derive(Clone, PartialEq, Eq)]
struct Pic {
    side: usize,
    cells: Vec<u8>,
}

impl Pic {
    fn one() -> Pic {
        Pic {
            side: 1,
            cells: vec![1],
        }
    }
    fn block(window: [u8; 4]) -> Pic {
        Pic {
            side: 2,
            cells: window.to_vec(),
        }
    }
    fn at(&self, row: usize, col: usize) -> u8 {
        self.cells[row * self.side + col]
    }
    fn pairs(&self) -> BTreeSet<[u8; 4]> {
        let mut out = BTreeSet::new();
        for row in 0..self.side - 1 {
            for col in 0..self.side - 1 {
                out.insert([
                    self.at(row, col),
                    self.at(row, col + 1),
                    self.at(row + 1, col),
                    self.at(row + 1, col + 1),
                ]);
            }
        }
        out
    }
}

#[derive(Clone)]
struct Design {
    base: usize,
    code: u32,
    tile: Vec<u8>,
}

impl Design {
    fn new(base: usize, code: u32) -> Design {
        let tile = (0..base * base)
            .map(|cell| ((code >> cell) & 1) as u8)
            .collect();
        Design { base, code, tile }
    }
    fn filled(&self) -> Vec<(usize, usize)> {
        let b = self.base;
        (0..b * b)
            .filter(|&cell| self.tile[cell] == 1)
            .map(|cell| (cell / b, cell % b))
            .collect()
    }
    fn empty(&self) -> bool {
        self.code == 0
    }
    fn full(&self) -> bool {
        self.tile.iter().all(|&cell| cell == 1)
    }
    fn degenerate(&self) -> bool {
        let cells = self.filled();
        let last = self.base - 1;
        !cells.is_empty()
            && (cells.iter().all(|c| c.0 == 0)
                || cells.iter().all(|c| c.0 == last)
                || cells.iter().all(|c| c.1 == 0)
                || cells.iter().all(|c| c.1 == last))
    }
    fn image(&self, g: usize) -> Design {
        let b = self.base;
        let mut code = 0u32;
        for (row, col) in self.filled() {
            let (mut r, mut c) = if g & 4 != 0 { (col, row) } else { (row, col) };
            for _ in 0..g & 3 {
                (r, c) = (c, b - 1 - r);
            }
            code |= 1 << (r * b + c);
        }
        Design::new(b, code)
    }
    fn canonical(&self) -> u32 {
        (0..8).map(|g| self.image(g).code).min().unwrap()
    }
    fn substitute(&self, pic: &Pic) -> Pic {
        let b = self.base;
        let side = pic.side * b;
        let mut cells = vec![0u8; side * side];
        for row in 0..side {
            for col in 0..side {
                cells[row * side + col] =
                    pic.at(row / b, col / b) & self.tile[(row % b) * b + col % b];
            }
        }
        Pic { side, cells }
    }
    fn render(&self, level: usize) -> Pic {
        let mut pic = Pic::one();
        for _ in 0..level {
            pic = self.substitute(&pic);
        }
        pic
    }
    fn language(&self) -> (Vec<[u8; 4]>, usize) {
        let mut set = self.render(1).pairs();
        let mut level = 1;
        loop {
            let mut next = BTreeSet::new();
            for window in &set {
                next.extend(self.substitute(&Pic::block(*window)).pairs());
            }
            assert!(
                next.is_superset(&set),
                "the 2 x 2 windows of a render shrank one level up"
            );
            if next == set {
                return (set.into_iter().collect(), level);
            }
            set = next;
            level += 1;
        }
    }
}

fn depth(base: usize, top: usize) -> usize {
    let mut n = 0;
    while base.pow(n as u32) + 1 < top {
        n += 1;
    }
    n
}

fn count(design: &Design, top: usize) -> Vec<usize> {
    let (pairs, _) = design.language();
    let n = depth(design.base, top);
    let pics: Vec<Pic> = pairs
        .iter()
        .map(|window| {
            let mut pic = Pic::block(*window);
            for _ in 0..n {
                pic = design.substitute(&pic);
            }
            pic
        })
        .collect();
    let side = 2 * design.base.pow(n as u32);
    let mut ids: Vec<Vec<u32>> = pics
        .iter()
        .map(|pic| pic.cells.iter().map(|&cell| cell as u32).collect())
        .collect();
    let values: BTreeSet<u8> = pics
        .iter()
        .flat_map(|pic| pic.cells.iter().copied())
        .collect();
    let mut out = vec![values.len()];
    for k in 1..top {
        let wide = side - k + 1;
        let narrow = wide - 1;
        let mut map: HashMap<u64, u32, Fast> = HashMap::default();
        for (pic, id) in pics.iter().zip(ids.iter_mut()) {
            let mut next = vec![0u32; narrow * narrow];
            for x in 0..narrow {
                for y in 0..narrow {
                    let key = (id[x * wide + y] as u64) << 34
                        | (id[(x + 1) * wide + y + 1] as u64) << 2
                        | (pic.at(x + k, y) as u64) << 1
                        | pic.at(x, y + k) as u64;
                    let fresh = map.len() as u32;
                    next[x * narrow + y] = *map.entry(key).or_insert(fresh);
                }
            }
            *id = next;
        }
        assert!(map.len() < 1 << 30, "window names overflow the key");
        out.push(map.len());
    }
    out
}

fn brute(pic: &Pic, top: usize) -> Vec<usize> {
    let side = pic.side;
    let mut out = Vec::new();
    for k in 1..=top.min(side).min(64) {
        let mask = if k == 64 { u64::MAX } else { (1u64 << k) - 1 };
        let wide = side - k + 1;
        let mut words = vec![0u64; side * wide];
        for row in 0..side {
            let mut word = 0u64;
            for col in 0..side {
                word = ((word << 1) | pic.at(row, col) as u64) & mask;
                if col + 1 >= k {
                    words[row * wide + col + 1 - k] = word;
                }
            }
        }
        let mut seen: HashSet<Vec<u64>, Fast> = HashSet::default();
        for x in 0..wide {
            for y in 0..wide {
                seen.insert((0..k).map(|t| words[(x + t) * wide + y]).collect());
            }
        }
        out.push(seen.len());
    }
    out
}

fn crate_render(design: &Design, level: usize) -> Pic {
    let b = design.base;
    let tensor = create(Code::from(design.code as u64), b, 2, b, level).unwrap();
    Pic {
        side: b.pow(level as u32),
        cells: tensor.bytes().unwrap().to_vec(),
    }
}

fn blocks(pairs: &BTreeSet<[u8; 4]>) -> Option<(usize, u32, u32)> {
    for s in 1..=2usize {
        let total = 1u32 << (s * s);
        for a in 0..total {
            for b in a + 1..total {
                let free = (0..16u32).all(|arrangement| {
                    let side = 2 * s;
                    let mut cells = vec![0u8; side * side];
                    for row in 0..side {
                        for col in 0..side {
                            let slot = (row / s) * 2 + col / s;
                            let block = if (arrangement >> slot) & 1 == 1 { b } else { a };
                            cells[row * side + col] =
                                ((block >> ((row % s) * s + col % s)) & 1) as u8;
                        }
                    }
                    Pic { side, cells }.pairs().is_subset(pairs)
                });
                if free {
                    return Some((s, a, b));
                }
            }
        }
    }
    None
}

fn branching(states: usize, edges: &[(usize, usize)]) -> bool {
    let mut reach = vec![vec![false; states]; states];
    for &(from, to) in edges {
        reach[from][to] = true;
    }
    for mid in 0..states {
        for from in 0..states {
            for to in 0..states {
                if reach[from][mid] && reach[mid][to] {
                    reach[from][to] = true;
                }
            }
        }
    }
    (0..states).any(|root| {
        let class: Vec<usize> = (0..states)
            .filter(|&v| reach[root][v] && reach[v][root])
            .collect();
        let inner = edges
            .iter()
            .filter(|(f, t)| class.contains(f) && class.contains(t))
            .count();
        !class.is_empty() && inner > class.len()
    })
}

fn lines(pairs: &BTreeSet<[u8; 4]>) -> Option<&'static str> {
    let has = |w: [u8; 4]| pairs.contains(&w);
    let mut two = Vec::new();
    let mut cols = Vec::new();
    for a in 0..2u8 {
        for b in 0..2u8 {
            if has([a, a, b, b]) {
                two.push((a as usize, b as usize));
            }
            if has([a, b, a, b]) {
                cols.push((a as usize, b as usize));
            }
        }
    }
    let mut diag = Vec::new();
    let mut anti = Vec::new();
    for p in 0..2u8 {
        for q in 0..2u8 {
            for s in 0..2u8 {
                let from = (p * 2 + q) as usize;
                let to = (q * 2 + s) as usize;
                if has([q, s, p, q]) {
                    diag.push((from, to));
                }
                if has([p, q, q, s]) {
                    anti.push((from, to));
                }
            }
        }
    }
    if branching(2, &two) {
        Some("rows")
    } else if branching(2, &cols) {
        Some("columns")
    } else if branching(4, &diag) {
        Some("diagonals")
    } else if branching(4, &anti) {
        Some("antidiagonals")
    } else {
        None
    }
}

fn xor(pairs: &BTreeSet<[u8; 4]>) -> bool {
    [(0, 1, 2), (1, 0, 3), (2, 0, 3), (3, 1, 2)]
        .iter()
        .any(|&(corner, left, right)| {
            let rule: BTreeSet<[u8; 4]> = (0..16u8)
                .map(|bits| [bits & 1, (bits >> 1) & 1, (bits >> 2) & 1, (bits >> 3) & 1])
                .filter(|w| w[corner] == w[left] ^ w[right])
                .collect();
            &rule == pairs
        })
}

fn shape(s: usize, block: u32) -> String {
    (0..s)
        .map(|row| {
            (0..s)
                .map(|col| ((block >> (row * s + col)) & 1).to_string())
                .collect::<String>()
        })
        .collect::<Vec<_>>()
        .join("/")
}

fn aligned(design: &Design, m: usize) -> bool {
    let (pairs, _) = design.language();
    let unit = design.base.pow(m as u32);
    let render = design.render(m);
    pairs.iter().all(|window| {
        let mut pic = Pic::block(*window);
        for _ in 0..m {
            pic = design.substitute(&pic);
        }
        (0..=unit).all(|u| {
            (0..=unit).all(|v| {
                (u % unit == 0 && v % unit == 0)
                    || (0..unit).any(|r| (0..unit).any(|c| pic.at(u + r, v + c) != render.at(r, c)))
            })
        })
    })
}

fn recognised(design: &Design) -> Option<usize> {
    (1..=3).find(|&m| aligned(design, m))
}

fn kind(design: &Design) -> String {
    if design.empty() {
        return "sft empty".into();
    }
    if design.full() {
        return "sft full".into();
    }
    if design.degenerate() {
        return "sft boundary".into();
    }
    let pairs: BTreeSet<[u8; 4]> = design.language().0.into_iter().collect();
    if let Some((s, a, b)) = blocks(&pairs) {
        return format!("not-sft blocks {} {}", shape(s, a), shape(s, b));
    }
    if xor(&pairs) {
        return "not-sft xor".into();
    }
    if let Some(direction) = lines(&pairs) {
        return format!("not-sft lines {direction}");
    }
    "open".into()
}

fn orbits(base: usize) -> Vec<(u32, usize)> {
    let mut sizes: BTreeMap<u32, usize> = BTreeMap::new();
    for code in 0..1u32 << (base * base) {
        *sizes
            .entry(Design::new(base, code).canonical())
            .or_insert(0) += 1;
    }
    sizes.into_iter().collect()
}

fn row(values: &[usize]) -> String {
    values
        .iter()
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn periodic(p: &[usize]) -> Option<usize> {
    p.iter()
        .enumerate()
        .map(|(i, &v)| (i + 1, v))
        .find(|&(k, v)| 2 * v <= k * k)
        .map(|(k, _)| k)
}

fn carpet(k: i64) -> i64 {
    let mut n = 1i64;
    while 3 * n + 1 < k {
        n *= 3;
    }
    if k <= 2 * n + 1 {
        6 * k * k + 12 * (n - 1) * k - 10 * n * n - 12 * n + 8
    } else {
        2 * k * k + (24 * n - 4) * k - 18 * n * n - 24 * n + 4
    }
}

fn gasket(k: i64) -> i64 {
    4 * k * k - 6 * k + 4
}

fn formula(base: usize, code: u32, p: &[usize]) -> Option<String> {
    let rule: fn(i64) -> i64 = match (base, code) {
        (3, 495) => carpet,
        (2, 7) => gasket,
        _ => return None,
    };
    let first = if base == 3 { 2 } else { 1 };
    let miss = (first..=p.len()).find(|&k| rule(k as i64) != p[k - 1] as i64);
    Some(match miss {
        Some(k) => format!("formula fails at {k}"),
        None => format!("formula holds for {first} <= k <= {}", p.len()),
    })
}

fn quadratic(p: &[usize]) -> Option<String> {
    let v: Vec<i64> = p.iter().map(|&x| x as i64).collect();
    (1..=3usize).find_map(|first| {
        let i = first - 1;
        if v.len() < i + 4 {
            return None;
        }
        let d2 = v[i + 2] - 2 * v[i + 1] + v[i];
        if (i..v.len() - 2).any(|j| v[j + 2] - 2 * v[j + 1] + v[j] != d2) {
            return None;
        }
        let k = first as i64;
        let b2 = 2 * (v[i + 1] - v[i]) - d2 * (2 * k + 1);
        let c2 = 2 * v[i] - d2 * k * k - b2 * k;
        if d2 % 2 != 0 || b2 % 2 != 0 || c2 % 2 != 0 {
            return Some(format!(
                "quadratic from {first} ({d2} k^2 + {b2} k + {c2})/2"
            ));
        }
        Some(format!(
            "quadratic from {first} {} k^2 + {} k + {}",
            d2 / 2,
            b2 / 2,
            c2 / 2
        ))
    })
}

fn verb_count(base: usize, top: usize, only: Option<u32>) {
    let codes: Vec<(u32, usize)> = match only {
        Some(code) => vec![(code, 1)],
        None => orbits(base),
    };
    let mut distinct: BTreeSet<Vec<usize>> = BTreeSet::new();
    let mut tally: BTreeMap<&str, (usize, usize)> = BTreeMap::new();
    for (code, size) in codes {
        let design = Design::new(base, code);
        let p = count(&design, top);
        if only.is_none() {
            for g in 1..8 {
                let image = design.image(g);
                assert_eq!(
                    count(&image, top.min(12)),
                    p[..top.min(12)],
                    "code {code} image {g}"
                );
            }
        }
        let (pairs, level) = design.language();
        let (tag, shown) = if design.degenerate() {
            (format!("boundary edge {}", row(&p)), vec![1; p.len()])
        } else {
            ("plane".to_string(), p)
        };
        let flat = match periodic(&shown) {
            Some(k) => format!("periodic@{k}"),
            None => "above".into(),
        };
        println!(
            "base {base} code {code} orbit {size} fill {} pairs {} settle {level} {flat} p {} {tag}",
            design.filled().len(),
            pairs.len(),
            row(&shown)
        );
        let class = if shown.iter().all(|&v| v == 1) {
            "constant"
        } else if periodic(&shown).is_some() {
            "periodic"
        } else {
            "above"
        };
        let entry = tally.entry(class).or_insert((0, 0));
        entry.0 += 1;
        entry.1 += size;
        distinct.insert(shown.clone());
        if let Some(line) = quadratic(&shown) {
            println!("base {base} code {code} {line}");
        }
        if let Some(line) = formula(base, code, &shown) {
            println!("base {base} code {code} {line}");
        }
    }
    if only.is_none() {
        summary(base, &distinct, &tally);
    }
}

fn summary(base: usize, distinct: &BTreeSet<Vec<usize>>, tally: &BTreeMap<&str, (usize, usize)>) {
    for (class, (orbits, codes)) in tally {
        println!("base {base} class {class} orbits {orbits} codes {codes}");
    }
    println!("base {base} distinct sequences {}", distinct.len());
}

fn verb_pairs(base: usize, code: u32) {
    let (pairs, level) = Design::new(base, code).language();
    let shown: Vec<String> = pairs
        .iter()
        .map(|w| format!("{}{}/{}{}", w[0], w[1], w[2], w[3]))
        .collect();
    println!(
        "base {base} code {code} settle {level} pairs {} {}",
        pairs.len(),
        shown.join(" ")
    );
}

fn verb_scan(base: usize, level: usize, top: usize, only: Option<u32>) {
    let carpet = crate_render(&Design::new(3, 495), level.min(6));
    let tile = create(Code::from(7u64), 3, 2, 2, level.min(6)).unwrap();
    assert_eq!(
        carpet.cells,
        tile.bytes().unwrap(),
        "code 7 side 3 is not base 3 code 495"
    );
    let mut worst = 0usize;
    let codes: Vec<(u32, usize)> = match only {
        Some(code) => vec![(code, 1)],
        None => orbits(base),
    };
    for (code, _) in codes {
        let design = Design::new(base, code);
        let pic = crate_render(&design, level);
        assert!(
            pic == design.render(level),
            "code {code} renders differently"
        );
        let seen = brute(&pic, top);
        let p = count(&design, seen.len());
        let gaps = seen.iter().zip(&p).filter(|(a, b)| a != b).count();
        worst = worst.max(gaps);
        println!(
            "base {base} code {code} level {level} scan {} gaps {gaps}",
            row(&seen)
        );
    }
    println!("base {base} level {level} worst {worst}");
}

fn verb_kind(base: usize) {
    let mut tally: BTreeMap<String, usize> = BTreeMap::new();
    for (code, size) in orbits(base) {
        let design = Design::new(base, code);
        let verdict = match (
            design.empty() || design.full() || design.degenerate(),
            recognised(&design),
        ) {
            (true, _) => kind(&design),
            (false, Some(m)) => format!("{} recognised {m}", kind(&design)),
            (false, None) => format!("{} unrecognised", kind(&design)),
        };
        let head = verdict.split(' ').take(2).collect::<Vec<_>>().join(" ");
        *tally.entry(head).or_insert(0) += size;
        println!("base {base} code {code} orbit {size} {verdict}");
    }
    for (head, total) in tally {
        println!("base {base} codes {total} {head}");
    }
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let arg =
        |i: usize, default: usize| args.get(i).and_then(|s| s.parse().ok()).unwrap_or(default);
    let clock = Instant::now();
    match args.get(1).map(String::as_str) {
        Some("count") => verb_count(arg(2, 3), arg(3, 28), args.get(4).and_then(|s| s.parse().ok())),
        Some("scan") => verb_scan(arg(2, 3), arg(3, 6), arg(4, 12), args.get(5).and_then(|s| s.parse().ok())),
        Some("kind") => verb_kind(arg(2, 3)),
        Some("pairs") => verb_pairs(arg(2, 3), arg(3, 495) as u32),
        _ => eprintln!("verbs: count <base> <top> [code], scan <base> <level> <top> [code], kind <base>, pairs <base> <code>"),
    }
    eprintln!("seconds {:.1}", clock.elapsed().as_secs_f64());
}
