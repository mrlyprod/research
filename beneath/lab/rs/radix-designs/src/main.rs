use mrly_math::bang::factory::create;
use mrly_math::bang::Code;
use mrly_num::gauss::Ring;
use mrly_math::radix::{flowsnake, gasket, koch, terdragon, tile, twindragon, Base, Radix};
use std::collections::{HashMap, HashSet};
use std::time::Instant;

type Point = (f64, f64);
type Map = (Point, Point);

const BUDGET: usize = 200_000;
const COMPARE: usize = 20_000;

fn ring_word(ring: Ring) -> &'static str {
    match ring {
        Ring::Gaussian => "Z[i]",
        Ring::Eisenstein => "Z[omega]",
    }
}

fn spell(ring: Ring, (a, b): (i64, i64)) -> String {
    let unit = match ring {
        Ring::Gaussian => "i",
        Ring::Eisenstein => "w",
    };
    match (a, b) {
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

fn complex_mul(z: (f64, f64), w: (f64, f64)) -> (f64, f64) {
    (z.0 * w.0 - z.1 * w.1, z.0 * w.1 + z.1 * w.0)
}

fn complex_div(z: (f64, f64), w: (f64, f64)) -> (f64, f64) {
    let den = w.0 * w.0 + w.1 * w.1;
    ((z.0 * w.0 + z.1 * w.1) / den, (z.1 * w.0 - z.0 * w.1) / den)
}

fn iterate(maps: &[Map], level: usize) -> Vec<(f64, f64)> {
    let mut out = vec![(0.0f64, 0.0f64)];
    for _ in 0..level {
        let mut next = Vec::with_capacity(out.len() * maps.len());
        for &(scale, shift) in maps {
            for &p in &out {
                let t = complex_mul(scale, p);
                next.push((t.0 + shift.0, t.1 + shift.1));
            }
        }
        out = next;
    }
    out
}

fn koch_maps() -> Vec<Map> {
    let root = 3f64.sqrt();
    vec![
        ((1.0 / 3.0, 0.0), (0.0, 0.0)),
        ((1.0 / 6.0, root / 6.0), (1.0 / 3.0, 0.0)),
        ((1.0 / 6.0, -root / 6.0), (0.5, root / 6.0)),
        ((1.0 / 3.0, 0.0), (2.0 / 3.0, 0.0)),
    ]
}

fn gasket_maps() -> Vec<Map> {
    let root = 3f64.sqrt();
    let corners = [(1.0, 1.0), (3.0, 1.0), (2.0, 1.0 + root)];
    corners
        .iter()
        .map(|&(x, y)| ((0.5, 0.0), (x / 2.0, y / 2.0)))
        .collect()
}

fn twindragon_maps() -> Vec<Map> {
    let scale = complex_div((1.0, 0.0), (1.0, 1.0));
    vec![(scale, (0.0, 0.0)), (scale, scale)]
}

fn tile_maps(base: Base) -> Vec<Map> {
    let ring = base.ring();
    let b = ring.place(base.value().0, base.value().1);
    let scale = complex_div((1.0, 0.0), b);
    base.residues()
        .unwrap()
        .into_iter()
        .map(|d| {
            let place = ring.place(d.0, d.1);
            (scale, complex_div(place, b))
        })
        .collect()
}

fn pin(got: &[(f64, f64)], want: &[(f64, f64)]) -> (f64, (f64, f64), f64) {
    let last = got.len() - 1;
    let dp = (got[last].0 - got[0].0, got[last].1 - got[0].1);
    let dt = (want[last].0 - want[0].0, want[last].1 - want[0].1);
    let den = dp.0 * dp.0 + dp.1 * dp.1;
    let scale = (dt.0 * dp.0 + dt.1 * dp.1) / den;
    let turn = (dt.1 * dp.0 - dt.0 * dp.1) / den;
    let shift = (want[0].0 - scale * got[0].0, want[0].1 - scale * got[0].1);
    (scale, shift, turn)
}

fn deviation(got: &[(f64, f64)], want: &[(f64, f64)], scale: f64, shift: (f64, f64)) -> f64 {
    got.iter()
        .zip(want.iter())
        .map(|(a, b)| {
            let x = scale * a.0 + shift.0 - b.0;
            let y = scale * a.1 + shift.1 - b.1;
            (x * x + y * y).sqrt()
        })
        .fold(0.0f64, f64::max)
}

fn top_level(size: usize, budget: usize) -> usize {
    let mut top = 1;
    while size.pow(top as u32 + 1) <= budget {
        top += 1;
    }
    top
}

fn face_off(name: &str, design: &Radix, maps: &[Map], free: bool, note: &str) {
    let top = top_level(design.size(), COMPARE);
    let got = design.plane(top);
    let want = iterate(maps, top);
    assert_eq!(got.len(), want.len(), "{name} point counts differ");
    let (scale, shift, turn) = if free {
        pin(&got, &want)
    } else {
        (1.0, (0.0, 0.0), 0.0)
    };
    let worst = deviation(&got, &want, scale, shift);
    let motion = if free {
        format!(
            "scale {scale:.6} shift ({:.6}, {:.6}) turn {turn:.1e}",
            shift.0, shift.1
        )
    } else {
        "none".to_string()
    };
    println!(
        "{name}  level {top}  points {}  worst {worst:.3e}",
        got.len()
    );
    println!("  motion allowed: {motion}");
    println!("  {note}");
}

fn lsystem(level: usize) -> Vec<u8> {
    let mut out = b"F".to_vec();
    for _ in 0..level {
        let mut next = Vec::with_capacity(out.len() * 5);
        for &c in &out {
            if c == b'F' {
                next.extend_from_slice(b"F+F-F");
            } else {
                next.push(c);
            }
        }
        out = next;
    }
    out
}

fn turtle(level: usize) -> Vec<(f64, f64)> {
    let root = 3f64.sqrt();
    let step = [(1.0, 0.0), (-0.5, root / 2.0), (-0.5, -root / 2.0)];
    let mut heading = 0usize;
    let mut at = (0.0f64, 0.0f64);
    let mut starts = Vec::new();
    for c in lsystem(level) {
        match c {
            b'F' => {
                starts.push(at);
                at = (at.0 + step[heading].0, at.1 + step[heading].1);
            }
            b'+' => heading = (heading + 1) % 3,
            _ => heading = (heading + 2) % 3,
        }
    }
    starts.into_iter().map(|p| complex_div(p, at)).collect()
}

fn terdragon_reading(name: &str, design: &Radix, level: usize) -> f64 {
    let got = design.plane(level);
    let want = turtle(level);
    assert_eq!(got.len(), want.len(), "{name} point counts differ");
    let worst = deviation(&got, &want, 1.0, (0.0, 0.0));
    println!(
        "{name}  twists {}  level {level}  points {}  worst {worst:.3e}",
        design
            .twists()
            .iter()
            .map(|&u| spell(design.ring(), u))
            .collect::<Vec<_>>()
            .join(" "),
        got.len()
    );
    worst
}

fn run_compare() {
    println!("COMPARE  each design against an independent f64 iterated function system, seed 0");
    face_off(
        "gasket    ",
        &gasket().unwrap(),
        &gasket_maps(),
        true,
        "independent: the three ratio 1/2 similarities fixing an equilateral triangle placed at (1,1), (3,1), (2,1+sqrt 3); a translation and a positive scaling are allowed because that statement fixes the gasket only up to similarity, and the turn residual printed above is the check that no rotation was needed",
    );
    face_off(
        "koch      ",
        &koch().unwrap(),
        &koch_maps(),
        false,
        "independent in f64 only: z/3, e^(i pi/3) z/3 + 1/3, e^(-i pi/3) z/3 + 1/2 + i sqrt(3)/6, z/3 + 2/3 are the crate maps coefficient for coefficient, so this is a float self-check of exact ring arithmetic and no motion is allowed",
    );
    face_off(
        "twindragon",
        &twindragon().unwrap(),
        &twindragon_maps(),
        false,
        "self-check by definition: z/(1+i) and (z+1)/(1+i) are the radix maps of base 1+i on its two residues, so no motion is allowed and the number is f64 round-off",
    );
    let seven = flowsnake().unwrap();
    face_off(
        "tile7     ",
        &seven,
        &tile_maps(seven.base()),
        false,
        "self-check by definition: the seven maps (z + d)/(3+w) over the canonical residues are the radix maps of that base, so no motion is allowed; the flowsnake name is not tested here",
    );
    println!("  terdragon against the L-system F -> F + F - F at 120 degrees, three segments, unsourced reading");
    let level = 8;
    let twisted = terdragon_reading("  twisted  ", &terdragon().unwrap(), level);
    let plain = terdragon_reading(
        "  untwisted",
        &terdragon().unwrap().with_twists(&[0, 0, 0]).unwrap(),
        level,
    );
    println!(
        "  verdict: the twisted reading {} and the untwisted code 7 {}",
        if twisted < 1e-9 { "matches" } else { "misses" },
        if plain < 1e-9 { "matches" } else { "misses" }
    );
    assert!(twisted < 1e-9, "the twisted reading misses by {twisted}");
    assert!(plain > 1e-3, "the untwisted code 7 matches after all");
}

fn run_koch() {
    let design = koch().unwrap();
    println!(
        "KOCH  base {} on {}",
        spell(design.ring(), design.base().value()),
        ring_word(design.ring())
    );
    println!(
        "  digits {}  twists {}",
        design
            .digits()
            .iter()
            .map(|&d| spell(design.ring(), d))
            .collect::<Vec<_>>()
            .join(" "),
        design
            .twists()
            .iter()
            .map(|&u| spell(design.ring(), u))
            .collect::<Vec<_>>()
            .join(" ")
    );
    println!(
        "  canonical digits {}  class code {}",
        design.canonical().unwrap(),
        design.code().unwrap()
    );
    println!("  the four maps are phi_d coefficient for coefficient, so the column below is a float self-check");
    println!("  level  words  fill  distinct  worst");
    let maps = koch_maps();
    for level in 1..=5 {
        let got = design.plane(level);
        let want = iterate(&maps, level);
        let worst = deviation(&got, &want, 1.0, (0.0, 0.0));
        println!(
            "  {level:>5}  {:>5}  {:>4}  {:>8}  {worst:.3e}",
            got.len(),
            design.fill(level),
            design.distinct(level)
        );
        assert_eq!(got.len() as u128, design.fill(level));
        assert_eq!(design.fill(level), 4u128.pow(level as u32));
        assert!(worst < 1e-9, "level {level} misses by {worst}");
    }
}

fn report(name: &str, design: &Radix) {
    let ring = design.ring();
    let q = design.base().norm();
    let size = design.size();
    let top = top_level(size, BUDGET);
    let counts: Vec<String> = (1..=top)
        .map(|level| design.distinct(level).to_string())
        .collect();
    let pairwise = (1..=top).all(|level| design.distinct(level) as u128 == design.fill(level));
    println!(
        "{name}  {}  base {}  q {q}  code {}  |F| {size}  fill |F|^L  dim {:.6}",
        ring_word(ring),
        spell(ring, design.base().value()),
        design.code().unwrap(),
        design.dimension()
    );
    println!(
        "  digits {}",
        design
            .digits()
            .iter()
            .map(|&d| spell(ring, d))
            .collect::<Vec<_>>()
            .join(" ")
    );
    println!("  distinct levels 1..{top}: {}", counts.join(" "));
    println!("  pairwise distinct {pairwise}");
}

fn run_named() {
    println!("NAMED  codes only; the names are tested by the verb compare and nowhere here");
    report("code 7 at 2        ", &gasket().unwrap());
    report("code 3 at 1+i      ", &twindragon().unwrap());
    report(
        "code 7 at 2+w      ",
        &terdragon().unwrap().with_twists(&[0, 0, 0]).unwrap(),
    );
    report("code 127 at 3+w    ", &flowsnake().unwrap());
    report("code 147 at 3      ", &koch().unwrap());
    let glue = Radix::from_code(Base::new(Ring::Gaussian, (2, 0)).unwrap(), 3)
        .unwrap()
        .with_twists(&[0, 2])
        .unwrap();
    println!(
        "twist glue  Z[i]  base 2  q 4  code 3  |F| 2  twists 1 -1  dim {:.6}",
        glue.dimension()
    );
    println!("  words scaled by b^2: {:?}", glue.words(2));
    let top = 16;
    println!(
        "  fill levels 1..{top}:     {}",
        (1..=top)
            .map(|l| glue.fill(l).to_string())
            .collect::<Vec<_>>()
            .join(" ")
    );
    let counts: Vec<String> = (1..=top)
        .map(|level| glue.distinct(level).to_string())
        .collect();
    println!("  distinct levels 1..{top}: {}", counts.join(" "));
    for level in 1..=top {
        assert_eq!(
            glue.distinct(level),
            (1usize << (level - 1)) + 1,
            "the distinct count leaves 2^(L-1) + 1 at level {level}"
        );
    }
    println!("  distinct equals 2^(L-1) + 1 at every level to {top}");
}

fn run_today() {
    let level = 2;
    let mut checked = 0usize;
    let mut missed = 0usize;
    let mut first: Option<String> = None;
    for m in [2u64, 3u64] {
        let cells = (m * m) as usize;
        for code in 0..(1u128 << cells) {
            let design = tile(m, code).unwrap();
            let got: HashSet<(i64, i64)> = design.words(level).into_iter().collect();
            let tensor = create(Code::from(code), m as usize, 2, m as usize, level).unwrap();
            let side = m.pow(level as u32) as usize;
            let mut want = HashSet::new();
            for row in 0..side {
                for col in 0..side {
                    if tensor.at(row * side + col) == 1 {
                        want.insert((col as i64, row as i64));
                    }
                }
            }
            checked += 1;
            assert_eq!(design.words(level).len() as u128, design.fill(level));
            if got != want {
                missed += 1;
                if first.is_none() {
                    first = Some(format!("m {m} code {code}"));
                }
            }
        }
    }
    println!("TODAY  level {level}, bases 2 and 3, every plane code");
    println!("  codes checked {checked}  mismatches {missed}");
    println!(
        "  the code is read in box row-major order, bit r m + c, not in canonical residue order"
    );
    for m in [2u64, 3u64] {
        let full = tile(m, (1u128 << (m * m)) - 1).unwrap();
        let base = Base::new(Ring::Gaussian, (m as i64, 0)).unwrap();
        println!(
            "  m {m}: box residues are the canonical system {}",
            full.canonical().unwrap()
        );
        println!(
            "    box {}",
            full.digits()
                .iter()
                .map(|&d| spell(Ring::Gaussian, d))
                .collect::<Vec<_>>()
                .join(" ")
        );
        println!(
            "    canonical {}",
            base.residues()
                .unwrap()
                .iter()
                .map(|&z| spell(Ring::Gaussian, z))
                .collect::<Vec<_>>()
                .join(" ")
        );
    }
    match first {
        None => println!("  first mismatch none"),
        Some(where_) => println!("  first mismatch {where_}"),
    }
}

fn cycles(map: &[usize]) -> u32 {
    let mut seen = vec![false; map.len()];
    let mut count = 0;
    for start in 0..map.len() {
        if seen[start] {
            continue;
        }
        count += 1;
        let mut at = start;
        while !seen[at] {
            seen[at] = true;
            at = map[at];
        }
    }
    count
}

fn burnside(group: &[Vec<usize>]) -> u128 {
    let total: u128 = group.iter().map(|map| 1u128 << cycles(map)).sum();
    assert_eq!(total % group.len() as u128, 0, "Burnside is not an integer");
    total / group.len() as u128
}

fn image(group: &[Vec<usize>]) -> usize {
    group.iter().collect::<HashSet<_>>().len()
}

fn act(map: &[usize], code: usize) -> usize {
    let mut out = 0usize;
    for (i, &j) in map.iter().enumerate() {
        if (code >> i) & 1 == 1 {
            out |= 1 << j;
        }
    }
    out
}

fn orbit(group: &[Vec<usize>], code: usize) -> Vec<usize> {
    let mut seen = HashSet::new();
    let mut stack = vec![code];
    seen.insert(code);
    while let Some(at) = stack.pop() {
        for map in group {
            let next = act(map, at);
            if seen.insert(next) {
                stack.push(next);
            }
        }
    }
    let mut out: Vec<usize> = seen.into_iter().collect();
    out.sort_unstable();
    out
}

fn walk(group: &[Vec<usize>], q: usize) -> u128 {
    let mut seen = vec![false; 1usize << q];
    let mut orbits = 0u128;
    for code in 0..(1usize << q) {
        if seen[code] {
            continue;
        }
        orbits += 1;
        for member in orbit(group, code) {
            seen[member] = true;
        }
    }
    orbits
}

fn run_census() {
    let bases = [
        (Ring::Gaussian, (2i64, 0i64)),
        (Ring::Gaussian, (1, 1)),
        (Ring::Gaussian, (2, 1)),
        (Ring::Eisenstein, (2, 0)),
        (Ring::Eisenstein, (2, 1)),
        (Ring::Eisenstein, (3, 0)),
        (Ring::Eisenstein, (3, 1)),
    ];
    println!(
        "CENSUS  classes of digit CODES under the residue action, never designs up to similarity"
    );
    println!("  ring  base  q  abstract  image  mirror  codes  classes  walk");
    for (ring, value) in bases {
        let base = Base::new(ring, value).unwrap();
        let q = base.norm() as usize;
        let group = base.group().unwrap();
        let count = burnside(&group);
        let seen = walk(&group, q);
        assert_eq!(count, seen, "Burnside and the orbit walk disagree");
        println!(
            "  {:>9}  {:>4}  {q}  {:>8}  {:>5}  {:>6}  {:>5}  {count:>7}  {seen:>4}",
            ring_word(ring),
            spell(ring, value),
            group.len(),
            image(&group),
            base.mirrored(),
            1u128 << q
        );
        println!(
            "    residues {}",
            base.residues()
                .unwrap()
                .iter()
                .map(|&z| spell(ring, z))
                .collect::<Vec<_>>()
                .join(" ")
        );
    }
    println!("  abstract is the order of R^* semidirect <conj>, image the order it acts through on the residues");
    let base = Base::new(Ring::Eisenstein, (3, 0)).unwrap();
    let q = base.norm() as usize;
    let twisted: u128 = (0..1u128 << q)
        .map(|code| 6u128.pow(code.count_ones()))
        .sum();
    println!(
        "  base 3 on Z[omega]: {} codes in {} classes of codes; over the codes, not the classes, the twist vectors number sum_k binom(9,k) 6^k = {twisted} = 7^{q}",
        1u128 << q,
        burnside(&base.group().unwrap())
    );
    assert_eq!(twisted, 7u128.pow(q as u32));
    let mut sizes: HashMap<usize, usize> = HashMap::new();
    for code in 0..1u128 << q {
        *sizes.entry(code.count_ones() as usize).or_default() += 1;
    }
    let mut keys: Vec<usize> = sizes.keys().copied().collect();
    keys.sort_unstable();
    println!(
        "  codes by |F|: {}",
        keys.iter()
            .map(|k| format!("{k}:{}", sizes[k]))
            .collect::<Vec<_>>()
            .join(" ")
    );
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Alg {
    p: i128,
    q: i128,
    r: i128,
    m: i128,
}

fn gcd(a: i128, b: i128) -> i128 {
    if b == 0 {
        a.abs().max(1)
    } else {
        gcd(b, a % b)
    }
}

fn twist_of(ring: Ring) -> i128 {
    match ring {
        Ring::Gaussian => 0,
        Ring::Eisenstein => 1,
    }
}

impl Alg {
    fn new(m: i128, p: i128, q: i128, r: i128) -> Alg {
        assert!(r != 0, "a rational needs a nonzero denominator");
        let (p, q, r) = if r < 0 { (-p, -q, -r) } else { (p, q, r) };
        if p == 0 && q == 0 {
            return Alg {
                p: 0,
                q: 0,
                r: 1,
                m,
            };
        }
        let g = gcd(gcd(p, q), r);
        Alg {
            p: p / g,
            q: q / g,
            r: r / g,
            m,
        }
    }
    fn whole(m: i128, p: i64, q: i64) -> Alg {
        Alg::new(m, p as i128, q as i128, 1)
    }
    fn sub(self, other: Alg) -> Alg {
        Alg::new(
            self.m,
            self.p * other.r - other.p * self.r,
            self.q * other.r - other.q * self.r,
            self.r * other.r,
        )
    }
    fn add(self, other: Alg) -> Alg {
        Alg::new(
            self.m,
            self.p * other.r + other.p * self.r,
            self.q * other.r + other.q * self.r,
            self.r * other.r,
        )
    }
    fn mul(self, other: Alg) -> Alg {
        Alg::new(
            self.m,
            self.p * other.p - self.q * other.q,
            self.p * other.q + self.q * other.p - self.m * self.q * other.q,
            self.r * other.r,
        )
    }
    fn conjugate(self) -> Alg {
        Alg::new(self.m, self.p - self.m * self.q, -self.q, self.r)
    }
    fn div(self, other: Alg) -> Alg {
        let n = other.p * other.p - self.m * other.p * other.q + other.q * other.q;
        assert!(n != 0, "no division by zero");
        self.mul(Alg::new(
            self.m,
            other.r * (other.p - self.m * other.q),
            -other.r * other.q,
            n,
        ))
    }
}

fn direct(left: &[Alg], right: &[Alg]) -> bool {
    if left.len() != right.len() {
        return false;
    }
    if left.len() < 2 {
        return true;
    }
    let span = left[1].sub(left[0]);
    for first in 0..right.len() {
        for second in 0..right.len() {
            if first == second {
                continue;
            }
            let v = right[second].sub(right[first]).div(span);
            let t = right[first].sub(v.mul(left[0]));
            if left.iter().all(|&d| right.contains(&v.mul(d).add(t))) {
                return true;
            }
        }
    }
    false
}

fn conjugate_sets(left: &[Alg], right: &[Alg], mirror: bool) -> bool {
    if direct(left, right) {
        return true;
    }
    if !mirror {
        return false;
    }
    let flipped: Vec<Alg> = left.iter().map(|d| d.conjugate()).collect();
    direct(&flipped, right)
}

fn digits_of(ring: Ring, residues: &[(i64, i64)], code: usize) -> Vec<Alg> {
    let m = twist_of(ring);
    residues
        .iter()
        .enumerate()
        .filter(|(i, _)| (code >> i) & 1 == 1)
        .map(|(_, &z)| Alg::whole(m, z.0, z.1))
        .collect()
}

fn similar_classes(
    ring: Ring,
    residues: &[(i64, i64)],
    codes: &[usize],
    mirror: bool,
) -> Vec<usize> {
    let sets: Vec<Vec<Alg>> = codes
        .iter()
        .map(|&c| digits_of(ring, residues, c))
        .collect();
    let mut label = vec![usize::MAX; codes.len()];
    let mut classes = 0;
    for i in 0..codes.len() {
        if label[i] != usize::MAX {
            continue;
        }
        label[i] = classes;
        for j in i + 1..codes.len() {
            if label[j] == usize::MAX && conjugate_sets(&sets[i], &sets[j], mirror) {
                label[j] = classes;
            }
        }
        classes += 1;
    }
    label
}

fn orbit_labels(group: &[Vec<usize>], codes: &[usize]) -> (Vec<usize>, usize) {
    let mut seat: HashMap<usize, usize> = HashMap::new();
    let mut orbits = 0usize;
    for &code in codes {
        if seat.contains_key(&code) {
            continue;
        }
        for member in orbit(group, code) {
            seat.insert(member, orbits);
        }
        orbits += 1;
    }
    (codes.iter().map(|c| seat[c]).collect(), orbits)
}

type Witness = Option<(usize, usize)>;

fn crossing(orbits: &[usize], classes: &[usize], codes: &[usize]) -> (Witness, Witness) {
    let mut split = None;
    let mut merge = None;
    for i in 0..codes.len() {
        for j in i + 1..codes.len() {
            let same_orbit = orbits[i] == orbits[j];
            let same_class = classes[i] == classes[j];
            if same_orbit && !same_class && split.is_none() {
                split = Some((codes[i], codes[j]));
            }
            if !same_orbit && same_class && merge.is_none() {
                merge = Some((codes[i], codes[j]));
            }
        }
    }
    (split, merge)
}

fn pair_word(pair: Option<(usize, usize)>) -> String {
    match pair {
        Some((a, b)) => format!("{a}/{b}"),
        None => "-".to_string(),
    }
}

fn spell_code(ring: Ring, residues: &[(i64, i64)], code: usize) -> String {
    format!(
        "{code} ({})",
        residues
            .iter()
            .enumerate()
            .filter(|(i, _)| (code >> i) & 1 == 1)
            .map(|(_, &z)| spell(ring, z))
            .collect::<Vec<_>>()
            .join(", ")
    )
}

fn grid(z: (f64, f64)) -> (i64, i64) {
    ((z.0 * 1e9).round() as i64, (z.1 * 1e9).round() as i64)
}

fn shape(points: &[(f64, f64)]) -> Vec<(i64, i64)> {
    if points.len() < 2 {
        return Vec::new();
    }
    let mut best: Option<Vec<(i64, i64)>> = None;
    for first in 0..points.len() {
        for second in 0..points.len() {
            if first == second {
                continue;
            }
            let span = (
                points[second].0 - points[first].0,
                points[second].1 - points[first].1,
            );
            let mut key: Vec<(i64, i64)> = points
                .iter()
                .map(|p| {
                    grid(complex_div(
                        (p.0 - points[first].0, p.1 - points[first].1),
                        span,
                    ))
                })
                .collect();
            key.sort_unstable();
            if best.as_ref().is_none_or(|held| key < *held) {
                best = Some(key);
            }
        }
    }
    best.expect("two digits fix a shape")
}

fn float_similar(ring: Ring, residues: &[(i64, i64)], codes: &[usize], mirror: bool) -> usize {
    let mut seen = HashSet::new();
    for &code in codes {
        let points: Vec<(f64, f64)> = residues
            .iter()
            .enumerate()
            .filter(|(i, _)| (code >> i) & 1 == 1)
            .map(|(_, &z)| ring.place(z.0, z.1))
            .collect();
        let mut key = shape(&points);
        if mirror {
            let flipped: Vec<(f64, f64)> = points.iter().map(|&(x, y)| (x, -y)).collect();
            key = key.min(shape(&flipped));
        }
        seen.insert(key);
    }
    seen.len()
}

#[derive(Clone, Copy, Debug, Hash, PartialEq, Eq, PartialOrd, Ord)]
struct Frac {
    n: i128,
    d: i128,
}

impl Frac {
    fn new(n: i128, d: i128) -> Frac {
        assert!(d != 0, "a fraction needs a nonzero denominator");
        let (n, d) = if d < 0 { (-n, -d) } else { (n, d) };
        if n == 0 {
            return Frac { n: 0, d: 1 };
        }
        let g = gcd(n, d);
        Frac { n: n / g, d: d / g }
    }
}

fn det(u: (i64, i64), v: (i64, i64)) -> i128 {
    u.0 as i128 * v.1 as i128 - u.1 as i128 * v.0 as i128
}

fn lattice_of(residues: &[(i64, i64)], code: usize) -> Vec<(i64, i64)> {
    residues
        .iter()
        .enumerate()
        .filter(|(i, _)| (code >> i) & 1 == 1)
        .map(|(_, &z)| z)
        .collect()
}

fn affine_key(points: &[(i64, i64)]) -> Vec<(Frac, Frac)> {
    let n = points.len();
    if n < 2 {
        return Vec::new();
    }
    let sub = |a: (i64, i64), b: (i64, i64)| (a.0 - b.0, a.1 - b.1);
    let mut best: Option<Vec<(Frac, Frac)>> = None;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                if i == j || i == k || j == k {
                    continue;
                }
                let e1 = sub(points[j], points[i]);
                let e2 = sub(points[k], points[i]);
                let base = det(e1, e2);
                if base == 0 {
                    continue;
                }
                let mut key: Vec<(Frac, Frac)> = points
                    .iter()
                    .map(|&p| {
                        let u = sub(p, points[i]);
                        (Frac::new(det(u, e2), base), Frac::new(det(e1, u), base))
                    })
                    .collect();
                key.sort_unstable();
                if best.as_ref().is_none_or(|held| key < *held) {
                    best = Some(key);
                }
            }
        }
    }
    if let Some(key) = best {
        return key;
    }
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let v = sub(points[j], points[i]);
            let mut key: Vec<(Frac, Frac)> = points
                .iter()
                .map(|&p| {
                    let u = sub(p, points[i]);
                    let along = if v.0 != 0 {
                        Frac::new(u.0 as i128, v.0 as i128)
                    } else {
                        Frac::new(u.1 as i128, v.1 as i128)
                    };
                    (along, Frac::new(0, 1))
                })
                .collect();
            key.sort_unstable();
            if best.as_ref().is_none_or(|held| key < *held) {
                best = Some(key);
            }
        }
    }
    best.expect("two distinct digits fix a line")
}

fn plane_classes(residues: &[(i64, i64)], codes: &[usize]) -> Vec<usize> {
    let mut seat: HashMap<Vec<(Frac, Frac)>, usize> = HashMap::new();
    let mut out = Vec::with_capacity(codes.len());
    for &code in codes {
        let key = affine_key(&lattice_of(residues, code));
        let next = seat.len();
        out.push(*seat.entry(key).or_insert(next));
    }
    out
}

fn carry(matrix: [[i64; 2]; 2], points: &[(i64, i64)]) -> Vec<(i64, i64)> {
    points
        .iter()
        .map(|&(a, b)| {
            (
                matrix[0][0] * a + matrix[0][1] * b,
                matrix[1][0] * a + matrix[1][1] * b,
            )
        })
        .collect()
}

fn same_set(left: &[(i64, i64)], right: &[(i64, i64)]) -> bool {
    left.len() == right.len() && left.iter().all(|z| right.contains(z))
}

fn run_affine() {
    let mut below = 0usize;
    let mut equal = 0usize;
    let mut above = 0usize;
    let mut identical = 0usize;
    let mut crossed = 0usize;
    let mut affine_below = 0usize;
    let mut affine_equal = 0usize;
    let mut affine_above = 0usize;
    let mut affine_crossed = 0usize;
    let bases = [
        (Ring::Gaussian, (2i64, 0i64)),
        (Ring::Gaussian, (1, 1)),
        (Ring::Gaussian, (2, 1)),
        (Ring::Eisenstein, (2, 0)),
        (Ring::Eisenstein, (2, 1)),
        (Ring::Eisenstein, (3, 0)),
        (Ring::Eisenstein, (3, 1)),
    ];
    println!("AFFINE  the untwisted canonical digit sets of every base up to SIMILARITY and up to AFFINE conjugacy, beside the code census, in exact arithmetic");
    println!("  conjugating the place maps by an invertible real affine h(x) = H x + s gives (y + H d + s (b - 1))/b, again an untwisted base-b place map exactly when H commutes with multiplication by 1/b, and s (b - 1) sweeps the plane since N(b) >= 2 forces b != 1");
    println!("  at a non-real base the centraliser of 1/b in the two by two real matrices is C, so the conjugacy group is the similarity group x -> v x + t; at a real base 1/b is the scalar (1/b) I, it commutes with every H, and the conjugacy group is the whole real affine group GL_2 semidirect R^2");
    println!("  the mirror x -> v conj(x) + t preserves the untwisted base-b family exactly when conj(b) = b, since a direct conjugacy keeps the derivative 1/b and a mirror one sends it to 1/conj(b); at a real base it is one element of the full affine group and not the only new one");
    println!("  simil is the class count under the similarity group, affine the class count under the conjugacy group, equal to simil at the four non-real bases by the centraliser lemma and computed over GL_2(Q) semidirect Q^2 at the three real bases");
    println!("  orbits counts digit CODES under the unit group joined by conjugation where conj(b) is an associate of b; ssplit and smerge witness the crossing of orbits with simil, asplit and amerge the crossing of orbits with affine");
    println!(
        "  ring  base  q  |F|  codes  orbits  simil  affine   ssplit   smerge   asplit   amerge"
    );
    for (ring, value) in bases {
        let base = Base::new(ring, value).unwrap();
        let residues = base.residues().unwrap();
        let q = residues.len();
        let group = base.group().unwrap();
        let real = value.1 == 0;
        let mut codes_seen = 0u128;
        let mut orbits_seen = 0usize;
        let mut simil_seen = 0usize;
        let mut affine_seen = 0usize;
        for size in 0..=q {
            let codes: Vec<usize> = (0..1usize << q)
                .filter(|c| c.count_ones() as usize == size)
                .collect();
            let (orbits, orbit_count) = orbit_labels(&group, &codes);
            let simil = similar_classes(ring, &residues, &codes, real);
            let simil_count = simil.iter().copied().max().unwrap() + 1;
            assert_eq!(
                float_similar(ring, &residues, &codes, real),
                simil_count,
                "the float rerun of the normal form disagrees with the exact similarity classes"
            );
            let affine = if real {
                plane_classes(&residues, &codes)
            } else {
                simil.clone()
            };
            let affine_count = affine.iter().copied().max().unwrap() + 1;
            for i in 0..codes.len() {
                for j in i + 1..codes.len() {
                    assert!(
                        simil[i] != simil[j] || affine[i] == affine[j],
                        "the similarity classes do not refine the affine classes"
                    );
                }
            }
            let (split, merge) = crossing(&orbits, &simil, &codes);
            let (asplit, amerge) = crossing(&orbits, &affine, &codes);
            println!(
                "  {:>9}  {:>4}  {q}  {size:>3}  {:>5}  {orbit_count:>6}  {simil_count:>5}  {affine_count:>6}  {:>7}  {:>7}  {:>7}  {:>7}",
                ring_word(ring),
                spell(ring, value),
                codes.len(),
                pair_word(split),
                pair_word(merge),
                pair_word(asplit),
                pair_word(amerge)
            );
            match simil_count.cmp(&orbit_count) {
                std::cmp::Ordering::Less => below += 1,
                std::cmp::Ordering::Equal => equal += 1,
                std::cmp::Ordering::Greater => above += 1,
            }
            match affine_count.cmp(&orbit_count) {
                std::cmp::Ordering::Less => affine_below += 1,
                std::cmp::Ordering::Equal => affine_equal += 1,
                std::cmp::Ordering::Greater => affine_above += 1,
            }
            if split.is_none() && merge.is_none() {
                identical += 1;
            }
            if split.is_some() && merge.is_some() {
                crossed += 1;
            }
            if asplit.is_some() && amerge.is_some() {
                affine_crossed += 1;
            }
            codes_seen += codes.len() as u128;
            orbits_seen += orbit_count;
            simil_seen += simil_count;
            affine_seen += affine_count;
        }
        println!(
            "    totals: codes {codes_seen} = 2^{q}, code classes {orbits_seen}, similarity classes {simil_seen}, affine classes {affine_seen}, real base {real}, mirror in the code group {}",
            base.mirrored()
        );
        assert_eq!(codes_seen, 1u128 << q, "the sizes do not exhaust the codes");
        assert_eq!(
            orbits_seen as u128,
            burnside(&group),
            "the per-size orbit counts do not sum to the census, which controls the code column alone"
        );
    }
    println!(
        "  over the {} cells the similarity count is below the code count in {below}, equal in {equal} and above it in {above}, the two partitions identical in {identical} and crossing in {crossed}",
        below + equal + above
    );
    println!(
        "  over the same cells the affine count is below the code count in {affine_below}, equal in {affine_equal} and above it in {affine_above}, and the two partitions cross in {affine_crossed}"
    );
    let base = Base::new(Ring::Eisenstein, (3, 0)).unwrap();
    let residues = base.residues().unwrap();
    let group = base.group().unwrap();
    println!(
        "  the fixed witnesses at base 3 on Z[omega], |F| = 3, untwisted canonical digit sets"
    );
    for (a, b, note) in [
        (
            131usize,
            137usize,
            "share a census orbit and are not similar",
        ),
        (7, 42, "are similar and sit in different census orbits"),
    ] {
        let codes = vec![a, b];
        let (orbits, _) = orbit_labels(&group, &codes);
        let simil = similar_classes(Ring::Eisenstein, &residues, &codes, true);
        println!(
            "    codes {} and {} {note}: same orbit {}, same similarity class {}",
            spell_code(Ring::Eisenstein, &residues, a),
            spell_code(Ring::Eisenstein, &residues, b),
            orbits[0] == orbits[1],
            simil[0] == simil[1]
        );
        assert_eq!(orbits[0] == orbits[1], note.starts_with("share"));
        assert_eq!(simil[0] == simil[1], !note.starts_with("share"));
    }
    for (a, b, matrix) in [
        (131usize, 137usize, [[0i64, 2], [1, -1]]),
        (7, 131, [[1, 1], [0, 1]]),
    ] {
        let left = lattice_of(&residues, a);
        let right = lattice_of(&residues, b);
        let moved = carry(matrix, &left);
        let turn = det((matrix[0][0], matrix[1][0]), (matrix[0][1], matrix[1][1]));
        println!(
            "    codes {} and {} are affinely conjugate by H = [[{}, {}], [{}, {}]] of determinant {turn}, same affine key {}",
            spell_code(Ring::Eisenstein, &residues, a),
            spell_code(Ring::Eisenstein, &residues, b),
            matrix[0][0],
            matrix[0][1],
            matrix[1][0],
            matrix[1][1],
            affine_key(&left) == affine_key(&right)
        );
        assert!(turn != 0, "a conjugacy needs an invertible H");
        assert!(same_set(&moved, &right), "H does not carry the digit set");
        assert_eq!(
            affine_key(&left),
            affine_key(&right),
            "the affine normal form misses a witnessed conjugacy"
        );
    }
    let three: Vec<usize> = (0..512usize).filter(|c| c.count_ones() == 3).collect();
    let (_, orbit_count) = orbit_labels(&group, &three);
    let simil = similar_classes(Ring::Eisenstein, &residues, &three, true);
    let affine = plane_classes(&residues, &three);
    assert_eq!(three.len(), 84);
    assert_eq!(orbit_count, 13);
    assert_eq!(simil.iter().copied().max().unwrap() + 1, 9);
    assert_eq!(affine.iter().copied().max().unwrap() + 1, 2);
    let bare = similar_classes(Ring::Eisenstein, &residues, &three, false);
    println!(
        "    the mirror is load-bearing for similarity: the 84 three-digit codes fall in {} direct classes and {} once the mirror joins, and in {} affine classes, the collinear triples against the rest",
        bare.iter().copied().max().unwrap() + 1,
        simil.iter().copied().max().unwrap() + 1,
        affine.iter().copied().max().unwrap() + 1
    );
}

fn main() {
    let verbs: Vec<String> = std::env::args().skip(1).collect();
    let verbs = if verbs.is_empty() {
        ["koch", "compare", "named", "today", "census", "affine"]
            .iter()
            .map(|s| s.to_string())
            .collect()
    } else {
        verbs
    };
    for verb in verbs {
        let clock = Instant::now();
        match verb.as_str() {
            "koch" => run_koch(),
            "compare" => run_compare(),
            "named" => run_named(),
            "today" => run_today(),
            "census" => run_census(),
            "affine" => run_affine(),
            other => panic!("unknown verb {other}"),
        }
        println!("  {verb} {:.2} s", clock.elapsed().as_secs_f64());
    }
}
