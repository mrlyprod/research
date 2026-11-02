use crate::cocycle::{morse, value, COLUMN, DIAGONAL, FULL, GASKET, ROW, UNIT};
use crate::series::{Frac, Rep};
use crate::word::{render, CODES};
use mrly_core::rng::Rng;

pub const DEEP: usize = 13;
pub const SEEN: usize = 7;
pub const RANGE: usize = 1 << 14;
pub const FAR: usize = 1 << 15;
pub const EXACT: usize = 60;
const SEED: u64 = 1618033988;

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Shape {
    ZeroContact,
    GasketDomino,
}

pub struct Open {
    pub name: &'static str,
    pub rule: &'static str,
    pub shape: Shape,
    pub pairs: Vec<(u8, u8)>,
}

fn fill_of(word: &[u8]) -> i128 {
    word.iter().map(|code| code.count_ones() as i128).product()
}

pub fn count(shape: Shape, word: &[u8], heavy: u8, light: u8) -> i128 {
    let last = word.iter().rposition(|code| *code == light);
    match (shape, last) {
        (_, None) => 1,
        (Shape::ZeroContact, Some(at)) => fill_of(&word[..=at]),
        (Shape::GasketDomino, Some(at)) => {
            let mut total = 1i128;
            for (i, code) in word[..=at].iter().enumerate() {
                if *code == heavy {
                    total += fill_of(&word[..i]);
                }
            }
            total
        }
    }
}

fn pair_word(heavy: u8, light: u8, mask: usize, length: usize) -> Vec<u8> {
    (0..length)
        .map(|i| if (mask >> i) & 1 == 1 { light } else { heavy })
        .collect()
}

fn cross(left: &[u8], right: &[u8]) -> Vec<(u8, u8)> {
    let mut out: Vec<(u8, u8)> = Vec::new();
    for &a in left.iter() {
        for &b in right.iter() {
            out.push((a, b));
        }
    }
    out
}

pub fn open_families() -> Vec<Open> {
    let dominoes: Vec<u8> = ROW.iter().chain(COLUMN.iter()).copied().collect();
    vec![
        Open {
            name: "gasket and unit",
            rule: "3^(g - j), g the gasket count, j the terminal gasket run",
            shape: Shape::ZeroContact,
            pairs: cross(&GASKET, &UNIT),
        },
        Open {
            name: "gasket and diagonal",
            rule: "2^d 3^(g - j), d the diagonal count",
            shape: Shape::ZeroContact,
            pairs: cross(&GASKET, &DIAGONAL),
        },
        Open {
            name: "full and unit",
            rule: "4^(F - j), F the full count, j the terminal full run",
            shape: Shape::ZeroContact,
            pairs: cross(&[FULL], &UNIT),
        },
        Open {
            name: "full and diagonal",
            rule: "2^d 4^(F - j)",
            shape: Shape::ZeroContact,
            pairs: cross(&[FULL], &DIAGONAL),
        },
        Open {
            name: "gasket and domino",
            rule:
                "1 + sum over the gasket places i <= m of fill(w_1..i-1), m the last domino place",
            shape: Shape::GasketDomino,
            pairs: cross(&GASKET, &dominoes),
        },
    ]
}

fn forms(rep: &Rep) {
    println!("CLOSED FORMS ON THE 46");
    let mut pairs = 0usize;
    let mut words = 0usize;
    let mut drawn = 0usize;
    for family in open_families().iter() {
        let mut checked = 0usize;
        let mut bad = 0usize;
        let mut seen = 0usize;
        let mut wrong = 0usize;
        for (heavy, light) in family.pairs.iter() {
            for length in 1..=DEEP {
                for mask in 0..(1usize << length) {
                    let word = pair_word(*heavy, *light, mask, length);
                    let want = count(family.shape, &word, *heavy, *light);
                    if value(rep, &word) != want {
                        bad += 1;
                    }
                    checked += 1;
                    if length <= SEEN {
                        if render(&word).components() as i128 != want {
                            wrong += 1;
                        }
                        seen += 1;
                    }
                }
            }
        }
        println!(
            "{}: {} pairs, comp = {}, {checked} words to L = {DEEP}, mismatches {bad}, {seen} words drawn to L = {SEEN}, mismatches {wrong}",
            family.name,
            family.pairs.len(),
            family.rule
        );
        assert_eq!(
            bad, 0,
            "the closed form is exact against the representation"
        );
        assert_eq!(wrong, 0, "the closed form is exact against the drawn cells");
        pairs += family.pairs.len();
        words += checked;
        drawn += seen;
    }
    println!("{pairs} pairs, {words} words against the representation and {drawn} against the drawn cells, mismatches 0; with the 59 above every one of the 105 letter pairs now carries an exact closed form");
    assert_eq!(pairs, 46, "the five families cover the 46 open pairs");
}

fn class_of(code: u8) -> usize {
    if UNIT.contains(&code) {
        0
    } else if ROW.contains(&code) {
        1
    } else if COLUMN.contains(&code) {
        2
    } else if DIAGONAL.contains(&code) {
        3
    } else if GASKET.contains(&code) {
        4
    } else {
        5
    }
}

fn class_name(which: usize) -> &'static str {
    [
        "unit",
        "row domino",
        "column domino",
        "diagonal",
        "gasket",
        "full",
    ][which]
}

fn weight(code: u8) -> (i64, i64) {
    match code.count_ones() {
        1 => (0, 0),
        2 => (1, 0),
        3 => (0, 1),
        _ => (2, 0),
    }
}

fn phi_weight(code: u8) -> (i64, i64) {
    if DIAGONAL.contains(&code) {
        (1, 0)
    } else {
        (0, 0)
    }
}

pub fn exponent_weights(a: u8, b: u8) -> ((i64, i64), (i64, i64)) {
    let (ca, cb) = (class_of(a), class_of(b));
    let light = (0i64, 0i64);
    let two = (1i64, 0i64);
    let three = (0i64, 1i64);
    let four = (2i64, 0i64);
    let (low, high) = if ca <= cb { (ca, cb) } else { (cb, ca) };
    let (first, second) = match (low, high) {
        (0, 0) => (light, light),
        (0, 1) | (0, 2) => (light, two),
        (0, 3) => (light, two),
        (0, 4) => (light, three),
        (0, 5) => (light, four),
        (1, 1) | (2, 2) => (light, light),
        (1, 2) => (two, two),
        (1, 3) | (2, 3) => (two, two),
        (1, 4) | (2, 4) => (two, three),
        (1, 5) | (2, 5) => (light, two),
        (3, 3) => (two, two),
        (3, 4) => (two, three),
        (3, 5) => (two, four),
        (4, 4) => (light, light),
        _ => (light, light),
    };
    if ca <= cb {
        (first, second)
    } else {
        (second, first)
    }
}

fn nats(pair: (i64, i64)) -> f64 {
    pair.0 as f64 * 2f64.ln() + pair.1 as f64 * 3f64.ln()
}

fn ledger() {
    println!();
    println!("THE FREQUENCY LEDGER");
    let mut saturating = 0usize;
    let mut short: Vec<(usize, usize)> = Vec::new();
    let mut exact = 0usize;
    let mut refuted: Vec<(usize, usize)> = Vec::new();
    let mut between: Vec<(usize, usize)> = Vec::new();
    for (i, &a) in CODES.iter().enumerate() {
        for &b in &CODES[i + 1..] {
            let (wa, wb) = exponent_weights(a, b);
            let key = (class_of(a).min(class_of(b)), class_of(a).max(class_of(b)));
            if wa == weight(a) && wb == weight(b) {
                saturating += 1;
            } else {
                short.push(key);
            }
            if wa == phi_weight(a) && wb == phi_weight(b) {
                exact += 1;
            } else {
                refuted.push(key);
            }
            let middle = (nats(wa) + nats(wb)) / 2.0;
            let prediction = (nats(phi_weight(a)) + nats(phi_weight(b))) / 2.0;
            let ceiling = (nats(weight(a)) + nats(weight(b))) / 2.0;
            if middle > prediction + 1e-12 && middle < ceiling - 1e-12 {
                between.push(key);
            }
        }
    }
    let tally = |list: &[(usize, usize)]| {
        let mut seen: Vec<((usize, usize), usize)> = Vec::new();
        for key in list.iter() {
            match seen.iter_mut().find(|(at, _)| at == key) {
                Some((_, n)) => *n += 1,
                None => seen.push((*key, 1)),
            }
        }
        seen.iter()
            .map(|((low, high), n)| {
                if low == high {
                    format!("{} {n}", class_name(*low))
                } else {
                    format!("{} and {} {n}", class_name(*low), class_name(*high))
                }
            })
            .collect::<Vec<_>>()
            .join(", ")
    };
    println!(
        "at interior frequency the exponent saturates the fill ceiling on {saturating} of the 105 pairs and falls short on {}: {}",
        short.len(),
        tally(&short)
    );
    println!(
        "Phi(f) = (f_6 + f_9) log 2 is exact on {exact} pairs and refuted on {}: {}",
        refuted.len(),
        tally(&refuted)
    );
    println!(
        "at equal frequencies exactly {} pairs sit strictly between the prediction and the fill ceiling: {}",
        between.len(),
        tally(&between)
    );
    assert_eq!(saturating + short.len(), 105, "every pair is counted once");
    assert_eq!(exact + refuted.len(), 105, "every pair is counted once");
    assert_eq!(between.len(), 4, "domino against the full tile is alone");
}

struct Track {
    two: Vec<i64>,
    three: Vec<i64>,
    gaskets: Vec<usize>,
    reached: Vec<usize>,
    last: Vec<usize>,
}

impl Track {
    fn new(word: &[u8], gasket: u8) -> Track {
        let n = word.len();
        let mut two = vec![0i64; n + 1];
        let mut three = vec![0i64; n + 1];
        let mut gaskets: Vec<usize> = Vec::new();
        let mut reached = vec![0usize; n + 1];
        let mut last = vec![0usize; n + 1];
        for i in 1..=n {
            let heavy = word[i - 1] == gasket;
            two[i] = two[i - 1] + if heavy { 0 } else { 1 };
            three[i] = three[i - 1] + if heavy { 1 } else { 0 };
            if heavy {
                gaskets.push(i);
                last[i] = last[i - 1];
            } else {
                last[i] = i;
            }
            reached[i] = gaskets.len();
        }
        Track {
            two,
            three,
            gaskets,
            reached,
            last,
        }
    }

    fn log2_fill(&self, at: usize, scale: f64) -> f64 {
        self.two[at] as f64 + self.three[at] as f64 * scale
    }

    fn top(&self, length: usize) -> Option<usize> {
        let cut = self.last[length];
        if cut == 0 {
            return None;
        }
        let reach = self.reached[cut];
        if reach == 0 {
            None
        } else {
            Some(self.gaskets[reach - 1])
        }
    }

    fn ratio(&self, length: usize, scale: f64) -> Option<f64> {
        let top = self.top(length)?;
        let base = self.log2_fill(top - 1, scale);
        let reach = self.reached[self.last[length]];
        let mut total = (-base).exp2();
        for index in (0..reach).rev() {
            let place = self.gaskets[index];
            let term = (self.log2_fill(place - 1, scale) - base).exp2();
            if term < 1e-25 {
                break;
            }
            total += term;
        }
        Some(total)
    }

    fn log2_comp(&self, length: usize, scale: f64) -> f64 {
        match self.top(length) {
            None => 0.0,
            Some(top) => {
                self.log2_fill(top - 1, scale)
                    + self
                        .ratio(length, scale)
                        .expect("a top place has a ratio")
                        .log2()
            }
        }
    }
}

fn gcd(a: u128, b: u128) -> u128 {
    let (mut a, mut b) = (a, b);
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    if a == 0 {
        1
    } else {
        a
    }
}

fn exact_saturation(word: &[u8], heavy: u8, light: u8, length: usize) -> (u128, u128) {
    let comp = count(Shape::GasketDomino, &word[..length], heavy, light) as u128;
    let fill = fill_of(&word[..length]) as u128;
    let g = gcd(comp, fill);
    (comp / g, fill / g)
}

fn zero_contact_log2(word: &[u8], light: u8, scale: f64) -> f64 {
    match word.iter().rposition(|code| *code == light) {
        None => 0.0,
        Some(at) => {
            let mut two = 0i64;
            let mut three = 0i64;
            let mut four = 0i64;
            for code in word[..=at].iter() {
                match code.count_ones() {
                    2 => two += 1,
                    3 => three += 1,
                    4 => four += 1,
                    _ => (),
                }
            }
            two as f64 + three as f64 * scale + four as f64 * 2.0
        }
    }
}

fn morse_word(gasket: u8, domino: u8, swap: usize, length: usize) -> Vec<u8> {
    morse(length)
        .iter()
        .map(|bit| {
            if (*bit as usize) == swap {
                gasket
            } else {
                domino
            }
        })
        .collect()
}

fn morse_value(rep: &Rep) {
    println!();
    println!("THE THUE-MORSE VALUE ON THE GASKET-DOMINO PAIRS");
    let scale = 3f64.log2();
    let target = (1.0 + scale) / 2.0;
    let bound = 108f64.ln() + 1.5f64.ln() / 2.0;
    let dominoes: Vec<u8> = ROW.iter().chain(COLUMN.iter()).copied().collect();
    let mut worst = 0.0f64;
    let mut short = 0usize;
    let mut empty = 0usize;
    for &gasket in GASKET.iter() {
        for &domino in dominoes.iter() {
            for swap in 0..2usize {
                let word = morse_word(gasket, domino, swap, RANGE);
                let track = Track::new(&word, gasket);
                for length in 1..=DEEP {
                    if value(rep, &word[..length])
                        != count(Shape::GasketDomino, &word[..length], gasket, domino)
                    {
                        short += 1;
                    }
                }
                for length in 4..=RANGE {
                    if track.top(length).is_none() {
                        empty += 1;
                    }
                    let gap = (track.log2_comp(length, scale)
                        - length as f64 * (1.0 + scale) / 2.0)
                        .abs()
                        * 2f64.ln();
                    worst = worst.max(gap);
                }
            }
        }
    }
    assert_eq!(short, 0, "the closed form reads every Thue-Morse prefix");
    assert_eq!(empty, 0, "a last gasket place exists from L = 4");
    println!("all 16 gasket-domino pairs, both letter readings, closed form against the representation to L = {DEEP}, mismatches 0, and the prefixes of length 4 or more with no gasket place at or before the last domino place number {empty}");
    let word = morse_word(7, 3, 0, RANGE);
    let track = Track::new(&word, 7);
    let other = morse_word(7, 3, 1, RANGE);
    let flip = Track::new(&other, 7);
    let marks = [256usize, 1024, 4096, RANGE];
    let line = |run: &Track| {
        marks
            .iter()
            .map(|length| {
                format!(
                    "{:.12} ({length})",
                    run.log2_comp(*length, scale) / *length as f64
                )
            })
            .collect::<Vec<_>>()
            .join(", ")
    };
    println!(
        "prefix rate log2 comp / L over (3,7), gasket at t = 0: {}",
        line(&track)
    );
    println!(
        "prefix rate log2 comp / L over (3,7), gasket at t = 1: {}",
        line(&flip)
    );
    println!("the limit is (1/2) log 6, which is {target:.15} in log 2 units and 0.895879734614027... in nats, both floats");
    println!(
        "over 4 <= L <= {RANGE}, all 16 pairs and both readings, the largest |log comp - (L/2) log 6| is {worst:.6} nats against the certificate log 108 + (1/2) log(3/2) = {bound:.6} nats, itself below 4.885"
    );
    assert!(worst < bound, "the certificate holds at every length");
    let mut deficit = 0.0f64;
    let mut ties = 0usize;
    for (name, run, letters) in [
        ("gasket at t = 0", &track, &word),
        ("gasket at t = 1", &flip, &other),
    ] {
        let sat =
            |length: usize| (run.log2_comp(length, scale) - run.log2_fill(length, scale)).exp2();
        let mut floor = (1.0f64, 0usize);
        let mut roof = (0.0f64, 0usize);
        let mut past = (0.0f64, 0usize);
        let mut over = 0usize;
        for length in 1..=RANGE {
            let here = sat(length);
            if here < floor.0 {
                floor = (here, length);
            }
            if length >= 4 && here > roof.0 {
                roof = (here, length);
            }
            if length >= 5 && here > past.0 {
                past = (here, length);
            }
            if here > 5.0 / 12.0 {
                over += 1;
            }
            deficit = deficit.max(-here.log2());
        }
        let show = |at: usize| {
            if at > EXACT {
                return format!("{:.10} at L = {at}, a float", sat(at));
            }
            let (n, d) = exact_saturation(letters, 7, 3, at);
            format!("{n}/{d} at L = {at}")
        };
        println!(
            "saturation comp/fill, {name}: over 1 <= L <= {RANGE} the minimum is {:.10} at L = {}, the value at L = 4096 is {:.10}, both floats, against the proved floor 1/108 = 0.009259259",
            floor.0,
            floor.1,
            sat(4096)
        );
        println!(
            "saturation comp/fill, {name}: the largest value at L >= 4 is {} and at L >= 5 is {}, exact rationals, and the number of lengths in the whole sweep above 5/12 is {over}, all of them below L = 4, which the certificate excludes",
            show(roof.1),
            show(past.1)
        );
    }
    for length in 1..=RANGE {
        let here = (track.log2_comp(length, scale) - track.log2_fill(length, scale)).exp2();
        if (-here.log2() - deficit).abs() < 1e-9 {
            ties += 1;
        }
    }
    println!(
        "the largest fill deficit log2 fill - log2 comp over both readings is {deficit:.4}, a float, and it is attained on a tie set, {ties} lengths in the reading gasket at t = 0 alone, so no single length may be named for it; the proved ceiling is log2 108 = {:.4}",
        108f64.log2()
    );
    assert!(
        deficit < 108f64.log2(),
        "the deficit stays under its ceiling"
    );
    let mut rates: Vec<String> = Vec::new();
    for big in [7u8, 15] {
        for swap in 0..2usize {
            let word: Vec<u8> = morse(RANGE)
                .iter()
                .map(|bit| if (*bit as usize) == swap { big } else { 6 })
                .collect();
            let exponent = zero_contact_log2(&word, 6, scale);
            let first = if swap == 0 { big } else { 6 };
            let second = if swap == 0 { 6 } else { big };
            rates.push(format!(
                "({first},{second}) {:.12}",
                exponent / RANGE as f64
            ));
        }
    }
    println!(
        "the same word over four zero-contact pairs at L = {RANGE}, log 2 units, floats: {}, against (1/2) log2 6 = {target:.12} and (1/2) log2 8 = 1.5",
        rates.join(", ")
    );
}

fn named_words() -> Vec<(String, Vec<u8>)> {
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    out.push(("Thue-Morse (3,7)".into(), morse_word(7, 3, 0, 4096)));
    out.push(("Thue-Morse swapped".into(), morse_word(7, 3, 1, 4096)));
    out.push((
        "(3,7)^2000".into(),
        (0..4000).map(|i| if i % 2 == 0 { 3 } else { 7 }).collect(),
    ));
    out.push((
        "(7,3)^2000".into(),
        (0..4000).map(|i| if i % 2 == 0 { 7 } else { 3 }).collect(),
    ));
    let mut rng = Rng::new(SEED);
    out.push((
        format!("Bernoulli(1/2) at seed {SEED}"),
        (0..4000)
            .map(|_| if rng.boolean() { 7 } else { 3 })
            .collect(),
    ));
    let mut block: Vec<u8> = vec![3; 2000];
    block.extend(vec![7; 2000]);
    block.extend(vec![3; 2000]);
    out.push(("3^2000 7^2000 3^2000".into(), block));
    out
}

fn sandwich() {
    println!();
    println!("THE SANDWICH");
    let scale = 3f64.log2();
    println!("T < comp <= 1 + (3/2) T at T = fill(w_1..i*-1), i* the last gasket place at or before the last domino place");
    for (name, word) in named_words().iter() {
        let track = Track::new(word, 7);
        let mut low = f64::MAX;
        let mut high = 0.0f64;
        let mut bad = 0usize;
        for length in 1..=word.len() {
            let ratio = match track.ratio(length, scale) {
                None => continue,
                Some(value) => value,
            };
            low = low.min(ratio);
            high = high.max(ratio);
            let top = track.top(length).expect("a ratio needs a top place");
            let size = track.log2_fill(top - 1, scale);
            if ratio > 1.5 + (-size).exp2() + 1e-12 || ratio <= 1.0 {
                bad += 1;
            }
        }
        println!("{name}: comp / T in [{low:.4}, {high:.4}], floats, violations {bad}");
        assert_eq!(bad, 0, "the sandwich holds at every length");
    }
}

fn text(value: &Frac) -> String {
    let (num, den) = value.parts();
    if den == 1 {
        format!("{num}")
    } else {
        format!("{num}/{den}")
    }
}

fn chart_step(state: &[Frac; 4], matrix: &[Vec<i64>], divisor: i64) -> [Frac; 4] {
    let mut out = [Frac::zero(); 4];
    for i in 0..4 {
        for j in 0..4 {
            out[j] = out[j].add(&state[i].mul(&Frac::int(matrix[i][j])));
        }
    }
    for slot in out.iter_mut() {
        *slot = slot.div(&Frac::int(divisor));
    }
    out
}

fn image(point: (Frac, Frac), gasket: bool) -> (Frac, Frac) {
    let one = Frac::int(1);
    let (b, c) = point;
    if gasket {
        (
            one.add(&b).div(&Frac::int(3)),
            one.add(&c).div(&Frac::int(3)),
        )
    } else {
        (one.add(&b).div(&Frac::int(2)), Frac::zero())
    }
}

fn inside(point: &(Frac, Frac)) -> bool {
    let (b, c) = point;
    let zero = Frac::zero();
    let one = Frac::int(1);
    let half = Frac::new(1, 2);
    !b.below(&zero) && !one.below(b) && !c.below(&zero) && !half.below(c) && !one.below(&b.add(c))
}

fn cone(rep: &Rep) {
    println!();
    println!("THE INVARIANT CONE AND WHY IT IS NOT THE MECHANISM");
    println!("phi = (1,2,2,4)^T is a common RIGHT eigenvector, M_c phi = popcount(c) phi, so it normalises the ROW orbit lambda M_(c_1) ... M_(c_k) alone");
    let phi = [1i64, 2, 2, 4];
    let mut bad = 0usize;
    for &code in CODES.iter() {
        let matrix = &rep.matrices[&code];
        for i in 0..4 {
            let row: i64 = (0..4).map(|j| matrix[i][j] * phi[j]).sum();
            if row != code.count_ones() as i64 * phi[i] {
                bad += 1;
            }
        }
    }
    assert_eq!(bad, 0, "phi is a right eigenvector at every code");
    println!("checked on all 15 code matrices, mismatches {bad}; M_gasket gamma = M_full gamma = gamma, M_unit gamma = phi, M_diagonal gamma = 2 phi, and M_diagonal = 2 phi lambda");
    let word: [u8; 7] = [7, 3, 7, 7, 3, 3, 7];
    let mut state = [Frac::int(1), Frac::zero(), Frac::zero(), Frac::zero()];
    let mut point = (Frac::zero(), Frac::zero());
    let mut orbit: Vec<String> = Vec::new();
    let mut drift = 0usize;
    for code in word.iter() {
        state = chart_step(&state, &rep.matrices[code], code.count_ones() as i64);
        point = image(point, *code == 7);
        if state[1] != point.0 || state[2] != point.1 || state[3] != Frac::zero() {
            drift += 1;
        }
        orbit.push(format!("({},{})", text(&point.0), text(&point.1)));
    }
    assert_eq!(drift, 0, "the chart is the normalised row orbit");
    println!("in the chart n_4 = 0, (b,c) = (n_2,n_3) and comp/fill = 1 - b - c, the two maps are N_gasket(b,c) = ((1+b)/3, (1+c)/3) and N_domino(b,c) = ((1+b)/2, 0)");
    println!("along the word 7,3,7,7,3,3,7 the chart reads {}, matching the raw matrices at every step, mismatches {drift}", orbit.join(" "));
    let vertices = [
        (Frac::zero(), Frac::zero()),
        (Frac::int(1), Frac::zero()),
        (Frac::zero(), Frac::new(1, 2)),
        (Frac::new(1, 2), Frac::new(1, 2)),
    ];
    let mut escapes = 0usize;
    for vertex in vertices.iter() {
        for gasket in [true, false] {
            if !inside(&image(*vertex, gasket)) {
                escapes += 1;
            }
        }
    }
    assert_eq!(escapes, 0, "S is invariant under both maps");
    println!("S = {{0 <= b <= 1, 0 <= c <= 1/2, b + c <= 1}} is invariant under both maps, checked at all four vertices, escapes {escapes}");
    let mut images: Vec<String> = Vec::new();
    let mut widest = Frac::zero();
    for vertex in vertices.iter() {
        let mut at = *vertex;
        for gasket in [true, false, true] {
            at = image(at, gasket);
        }
        let sum = at.0.add(&at.1);
        if widest.below(&sum) {
            widest = sum;
        }
        images.push(format!("({},{})", text(&at.0), text(&at.1)));
    }
    println!("gasket-domino-gasket sends the vertices to {} with largest b + c = {}, strictly inside S, so the pair semigroup is primitive in this chart", images.join(" "), text(&widest));
    let edge = image(image((Frac::int(1), Frac::zero()), false), true);
    println!(
        "domino-gasket sends (1,0) to ({},{}) with b + c = 1, on the face, so length 3 is minimal",
        text(&edge.0),
        text(&edge.1)
    );
    let mut positive = 0usize;
    let mut products = 0usize;
    for length in 1..=12usize {
        for mask in 0..(1usize << length) {
            let mut matrix: Vec<Vec<i128>> = (0..4)
                .map(|i| (0..4).map(|j| if i == j { 1i128 } else { 0 }).collect())
                .collect();
            for i in 0..length {
                let code = if (mask >> i) & 1 == 1 { 7u8 } else { 3u8 };
                let next = &rep.matrices[&code];
                let mut out = vec![vec![0i128; 4]; 4];
                for r in 0..4 {
                    for k in 0..4 {
                        if matrix[r][k] == 0 {
                            continue;
                        }
                        for c in 0..4 {
                            out[r][c] += matrix[r][k] * next[k][c] as i128;
                        }
                    }
                }
                matrix = out;
            }
            products += 1;
            if matrix.iter().all(|row| row.iter().all(|entry| *entry > 0)) {
                positive += 1;
            }
        }
    }
    let negatives = |code: u8| {
        rep.matrices[&code]
            .iter()
            .flat_map(|row| row.iter())
            .filter(|entry| **entry < 0)
            .count()
    };
    println!(
        "entrywise positivity in the standard basis is the wrong test and fails: M_3 carries {} negative entries and M_7 carries {}, and {positive} of the {products} products of {{M_3, M_7}} of length 1 to 12 are entrywise positive",
        negatives(3),
        negatives(7)
    );
    assert_eq!(positive, 0, "no product is entrywise positive");
    let mut power: Vec<Vec<i128>> = (0..4)
        .map(|i| (0..4).map(|j| if i == j { 1i128 } else { 0 }).collect())
        .collect();
    let mut norms: Vec<String> = Vec::new();
    let mut wrong = 0usize;
    for length in 1..=32usize {
        let next = &rep.matrices[&3];
        let mut out = vec![vec![0i128; 4]; 4];
        for r in 0..4 {
            for k in 0..4 {
                if power[r][k] == 0 {
                    continue;
                }
                for c in 0..4 {
                    out[r][c] += power[r][k] * next[k][c] as i128;
                }
            }
        }
        power = out;
        let peak = power
            .iter()
            .flat_map(|row| row.iter())
            .map(|entry| entry.abs())
            .max()
            .expect("a matrix has entries");
        if peak != (1i128 << (length + 2)) - 2 {
            wrong += 1;
        }
        if [1usize, 2, 4, 8, 16, 32].contains(&length) {
            norms.push(format!("{peak}"));
        }
        if value(rep, &vec![3u8; length]) != 1 {
            wrong += 1;
        }
    }
    println!("the norm exponent is a different number from the component exponent: max |entry(M_3^L)| = 2^(L+2) - 2 reads {} at L = 1, 2, 4, 8, 16, 32, so the norm rate along 3^inf is log 2, while comp(A_(3^L)) = 1 at every L to 32 and the component rate is 0; mismatches {wrong}", norms.join(", "));
    assert_eq!(wrong, 0, "the norm witness is exact");
}

fn boundary(rep: &Rep) {
    println!();
    println!("THE BOUNDARY FREQUENCY (1,0) OVER (3,7)");
    let scale = 3f64.log2();
    let squares = |i: usize| {
        let root = (i as f64).sqrt() as usize;
        (root.saturating_sub(1)..root + 2).any(|r| r * r == i)
    };
    let powers: Vec<u8> = (1..=FAR)
        .map(|i| if i.is_power_of_two() { 7 } else { 3 })
        .collect();
    let boxes: Vec<u8> = (1..=FAR).map(|i| if squares(i) { 7 } else { 3 }).collect();
    let flats = vec![3u8; FAR];
    let mut bad = 0usize;
    for word in [&powers, &boxes, &flats] {
        for length in 1..=DEEP {
            if value(rep, &word[..length]) != count(Shape::GasketDomino, &word[..length], 7, 3) {
                bad += 1;
            }
        }
    }
    assert_eq!(bad, 0, "the closed form covers the three boundary words");
    println!("three words whose gasket frequency is 0, closed form against the representation to L = {DEEP}, mismatches {bad}");
    let run = Track::new(&powers, 7);
    let marks = [2048usize, 2049, 4096, 4097, 8192, 8193, 16384, 16385, 32768];
    let rates = marks
        .iter()
        .map(|length| {
            format!(
                "{:.9} ({length})",
                run.log2_comp(*length, scale) / *length as f64
            )
        })
        .collect::<Vec<_>>()
        .join(", ");
    println!("the gasket at the powers of 2, prefix rate in log 2 units, floats: {rates}");
    let mut low = (f64::MAX, 0usize);
    let mut high = (0.0f64, 0usize);
    for length in 4097..=8192usize {
        let rate = run.log2_comp(length, scale) / length as f64;
        if rate < low.0 {
            low = (rate, length);
        }
        if rate > high.0 {
            high = (rate, length);
        }
    }
    println!(
        "over the single block 4097 <= L <= 8192 the rate sweeps from {:.5} at L = {} down to {:.5} at L = {}, floats, so the accumulation set is the whole interval and not two points",
        high.0, high.1, low.0, low.1
    );
    let boxed = Track::new(&boxes, 7);
    let squared = [4096usize, 8192, 16384, 32768]
        .iter()
        .map(|length| {
            format!(
                "{:.9} ({length})",
                boxed.log2_comp(*length, scale) / *length as f64
            )
        })
        .collect::<Vec<_>>()
        .join(", ");
    println!(
        "the gasket at the squares, prefix rate in log 2 units, floats: {squared}, closing on 1"
    );
    let flat = Track::new(&flats, 7);
    let steady = (1..=FAR)
        .map(|length| flat.log2_comp(length, scale))
        .fold(0.0f64, f64::max);
    println!("the constant word 3^L has comp 1 and rate {steady} at every L to {FAR}");
    assert_eq!(steady, 0.0, "the constant word never grows");
}

fn by_products(rep: &Rep) {
    println!();
    println!("BY-PRODUCTS");
    let mut terms: Vec<String> = Vec::new();
    let mut bad = 0usize;
    for k in 1..=8usize {
        let word: Vec<u8> = (0..2 * k)
            .map(|i| if i % 2 == 0 { 7u8 } else { 3u8 })
            .collect();
        let want = (6i128.pow(k as u32) + 4) / 5;
        let seen = count(Shape::GasketDomino, &word, 7, 3);
        if seen != want || value(rep, &word) != want {
            bad += 1;
        }
        if k <= 6 {
            terms.push(format!("{seen}"));
        }
    }
    println!("comp(A_((7,3)^k)) = (6^k + 4)/5, reading {} at k = 1 to 6, checked to k = 8 against the closed form and the representation, mismatches {bad}", terms.join(", "));
    assert_eq!(bad, 0, "the periodic word has that closed form");
    let mut best = 0i128;
    let mut witness: Vec<u8> = Vec::new();
    for mask in 0..(1usize << 8) {
        let word = pair_word(7, 3, mask, 8);
        let seen = count(Shape::GasketDomino, &word, 7, 3);
        if seen > best {
            best = seen;
            witness = word.clone();
        }
    }
    println!(
        "the largest component count at L = 8 over (3,7) is {best} at the word {}, and 1094 = 2 x 547, while every closed form on the other 89 pairs gives a count of the form 2^a 3^b",
        witness
            .iter()
            .map(|code| format!("{code}"))
            .collect::<Vec<_>>()
            .join(",")
    );
    let scale = 3f64.log2();
    let mut block: Vec<u8> = vec![3, 7];
    while block.len() < 4096 {
        let size = block.len();
        let mut next = block.clone();
        next.extend(vec![7u8; size]);
        next.extend(vec![3u8; size]);
        block = next;
    }
    let track = Track::new(&block, 7);
    let mut low = f64::MAX;
    let mut high = 0.0f64;
    let mut thin = 1.0f64;
    for length in 1024..=4096usize {
        let rate = track.log2_comp(length, scale) / length as f64;
        low = low.min(rate);
        high = high.max(rate);
        let heavy = track.three[length] as f64 / length as f64;
        thin = thin.min(heavy.min(1.0 - heavy));
    }
    println!("the tripling word W_(k+1) = W_k 7^|W_k| 3^|W_k| at seed 3,7 has both letters at density at least {thin:.4} over 1024 <= L <= 4096 yet its prefix rate ranges over [{low:.4}, {high:.4}] in log 2 units with no narrowing, all floats");
}

pub fn study(rep: &Rep) {
    println!();
    println!("THE FORTY-SIX");
    forms(rep);
    ledger();
    morse_value(rep);
    sandwich();
    cone(rep);
    boundary(rep);
    by_products(rep);
}
