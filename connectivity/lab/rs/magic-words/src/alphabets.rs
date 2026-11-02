use crate::unsolved::exponent_weights;
use crate::word::{render, CODES};
use std::collections::HashMap;

pub const REACH: [usize; 2] = [7, 6];
pub const SHORT: usize = 5;
pub const READ: usize = 1 << 14;

#[derive(Clone, Copy, PartialEq)]
enum Class {
    Unit,
    Row,
    Column,
    Diagonal,
    Gasket,
    Full,
}

fn contacts(code: u8) -> (u32, u32) {
    let h = (code & 3 == 3) as u32 + (code & 12 == 12) as u32;
    let v = (code & 5 == 5) as u32 + (code & 10 == 10) as u32;
    (h, v)
}

fn free(code: u8) -> bool {
    contacts(code) == (0, 0)
}

fn class(code: u8) -> Class {
    let (h, v) = contacts(code);
    match code.count_ones() {
        1 => Class::Unit,
        2 if h == 1 => Class::Row,
        2 if v == 1 => Class::Column,
        2 => Class::Diagonal,
        3 => Class::Gasket,
        _ => Class::Full,
    }
}

fn step(state: [u64; 4], code: u8) -> [u64; 4] {
    let [k, h, v, f] = state;
    match class(code) {
        Class::Unit => [f, f, f, f],
        Class::Diagonal => [2 * f, 2 * f, 2 * f, 2 * f],
        Class::Row => [h, h, 2 * f, 2 * f],
        Class::Column => [v, 2 * f, v, 2 * f],
        Class::Gasket => [k, h + f, v + f, 3 * f],
        Class::Full => [k, 2 * h, 2 * v, 4 * f],
    }
}

fn add(a: f64, b: f64) -> f64 {
    let (high, low) = if a > b { (a, b) } else { (b, a) };
    high + (low - high).exp().ln_1p()
}

fn log_step(state: [f64; 4], code: u8) -> [f64; 4] {
    let [k, h, v, f] = state;
    let two = 2f64.ln();
    match class(code) {
        Class::Unit => [f, f, f, f],
        Class::Diagonal => [f + two; 4],
        Class::Row => [h, h, f + two, f + two],
        Class::Column => [v, f + two, v, f + two],
        Class::Gasket => [k, add(h, f), add(v, f), f + 3f64.ln()],
        Class::Full => [k, h + two, v + two, f + 2.0 * two],
    }
}

fn spread(x: u64) -> u128 {
    let half = |y: u32| -> u128 {
        let mut z = y as u64;
        z = (z | (z << 16)) & 0x0000_FFFF_0000_FFFF;
        z = (z | (z << 8)) & 0x00FF_00FF_00FF_00FF;
        z = (z | (z << 4)) & 0x0F0F_0F0F_0F0F_0F0F;
        z = (z | (z << 2)) & 0x3333_3333_3333_3333;
        z = (z | (z << 1)) & 0x5555_5555_5555_5555;
        z as u128
    };
    half(x as u32) | (half((x >> 32) as u32) << 64)
}

fn grow(rows: &[u128], code: u8) -> Vec<u128> {
    let lay = |s: u128, pattern: u8| -> u128 {
        (if pattern & 1 == 1 { s } else { 0 }) | (if pattern & 2 == 2 { s << 1 } else { 0 })
    };
    let mut out = Vec::with_capacity(rows.len() * 2);
    for row in rows {
        assert!(*row >> 64 == 0, "a parent row fits 64 cells");
        let s = spread(*row as u64);
        out.push(lay(s, code & 3));
        out.push(lay(s, (code >> 2) & 3));
    }
    out
}

fn draw(word: &[u8]) -> Vec<u128> {
    let mut rows = vec![1u128];
    for code in word {
        rows = grow(&rows, *code);
    }
    rows
}

fn find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

fn observe(rows: &[u128]) -> [u64; 4] {
    let mut parent: Vec<usize> = Vec::with_capacity(rows.len() * 8);
    let mut pieces = 0u64;
    let mut previous: Vec<(u32, u32, usize)> = Vec::with_capacity(64);
    let mut current: Vec<(u32, u32, usize)> = Vec::with_capacity(64);
    let mut above = 0u128;
    let (mut h, mut v, mut f) = (0u64, 0u64, 0u64);
    for row in rows {
        f += row.count_ones() as u64;
        h += (row & !(row << 1)).count_ones() as u64;
        v += (row & !above).count_ones() as u64;
        above = *row;
        current.clear();
        let mut x = *row;
        while x != 0 {
            let start = x.trailing_zeros();
            let length = (!(x >> start)).trailing_zeros();
            let id = parent.len();
            parent.push(id);
            pieces += 1;
            current.push((start, start + length - 1, id));
            x = if start + length >= 128 {
                0
            } else {
                x & !((1u128 << (start + length)) - 1)
            };
        }
        let (mut i, mut j) = (0usize, 0usize);
        while i < previous.len() && j < current.len() {
            let (pa, pb, pid) = previous[i];
            let (ca, cb, cid) = current[j];
            if pa <= cb && ca <= pb {
                let (x, y) = (find(&mut parent, pid), find(&mut parent, cid));
                if x != y {
                    parent[x] = y;
                    pieces -= 1;
                }
            }
            if pb < cb {
                i += 1;
            } else {
                j += 1;
            }
        }
        std::mem::swap(&mut previous, &mut current);
    }
    [pieces, h, v, f]
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
    if class(code) == Class::Diagonal {
        (1, 0)
    } else {
        (0, 0)
    }
}

fn nats(pair: (i64, i64)) -> f64 {
    pair.0 as f64 * 2f64.ln() + pair.1 as f64 * 3f64.ln()
}

fn has(set: &[u8], which: Class) -> bool {
    set.iter().any(|code| class(*code) == which)
}

fn saturates(set: &[u8]) -> bool {
    let (row, column, gasket) = (
        has(set, Class::Row),
        has(set, Class::Column),
        has(set, Class::Gasket),
    );
    set.iter().any(|code| free(*code)) || (row && column) || ((row || column) && gasket)
}

fn law(set: &[u8]) -> Vec<(i64, i64)> {
    let domino = has(set, Class::Row) || has(set, Class::Column);
    set.iter()
        .map(|code| {
            if saturates(set) {
                weight(*code)
            } else if domino && class(*code) == Class::Full {
                (1, 0)
            } else {
                (0, 0)
            }
        })
        .collect()
}

fn law_text(set: &[u8]) -> &'static str {
    if saturates(set) {
        "the fill exponent"
    } else if has(set, Class::Full) && (has(set, Class::Row) || has(set, Class::Column)) {
        "f_15 log 2"
    } else {
        "0"
    }
}

fn kind(set: &[u8]) -> String {
    if set.iter().any(|code| free(*code)) {
        return "a contact-free letter".to_string();
    }
    let count = |which: Class| set.iter().filter(|code| class(**code) == which).count();
    let (rows, columns, gaskets) = (
        count(Class::Row),
        count(Class::Column),
        count(Class::Gasket),
    );
    let dominoes = rows + columns;
    let mut parts: Vec<String> = Vec::new();
    if rows > 0 && columns > 0 {
        parts.push(format!("{dominoes} dominoes of both orientations"));
    } else if dominoes == 1 {
        parts.push("1 domino".to_string());
    } else if dominoes > 1 {
        parts.push(format!("{dominoes} like dominoes"));
    }
    if gaskets > 0 {
        parts.push(format!(
            "{gaskets} {}",
            if gaskets == 1 { "gasket" } else { "gaskets" }
        ));
    }
    if has(set, Class::Full) {
        parts.push("the full tile".to_string());
    }
    parts.join(", ")
}

fn triples() -> Vec<[u8; 3]> {
    let mut out = Vec::new();
    for (i, &a) in CODES.iter().enumerate() {
        for (j, &b) in CODES.iter().enumerate().skip(i + 1) {
            for &c in &CODES[j + 1..] {
                out.push([a, b, c]);
            }
        }
    }
    out
}

fn letters() {
    let loose: Vec<u8> = CODES.iter().copied().filter(|code| free(*code)).collect();
    let bound: Vec<u8> = CODES.iter().copied().filter(|code| !free(*code)).collect();
    println!(
        "contact-free codes, h = v = 0: {} ({}); the other {}: {}",
        loose
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join(", "),
        loose.len(),
        bound.len(),
        bound
            .iter()
            .map(|c| c.to_string())
            .collect::<Vec<_>>()
            .join(", ")
    );
    assert!(
        loose
            .iter()
            .all(|code| matches!(class(*code), Class::Unit | Class::Diagonal)),
        "contact-free is unit or diagonal"
    );
    assert_eq!(loose.len(), 6);
    let all = triples();
    let mut by_free = [0usize; 4];
    for set in all.iter() {
        by_free[set.iter().filter(|code| free(**code)).count()] += 1;
    }
    println!(
        "three-letter alphabets {}: with 0, 1, 2, 3 contact-free letters {}, {}, {}, {}; with at least one {}",
        all.len(),
        by_free[0],
        by_free[1],
        by_free[2],
        by_free[3],
        by_free[1] + by_free[2] + by_free[3]
    );
    assert_eq!((all.len(), by_free[0]), (455, 84));
}

fn ledgers() {
    let mut pairs = 0usize;
    for (i, &a) in CODES.iter().enumerate() {
        for &b in &CODES[i + 1..] {
            let (wa, wb) = exponent_weights(a, b);
            assert_eq!(
                law(&[a, b]),
                vec![wa, wb],
                "the law reads the two-letter table"
            );
            pairs += 1;
        }
    }
    println!("the law reproduces the two-letter exponent weights on all {pairs} pairs");
    let mut kinds: Vec<(String, &'static str, usize, usize)> = Vec::new();
    let (mut saturating, mut exact, mut between) = (0usize, 0usize, 0usize);
    for set in triples() {
        let weights = law(&set);
        if weights
            .iter()
            .zip(set.iter())
            .all(|(w, code)| *w == weight(*code))
        {
            saturating += 1;
        }
        if weights
            .iter()
            .zip(set.iter())
            .all(|(w, code)| *w == phi_weight(*code))
        {
            exact += 1;
        }
        let rate: f64 = weights.iter().map(|w| nats(*w)).sum::<f64>() / 3.0;
        let phi: f64 = set.iter().map(|c| nats(phi_weight(*c))).sum::<f64>() / 3.0;
        let ceiling: f64 = set.iter().map(|c| nats(weight(*c))).sum::<f64>() / 3.0;
        if rate > phi + 1e-12 && rate < ceiling - 1e-12 {
            between += 1;
        }
        let name = kind(&set);
        let classes = {
            let mut seen: Vec<u8> = Vec::new();
            for code in set.iter() {
                let c = class(*code) as u8;
                if !seen.contains(&c) {
                    seen.push(c);
                }
            }
            if set.iter().any(|code| free(*code)) {
                0
            } else {
                seen.len()
            }
        };
        match kinds.iter_mut().find(|(at, _, _, _)| *at == name) {
            Some((_, _, n, _)) => *n += 1,
            None => kinds.push((name, law_text(&set), 1, classes)),
        }
    }
    println!(
        "{:<46} {:>9} {:>8}  exponent at interior frequency",
        "alphabet kind", "alphabets", "classes"
    );
    for (name, text, n, classes) in kinds.iter() {
        let shown = if *classes == 0 {
            "any".to_string()
        } else {
            classes.to_string()
        };
        println!("{:<46} {:>9} {:>8}  {}", name, n, shown, text);
    }
    println!(
        "three letters: the exponent saturates the fill ceiling on {saturating} of 455, Phi is exact on {exact} and refuted on {}, and {between} sit strictly between Phi and the ceiling at equal frequencies",
        455 - exact
    );
    assert_eq!((saturating, exact, between), (443, 30, 2));
    let (mut every, mut full, mut flat, mut middle) = (0usize, 0usize, 0usize, 0usize);
    for mask in 1u32..(1 << 15) {
        let set: Vec<u8> = (0..15)
            .filter(|i| (mask >> i) & 1 == 1)
            .map(|i| CODES[i])
            .collect();
        every += 1;
        let weights = law(&set);
        if weights
            .iter()
            .zip(set.iter())
            .all(|(w, code)| *w == weight(*code))
        {
            full += 1;
        }
        if weights
            .iter()
            .zip(set.iter())
            .all(|(w, code)| *w == phi_weight(*code))
        {
            flat += 1;
        }
        let n = set.len() as f64;
        let rate: f64 = weights.iter().map(|w| nats(*w)).sum::<f64>() / n;
        let phi: f64 = set.iter().map(|c| nats(phi_weight(*c))).sum::<f64>() / n;
        let ceiling: f64 = set.iter().map(|c| nats(weight(*c))).sum::<f64>() / n;
        if rate > phi + 1e-12 && rate < ceiling - 1e-12 {
            middle += 1;
        }
    }
    println!(
        "every alphabet: of the {every} non-empty sets of codes the exponent saturates on {full} and falls short on {}, Phi is exact on {flat}, and {middle} sit strictly between at equal frequencies",
        every - full
    );
    assert_eq!((every, full), (32767, 32724));
}

#[derive(Default)]
struct Tally {
    words: usize,
    frame: usize,
    reading: usize,
    cuts: usize,
    floor: usize,
    crossed: usize,
    crossing: usize,
    witness: Option<Vec<u8>>,
    memo: HashMap<Vec<u8>, u64>,
}

fn terminal(word: &[u8], kinds: &[Class]) -> u32 {
    word.iter()
        .rev()
        .take_while(|code| kinds.contains(&class(**code)))
        .count() as u32
}

fn check(
    set: &[u8; 3],
    word: &[u8],
    frame: [u64; 4],
    seen: [u64; 4],
    drawn: &[[u64; 4]],
    tally: &mut Tally,
) {
    tally.words += 1;
    let mut bad = false;
    if frame != seen {
        tally.frame += 1;
        bad = true;
    }
    if set.iter().any(|code| free(*code)) {
        let mut fill = 1u64;
        let mut wrong = false;
        for (p, code) in word.iter().enumerate() {
            fill *= code.count_ones() as u64;
            if free(*code) {
                tally.cuts += 1;
                let suffix = &word[p + 1..];
                let pieces = match tally.memo.get(suffix) {
                    Some(n) => *n,
                    None => {
                        let n = observe(&draw(suffix))[0];
                        tally.memo.insert(suffix.to_vec(), n);
                        n
                    }
                };
                if fill * pieces != seen[0] {
                    wrong = true;
                }
            }
        }
        if wrong {
            tally.reading += 1;
            bad = true;
        }
    } else {
        let last = word
            .iter()
            .rposition(|code| matches!(class(*code), Class::Row | Class::Column));
        let want = match last {
            None => 1,
            Some(m) if class(word[m]) == Class::Row => drawn[m][1],
            Some(m) => drawn[m][2],
        };
        if want != seen[0] {
            tally.reading += 1;
            bad = true;
        }
    }
    let [_, h, v, f] = seen;
    let t = terminal(word, &[Class::Row, Class::Full]);
    let s = terminal(word, &[Class::Column, Class::Full]);
    if h > f || v > f || f > 3 * (1u64 << t) * h || f > 3 * (1u64 << s) * v {
        tally.floor += 1;
        bad = true;
    }
    if word.len() <= SHORT {
        tally.crossing += 1;
        let grid = render(word);
        if grid.components() != seen[0] || grid.fill() != seen[3] {
            tally.crossed += 1;
            bad = true;
        }
    }
    if bad && tally.witness.is_none() {
        tally.witness = Some(word.to_vec());
    }
}

fn descend(
    set: &[u8; 3],
    reach: usize,
    word: &mut Vec<u8>,
    pics: &mut Vec<Vec<u128>>,
    frames: &mut Vec<[u64; 4]>,
    drawn: &mut Vec<[u64; 4]>,
    tally: &mut Tally,
) {
    if word.len() == reach {
        return;
    }
    for &code in set.iter() {
        let pic = grow(pics.last().expect("a picture"), code);
        let frame = step(*frames.last().expect("a frame"), code);
        let seen = observe(&pic);
        word.push(code);
        if set.iter().all(|letter| word.contains(letter)) {
            check(set, word, frame, seen, drawn, tally);
        }
        pics.push(pic);
        frames.push(frame);
        drawn.push(seen);
        descend(set, reach, word, pics, frames, drawn, tally);
        pics.pop();
        frames.pop();
        drawn.pop();
        word.pop();
    }
}

fn census() {
    let mut groups = [Tally::default(), Tally::default()];
    let mut alphabets = [0usize; 2];
    for set in triples() {
        let at = if set.iter().any(|code| free(*code)) {
            1
        } else {
            0
        };
        alphabets[at] += 1;
        let mut word: Vec<u8> = Vec::new();
        let mut pics = vec![vec![1u128]];
        let mut frames = vec![[1u64; 4]];
        let mut drawn = vec![[1u64; 4]];
        descend(
            &set,
            REACH[at],
            &mut word,
            &mut pics,
            &mut frames,
            &mut drawn,
            &mut groups[at],
        );
    }
    let labels = ["no contact-free letter", "a contact-free letter"];
    let readings = [
        "comp = H or V of the prefix before the last domino place",
        "comp = fill(prefix to p) * comp(suffix) at every contact-free place p",
    ];
    for at in 0..2 {
        let t = &groups[at];
        println!(
            "{} alphabets with {}: {} words of length 3 to {} carrying all three letters, drawn (comp, H, V, fill) against the frame product {} mismatches; {}: {} mismatches{}; run floor fill/(3 * 2^t) <= H <= fill and its transpose: {} violations",
            alphabets[at],
            labels[at],
            t.words,
            REACH[at],
            t.frame,
            readings[at],
            t.reading,
            if at == 1 { format!(" over {} cut places", t.cuts) } else { String::new() },
            t.floor
        );
        println!(
            "  the bit-row drawing against the cell renderer on the {} of those words of length at most {SHORT}: {} mismatches; first failing word {}",
            t.crossing,
            t.crossed,
            t.witness.as_ref().map(|w| format!("{w:?}")).unwrap_or_else(|| "none".to_string())
        );
        assert_eq!(
            (t.frame, t.reading, t.floor, t.crossed),
            (0, 0, 0, 0),
            "the frame, the reading and the floor hold on every drawn word"
        );
    }
    assert_eq!(
        (groups[0].words, groups[1].words),
        (84 * (6 + 36 + 150 + 540 + 1806), 371 * (6 + 36 + 150 + 540))
    );
}

fn fixed(images: [&[u8]; 3], length: usize) -> Vec<u8> {
    let mut word = vec![0u8];
    while word.len() < length {
        word = word
            .iter()
            .flat_map(|letter| images[*letter as usize].iter().copied())
            .collect();
    }
    word.truncate(length);
    word
}

fn readings() {
    let words: [(&str, Vec<u8>); 3] = [
        ("(0 1 2)^inf", (0..READ).map(|i| (i % 3) as u8).collect()),
        (
            "0 -> 012, 1 -> 02, 2 -> 1",
            fixed([&[0, 1, 2], &[0, 2], &[1]], READ),
        ),
        (
            "0 -> 01, 1 -> 02, 2 -> 0",
            fixed([&[0, 1], &[0, 2], &[0]], READ),
        ),
    ];
    for (name, word) in words.iter() {
        let counts: Vec<usize> = (0..3u8)
            .map(|l| word.iter().filter(|x| **x == l).count())
            .collect();
        println!(
            "word {name}: letter frequencies at level {READ} {:.6}, {:.6}, {:.6}",
            counts[0] as f64 / READ as f64,
            counts[1] as f64 / READ as f64,
            counts[2] as f64 / READ as f64
        );
    }
    let mut table: Vec<(String, usize, [f64; 3], f64)> = Vec::new();
    for set in triples() {
        let weights = law(&set);
        let mut deviation = [0f64; 3];
        let mut rates = [0f64; 3];
        for (at, (_, word)) in words.iter().enumerate() {
            let mut state = [0f64; 4];
            let mut counts = [0usize; 3];
            for letter in word.iter() {
                state = log_step(state, set[*letter as usize]);
                counts[*letter as usize] += 1;
            }
            let rate = state[0] / READ as f64;
            let want: f64 = (0..3)
                .map(|i| counts[i] as f64 / READ as f64 * nats(weights[i]))
                .sum();
            rates[at] = rate;
            deviation[at] = (rate - want).abs();
        }
        let blind = (rates[0] - rates[1]).abs();
        let name = kind(&set);
        match table.iter_mut().find(|(at, _, _, _)| *at == name) {
            Some((_, n, worst, gap)) => {
                *n += 1;
                for i in 0..3 {
                    worst[i] = worst[i].max(deviation[i]);
                }
                *gap = gap.max(blind);
            }
            None => table.push((name, 1, deviation, blind)),
        }
    }
    println!(
        "{:<46} {:>9} {:>12} {:>12} {:>12} {:>12}",
        "alphabet kind", "alphabets", "periodic", "ternary TM", "tribonacci", "|per - TM|"
    );
    let mut worst = 0f64;
    for (name, n, dev, gap) in table.iter() {
        println!(
            "{:<46} {:>9} {:>12.3e} {:>12.3e} {:>12.3e} {:>12.3e}",
            name, n, dev[0], dev[1], dev[2], gap
        );
        worst = worst.max(dev[0]).max(dev[1]).max(dev[2]);
    }
    println!("largest |log comp / level - law at the prefix frequencies| over all 455 alphabets and the three words at level {READ}: {worst:.3e}, floats");
    assert!(worst < 1e-3, "every reading sits on its law");
}

pub fn study() {
    println!();
    println!("EVERY ALPHABET");
    letters();
    ledgers();
    census();
    readings();
}
