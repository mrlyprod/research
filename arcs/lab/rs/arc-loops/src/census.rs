use crate::count::series;
use crate::design::Design;
use crate::fit::{fit, roots, show};
use std::collections::BTreeMap;

pub fn klein(code: u128, base: usize) -> u128 {
    let d = Design::full(code, base);
    (1..4).map(|g| d.image(g).full_code()).fold(code, u128::min)
}

pub fn describe(s: &[u128]) -> String {
    let best = (0..s.len().min(8))
        .filter_map(|k| fit(&s[k..]).map(|f| (k, f)))
        .filter(|(_, f)| f.margin >= 2 && f.poly.iter().all(|&c| c.abs() < 1 << 40))
        .min_by_key(|(k, f)| (f.poly.len(), *k));
    match best.map(|(k, f)| (k, Some(f))).unwrap_or((0, None)) {
        (_, None) => "no integer recurrence in range".to_string(),
        (k, Some(f)) => {
            let (found, rest) = roots(&f.poly);
            let tail = if rest.len() > 1 {
                format!(", factor with no integer root {}", show(&rest))
            } else {
                String::new()
            };
            let list: Vec<String> = found.iter().map(|r| r.to_string()).collect();
            format!(
                "from level {k}, order {}, margin {}, roots [{}]{tail}",
                f.poly.len() - 1,
                f.margin,
                list.join(" ")
            )
        }
    }
}

pub fn census(base: usize, top: usize, deep: usize, codes: &[u128]) {
    let mut classes: BTreeMap<u128, Vec<u128>> = BTreeMap::new();
    for &code in codes {
        classes.entry(klein(code, base)).or_default().push(code);
    }
    let mut groups: BTreeMap<Vec<u128>, Vec<u128>> = BTreeMap::new();
    for &rep in classes.keys() {
        let s = series(base, &Design::full(rep, base).tile, top);
        groups.entry(s).or_default().push(rep);
    }
    let zero = groups
        .iter()
        .filter(|(s, _)| s.iter().all(|&x| x == 0))
        .map(|(_, r)| r.len())
        .sum::<usize>();
    let silent = groups
        .iter()
        .filter(|(s, _)| s.iter().all(|&x| x == 0))
        .flat_map(|(_, r)| r.iter().map(|c| classes[c].len()))
        .sum::<usize>();
    println!("BASE {base}: {} codes, {} Klein classes, {} nonzero loop sequences, {} classes and {} codes never loop, levels 0..{top}", codes.len(), classes.len(), groups.len() - usize::from(zero > 0), zero, silent);
    let mut rows: Vec<(&Vec<u128>, &Vec<u128>)> = groups
        .iter()
        .filter(|(s, _)| s.iter().any(|&x| x > 0))
        .collect();
    rows.sort_by_key(|(s, _)| (*s).clone());
    let mut tally = [0usize; 3];
    let mut known = 0;
    for (s, reps) in rows {
        let shown: Vec<String> = s.iter().take(9).map(|x| x.to_string()).collect();
        let names: Vec<String> = reps
            .iter()
            .map(|&c| {
                let d = Design::full(c, base);
                match d.parity_code() {
                    Some(p) if base > 2 => format!("{c} (parity {p})"),
                    _ => c.to_string(),
                }
            })
            .collect();
        let said = describe(s);
        tally[usize::from(said.contains("no integer root"))
            + 2 * usize::from(said.starts_with("no"))] += 1;
        let size: usize = reps.iter().map(|r| classes[r].len()).sum();
        let found = oeis(s);
        known += usize::from(found.starts_with("OEIS A"));
        println!(
            "  {} | {} | codes {} | {size} codes | {found}",
            shown.join(" "),
            said,
            names.join(" ")
        );
    }
    println!("  fits of margin at least 2: {} with integer roots only, {} with a factor of no integer root, {} with no fit", tally[0], tally[1], tally[2]);
    println!(
        "  in the OEIS by the first seven terms from the first nonzero one: {known} of {}",
        tally.iter().sum::<usize>()
    );
    let open: Vec<u128> = groups
        .iter()
        .filter(|(s, _)| describe(s).starts_with("no"))
        .map(|(_, r)| r[0])
        .collect();
    if deep > top && !open.is_empty() {
        census_deep(base, deep, &open);
    }
}

fn census_deep(base: usize, top: usize, codes: &[u128]) {
    println!("  the sequences with no fit, rerun to level {top}");
    for &code in codes {
        let d = Design::full(code, base);
        let s = series(base, &d.tile, top);
        let shown: Vec<String> = s.iter().map(|x| x.to_string()).collect();
        println!(
            "    {} | kept {} | {} | {}",
            d.name(),
            d.kept(),
            shown.join(" "),
            describe(&s)
        );
    }
}

pub fn bases() {
    for (base, top) in [(4usize, 12usize), (5, 10)] {
        let mut codes: Vec<u128> = (1..16u128)
            .map(|p| Design::new(p, base, 2).full_code())
            .collect();
        let all = (1u128 << (base * base)) - 1;
        codes.extend((0..base * base).map(|c| all ^ (1 << c)));
        let mut reps: Vec<u128> = codes.iter().map(|&c| klein(c, base)).collect();
        reps.sort();
        reps.dedup();
        let mut distinct: BTreeMap<Vec<u128>, bool> = BTreeMap::new();
        println!("BASE {base}: the fifteen parity codes and every one-cell deletion, {} Klein classes, levels 0..{top}", reps.len());
        for rep in reps {
            let d = Design::full(rep, base);
            let s = series(base, &d.tile, top);
            if s.iter().all(|&x| x == 0) {
                println!("  {} | kept {} | no loop", d.name(), d.kept());
                continue;
            }
            let shown: Vec<String> = s.iter().take(8).map(|x| x.to_string()).collect();
            let found = oeis(&s);
            distinct.insert(s.clone(), found.starts_with("OEIS A"));
            println!(
                "  {} | kept {} | {} | {} | {found}",
                d.name(),
                d.kept(),
                shown.join(" "),
                describe(&s)
            );
        }
        let known = distinct.values().filter(|&&hit| hit).count();
        println!("  {} distinct nonzero loop sequences, in the OEIS by the first seven terms from the first nonzero one: {known} of {}", distinct.len(), distinct.len());
    }
}

pub fn lengths(base: usize, top: usize, codes: &[u128]) {
    for &code in codes {
        let d = Design::full(code, base);
        println!(
            "NEW LOOPS BY LENGTH in block strands, void strands counted, {}",
            d.name()
        );
        let mut kept = crate::count::unit(true);
        let mut longest = 0;
        for level in 0..top {
            kept = crate::count::glue(base, &d.tile, &kept, level + 1 < top);
            let row: Vec<String> = kept
                .lengths
                .iter()
                .map(|(l, c)| format!("{l}:{c}"))
                .collect();
            longest = longest.max(kept.lengths.keys().last().copied().unwrap_or(0));
            println!("  gluing {level} to {}: {}", level + 1, row.join(" "));
        }
        println!(
            "  longest new loop over gluings 0 to {}: {longest} strands",
            top - 1
        );
    }
}

fn fibonacci(n: i64) -> i128 {
    let (mut a, mut b) = (1i128, 0i128);
    for _ in 0..n + 1 {
        (a, b) = (b, a + b);
    }
    a
}

fn gain_law(code: u128, number: usize, rule: usize, n: i64) -> Option<i128> {
    let p = |b: i128, e: i64| b.pow(e as u32);
    let first = |v: i128, f: i128| Some(if n == 0 { v } else { f });
    match (code, number, rule) {
        (7, 2, 2) => first(0, p(2, n) - 2),
        (11, 2, 2) => first(0, if n > 0 { p(2, n - 1) } else { 0 }),
        (9, 2, 2) => Some(1),
        (7, 3, 2) => Some(5 * p(3, n) - 7 * n as i128 - 5),
        (14, 3, 2) => first(1, if n > 0 { 2 * p(3, n - 1) } else { 0 }),
        (9, 3, 2) => Some(2 * p(3, n)),
        (6, 3, 2) => Some(p(2, n + 1)),
        (11, 3, 2) | (13, 3, 2) => Some(p(2, n + 1) - 2),
        (13, 3, 3) => first(0, fibonacci(2 * n - 3)),
        _ => None,
    }
}

fn loop_law(code: u128, number: usize, rule: usize, n: i64) -> Option<i128> {
    let p = |b: i128, e: i64| b.pow(e as u32);
    let first = |f: &dyn Fn() -> i128| Some(if n == 0 { 0 } else { f() });
    match (code, number, rule) {
        (7, 2, 2) => first(&|| p(3, n - 1) - p(2, n) + 1),
        (11, 2, 2) => first(&|| p(3, n - 1) - p(2, n - 1)),
        (9, 2, 2) => Some(p(2, n) - 1),
        (7, 3, 2) => Some((p(8, n) - 1) / 7 - p(3, n) + n as i128 + 1),
        (14, 3, 2) => first(&|| 2 * p(5, n - 1) - p(3, n - 1)),
        (9, 3, 2) => Some(p(5, n) - p(3, n)),
        (6, 3, 2) => Some(p(4, n) - p(2, n)),
        (11, 3, 2) | (13, 3, 2) => Some((p(7, n) - 6 * p(2, n) + 5) / 15),
        _ => None,
    }
}

pub fn gains() {
    println!(
        "NAMED DESIGNS: loops L by level and new loops J = L(level + 1) - kept L(level) per gluing"
    );
    let named = [
        (7u128, 2usize, 2usize),
        (11, 2, 2),
        (9, 2, 2),
        (7, 3, 2),
        (14, 3, 2),
        (9, 3, 2),
        (6, 3, 2),
        (11, 3, 2),
        (13, 3, 2),
        (13, 3, 3),
        (287, 3, 3),
    ];
    for (code, number, rule) in named {
        let d = Design::new(code, number, rule);
        let top = if number == 2 { 20 } else { 14 };
        let s = series(number, &d.tile, top);
        let k = d.kept() as i128;
        let j: Vec<i128> = (0..top)
            .map(|n| s[n + 1] as i128 - k * s[n] as i128)
            .collect();
        let gain: Vec<String> = j.iter().map(|x| x.to_string()).collect();
        let shown: Vec<String> = s.iter().take(12).map(|x| x.to_string()).collect();
        println!(
            "  {} | kept {k}\n    L {}\n    J {}\n    {} | {}",
            d.name(),
            shown.join(" "),
            gain.join(" "),
            describe(&s),
            oeis(&s)
        );
        if gain_law(code, number, rule, 0).is_some() {
            let held = (0..top)
                .filter(|&n| gain_law(code, number, rule, n as i64) == Some(j[n]))
                .count();
            println!(
                "    stated J law holds at {held} of {top} gluings 0 to {}",
                top - 1
            );
        }
        if loop_law(code, number, rule, 0).is_some() {
            let held = (0..=top)
                .filter(|&n| loop_law(code, number, rule, n as i64) == Some(s[n] as i128))
                .count();
            println!(
                "    stated L law holds at {held} of {} levels 0 to {top}",
                top + 1
            );
        }
    }
}

static DUMP: std::sync::OnceLock<Option<String>> = std::sync::OnceLock::new();

pub fn oeis(s: &[u128]) -> String {
    let dump = DUMP.get_or_init(|| {
        std::env::var("OEIS_STRIPPED")
            .ok()
            .and_then(|p| std::fs::read_to_string(p).ok())
    });
    let Some(dump) = dump else {
        return "OEIS not read".to_string();
    };
    let Some(k) = s.iter().position(|&x| x > 0) else {
        return String::new();
    };
    let window: Vec<String> = s[k..].iter().take(7).map(|x| x.to_string()).collect();
    let key = format!(",{},", window.join(","));
    let hits: Vec<&str> = dump
        .lines()
        .filter(|l| l.contains(&key))
        .filter_map(|l| l.split(' ').next())
        .take(3)
        .collect();
    if hits.is_empty() {
        "OEIS absent".to_string()
    } else {
        format!("OEIS {}", hits.join(" "))
    }
}
