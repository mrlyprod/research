mod laws;
mod machine;
mod span;

use laws::{conditions, depths, law, Conditions};
use machine::{
    affine, commutes, full, image, orbits, reps, rule_of, symmetries, word, Machine, Split,
    NEIGHBOURHOODS, SOURCES, SPREADS,
};
use span::{closure, first_block, lone_prime, long_prime, wernsdorf, Rng};
use std::collections::BTreeMap;
use std::env;
use std::time::Instant;

const TRIES: usize = 512;

// LEMMA

fn lemma() {
    for m in 2..=10usize {
        let size = 1u32 << m;
        let (mut nonaffine, mut ones_structure, mut faults) = (0, 0, 0);
        for code in 0..1u32 << (m + 1) {
            let g: Vec<u8> = (0..=m).map(|t| (code >> t & 1) as u8).collect();
            let f = |x: u32| g[x.count_ones() as usize];
            let structures: Vec<u32> = (1..size)
                .filter(|&v| {
                    let c = f(v) ^ f(0);
                    (0..size).all(|x| f(x ^ v) ^ f(x) == c)
                })
                .collect();
            if affine(&g) {
                continue;
            }
            nonaffine += 1;
            let all = size - 1;
            let predicted = (0..=m).all(|t| g[t] ^ g[m - t] == g[0] ^ g[m]);
            if structures.iter().any(|&v| v != all) || structures.contains(&all) != predicted {
                faults += 1;
            }
            ones_structure += predicted as usize;
        }
        println!("  m {m}: {nonaffine} non-affine symmetric functions, {ones_structure} with the all-ones linear structure, {faults} with any other or mispredicted");
    }
}

// TRAP

fn trap() {
    for n in 3..=5usize {
        let cells = n * n;
        for spread in SPREADS {
            let perm = Split::new(n, 0, 0, 0, spread).perm;
            let orbs = orbits(&perm);
            let product: u128 = orbs.iter().map(|o| (1u128 << o.count_ones()) - 1).product();
            let formula = (1u128 << cells) - 1 - product;
            let small: Vec<u32> = orbs.iter().map(|&o| o as u32).collect();
            let brute = (1u32..(1u32 << cells) - 1)
                .filter(|&a| small.iter().any(|&o| a & o == o))
                .count() as u128;
            println!(
                "  n {n} {spread:<5} orbits {:>2}, sizes {:?}: trapped splits {brute} of {}, formula {formula}",
                orbs.len(),
                orbs.iter().map(|o| o.count_ones()).collect::<Vec<_>>(),
                (1u128 << cells) - 2
            );
        }
    }
    for n in [3usize, 4, 5, 6, 8, 12, 16, 32] {
        let row: Vec<String> = SPREADS
            .iter()
            .map(|spread| {
                let perm = Split::new(n, 0, 0, 0, spread).perm;
                let free: f64 = orbits_of(&perm)
                    .iter()
                    .map(|&s| 1.0 - 0.5f64.powi(s as i32))
                    .product();
                format!("{spread} {:.6}", 1.0 - free)
            })
            .collect();
        println!("  n {n}: trapped fraction of all masks: {}", row.join(", "));
    }
    let mut faults = 0;
    let mut checked = 0;
    for spread in SPREADS {
        for hood in NEIGHBOURHOODS {
            for source in &SOURCES[..3] {
                for fill in 1u128..511 {
                    let s = Split::new(3, fill, hood.1, rule_of(source.1), spread);
                    let h = Machine::new(&s).table();
                    let mut rng = Rng(0x7472_6170);
                    let none = vec![false; 512];
                    for o in orbits(&s.perm) {
                        let p = o.trailing_zeros();
                        let inside = o & fill == o;
                        let u = closure(&[h.as_slice()], 1 << p, 9, &none, &mut rng);
                        let within = u
                            .as_ref()
                            .is_some_and(|u| u.gens.iter().all(|&g| g as u128 & !fill == 0));
                        let exact = u.as_ref().is_some_and(|u| {
                            u.dim == o.count_ones() as usize
                                && u.gens.iter().all(|&g| g as u128 & !o == 0)
                        });
                        if inside && !exact || !inside && within {
                            faults += 1;
                        }
                        checked += 1;
                    }
                }
            }
        }
    }
    println!("  3 x 3: {checked} orbit closures, an orbit inside the written side is exactly its own block and no other unit closure stays inside it: {faults} faults");
}

fn orbits_of(perm: &[usize]) -> Vec<usize> {
    let mut seen = vec![false; perm.len()];
    let mut out = Vec::new();
    for s in 0..perm.len() {
        let (mut p, mut len) = (s, 0);
        while !seen[p] {
            seen[p] = true;
            p = perm[p];
            len += 1;
        }
        if len > 0 {
            out.push(len);
        }
    }
    out
}

// CELLS

struct Cell {
    name: String,
    nonlinear: bool,
    c: Conditions,
    units: usize,
    ones: usize,
    group: &'static str,
    law: String,
    block: Vec<u32>,
    prime: usize,
}

fn judge(schedule: &[Split], name: String, nonlinear: bool, verdict_always: bool) -> Cell {
    let bits = schedule[0].cells();
    let tables: Vec<Vec<u32>> = schedule.iter().map(|s| Machine::new(s).table()).collect();
    let gens: Vec<&[u32]> = tables.iter().map(|t| t.as_slice()).collect();
    let c = conditions(schedule);
    let none = vec![false; 1 << bits];
    let mut rng = Rng(0x756e_6974);
    let mut proper = None;
    let mut units = 0;
    for p in 0..bits {
        match closure(&gens, 1 << p, bits, &none, &mut rng) {
            None => units += 1,
            Some(u) => {
                if proper.is_none() {
                    proper = Some(u)
                }
            }
        }
    }
    let all = (1u64 << bits) as u32 - 1;
    let one = closure(&gens, all, bits, &none, &mut rng);
    let ones = one.as_ref().map_or(bits, |u| u.dim);
    let block = match (proper, one) {
        (Some(u), _) => Some(u),
        (None, Some(u)) => Some(u),
        (None, None) if verdict_always || units == bits => first_block(&gens, bits),
        _ => None,
    };
    let (group, law_name, kept, prime) = match block {
        Some(u) => {
            assert!(span::is_block(&gens, &u), "a closure that is not a block");
            let l = if nonlinear {
                law(schedule, &c, &gens)
            } else {
                "affine".to_string()
            };
            ("imprimitive", l, u.perp(), 0)
        }
        None => {
            let w = wernsdorf(&gens, bits, TRIES, 0x7765_726e);
            match w {
                Some(p) => ("alternating", String::new(), Vec::new(), p),
                None => ("primitive", String::new(), Vec::new(), 0),
            }
        }
    };
    Cell {
        name,
        nonlinear,
        c,
        units,
        ones,
        group,
        law: law_name,
        block: kept,
        prime,
    }
}

fn flags(c: &Conditions) -> String {
    let f = |b: bool, s: &str| {
        if b {
            s.to_string()
        } else {
            "-".repeat(s.len())
        }
    };
    format!(
        "{}{}{}{}{}",
        f(c.reads, "R"),
        f(c.spans, "S"),
        f(c.windows, "W"),
        f(c.isolated, "I"),
        f(c.gamma, "G")
    )
}

fn tabulate(cells: &[Cell]) {
    let mut table: BTreeMap<(String, &str), usize> = BTreeMap::new();
    for cell in cells.iter().filter(|c| c.nonlinear) {
        *table.entry((flags(&cell.c), cell.group)).or_default() += 1;
    }
    for ((f, g), k) in &table {
        println!("    {f} {g:<12} {k}");
    }
    let mut laws: BTreeMap<String, usize> = BTreeMap::new();
    for cell in cells
        .iter()
        .filter(|c| c.nonlinear && c.group == "imprimitive")
    {
        *laws
            .entry(cell.law.split(" of ").next().unwrap().to_string())
            .or_default() += 1;
    }
    println!("    imprimitive nonlinear by first law: {laws:?}");
    let theorem = cells.iter().filter(|c| c.c.theorem()).count();
    let theorem_bad = cells
        .iter()
        .filter(|c| c.c.theorem() && c.group == "imprimitive")
        .count();
    let gamma_unit = cells
        .iter()
        .filter(|c| c.c.gamma && c.units < c.c.edges.len())
        .count();
    println!("    theorem hypotheses met {theorem}, of them imprimitive {theorem_bad}; strongly connected second-difference graph with a proper unit closure {gamma_unit}");
    let alt = cells.iter().filter(|c| c.group == "alternating").count();
    let prim = cells.iter().filter(|c| c.group == "primitive").count();
    let mut primes: BTreeMap<usize, usize> = BTreeMap::new();
    for c in cells.iter().filter(|c| c.prime > 0) {
        *primes.entry(c.prime).or_default() += 1;
    }
    println!("    primitive {}, alternating by a prime cycle longer than half the states {alt}, primitive without one {prim}", alt + prim);
    let unitfree: Vec<&Cell> = cells
        .iter()
        .filter(|c| c.nonlinear && c.group == "imprimitive" && c.units == c.c.edges.len())
        .collect();
    println!(
        "    imprimitive with every unit closure full {}, of them all-ones closure full {}",
        unitfree.len(),
        unitfree
            .iter()
            .filter(|c| c.ones == c.c.edges.len())
            .count()
    );
    for c in cells.iter().filter(|c| {
        c.nonlinear
            && c.group == "imprimitive"
            && (c.c.brief() || c.c.theorem() || c.law == "other")
    }) {
        println!(
            "    falsifier {}: {} {}, units full {}, all-ones closure {}, block codim {}, kept parities {:x?}",
            c.name,
            flags(&c.c),
            c.law,
            c.units,
            c.ones,
            c.block.len(),
            c.block
        );
    }
}

// PAIRS

fn catalog(n: usize) -> Vec<u128> {
    if n == 3 {
        (1u128..511).collect()
    } else {
        (1..16)
            .flat_map(|a| (1..16).map(move |b| word(a, b)))
            .filter(|&m| m != 0xffff)
            .collect()
    }
}

fn pairs(n: usize, verbose: bool) -> Vec<Cell> {
    let cells = n * n;
    let all = full(cells);
    let syms = symmetries(n);
    let group: Vec<&Vec<usize>> = syms.iter().collect();
    let splits = reps(&catalog(n), &group, Some(all));
    println!(
        "  n {n}: {} splits up to torus symmetry and complement",
        splits.len()
    );
    let mut out = Vec::new();
    for hood in NEIGHBOURHOODS {
        for source in SOURCES {
            for &fill in &splits {
                let rule = rule_of(source.1);
                let schedule = [
                    Split::new(n, fill, hood.1, rule, "id"),
                    Split::new(n, all & !fill, hood.1, rule, "id"),
                ];
                let name = format!("pair n {n} fill {fill:#06x} {} {}", hood.0, source.0);
                let cell = judge(&schedule, name, source.0 != "odds", true);
                if verbose && cell.nonlinear {
                    println!(
                        "    {} {} {} {} units {} ones {} prime {}",
                        cell.name,
                        flags(&cell.c),
                        cell.group,
                        cell.law,
                        cell.units,
                        cell.ones,
                        cell.prime
                    );
                }
                out.push(cell);
            }
        }
    }
    let nonlinear = out.iter().filter(|c| c.nonlinear).count();
    println!(
        "  n {n}: {} pair cells, {} splits x {} neighbourhoods x {} Sources, {nonlinear} nonlinear",
        out.len(),
        splits.len(),
        NEIGHBOURHOODS.len(),
        SOURCES.len()
    );
    tabulate(&out);
    out
}

// SINGLES

fn singles(only: bool) -> Vec<Cell> {
    let n = 4;
    let all = full(16);
    let syms = symmetries(n);
    let mut out = Vec::new();
    let mut carried = 0;
    for spread in SPREADS {
        let probe = Split::new(n, 0, 0, 0, spread);
        let central: Vec<&Vec<usize>> = syms.iter().filter(|s| commutes(s, &probe.perm)).collect();
        for hood in NEIGHBOURHOODS {
            for source in SOURCES {
                for fill in reps(&catalog(4), &central, None) {
                    let s = Split::new(n, fill, hood.1, rule_of(source.1), spread);
                    let onto = image(&s.perm, fill) == all & !fill;
                    carried += onto as usize;
                    if only && !onto {
                        continue;
                    }
                    let name = format!(
                        "single fill {fill:#06x} {} {} {spread}{}",
                        hood.0,
                        source.0,
                        if onto { " onto" } else { "" }
                    );
                    out.push(judge(&[s], name, source.0 != "odds", false));
                }
            }
        }
    }
    println!(
        "  4 x 4 single rounds: {} cells, {carried} whose spread carries the split onto its complement",
        out.len()
    );
    tabulate(&out);
    let mut by: BTreeMap<String, [usize; 3]> = BTreeMap::new();
    for c in out.iter().filter(|c| c.nonlinear) {
        let spread = c.name.split(' ').nth(5).unwrap().to_string();
        let e = by.entry(spread).or_default();
        e[0] += 1;
        e[1] += (c.units == 16) as usize;
        e[2] += (c.group != "imprimitive") as usize;
    }
    for (s, [k, u, p]) in by {
        println!("    {s:<5} nonlinear {k}, every unit closure full {u}, primitive {p}");
    }
    out
}

// FAMILIES

fn checker(n: usize) -> u128 {
    (0..n * n)
        .filter(|&p| (p / n + p % n) % 2 == 0)
        .fold(0, |m, p| m | 1u128 << p)
}

fn stripes(n: usize) -> u128 {
    (0..n * n)
        .filter(|&p| p % n % 2 == 0)
        .fold(0, |m, p| m | 1u128 << p)
}

fn family() {
    for n in [4usize, 6, 8, 10] {
        for (name, fill) in [("checkerboard", checker(n)), ("stripes", stripes(n))] {
            let all = full(n * n);
            let row: Vec<String> = NEIGHBOURHOODS
                .iter()
                .flat_map(|hood| {
                    SOURCES[..3].iter().map(move |source| {
                        let rule = rule_of(source.1);
                        let schedule = [
                            Split::new(n, fill, hood.1, rule, "id"),
                            Split::new(n, all & !fill, hood.1, rule, "id"),
                        ];
                        let c = conditions(&schedule);
                        format!(
                            "{} {} {}{}",
                            hood.0,
                            source.0,
                            flags(&c),
                            if c.theorem() { " primitive" } else { "" }
                        )
                    })
                })
                .collect();
            println!("  n {n} {name}: {}", row.join("; "));
        }
    }
    let n = 10;
    let fill = checker(n);
    let rule = rule_of(SOURCES[1].1);
    let schedule = [
        Split::new(n, fill, NEIGHBOURHOODS[0].1, rule, "id"),
        Split::new(n, full(n * n) & !fill, NEIGHBOURHOODS[0].1, rule, "id"),
    ];
    let c = conditions(&schedule);
    let d = depths(&c.edges, 0);
    let mut hist: BTreeMap<usize, usize> = BTreeMap::new();
    for &k in &d {
        *hist.entry(k).or_default() += 1;
    }
    println!("  n 10 checkerboard carpet binary: {}, second-difference graph depths from cell 0 {hist:?}", flags(&c));
}

// FIVE

fn lines() -> Vec<(String, u128)> {
    let mut out = Vec::new();
    for (dir, slope) in [("axis", 0usize), ("diagonal", 1), ("knight", 2)] {
        for (label, set) in [
            ("{0}", vec![0usize]),
            ("{0,1}", vec![0, 1]),
            ("{0,2}", vec![0, 2]),
        ] {
            let fill = (0..25)
                .filter(|&p| set.contains(&((p % 5 + 5 * 5 - slope * (p / 5)) % 5)))
                .fold(0u128, |m, p| m | 1u128 << p);
            out.push((format!("{dir} {label}"), fill));
        }
    }
    out
}

fn decompose(step: &dyn Fn(u32) -> u32, bits: usize) -> Vec<usize> {
    let size = 1usize << bits;
    let mut seen = vec![0u64; size / 64];
    let mut out = Vec::new();
    for s in 0..size {
        if seen[s >> 6] >> (s & 63) & 1 == 1 {
            continue;
        }
        let (mut x, mut len) = (s as u32, 0usize);
        while seen[x as usize >> 6] >> (x & 63) & 1 == 0 {
            seen[x as usize >> 6] |= 1 << (x & 63);
            x = step(x);
            len += 1;
        }
        out.push(len);
    }
    out
}

fn five(budget: f64) {
    let start = Instant::now();
    let all = full(25);
    let mut chosen = None;
    for hood in NEIGHBOURHOODS {
        for source in &SOURCES[..3] {
            let row: Vec<String> = lines()
                .into_iter()
                .map(|(name, fill)| {
                    let rule = rule_of(source.1);
                    let schedule = [
                        Split::new(5, fill, hood.1, rule, "id"),
                        Split::new(5, all & !fill, hood.1, rule, "id"),
                    ];
                    let c = conditions(&schedule);
                    if c.theorem() && chosen.is_none() {
                        chosen = Some((format!("{name} {} {}", hood.0, source.0), schedule));
                    }
                    format!("{name} {}", flags(&c))
                })
                .collect();
            println!("  5 x 5 {} {}: {}", hood.0, source.0, row.join("; "));
        }
    }
    let Some((name, schedule)) = chosen else {
        println!("  5 x 5: no line pair meets the hypotheses");
        return;
    };
    let (ma, mb) = (Machine::new(&schedule[0]), Machine::new(&schedule[1]));
    let size = 1usize << 25;
    let mut rng = Rng(0x6669_7665);
    let (mut jordan, mut long) = (None, None);
    let mut tries = 0;
    while start.elapsed().as_secs_f64() < budget && long.is_none() {
        tries += 1;
        let len = 2 + (tries - 1) % 2;
        let stamps: Vec<u32> = (0..len)
            .map(|_| rng.next() as u32 & (size as u32 - 1))
            .collect();
        let step = |x: u32| stamps.iter().fold(x, |y, &c| mb.apply(ma.apply(y)) ^ c);
        let lengths = decompose(&step, 25);
        if jordan.is_none() {
            jordan = lone_prime(&lengths, size).map(|p| (p, tries, len));
        }
        long = long_prime(&lengths, size).map(|p| (p, tries, len));
        println!(
            "    word {tries}: {len} stamped pair rounds, {} cycles, longest {}, lone prime {:?}, long prime {:?}",
            lengths.len(),
            lengths.iter().max().unwrap(),
            lone_prime(&lengths, size),
            long.map(|l| l.0)
        );
    }
    println!("  5 x 5 {name}: hypotheses met, so primitive; Jordan by a lone prime {jordan:?}; a prime cycle longer than half the states {long:?} (prime, word, pair rounds)");
}

fn main() {
    let verb = env::args().nth(1).unwrap_or_else(|| "all".to_string());
    let start = Instant::now();
    if verb == "all" || verb == "lemma" {
        println!("# the linear structures of symmetric functions");
        lemma();
    }
    if verb == "all" || verb == "trap" {
        println!("# law (a): a spread orbit inside the written side");
        trap();
    }
    if verb == "all" || verb == "pairs" || verb == "verbose" {
        println!("# a split and its complement in turn, identity spread");
        let three = pairs(3, verb == "verbose");
        let four = pairs(4, verb == "verbose");
        let nonlinear = three.iter().chain(&four).filter(|c| c.nonlinear).count();
        println!(
            "  both catalogs: {} pair cells, {nonlinear} nonlinear",
            three.len() + four.len()
        );
    }
    if verb == "all" || verb == "singles" || verb == "onto" {
        println!("# one split, every spread, 4 x 4");
        singles(verb == "onto");
    }
    if verb == "all" || verb == "family" {
        println!("# periodic pairs at every even n");
        family();
    }
    if verb == "five" {
        let budget = env::args().nth(2).map_or(150.0, |b| b.parse().unwrap());
        println!("# 5 x 5 line pairs");
        five(budget);
    }
    println!("# {:.1} s", start.elapsed().as_secs_f64());
}
