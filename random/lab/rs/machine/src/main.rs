mod battery;
mod design;
mod gen;

use battery::{corrected, extra, last, screen, window, Corrected, Extra};
use design::{
    bits, census_check, dense_at, obstruction, offsets, parity_agreement, rules_at, sources_at,
    splits, triples, Mess, Rule, Split, Spread, HOODS, SPREADS,
};
use gen::{counter_words, ChaCha, Counter, Gen, Grid, Load, Machine, Rule30, CTR, OUT};
use std::env;
use std::time::Instant;

const TOP: usize = 6;
const SLACK: f64 = 1.5;
const FINAL: usize = 4;

// CELLS

#[derive(Clone)]
struct Cell {
    split: Split,
    hood: (&'static str, u32),
    spread: Spread,
    rule: Rule,
}

impl Cell {
    fn name(&self) -> String {
        format!(
            "{} {} {} {}",
            self.split.name,
            self.hood.0,
            self.spread.name(),
            self.rule.name
        )
    }
    fn build(&self, level: u32, load: Load) -> Machine {
        let mut m = Machine::new(
            self.split.cells(level),
            self.hood.1,
            &self.rule.mess,
            &self.spread,
            load,
        );
        m.place(m.loads()[0].0);
        m
    }
}

fn cost(r: usize, cells: usize, k: usize) -> f64 {
    (r * cells * k) as f64 / OUT as f64
}

fn cap_for(n: usize) -> usize {
    4 * n
}

fn data_dir() -> std::path::PathBuf {
    let here = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let repo = here.ancestors().nth(5).unwrap().to_path_buf();
    repo.join("data")
        .join(here.strip_prefix(&repo).unwrap_or(&here))
}

// REPORT

fn extra_text(e: &Option<Extra>) -> String {
    e.as_ref().map_or(String::new(), |e| {
        format!(
            "; full walsh {:.2}/2^10, rank 512 {:?} chi2 {:.2} min {}, rank 1024 {:?} chi2 {:.2} min {}",
            e.full, e.big[0].0, e.big[0].1, e.big[0].2, e.big[1].0, e.big[1].1, e.big[1].2
        )
    })
}

fn line(tag: &str, c: &Corrected, cells: usize, k: usize) -> String {
    let v = &c.v;
    let six = c.six.map_or("-".into(), |r| r.to_string());
    let moved = if c.blocked.is_empty() {
        String::new()
    } else {
        format!(
            ", moved by {}",
            c.blocked
                .iter()
                .map(|(t, n)| format!("{n} at {t}"))
                .collect::<Vec<_>>()
                .join(", ")
        )
    };
    let cone = v.cone.map_or("-".into(), |c| c.to_string());
    let sac = v.sac_from.map_or("-".into(), |c| c.to_string());
    let stats = format!(
        "monobit |z| {:.2}, sac {:.4}, walsh {:.2}/2^10, cube |z| {:.2}, complexity L-2^15 {:+}, rank chi2 {:.2}",
        v.stats[0],
        v.stats[1],
        v.stats[2] * 1024.0,
        v.stats[3],
        v.stats[4],
        v.stats[5]
    );
    match v.rstar {
        Some(r) => format!(
            "{tag}: R* {r} (six-test {six}{moved}), C {:.0}, cone {cone}, sac from {sac}, last of the six {}; at R*: {stats}{}",
            cost(r, cells, k),
            last(v),
            extra_text(&c.extra)
        ),
        None => format!(
            "{tag}: no R* through round {} (six-test {six}{moved}, cone {cone}, sac from {sac}); {stats}",
            v.horizon
        ),
    }
}

fn speed<G: Gen>(g: &G, rounds: usize) -> f64 {
    let batches = 1024u64;
    let mut o = [0u64; OUT];
    let mut sink = 0u64;
    let t = Instant::now();
    for m in 0..batches {
        let mut s = g.start(&counter_words(64 * m, None));
        for r in 1..=rounds {
            g.step(&mut s, r);
        }
        g.out(&s, &mut o);
        sink ^= o[0];
    }
    let ns = t.elapsed().as_nanos() as f64;
    std::hint::black_box(sink);
    ns / (batches * 64 * 8) as f64
}

fn chacha_speed(c: &ChaCha) -> f64 {
    let count = 1u32 << 16;
    let mut sink = 0u64;
    let t = Instant::now();
    for i in 0..count {
        sink ^= c.value(std::hint::black_box(i));
    }
    let ns = t.elapsed().as_nanos() as f64;
    std::hint::black_box(sink);
    ns / (count as f64 * 8.0)
}

fn symmetries(n: usize, m: &Machine, rounds: usize) -> [usize; 3] {
    let mats: [[i64; 4]; 8] = [
        [1, 0, 0, 1],
        [0, -1, 1, 0],
        [-1, 0, 0, -1],
        [0, 1, -1, 0],
        [0, 1, 1, 0],
        [0, -1, -1, 0],
        [-1, 0, 0, 1],
        [1, 0, 0, -1],
    ];
    let cells = &m.fill;
    let w = n as i64;
    let stamp = |t: usize| -> Vec<bool> {
        let mut s = vec![false; n * n];
        for p in (0..n * n).filter(|&p| cells[p]) {
            s[((p / n + t % n) % n) * n + (p % n + (t / n) % n) % n] = true;
        }
        s
    };
    let stamps: Vec<Vec<bool>> = (1..=rounds).map(stamp).collect();
    let mut counts = [0usize; 3];
    for mat in mats {
        for ti in 0..w {
            for tj in 0..w {
                let g: Vec<usize> = (0..n * n)
                    .map(|p| {
                        let (i, j) = ((p / n) as i64, (p % n) as i64);
                        let a = (mat[0] * i + mat[1] * j + ti).rem_euclid(w);
                        let b = (mat[2] * i + mat[3] * j + tj).rem_euclid(w);
                        (a * w + b) as usize
                    })
                    .collect();
                let keeps = |set: &[bool]| (0..n * n).all(|p| set[p] == set[g[p]]);
                if !keeps(cells) {
                    continue;
                }
                counts[0] += 1;
                let commutes = (0..n * n).all(|q| {
                    let mut a: Vec<usize> =
                        m.src[q][..m.terms].iter().map(|&p| g[p as usize]).collect();
                    let mut b: Vec<usize> =
                        m.src[g[q]][..m.terms].iter().map(|&p| p as usize).collect();
                    a.sort_unstable();
                    b.sort_unstable();
                    a == b
                });
                if !commutes {
                    continue;
                }
                counts[1] += 1;
                if stamps.iter().all(|s| keeps(s)) {
                    counts[2] += 1;
                }
            }
        }
    }
    counts
}

// VERBS

fn tail(x: f64, count: u64) -> f64 {
    let z = x / std::f64::consts::SQRT_2;
    let steps = 20000;
    let h = 8.0 / steps as f64;
    let f = |t: f64| (-t * t).exp();
    let mut s = f(z) + f(z + 8.0);
    for i in 1..steps {
        s += f(z + i as f64 * h) * if i % 2 == 1 { 4.0 } else { 2.0 };
    }
    let q = s * h / 3.0 * 2.0 / std::f64::consts::PI.sqrt();
    -(count as f64 * (-q).ln_1p()).exp_m1()
}

fn profile(g: &[u8]) {
    let k = g.len() - 1;
    let binom = |n: usize, r: usize| -> i64 {
        if r > n {
            return 0;
        }
        (0..r).fold(1i64, |acc, i| acc * (n - i) as i64 / (i + 1) as i64)
    };
    let walsh = (0..=k).map(|v| {
        (0..=k)
            .map(|w| {
                let kr: i64 = (0..=w)
                    .map(|i| if i % 2 == 0 { 1 } else { -1 } * binom(v, i) * binom(k - v, w - i))
                    .sum();
                if g[w] == 1 {
                    -kr
                } else {
                    kr
                }
            })
            .sum::<i64>()
            .abs()
    });
    let nl = (1i64 << (k - 1)) - walsh.max().unwrap() / 2;
    let ceiling = if k % 2 == 0 {
        (1i64 << (k - 1)) - (1i64 << (k / 2 - 1))
    } else {
        (1i64 << (k - 1)) - (1i64 << ((k - 1) / 2))
    };
    let ones: i64 = (0..=k).filter(|&j| g[j] == 1).map(|j| binom(k, j)).sum();
    let sel = design::anf(g);
    let degree = (0..=k.min(8)).rev().find(|&w| sel[w] != 0);
    println!(
        "g {} at k {k}: degree {degree:?}, ones {ones}, nonlinearity {nl} of the symmetric ceiling {ceiling} = {:.4}, parity distance {}",
        bits(g),
        nl as f64 / ceiling as f64,
        design::parity_distance(g)
    );
}

fn sources() {
    for k in [8usize, 4] {
        let named = sources_at(k);
        let (rules, _, dense) = rules_at(k);
        println!(
            "k {k}: {} distinct Source sets, {} rules at parity distance {} or {} of {} counts, {} rules in the union",
            named.len(),
            dense,
            (k + 1) / 2,
            (k + 2) / 2,
            k + 1,
            rules.len()
        );
        for (name, g) in &named {
            println!(
                "  g {} parity distance {} {}",
                bits(g),
                design::parity_distance(g),
                name
            );
        }
        assert_eq!(dense, dense_at(k).len());
    }
    match census_check(8) {
        Some((read, agree)) => {
            println!("the mess census: {read} Sources read, {agree} agree with the rebuild on counts 0..8")
        }
        None => println!("the mess census JSON is absent; the rebuild stands alone"),
    }
    if let Mess::Life { birth, survive } = obstruction(8).mess {
        println!(
            "primes birth, odds survive at k 8: agrees with the XOR of the 8 neighbours on {:.4} of (centre, neighbourhood) pairs",
            parity_agreement(&birth, &survive)
        );
    }
}

fn reference<G: Gen>(name: &str, g: &G) {
    let c = corrected(g, 1, 1);
    println!("{}", line(name, &c, 0, 0));
    let e = c
        .extra
        .clone()
        .unwrap_or_else(|| extra(&window(g, 1, 1)[0]));
    println!(
        "  {name}: 64x64 rank classes 64/63/62/<=61 over 2^14 matrices {:?}{}",
        probe_rank(g),
        extra_text(&Some(e))
    );
}

fn calibrate() {
    for c in [Counter::splitmix(), Counter::xoshiro(battery::MAIN_LOG)] {
        reference(c.name, &c);
    }
    reference("chacha8", &ChaCha::new(8));
    let sm = Counter::splitmix();
    let (mean, sd, max, over) = battery::screen_null(&sm, 1024);
    println!(
        "screen statistic on splitmix64 over 1024 batches: mean {mean:.1}, sd {sd:.1}, max {max:.1}, {over} at or over the bar {:.0}",
        battery::SCREEN_BAR
    );
    let (mut lo, mut hi) = (4.0, 8.0);
    for _ in 0..50 {
        let mid = (lo + hi) / 2.0;
        if tail(mid, 1 << 26) > 0.5 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    println!(
        "full walsh bar sqrt(2 ln 2^27) + 1 = {:.3}/2^10, exceeded by the largest of 2^26 ideal coefficients with probability {:.2e}; that largest has median {lo:.2}/2^10",
        battery::full_bar(),
        tail(battery::full_bar(), 1 << 26)
    );
}

fn probe_rank<G: Gen>(g: &G) -> [u64; 4] {
    let mut h = [0u64; 4];
    let mut o = [0u64; OUT];
    for m in 0..1u64 << (battery::MAIN_LOG - 6) {
        let mut s = g.start(&counter_words(64 * m, None));
        g.step(&mut s, 1);
        g.out(&s, &mut o);
        h[match battery::rank(&o) {
            64 => 0,
            63 => 1,
            62 => 2,
            _ => 3,
        }] += 1;
    }
    h
}

fn chacha() {
    let c = ChaCha::new(8);
    println!(
        "chacha8 key {}, nonce {}",
        c.key.map(|w| format!("{w:08x}")).join(" "),
        c.nonce.map(|w| format!("{w:08x}")).join(" ")
    );
    for ctr in [0u32, 1, 2, 63, 64, 1 << 20, 0x8000_0005, u32::MAX] {
        println!("chacha8 counter {ctr:08x} out {:016x}", c.value(ctr));
    }
    let key: [u32; 8] =
        std::array::from_fn(|i| u32::from_le_bytes(std::array::from_fn(|b| (4 * i + b) as u8)));
    let b = gen::chacha_block(&key, 1, &[0x0900_0000, 0x4a00_0000, 0], 20);
    println!(
        "chacha20 rfc 8439 2.3.2 block: {}",
        b.map(|w| format!("{w:08x}")).join(" ")
    );
    println!("chacha8 scalar: {:.2} ns/byte", chacha_speed(&c));
}

struct Screened {
    cell: Cell,
    r: usize,
    c: f64,
}

fn second(out: Vec<Screened>, level: u32, tag: &str) -> Vec<Screened> {
    let n = 1usize << level;
    let cells = n * n;
    let rmax = cap_for(n);
    let t = Instant::now();
    let mut best2 = f64::INFINITY;
    let mut kept = Vec::new();
    let total = out.len();
    for s in out {
        let k = offsets(s.cell.hood.1).len();
        let cap = if best2.is_finite() {
            ((SLACK * best2 * OUT as f64 / (cells * k) as f64) as usize).min(rmax)
        } else {
            rmax
        };
        if s.r > cap {
            continue;
        }
        let m = s.cell.build(level, Load::Fill);
        if let Some(r) = battery::lite(&m, s.r, cap) {
            let c = cost(r, cells, k);
            best2 = best2.min(c);
            kept.push(Screened { cell: s.cell, r, c });
        }
    }
    kept.retain(|s| s.c <= SLACK * best2);
    kept.sort_by(|a, b| a.c.partial_cmp(&b.c).unwrap().then(a.r.cmp(&b.r)));
    println!(
        "{tag} second screen n {n} (walsh, monobit, rank over 2^14 counters from the screen round): {total} cells in, {} within {SLACK} of the best C {best2:.0}, {:.1} s",
        kept.len(),
        t.elapsed().as_secs_f64()
    );
    for s in kept.iter().take(12) {
        println!("  second R {} C {:.0}: {}", s.r, s.c, s.cell.name());
    }
    kept
}

fn catalog(level: u32) -> Vec<Screened> {
    let n = 1usize << level;
    let cells = n * n;
    let rmax = cap_for(n);
    let mut best = f64::INFINITY;
    let mut out = Vec::new();
    let mut tried = 0usize;
    let mut passed = 0usize;
    let mut blind = 0usize;
    let mut geometries = 0usize;
    let t = Instant::now();
    for split in splits() {
        let mask = split.cells(level);
        if mask.iter().all(|&a| a) || !mask.iter().any(|&a| a) {
            continue;
        }
        for hood in HOODS {
            let k = offsets(hood.1).len();
            let (rules, _, _) = rules_at(k);
            for spread in SPREADS {
                geometries += 1;
                let probe = Cell {
                    split: split.clone(),
                    hood,
                    spread: Spread::Map(spread),
                    rule: rules[0].clone(),
                }
                .build(level, Load::Fill);
                let Some(reach) = probe.structural(rmax) else {
                    blind += 1;
                    tried += rules.len();
                    continue;
                };
                for rule in &rules {
                    let cell = Cell {
                        split: split.clone(),
                        hood,
                        spread: Spread::Map(spread),
                        rule: rule.clone(),
                    };
                    let cap = if best.is_finite() {
                        ((SLACK * best * OUT as f64 / (cells * k) as f64) as usize).min(rmax)
                    } else {
                        rmax
                    };
                    tried += 1;
                    if reach > cap {
                        continue;
                    }
                    let m = cell.build(level, Load::Fill);
                    if let Some(r) = screen(&m, cap) {
                        passed += 1;
                        let c = cost(r, cells, k);
                        best = best.min(c);
                        out.push(Screened { cell, r, c });
                    }
                }
            }
        }
    }
    out.retain(|s| s.c <= SLACK * best);
    out.sort_by(|a, b| a.c.partial_cmp(&b.c).unwrap().then(a.r.cmp(&b.r)));
    println!(
        "screen n {n}: {geometries} geometries, {blind} whose structural cone misses an output pair through round {rmax}; {tried} cells, {passed} passed under their cap, {} within {SLACK} of the best screen C {best:.0}, {:.1} s",
        out.len(),
        t.elapsed().as_secs_f64()
    );
    let mut by = std::collections::BTreeMap::new();
    for s in &out {
        *by.entry((s.cell.hood.0, s.cell.spread.name()))
            .or_insert(0usize) += 1;
    }
    println!("  within the slack by hood and spread: {by:?}");
    for s in out.iter().take(6) {
        println!("  screen R {} C {:.0}: {}", s.r, s.c, s.cell.name());
    }
    second(out, level, "catalog")
}

fn full(cell: &Cell, level: u32, guess: Option<usize>) -> (Machine, Corrected) {
    let n = 1usize << level;
    let m = cell.build(level, Load::Fill);
    let g = guess
        .or_else(|| screen(&m, cap_for(n)))
        .unwrap_or(cap_for(n));
    let v = corrected(&m, g, cap_for(n));
    (m, v)
}

fn rule30(n: usize, near: bool) -> Corrected {
    let r = Rule30::new(n * n, near);
    let cap = 4 * n * n;
    let g = screen(&r, cap).unwrap_or(cap);
    let v = corrected(&r, g, cap);
    let tag = if near { " near load" } else { "" };
    println!(
        "{}",
        line(&format!("rule 30 ring {}{tag}", n * n), &v, n * n, 3)
    );
    v
}

fn r30slide(n: usize, at: (usize, usize)) -> (Rule30, Corrected) {
    let r = Rule30::slide(n, at);
    let cap = n * n;
    let g = screen(&r, cap).unwrap_or(cap);
    let v = corrected(&r, g, cap);
    (r, v)
}

fn shortlist(screened: &[Screened]) -> Vec<&Screened> {
    let mut list: Vec<&Screened> = screened.iter().take(TOP).collect();
    for s in screened {
        let geo = |c: &Cell| (c.split.name.clone(), c.hood.0, c.spread.name());
        if !list.iter().any(|t| geo(&t.cell) == geo(&s.cell)) {
            list.push(s);
        }
    }
    list
}

fn positions(m: &Machine) -> [(&'static str, (usize, usize)); 3] {
    let n = m.n;
    [
        ("antipode", m.loads()[0].0),
        ("upstream", (n - 4, n - 8)),
        ("below", (n / 2 + 2, 0)),
    ]
}

fn named(letter: u32, hood: &str, spread: Spread, g: &str) -> Cell {
    Cell {
        split: Split {
            name: format!("{letter}^L"),
            letter,
            anti: false,
        },
        hood: HOODS.into_iter().find(|h| h.0 == hood).unwrap(),
        spread,
        rule: Rule {
            name: format!("g={g}"),
            mess: Mess::Feistel(g.bytes().map(|b| b - b'0').collect()),
        },
    }
}

fn winner() -> Cell {
    named(14, "carpet", Spread::Map("slide"), "010100101")
}

fn branch() -> Cell {
    let text = std::fs::read_to_string(data_dir().join("branch.txt")).unwrap_or_default();
    let v: Vec<i64> = text
        .split_whitespace()
        .filter_map(|w| w.parse().ok())
        .collect();
    if v.len() != 6 {
        panic!(
            "run the family verb first: no branch.txt in {:?}",
            data_dir()
        );
    }
    Cell {
        spread: Spread::Sum([(v[0], v[1]), (v[2], v[3]), (v[4], v[5])]),
        ..winner()
    }
}

fn best() {
    let t = Instant::now();
    let screened = catalog(4);
    let mut fulls: Vec<(Cell, Corrected, f64)> = Vec::new();
    for s in shortlist(&screened) {
        let (m, v) = full(&s.cell, 4, Some(s.r));
        println!(
            "{}",
            line(
                &format!("n 16 {} (screen R {})", s.cell.name(), s.r),
                &v,
                256,
                m.k
            )
        );
        let c = v.v.rstar.map_or(f64::INFINITY, |r| cost(r, 256, m.k));
        fulls.push((s.cell.clone(), v, c));
    }
    fulls.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
    for (cell, v, c) in &fulls {
        println!("  ranked: C {c:.0}, R* {:?}: {}", v.v.rstar, cell.name());
    }
    println!("best: {:.1} s", t.elapsed().as_secs_f64());
}

fn bad(cell: &Cell) {
    for level in [4u32, 5] {
        let n = 1usize << level;
        let base = cell.build(level, Load::Fill);
        let sym = symmetries(n, &base, 4 * n);
        println!(
            "n {n} {}: |A| {}, written cells {}, B reads per round {}; {} torus symmetries fix the split, {} also commute with the spread, {} also fix every stamp through round {}",
            cell.name(),
            base.fill.iter().filter(|&&a| a).count(),
            base.write_count(),
            base.read_count(),
            sym[0],
            sym[1],
            sym[2],
            4 * n
        );
    }
    let n = 16usize;
    let base = cell.build(4, Load::Fill);
    let guess = screen(&base, cap_for(n)).unwrap_or(cap_for(n));
    let mut zero = cell.build(4, Load::Zero);
    zero.place(base.loads()[0].0);
    let v = corrected(&zero, guess, cap_for(n));
    println!("{}", line("  bad load all-zero fill", &v, n * n, zero.k));
    let mut inside = cell.build(4, Load::Fill);
    inside.place((2, 0));
    let v = corrected(&inside, 1, cap_for(n));
    println!(
        "{}",
        line(
            "  bad load counter inside the output block",
            &v,
            n * n,
            inside.k
        )
    );
    for off in [(n / 2, 0), (0, n / 2)] {
        let shifted = |p: usize| ((p / n + off.0) % n) * n + (p % n + off.1) % n;
        let symmetric = (0..n * n).all(|p| base.fill[p] == base.fill[shifted(p)]);
        let partner: [u32; OUT] = std::array::from_fn(|j| shifted(base.out[j] as usize) as u32);
        let touches = partner.iter().any(|p| base.ctr.contains(p));
        if symmetric || touches {
            println!("  fold partner offset {off:?}: skipped, symmetric {symmetric}, touches the counter {touches}");
            continue;
        }
        let mut m = cell.build(4, Load::Fill);
        m.fold = Some(partner);
        let v = corrected(&m, guess.saturating_sub(4).max(1), cap_for(n));
        println!(
            "{}",
            line(
                &format!("  fold width 2, partner offset {off:?}"),
                &v,
                n * n,
                m.k
            )
        );
        break;
    }
}

fn controls(best: &Cell) {
    let k = offsets(best.hood.1).len();
    let mut list = vec![obstruction(k)];
    let (rules, _, _) = rules_at(k);
    list.extend(rules.into_iter().filter(|r| r.name.contains("odds")));
    for rule in list {
        let cell = Cell {
            rule: rule.clone(),
            ..best.clone()
        };
        let m = cell.build(4, Load::Fill);
        let r = screen(&m, cap_for(16));
        println!(
            "control n 16 {}: screen R {}",
            cell.name(),
            r.map_or(format!("none through {}", cap_for(16)), |r| r.to_string())
        );
    }
}

// HONEST

struct Row {
    tag: String,
    rs: Vec<Option<usize>>,
    k: usize,
    s: usize,
    ns: f64,
}

fn spread_terms(cell: &Cell) -> usize {
    match cell.spread {
        Spread::Map(_) => 1,
        Spread::Sum(_) => 3,
    }
}

fn machine_row(cell: &Cell, level: u32) -> Row {
    let n = 1usize << level;
    let base = cell.build(level, Load::Fill);
    let guess = screen(&base, cap_for(n)).unwrap_or(cap_for(n));
    let mut rs = Vec::new();
    for (tag, at) in positions(&base) {
        let mut m = cell.build(level, Load::Fill);
        m.place(at);
        let c = corrected(&m, guess, cap_for(n));
        println!(
            "{}",
            line(
                &format!(
                    "n {n} {} load {tag} at {at:?}, {} of 32 in A",
                    cell.name(),
                    m.overlap(at)
                ),
                &c,
                n * n,
                m.k
            )
        );
        rs.push(c.v.rstar);
    }
    let ns = rs[0].map_or(f64::NAN, |r| speed(&base, r));
    Row {
        tag: cell.name(),
        rs,
        k: base.k,
        s: spread_terms(cell),
        ns,
    }
}

fn slide_cache(n: usize) -> Option<Vec<Option<usize>>> {
    let text = std::fs::read_to_string(data_dir().join("slide30.txt")).ok()?;
    let rs: Vec<Option<usize>> = ["antipode", "upstream", "below"]
        .iter()
        .map(|tag| {
            text.lines()
                .map(|l| l.split_whitespace().collect::<Vec<_>>())
                .find(|w| w.len() == 3 && w[0] == n.to_string() && w[1] == *tag)
                .and_then(|w| w[2].parse().ok())
        })
        .collect();
    rs.iter().all(|r| r.is_some()).then_some(rs)
}

fn slide_row(level: u32) -> Row {
    let n = 1usize << level;
    let probe = winner().build(level, Load::Fill);
    let at = positions(&probe);
    let rs = match slide_cache(n) {
        Some(rs) => {
            println!("n {n} rule 30 + slide T* at the loads antipode, upstream, below, from the slide30 cache: {rs:?}");
            rs
        }
        None => {
            let mut rs = Vec::new();
            let mut text =
                std::fs::read_to_string(data_dir().join("slide30.txt")).unwrap_or_default();
            for (tag, at) in at {
                let (_, c) = r30slide(n, at);
                println!(
                    "{}",
                    line(
                        &format!("n {n} rule 30 + slide load {tag} at {at:?}"),
                        &c,
                        n * n,
                        3
                    )
                );
                if let Some(t) = c.v.rstar {
                    text.push_str(&format!("{n} {tag} {t}\n"));
                }
                rs.push(c.v.rstar);
            }
            std::fs::create_dir_all(data_dir()).unwrap();
            std::fs::write(data_dir().join("slide30.txt"), text).unwrap();
            rs
        }
    };
    let ns = rs[0].map_or(f64::NAN, |t| speed(&Rule30::slide(n, at[0].1), t));
    Row {
        tag: "rule 30 + slide".into(),
        rs,
        k: 3,
        s: 1,
        ns,
    }
}

fn c_of(row: &Row, i: usize, n: usize) -> f64 {
    row.rs[i].map_or(f64::NAN, |r| cost(r, n * n, row.k))
}

fn c_all(row: &Row, i: usize, n: usize) -> f64 {
    row.rs[i].map_or(f64::NAN, |r| cost(r, n * n, row.k + row.s))
}

fn headline() {
    let t = Instant::now();
    let cha = chacha_speed(&ChaCha::new(8));
    let cells = [winner(), branch()];
    let mut a = Vec::new();
    let mut b = Vec::new();
    let mut c = Vec::new();
    for level in [4u32, 5] {
        let n = 1usize << level;
        let r30 = slide_row(level);
        let ring = rule30(n, false);
        let near = rule30(n, true);
        let rows: Vec<Row> = cells.iter().map(|cell| machine_row(cell, level)).collect();
        for row in &rows {
            let ratios: Vec<String> = (0..3)
                .map(|i| format!("{:.2}", c_of(&r30, i, n) / c_of(row, i, n)))
                .collect();
            let strict: Vec<String> = (0..3)
                .map(|i| format!("{:.2}", c_all(&r30, i, n) / c_all(row, i, n)))
                .collect();
            println!(
                "n {n} {}: R* {:?} at the loads antipode, upstream, below; rule 30 + slide T* {:?}; C(rule 30 + slide)/C(machine) {} by rule reads, {} counting the spread reads",
                row.tag,
                row.rs,
                r30.rs,
                ratios.join(" "),
                strict.join(" ")
            );
        }
        a.push(format!(
            "n {n}: {} rule 30 + slide {:.2}",
            rows.iter()
                .map(|r| format!("{} {:.2},", r.tag, r.ns))
                .collect::<Vec<_>>()
                .join(" "),
            r30.ns
        ));
        b.push(format!(
            "n {n}: {}",
            rows.iter()
                .map(|r| format!(
                    "{} {:.2} (loads {:.2} to {:.2})",
                    r.tag,
                    c_of(&r30, 0, n) / c_of(r, 0, n),
                    (0..3)
                        .map(|i| c_of(&r30, i, n) / c_of(r, i, n))
                        .fold(f64::INFINITY, f64::min),
                    (0..3)
                        .map(|i| c_of(&r30, i, n) / c_of(r, i, n))
                        .fold(0.0, f64::max)
                ))
                .collect::<Vec<_>>()
                .join(", ")
        ));
        let cr = |v: &Corrected| v.v.rstar.map_or(f64::NAN, |t| cost(t, n * n, 3));
        c.push(format!(
            "n {n}: ring antipodal {:.2}, ring near {:.2}",
            cr(&ring) / c_of(&rows[0], 0, n),
            cr(&near) / c_of(&rows[0], 0, n)
        ));
    }
    println!("headline: {:.1} s", t.elapsed().as_secs_f64());
    println!(
        "(a) ns/byte at R*, identical harness, 32-bit counter in, 64-bit block out: {}; chacha8 scalar {cha:.2}",
        a.join("; ")
    );
    println!(
        "(b) C(rule 30 + slide)/C(machine), antipodal load: {}",
        b.join("; ")
    );
    println!(
        "(c) cone advantage, dimensional: C(rule 30 untransported)/C(slide machine): {}",
        c.join("; ")
    );
}

// SPEED

const REPS: usize = 7;
const RSTAR: [(u32, [usize; 3]); 2] = [(4, [24, 15, 68]), (5, [45, 22, 141])];

fn spread_of(mut v: Vec<f64>) -> String {
    v.sort_by(|a, b| a.total_cmp(b));
    format!(
        "{:.1} ({:.1} to {:.1})",
        v[v.len() / 2],
        v[0],
        v[v.len() - 1]
    )
}

fn timed(f: impl Fn() -> f64) -> String {
    spread_of((0..REPS).map(|_| f()).collect())
}

fn speeds() {
    let t = Instant::now();
    for (level, [slide, branch3, rule30]) in RSTAR {
        let n = 1usize << level;
        let sm = winner().build(level, Load::Fill);
        let bm = branch().build(level, Load::Fill);
        let r30 = Rule30::slide(n, positions(&sm)[0].1);
        println!(
            "side {n} slide machine R* {slide}: {} ns/byte",
            timed(|| speed(&sm, slide))
        );
        println!(
            "side {n} branch-3 machine R* {branch3}: {} ns/byte",
            timed(|| speed(&bm, branch3))
        );
        println!(
            "side {n} rule 30 + slide R* {rule30}: {} ns/byte",
            timed(|| speed(&r30, rule30))
        );
    }
    let cha = ChaCha::new(8);
    println!("chacha8 scalar: {} ns/byte", timed(|| chacha_speed(&cha)));
    println!(
        "speed: median and range of {REPS} runs, {:.1} s",
        t.elapsed().as_secs_f64()
    );
}

// DIAGNOSIS

fn active<G: Gen<State = Grid>>(g: &G, rounds: usize, batches: u64) -> Vec<f64> {
    let mut sum = vec![0f64; rounds + 1];
    for m in 0..batches {
        let mut base = g.start(&counter_words(64 * m, None));
        let mut flips: Vec<Grid> = (0..CTR)
            .map(|b| g.start(&counter_words(64 * m, Some(b))))
            .collect();
        for (r, total) in sum.iter_mut().enumerate() {
            if r > 0 {
                g.step(&mut base, r);
                for f in flips.iter_mut() {
                    g.step(f, r);
                }
            }
            for f in &flips {
                *total += base
                    .cells()
                    .iter()
                    .zip(f.cells())
                    .map(|(a, b)| (a ^ b).count_ones() as f64)
                    .sum::<f64>();
            }
        }
    }
    sum.iter()
        .map(|s| s / (batches as f64 * 64.0 * CTR as f64))
        .collect()
}

fn branching(m: &Machine, g: &[u8]) -> f64 {
    let flip = |reads: usize| -> f64 {
        if reads == 0 {
            return 0.0;
        }
        let mut p = 0.0;
        let mut binom = 1.0;
        for c in 0..reads {
            if g[c] != g[c + 1] {
                p += binom;
            }
            binom = binom * (reads - 1 - c) as f64 / (c + 1) as f64;
        }
        p / 2f64.powi(reads as i32 - 1)
    };
    let total: f64 = m
        .written()
        .iter()
        .map(|&(_, reads)| reads as f64 * flip(reads))
        .sum();
    total / (m.n * m.n) as f64
}

fn diagnose() {
    let n = 32usize;
    let rounds = 24;
    let gs = [winner(), branch()];
    let mut lines = Vec::new();
    for cell in &gs {
        let m = cell.build(5, Load::Fill);
        let Mess::Feistel(g) = &cell.rule.mess else {
            continue;
        };
        let beta = branching(&m, g);
        let a = active(&m, rounds, 16);
        lines.push((cell.name(), a));
        println!(
            "n 32 {}: written cells {}, reads {}, expected new active cells per active cell per Feistel step {beta:.3}, predicted sparse multiplier {:.3}",
            cell.name(),
            m.write_count(),
            m.read_count(),
            m.terms as f64 * (1.0 + beta)
        );
    }
    let r = Rule30::slide(n, winner().build(5, Load::Fill).loads()[0].0);
    lines.push(("rule 30 + slide".into(), active(&r, rounds, 16)));
    for (name, a) in &lines {
        println!(
            "n 32 active cells per round, {name}: {}",
            a.iter()
                .map(|v| format!("{v:.1}"))
                .collect::<Vec<_>>()
                .join(" ")
        );
        println!(
            "  multiplier a(r)/a(r-1): {}",
            a.windows(2)
                .map(|w| format!("{:.2}", w[1] / w[0]))
                .collect::<Vec<_>>()
                .join(" ")
        );
    }
}

// BRANCH 3

fn translation_symmetric(fill: &[bool], n: usize, d: (i64, i64)) -> bool {
    let w = n as i64;
    (0..n * n).all(|p| {
        let (i, j) = ((p / n) as i64, (p % n) as i64);
        let q = ((i + d.0).rem_euclid(w) * w + (j + d.1).rem_euclid(w)) as usize;
        fill[p] == fill[q]
    })
}

fn spread_rank(m: &Machine) -> usize {
    let cells = m.n * m.n;
    let w = cells / 64;
    let mut rows = vec![0u64; cells * w];
    for (q, s) in m.src.iter().enumerate() {
        for &p in &s[..m.terms] {
            rows[q * w + p as usize / 64] ^= 1 << (p % 64);
        }
    }
    battery::rank_wide(&mut rows, w)
}

fn bijective(m: &Machine, rounds: usize) -> (bool, bool) {
    let cells = m.n * m.n;
    let mut st = 0x6269_6a65_6374u64;
    let x0: Vec<u64> = (0..cells).map(|_| gen::splitmix(&mut st)).collect();
    let mut g = m.grid(x0.clone());
    for _ in 0..m.n {
        m.spread(&mut g);
    }
    let order = g.cells() == x0.as_slice();
    let mut g = m.grid(x0.clone());
    for r in 1..=rounds {
        m.step(&mut g, r);
    }
    for r in (1..=rounds).rev() {
        m.undo(&mut g, r);
    }
    (order, g.cells() == x0.as_slice())
}

fn impulse(m: &Machine) -> Vec<usize> {
    let cells = m.n * m.n;
    let mut x = vec![0u64; cells];
    x[0] = 1;
    let mut g = m.grid(x);
    (1..=m.n)
        .map(|_| {
            m.spread(&mut g);
            g.cells().iter().filter(|&&v| v != 0).count()
        })
        .collect()
}

fn family_screen(level: u32, radius: i64) -> Vec<Screened> {
    let n = 1usize << level;
    let cells = n * n;
    let rmax = cap_for(n);
    let base = winner();
    let fill = base.split.cells(level);
    let k = offsets(base.hood.1).len();
    let t = Instant::now();
    let (mut excluded, mut blind, mut tried) = (0, 0, 0);
    let mut best = f64::INFINITY;
    let mut out = Vec::new();
    for v in triples(radius) {
        tried += 1;
        let pairs = [(0, 1), (0, 2), (1, 2)];
        if pairs
            .iter()
            .any(|&(x, y)| translation_symmetric(&fill, n, (v[x].0 - v[y].0, v[x].1 - v[y].1)))
        {
            excluded += 1;
            continue;
        }
        let cell = Cell {
            spread: Spread::Sum(v),
            ..base.clone()
        };
        let m = cell.build(level, Load::Fill);
        let Some(reach) = m.structural(rmax) else {
            blind += 1;
            continue;
        };
        let cap = if best.is_finite() {
            ((SLACK * best * OUT as f64 / (cells * k) as f64) as usize).min(rmax)
        } else {
            rmax
        };
        if reach > cap {
            continue;
        }
        if let Some(r) = screen(&m, cap) {
            let c = cost(r, cells, k);
            best = best.min(c);
            out.push(Screened { cell, r, c });
        }
    }
    out.retain(|s| s.c <= SLACK * best);
    out.sort_by(|a, b| a.c.partial_cmp(&b.c).unwrap().then(a.r.cmp(&b.r)));
    println!(
        "family n {n}, offsets in the box of radius {radius}: {tried} triples, {excluded} with an offset difference that is a translation symmetry of the split, {blind} structurally blind through round {rmax}, {} within {SLACK} of the best screen C {best:.0}, {:.1} s",
        out.len(),
        t.elapsed().as_secs_f64()
    );
    for s in out.iter().take(6) {
        println!("  screen R {} C {:.0}: {}", s.r, s.c, s.cell.name());
    }
    second(out, level, "family")
}

fn family(radius: i64) {
    let t = Instant::now();
    let r30: Vec<f64> = [16usize, 32]
        .iter()
        .map(|&n| {
            slide_cache(n)
                .and_then(|rs| rs[0])
                .map_or(f64::NAN, |t| cost(t, n * n, 3))
        })
        .collect();
    let lite: Vec<Vec<Screened>> = [4u32, 5]
        .iter()
        .map(|&l| family_screen(l, radius))
        .collect();
    let mut both: Vec<(Cell, f64, f64)> = lite[0]
        .iter()
        .filter_map(|a| {
            lite[1]
                .iter()
                .find(|b| b.cell.spread == a.cell.spread)
                .map(|b| (a.cell.clone(), a.c, b.c))
        })
        .collect();
    let lite_worse = |x: &(Cell, f64, f64)| (r30[0] / x.1).min(r30[1] / x.2);
    both.sort_by(|x, y| lite_worse(y).partial_cmp(&lite_worse(x)).unwrap());
    println!(
        "family: {} triples pass the second screen at both sides; the {FINAL} best by the worse side of C(rule 30 + slide)/C on the second screen go to the full battery",
        both.len()
    );
    let found: Vec<Cell> = both.iter().take(FINAL).map(|x| x.0.clone()).collect();
    let mut scored = Vec::new();
    for cell in &found {
        let mut cs = Vec::new();
        for level in [4u32, 5] {
            let n = 1usize << level;
            let (m, v) = full(cell, level, None);
            println!(
                "{}",
                line(&format!("n {n} {}", cell.name()), &v, n * n, m.k)
            );
            cs.push(v.v.rstar.map_or(f64::INFINITY, |r| cost(r, n * n, m.k)));
        }
        scored.push((cell.clone(), cs));
    }
    let worse = |c: &Vec<f64>| (r30[0] / c[0]).min(r30[1] / c[1]);
    scored.sort_by(|a, b| worse(&b.1).partial_cmp(&worse(&a.1)).unwrap());
    for (cell, cs) in &scored {
        println!(
            "  {}: C {:.0} at n 16 and {:.0} at n 32, C(rule 30 + slide)/C {:.2} and {:.2}",
            cell.name(),
            cs[0],
            cs[1],
            r30[0] / cs[0],
            r30[1] / cs[1]
        );
    }
    let (best, _) = &scored[0];
    println!("family best by the worse side: {}", best.name());
    if let Spread::Sum(v) = best.spread {
        let dir = data_dir();
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(
            dir.join("branch.txt"),
            v.iter()
                .map(|(a, b)| format!("{a} {b}"))
                .collect::<Vec<_>>()
                .join(" "),
        )
        .unwrap();
    }
    for level in [4u32, 5] {
        let n = 1usize << level;
        let m = best.build(level, Load::Fill);
        let (order, inverse) = bijective(&m, 8);
        let sym = symmetries(n, &m, 4 * n);
        println!(
            "n {n} {}: spread matrix rank {} of {}, L^n = I {order}, R^-1 R = I over 8 rounds on 64 states {inverse}; {} torus symmetries fix the split, {} also commute with the spread, {} also fix every stamp through round {}; impulse |L^t e| for t = 1..n: {:?}",
            best.name(),
            spread_rank(&m),
            n * n,
            sym[0],
            sym[1],
            sym[2],
            4 * n,
            impulse(&m)
        );
    }
    println!("family: {:.1} s", t.elapsed().as_secs_f64());
}

fn all() {
    let t = Instant::now();
    sources();
    calibrate();
    chacha();
    best();
    bad(&winner());
    controls(&winner());
    family(1);
    diagnose();
    headline();
    println!("all: {:.1} s", t.elapsed().as_secs_f64());
}

fn main() {
    let verb = env::args().nth(1).unwrap_or_else(|| "all".into());
    match verb.as_str() {
        "sources" => sources(),
        "tail" => {
            for word in env::args().skip(2) {
                let x: f64 = word.parse().unwrap();
                println!(
                    "the maximum of 2^26 independent |N(0, 1)| exceeds {x} with probability {:.4}",
                    tail(x, 1u64 << 26)
                );
            }
        }
        "mess" => {
            for word in env::args().skip(2) {
                profile(&word.bytes().map(|b| b - b'0').collect::<Vec<u8>>());
            }
        }
        "calibrate" => calibrate(),
        "chacha" => chacha(),
        "screen" => {
            catalog(4);
        }
        "best" => best(),
        "bad" => bad(&winner()),
        "controls" => controls(&winner()),
        "rule30" => {
            for near in [false, true] {
                rule30(16, near);
                rule30(32, near);
            }
        }
        "slide30" => {
            slide_row(4);
            slide_row(5);
        }
        "fscreen" => {
            let a: Vec<i64> = env::args().skip(2).filter_map(|w| w.parse().ok()).collect();
            family_screen(a[0] as u32, a[1]);
        }
        "family" => family(
            env::args()
                .nth(2)
                .and_then(|w| w.parse().ok())
                .unwrap_or(1),
        ),
        "diagnose" => diagnose(),
        "headline" => headline(),
        "speed" => speeds(),
        "all" => all(),
        "trace" => {
            let a: Vec<String> = env::args().skip(2).collect();
            let letter: u32 = a[0].parse().unwrap();
            let anti = a[1] == "anti";
            let hood = HOODS.into_iter().find(|h| h.0 == a[2]).unwrap();
            let spread = SPREADS.into_iter().find(|s| *s == a[3]).unwrap();
            let g: Vec<u8> = a[4].bytes().map(|b| b - b'0').collect();
            let level: u32 = a[5].parse().unwrap();
            let cell = Cell {
                split: Split {
                    name: format!("{}{letter}^L", if anti { "anti " } else { "" }),
                    letter,
                    anti,
                },
                hood,
                spread: Spread::Map(spread),
                rule: Rule {
                    name: format!("g={}", a[4]),
                    mess: Mess::Feistel(g),
                },
            };
            let m = cell.build(level, Load::Fill);
            println!(
                "structural reach {:?}, |A| {}, reads {}",
                m.structural(4 << level),
                m.write_count(),
                m.read_count()
            );
            battery::trace(&m, 4 * (1 << level));
        }
        other => eprintln!("unknown verb {other}; verbs: sources tail mess calibrate chacha screen best bad controls rule30 slide30 family diagnose headline speed all trace"),
    }
}
