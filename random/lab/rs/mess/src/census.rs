use crate::mess::{
    anf, balanced, ceiling, parseval, profile, symmetric_table, table, Profile, Tables,
};
use crate::sources::{designs, prefix, variant, Named, KMAX};
use std::collections::BTreeMap;

const PARSEVAL: usize = 40;

pub struct Entry {
    pub name: String,
    pub oeis: Option<String>,
    pub g: Vec<u8>,
    pub c: Vec<u8>,
    pub rows: Vec<Profile>,
}

pub fn build(sources: &[Named]) -> Vec<Entry> {
    let tables: Vec<Tables> = (1..=KMAX).map(Tables::new).collect();
    sources
        .iter()
        .map(|named| {
            let rows = (1..=KMAX)
                .map(|k| {
                    let g = prefix(named, k);
                    assert!(
                        k > PARSEVAL || parseval(&tables[k - 1], &g),
                        "Parseval fails for {} at {k}",
                        named.name
                    );
                    profile(&tables[k - 1], &g)
                })
                .collect();
            Entry {
                name: named.name.clone(),
                oeis: named.oeis.clone(),
                g: named.g.clone(),
                c: anf(&named.g),
                rows,
            }
        })
        .collect()
}

pub fn reps(entries: &[Entry]) -> Vec<&Entry> {
    let mut seen: Vec<&[u8]> = Vec::new();
    let mut out = Vec::new();
    for entry in entries {
        if !seen.contains(&&entry.g[..]) {
            seen.push(&entry.g[..]);
            out.push(entry);
        }
    }
    out
}

pub fn affine(entry: &Entry) -> bool {
    entry.c[2..].iter().all(|&bit| bit == 0)
}

pub fn brute(entries: &[Entry]) -> usize {
    let mut cells = 0;
    for entry in reps(entries) {
        for k in 1..=14 {
            let t = table(&symmetric_table(&entry.g[..=k], k), k);
            let p = &entry.rows[k - 1];
            assert_eq!(t.degree, p.degree, "{} degree at {k}", entry.name);
            assert_eq!(t.ones, p.ones, "{} ones at {k}", entry.name);
            assert_eq!(
                t.nonlinearity, p.nonlinearity,
                "{} nonlinearity at {k}",
                entry.name
            );
            assert_eq!(
                t.autocorrelation_max, p.autocorrelation_max,
                "{} autocorrelation at {k}",
                entry.name
            );
            assert_eq!(
                t.differential, p.differential,
                "{} differential at {k}",
                entry.name
            );
            cells += 5;
        }
    }
    cells
}

pub fn classes(entries: &[Entry]) -> Vec<Vec<&str>> {
    let mut by: BTreeMap<&[u8], Vec<&str>> = BTreeMap::new();
    for entry in entries {
        by.entry(&entry.g[..]).or_default().push(&entry.name);
    }
    by.into_values().collect()
}

pub fn terms(entry: &Entry, limit: usize) -> String {
    let list: Vec<String> = (0..=limit.min(KMAX))
        .filter(|&j| entry.g[j] == 1)
        .map(|j| j.to_string())
        .collect();
    list.join(",")
}

pub fn report_classes(entries: &[Entry]) {
    println!("SOURCES");
    for class in classes(entries) {
        let entry = entries.iter().find(|e| e.name == class[0]).unwrap();
        println!("class {} : {}", class.join(" = "), terms(entry, 40));
    }
}

pub fn ratio(p: &Profile) -> f64 {
    p.nonlinearity as f64 / ceiling(p.k) as f64
}

pub fn report_designs(entries: &[Entry]) {
    println!("DESIGN FILLS");
    println!("name k degree nl nl/ceil balanced diff/2^k");
    for k in designs().into_iter().filter(|&k| k >= 2) {
        println!("k={k} ceiling={}", ceiling(k));
        for entry in entries {
            let p = &entry.rows[k - 1];
            println!(
                "  {:14} deg={:2} nl={} r={:.4} bal={} diff={:.4} ci={}",
                entry.name,
                p.degree,
                p.nonlinearity,
                ratio(p),
                balanced(p) as u8,
                p.differential as f64 / (1u128 << k) as f64,
                p.immunity
            );
        }
    }
}

pub fn report_degrees(entries: &[Entry]) {
    println!("DEGREE LAWS");
    for entry in entries {
        let line: Vec<String> = entry.rows.iter().map(|p| p.degree.to_string()).collect();
        println!("{} deg(k) k=1..{KMAX}: {}", entry.name, line.join(","));
    }
    for entry in entries {
        let line: Vec<String> = (0..=127.min(KMAX))
            .map(|w| entry.c[w].to_string())
            .collect();
        println!(
            "{} c(w) w=0..{}: {}",
            entry.name,
            line.len() - 1,
            line.join("")
        );
    }
}

pub fn report_trend(entries: &[Entry]) {
    println!("TREND nl/ceiling");
    for entry in entries {
        let line: Vec<String> = [8usize, 16, 24, 32, 48, 64]
            .iter()
            .map(|&k| format!("{k}:{:.4}", ratio(&entry.rows[k - 1])))
            .collect();
        println!("{:14} {}", entry.name, line.join(" "));
    }
}

pub fn thin_bound(entries: &[Entry]) -> (usize, usize, usize) {
    let (mut cells, mut tight, mut total) = (0, 0, 0);
    for entry in reps(entries) {
        for p in &entry.rows {
            let side = p.ones.min((1u128 << p.k) - p.ones);
            assert!(
                p.nonlinearity <= side,
                "{} breaks NL <= min(wt, 2^k - wt) at {}",
                entry.name,
                p.k
            );
            assert!(
                p.nonlinearity <= ceiling(p.k),
                "{} breaks the ceiling at {}",
                entry.name,
                p.k
            );
            cells += 2;
            if p.k >= 8 && !affine(entry) {
                total += 1;
                tight += (p.nonlinearity == side) as usize;
            }
        }
    }
    (cells, tight, total)
}

pub fn report_tight(entries: &[Entry]) {
    println!("NL = min(wt, 2^k - wt) per source, k = 8..{KMAX}");
    for entry in reps(entries) {
        let tight: Vec<usize> = entry
            .rows
            .iter()
            .filter(|p| p.k >= 8 && p.nonlinearity == p.ones.min((1u128 << p.k) - p.ones))
            .map(|p| p.k)
            .collect();
        let loose: Vec<String> = entry
            .rows
            .iter()
            .filter(|p| p.k >= 8 && p.nonlinearity != p.ones.min((1u128 << p.k) - p.ones))
            .map(|p| p.k.to_string())
            .collect();
        println!(
            "{:14} tight at {} of {}; loose at k = {}",
            entry.name,
            tight.len(),
            KMAX - 7,
            loose.join(",")
        );
    }
}

pub fn report_density(entries: &[Entry]) {
    println!("P(f = 1) = ones / 2^k at k = 8, 16, 32, 64");
    for entry in entries {
        let line: Vec<String> = [8usize, 16, 32, 64]
            .iter()
            .map(|&k| {
                format!(
                    "{k}:{:.4}",
                    entry.rows[k - 1].ones as f64 / (1u128 << k) as f64
                )
            })
            .collect();
        println!("{:14} {}", entry.name, line.join(" "));
    }
}

pub fn best_by_k(entries: &[Entry]) {
    println!("BEST SOURCE PER k (nl/ceiling, degree>=3 only second)");
    for k in 1..=KMAX {
        let top = entries
            .iter()
            .max_by(|a, b| {
                a.rows[k - 1]
                    .nonlinearity
                    .cmp(&b.rows[k - 1].nonlinearity)
                    .then(a.rows[k - 1].degree.cmp(&b.rows[k - 1].degree))
            })
            .unwrap();
        let cubic = entries
            .iter()
            .filter(|e| e.rows[k - 1].degree >= 3)
            .max_by_key(|e| e.rows[k - 1].nonlinearity);
        let first = &top.rows[k - 1];
        let second = cubic.map_or("-".to_string(), |e| {
            format!(
                "{} {:.4} deg={}",
                e.name,
                ratio(&e.rows[k - 1]),
                e.rows[k - 1].degree
            )
        });
        println!(
            "k={k:2} {} {:.4} deg={} | {second}",
            top.name,
            ratio(first),
            first.degree
        );
    }
}

pub fn best_ratio(entries: &[Entry]) -> (f64, String, usize) {
    let mut top = (0.0f64, String::new(), 0usize);
    for e in reps(entries).into_iter().filter(|e| !affine(e)) {
        for p in e.rows.iter().filter(|p| p.k >= 6) {
            if ratio(p) > top.0 {
                top = (ratio(p), e.name.clone(), p.k);
            }
        }
    }
    top
}

pub fn full_degree(entries: &[Entry]) -> (usize, usize) {
    let (mut full, mut all) = (0, 0);
    for e in reps(entries).into_iter().filter(|e| !affine(e)) {
        for p in e.rows.iter().filter(|p| p.k >= 3) {
            all += 1;
            full += (p.degree == p.k as i32) as usize;
        }
    }
    (full, all)
}

pub fn summary(entries: &[Entry]) {
    let reps = reps(entries);
    let nonaffine = |e: &&&Entry| !affine(e);
    println!(
        "SUMMARY over {} classes (twins merged), the affine ones set aside",
        reps.len()
    );
    let top = best_ratio(entries);
    println!(
        "best nl/ceiling over k >= 6: {:.4} by {} at k = {}",
        top.0, top.1, top.2
    );
    let mut low = (2.0f64, "", 0usize);
    for e in reps.iter().filter(nonaffine) {
        for p in e.rows.iter().filter(|p| p.k >= 8) {
            let share = p.differential as f64 / (1u128 << p.k) as f64;
            if share < low.0 {
                low = (share, &e.name, p.k);
            }
        }
    }
    println!(
        "lowest differential/2^k over k >= 8: {:.4} by {} at k = {}; bent would be 0.5",
        low.0, low.1, low.2
    );
    let (full, all) = full_degree(entries);
    println!("cells of full degree k: {full} of {all}, k = 3..{KMAX}");
    let bal: Vec<String> = reps
        .iter()
        .filter(nonaffine)
        .flat_map(|e| {
            e.rows
                .iter()
                .filter(|p| p.k >= 2 && balanced(p))
                .map(move |p| format!("{}@{}", e.name, p.k))
        })
        .collect();
    println!(
        "balanced cells beyond the affine pair: {}",
        if bal.is_empty() {
            "none".to_string()
        } else {
            bal.join(" ")
        }
    );
    let primes = entries.iter().find(|e| e.name == "primes").unwrap();
    let sigma: Vec<u8> = (0..=KMAX).map(|w| ((w % 4) >= 2) as u8).collect();
    let first = (1..=KMAX).find(|&k| primes.g[..=k] != sigma[..=k]).unwrap();
    println!("primes g equals the quadratic sigma_2 pattern [w mod 4 in 2,3] for k < {first}");
    let mut runs: Vec<(String, usize, usize)> = Vec::new();
    for k in 6..=KMAX {
        let e = reps
            .iter()
            .filter(nonaffine)
            .max_by(|a, b| {
                a.rows[k - 1]
                    .nonlinearity
                    .cmp(&b.rows[k - 1].nonlinearity)
                    .then(a.rows[k - 1].degree.cmp(&b.rows[k - 1].degree))
            })
            .unwrap();
        match runs.last_mut() {
            Some(last) if last.0 == e.name => last.2 = k,
            _ => runs.push((e.name.clone(), k, k)),
        }
    }
    let line: Vec<String> = runs
        .iter()
        .map(|(n, a, b)| format!("{n} {a}..{b}"))
        .collect();
    println!("best class by nl, k = 6..{KMAX}: {}", line.join("; "));
}

pub fn report_dropped(entries: &[Entry]) {
    println!("DROPPED LOW TERMS: the Sources with 0 and 1 kept or dropped, as counts() drops them");
    let holds = |bit: usize| -> Vec<String> {
        classes(entries)
            .into_iter()
            .filter(|class| {
                entries
                    .iter()
                    .find(|e| e.name == class[0])
                    .is_some_and(|e| e.g[bit] == 1)
            })
            .map(|class| class.join(" = "))
            .collect()
    };
    for bit in [0, 1] {
        let sets = holds(bit);
        let sources: usize = sets.iter().map(|s| s.split(" = ").count()).sum();
        println!(
            "contain {bit}: {sources} of {} Sources in {} sets: {}",
            entries.len(),
            sets.len(),
            sets.join("; ")
        );
    }
    for (label, zeros, ones) in [
        ("keep 0 and 1", true, true),
        ("drop 0", false, true),
        ("drop 0 and 1", false, false),
    ] {
        let dropped = build(&variant(zeros, ones));
        let sets = reps(&dropped).len();
        let affine_sets = reps(&dropped).into_iter().filter(|e| affine(e)).count();
        let (_, tight, total) = thin_bound(&dropped);
        let (full, all) = full_degree(&dropped);
        let top = best_ratio(&dropped);
        println!(
            "{label}: {sets} sets, {affine_sets} affine, NL = min(wt, 2^k - wt) on {tight} of {total} cells with k >= 8, full degree on {full} of {all} cells, best nl/ceiling {:.4} by {} at k = {}, brute-force table agrees on {} cells",
            top.0, top.1, top.2, brute(&dropped)
        );
    }
}
