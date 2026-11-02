use crate::cache::data_dir;
use crate::design::classes;
use crate::grammar::Grammar;
use crate::mask::{Kernel, MaskSpec};
use crate::recipe::{Mode, Recipe, Side};
use crate::rows::{each_row, Row};
use mrly_core::{json, Json};
use mrly_life::Source;
use mrly_core::Result;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs;
use std::path::Path;

pub const STOPS: [&str; 6] = ["dead", "still", "loop", "mover", "wall", "limit"];

pub fn class(recipe: &Recipe) -> String {
    match recipe.mode {
        Mode::Island => format!("island x{}", recipe.tiling),
        Mode::Wallpaper => "wallpaper".to_string(),
    }
}

pub fn mask_kind(mask: MaskSpec) -> String {
    match mask {
        MaskSpec::Design { side, level, .. } => format!("{side}.{level}"),
        MaskSpec::Copy { inverted: false } => "copy".to_string(),
        MaskSpec::Copy { inverted: true } => "copyinv".to_string(),
    }
}

pub fn rule_class(birth: &Side, survive: &Side) -> String {
    match (birth, survive) {
        (Side::Drawn { seq: b, .. }, Side::Drawn { seq: s, .. }) if b == s => b.name(),
        (Side::Drawn { seq: b, .. }, Side::Drawn { seq: s, .. }) => {
            format!("{}/{}", b.name(), s.name())
        }
        _ => "life".to_string(),
    }
}

fn stop_index(stop: &str) -> usize {
    STOPS.iter().position(|s| *s == stop).unwrap_or(5)
}

fn mask_side(recipe: &Recipe) -> usize {
    match recipe.mask {
        MaskSpec::Design { side, level, .. } => side.pow(level as u32),
        MaskSpec::Copy { .. } => recipe.board_side(),
    }
}

fn light_bound(recipe: &Recipe) -> usize {
    let side = recipe.board_side();
    let cap = recipe.cap.max(side);
    let reach = (mask_side(recipe) - 1) / 2;
    let distance = ((cap - side) / 2 + 1).saturating_sub(reach);
    distance.div_ceil(reach.max(1)).max(1)
}

fn seed_key(name: &str) -> String {
    let words: Vec<&str> = name.split(' ').collect();
    let mut out: Vec<&str> = vec![words[0]];
    out.extend(words.iter().skip(3).filter(|w| !w.starts_with("cap:")));
    out.join(" ")
}

fn median(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    values[(values.len() - 1) / 2]
}

fn near(fraction: f64) -> &'static str {
    for (mark, word) in [
        (0.0, " (within 0.02 of 0)"),
        (0.5, " (within 0.02 of 1/2)"),
        (1.0, " (within 0.02 of 1)"),
    ] {
        if (fraction - mark).abs() <= 0.02 {
            return word;
        }
    }
    ""
}

fn is_linear(side: &Side) -> bool {
    matches!(
        side,
        Side::Drawn {
            seq: Source::Odds,
            ones: true,
            ..
        }
    )
}

fn order_of_two(n: usize) -> usize {
    let mut k = 1;
    let mut p = 2 % n;
    while p != 1 {
        p = p * 2 % n;
        k += 1;
    }
    k
}

#[derive(Default)]
struct Wall {
    islands: usize,
    walls: usize,
    movers_to_wall: usize,
    ratios: Vec<f64>,
}

#[derive(Default)]
struct Tally {
    rows: usize,
    cut: Option<(usize, String)>,
    by_class: BTreeMap<String, [usize; 6]>,
    by_seed: BTreeMap<String, usize>,
    grid: BTreeMap<(String, String), [usize; 6]>,
    wall: BTreeMap<String, Wall>,
    limit_growth: Vec<f64>,
    fates: HashMap<String, [u8; 4]>,
    linear: BTreeMap<(usize, Option<usize>), usize>,
    linear_islands: BTreeMap<String, usize>,
    affine: BTreeMap<(usize, Option<usize>), usize>,
    sequences: BTreeMap<Vec<u32>, (usize, String)>,
}

fn plane_class(recipe: &Recipe) -> bool {
    matches!(recipe.seed, crate::Shape::Code(code) if classes().contains(&code))
}

fn gather() -> Result<Tally> {
    let mut t = Tally::default();
    let fate_slot = |recipe: &Recipe| match (recipe.mode, recipe.tiling) {
        (Mode::Wallpaper, _) => 3,
        (_, 1) => 0,
        (_, 3) => 1,
        _ => 2,
    };
    each_row(|row: Row| {
        let Ok(recipe) = row.recipe() else {
            return;
        };
        t.rows += 1;
        if t.cut.as_ref().is_none_or(|(i, _)| row.index > *i) {
            t.cut = Some((row.index, row.name.clone()));
        }
        let stop = stop_index(&row.stop);
        let class = class(&recipe);
        t.by_class.entry(class.clone()).or_default()[stop] += 1;
        *t.by_seed
            .entry(format!("{}.{}", recipe.side, recipe.level))
            .or_default() += 1;
        t.grid
            .entry((
                mask_kind(recipe.mask),
                rule_class(&recipe.birth, &recipe.survive),
            ))
            .or_default()[stop] += 1;
        if recipe.mode == Mode::Island {
            let wall = t.wall.entry(mask_kind(recipe.mask)).or_default();
            wall.islands += 1;
            if row.stop == "wall" {
                wall.walls += 1;
                wall.ratios
                    .push(row.gen as f64 / light_bound(&recipe) as f64);
            }
            if row.wall_by.is_some() {
                wall.movers_to_wall += 1;
            }
        }
        if row.stop == "limit" && recipe.mode == Mode::Island {
            t.limit_growth.push(row.growth);
        }
        let copy = matches!(recipe.mask, MaskSpec::Copy { .. });
        if !(copy && recipe.tiling > 1) {
            t.fates.entry(seed_key(&row.name)).or_insert([255; 4])[fate_slot(&recipe)] = stop as u8;
        }
        if is_linear(&recipe.birth) && is_linear(&recipe.survive) {
            if recipe.mode == Mode::Wallpaper {
                *t.linear
                    .entry((recipe.board_side(), row.period))
                    .or_default() += 1;
            } else {
                *t.linear_islands.entry(row.stop.clone()).or_default() += 1;
            }
        }
        if recipe.mode == Mode::Wallpaper
            && recipe.birth == Side::drawn(Source::Evens, true, false)
            && recipe.survive == recipe.birth
        {
            *t.affine
                .entry((recipe.board_side(), row.period))
                .or_default() += 1;
        }
        if recipe.mode == Mode::Wallpaper && plane_class(&recipe) && !row.flips {
            let len = match (row.stop.as_str(), row.period) {
                ("limit", _) | ("wall", _) => row.pop.first.len(),
                (_, Some(_)) => row.gen.min(row.pop.first.len()),
                _ => row.pop.first.len(),
            };
            let seq: Vec<u32> = row.pop.first[..len].to_vec();
            let entry = t
                .sequences
                .entry(seq)
                .or_insert((row.index, row.name.clone()));
            if row.index < entry.0 {
                *entry = (row.index, row.name.clone());
            }
        }
    })?;
    Ok(t)
}

fn totals(grammar: &Grammar) -> Result<(usize, BTreeMap<String, usize>, BTreeMap<String, usize>)> {
    let mut by_class = BTreeMap::new();
    let mut by_seed = BTreeMap::new();
    let total = grammar.walk(0, usize::MAX, |_, recipe| {
        *by_class.entry(class(&recipe)).or_default() += 1;
        *by_seed
            .entry(format!("{}.{}", recipe.side, recipe.level))
            .or_default() += 1;
    })?;
    Ok((total, by_class, by_seed))
}

fn budgets(grammar: &Grammar) -> Result<BTreeSet<usize>> {
    let mut out = BTreeSet::new();
    for (spec, _) in grammar.masks()? {
        if let MaskSpec::Design { shape, side, level } = spec {
            out.insert(Kernel::design(shape, side, level)?.budget);
        }
    }
    for (shape, side, level, _) in grammar.seeds()? {
        for &tiling in &grammar.tilings {
            for inverted in [false, true] {
                if let Ok(kernel) = Kernel::copy(shape, side, level, tiling, inverted) {
                    if kernel.budget > 0 {
                        out.insert(kernel.budget);
                    }
                }
            }
        }
    }
    Ok(out)
}

fn budget_law(grammar: &Grammar) -> Result<()> {
    let budgets = budgets(grammar)?;
    println!(
        "budget law: {} distinct mask budgets in the grammar, {:?}",
        budgets.len(),
        budgets
    );
    println!("source: values inside budget 8 | 24 | 80, budgets of the grammar with no value (rolls off)");
    for &seq in &grammar.sources {
        let side = Side::drawn(seq, false, false);
        let empty = budgets
            .iter()
            .filter(|&&b| side.values(b).map(|v| v.is_empty()).unwrap_or(true))
            .count();
        println!(
            "{}: {:?} | {:?} | {:?}, empty on {empty} of {} budgets",
            seq.name(),
            side.values(8)?,
            side.values(24)?,
            side.values(80)?,
            budgets.len()
        );
    }
    let mut dead_rules = 0;
    let mut live_rules = 0;
    let rules = grammar.rules();
    for (birth, survive) in &rules {
        let any = budgets
            .iter()
            .any(|&b| !birth.values(b).unwrap_or_default().is_empty());
        let both = budgets.iter().any(|&b| {
            !birth.values(b).unwrap_or_default().is_empty()
                && !survive.values(b).unwrap_or_default().is_empty()
        });
        dead_rules += usize::from(!any);
        live_rules += usize::from(both);
    }
    println!(
        "rules: {} in the grammar, {dead_rules} with an empty birth side on every budget, {live_rules} with both sides live on some budget",
        rules.len()
    );
    Ok(())
}

fn fits_polynomial(seq: &[u32], degree: usize) -> bool {
    let mut d: Vec<i64> = seq.iter().map(|&v| i64::from(v)).collect();
    for _ in 0..degree {
        d = d.windows(2).map(|w| w[1] - w[0]).collect();
    }
    d.len() >= 2 && d.iter().all(|&v| v == d[0])
}

fn has_short_period(seq: &[u32]) -> bool {
    (1..=4).any(|p| seq.len() > 2 * p && seq.iter().skip(p).zip(seq.iter()).all(|(a, b)| a == b))
}

fn oeis_material(seq: &[u32]) -> bool {
    let distinct: BTreeSet<u32> = seq.iter().copied().collect();
    seq.len() >= 12
        && seq.windows(2).all(|w| w[1] >= w[0])
        && distinct.len() >= 8
        && !has_short_period(seq)
        && !fits_polynomial(seq, 1)
        && !fits_polynomial(seq, 2)
}

fn oeis_hits(dump: &Path, sequences: &[Vec<u32>]) -> (Vec<Vec<String>>, Vec<Vec<String>>) {
    let text = fs::read_to_string(dump).unwrap_or_default();
    let runs: Vec<String> = sequences
        .iter()
        .map(|seq| {
            format!(
                ",{},",
                seq.iter()
                    .take(12)
                    .map(|v| v.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )
        })
        .collect();
    let mut full: Vec<Vec<String>> = vec![Vec::new(); sequences.len()];
    let windows: Vec<Vec<String>> = sequences
        .iter()
        .map(|seq| {
            seq.windows(5)
                .take(8)
                .map(|w| {
                    format!(
                        ",{},",
                        w.iter()
                            .map(|v| v.to_string())
                            .collect::<Vec<_>>()
                            .join(",")
                    )
                })
                .collect()
        })
        .collect();
    let mut hits: Vec<Vec<String>> = vec![Vec::new(); sequences.len()];
    for line in text.lines() {
        let Some((id, data)) = line.split_once(' ') else {
            continue;
        };
        for (i, list) in windows.iter().enumerate() {
            if hits[i].len() < 3 && list.iter().any(|w| data.contains(w.as_str())) {
                hits[i].push(id.to_string());
            }
            if data.contains(runs[i].as_str()) {
                full[i].push(id.to_string());
            }
        }
    }
    (hits, full)
}

pub fn report(grammar: &Grammar, oeis: Option<&Path>) -> Result<()> {
    let (total, class_totals, seed_totals) = totals(grammar)?;
    let t = gather()?;
    let (index, name) = t.cut.clone().unwrap_or((0, String::new()));
    println!(
        "cut: {} rows read of {total} recipes, last index {index}, last recipe {name:?}",
        t.rows
    );
    for (class, n) in &class_totals {
        let done: usize = t.by_class.get(class).map(|c| c.iter().sum()).unwrap_or(0);
        println!(
            "class {class}: {done} of {n} done{}",
            if done == *n { ", complete" } else { ", cut" }
        );
    }
    for (seed, n) in &seed_totals {
        let done = t.by_seed.get(seed).copied().unwrap_or(0);
        println!(
            "seeds {seed}: {done} of {n} done{}",
            if done == *n { ", complete" } else { ", cut" }
        );
    }
    println!("stop tally by class: {:?}", STOPS);
    for (class, counts) in &t.by_class {
        println!("{class}: {counts:?}");
    }
    let mut growth = t.limit_growth.clone();
    let growing = growth.iter().filter(|&&g| g >= 0.5).count();
    println!(
        "island limit rows: {}, box growth slope median {:.3}, {growing} with slope at least 0.5",
        growth.len(),
        median(&mut growth)
    );
    println!("wall law by mask kind: islands, wall stops, fraction, movers with a wall_by, median wall generation over the light-speed bound");
    for (kind, w) in &t.wall {
        let fraction = w.walls as f64 / w.islands.max(1) as f64;
        let mut ratios = w.ratios.clone();
        let light = ratios.iter().filter(|&&r| r == 1.0).count();
        let share = light as f64 / w.walls.max(1) as f64;
        println!(
            "mask {kind}: {} islands, {} walls, {fraction:.3}{}, {} movers to the wall, ratio median {:.2}, {light} at light speed, {share:.3}{}",
            w.islands,
            w.walls,
            near(fraction),
            w.movers_to_wall,
            median(&mut ratios),
            near(share)
        );
    }
    budget_law(grammar)?;
    println!("island against wallpaper on one seed, mask and rule (copy masks only at tiling 1): pairs, equal fates, fraction");
    for (slot, label) in [(0usize, "x1"), (1, "x3"), (2, "x5")] {
        let (mut pairs, mut equal) = (0usize, 0usize);
        for fates in t.fates.values() {
            if fates[slot] != 255 && fates[3] != 255 {
                pairs += 1;
                equal += usize::from(fates[slot] == fates[3]);
            }
        }
        let fraction = equal as f64 / pairs.max(1) as f64;
        println!(
            "island {label} against wallpaper: {pairs} pairs, {equal} equal, {fraction:.3}{}",
            near(fraction)
        );
    }
    println!("odd-rule law: linear rules are odds_ones and odds_zeros_ones (odd counts on both sides); wallpaper periods by torus side, with 2^ord_n(2) - 1");
    for ((side, period), n) in &t.linear {
        let k = order_of_two(*side);
        let bound = (1u128 << k) - 1;
        let divides = period.map(|p| bound % p as u128 == 0);
        println!(
            "torus {side}: period {period:?} x {n}, ord_{side}(2) = {k}, bound {bound}, divides {divides:?}"
        );
    }
    println!("linear rules on islands: {:?}", t.linear_islands);
    println!("affine rule evens_zeros on wallpaper: periods by torus side");
    for ((side, period), n) in &t.affine {
        println!("torus {side}: period {period:?} x {n}");
    }
    let mut candidates: Vec<(Vec<u32>, usize, String)> = t
        .sequences
        .iter()
        .filter(|(seq, _)| oeis_material(seq))
        .map(|(seq, (i, name))| (seq.clone(), *i, name.clone()))
        .collect();
    candidates.sort_by_key(|(seq, i, _)| {
        let distinct: BTreeSet<u32> = seq.iter().copied().collect();
        (std::cmp::Reverse(distinct.len()), *i)
    });
    candidates.truncate(10);
    println!(
        "population sequences: {} distinct on wallpaper plane-class seeds without a flipping background, {} OEIS material (12 terms or more, non-decreasing, 8 distinct values, no fit of degree 1 or 2), ten listed",
        t.sequences.len(),
        t.sequences.keys().filter(|s| oeis_material(s)).count()
    );
    let (hits, full) = match oeis {
        Some(path) => oeis_hits(
            path,
            &candidates
                .iter()
                .map(|(s, _, _)| s.clone())
                .collect::<Vec<_>>(),
        ),
        None => (
            vec![Vec::new(); candidates.len()],
            vec![Vec::new(); candidates.len()],
        ),
    };
    println!(
        "sequence: first 12 terms | 5-term window hits in the dump (at most 3) | 12-term run hits"
    );
    for (((seq, _, name), hit), run) in candidates.iter().zip(&hits).zip(&full) {
        let terms: Vec<String> = seq.iter().take(12).map(|v| v.to_string()).collect();
        println!(
            "{name}: {} | {} | {}",
            terms.join(","),
            if hit.is_empty() {
                "absent".to_string()
            } else {
                hit.join(" ")
            },
            if run.is_empty() {
                "absent".to_string()
            } else {
                run.join(" ")
            }
        );
    }
    write_tally(&t, total)?;
    Ok(())
}

fn write_tally(t: &Tally, total: usize) -> Result<()> {
    let masks: Vec<String> = t
        .grid
        .keys()
        .map(|(m, _)| m.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let rules: Vec<String> = t
        .grid
        .keys()
        .map(|(_, r)| r.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let grid: Vec<Vec<Vec<usize>>> = rules
        .iter()
        .map(|r| {
            masks
                .iter()
                .map(|m| {
                    t.grid
                        .get(&(m.clone(), r.clone()))
                        .map(|c| c.to_vec())
                        .unwrap_or_else(|| vec![0; 6])
                })
                .collect()
        })
        .collect();
    let (index, name) = t.cut.clone().unwrap_or((0, String::new()));
    let value: Json = json!({
        "rows": t.rows,
        "total": total,
        "cut": { "index": index, "name": name },
        "stops": STOPS,
        "masks": masks,
        "rules": rules,
        "grid": grid,
    });
    let path = data_dir().join("tally.json");
    fs::create_dir_all(data_dir()).map_err(|e| mrly_core::Error::Value(e.to_string()))?;
    fs::write(&path, format!("{value}\n")).map_err(|e| mrly_core::Error::Value(e.to_string()))?;
    println!("tally written to {}", path.display());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_light_speed_bound_and_the_frobenius_order_read_as_expected() {
        let recipe = Recipe::parse(
            "code1.3.1 x1 island mask:carpet.3.1 birth:evens survive:evens x:256 cap:151",
        )
        .unwrap();
        assert_eq!(light_bound(&recipe), 74);
        assert_eq!(order_of_two(3), 2);
        assert_eq!(order_of_two(7), 3);
        assert_eq!(order_of_two(9), 6);
        assert!(oeis_material(&[
            1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377
        ]));
        assert!(!oeis_material(&[
            1, 2, 4, 7, 11, 16, 22, 29, 37, 46, 56, 67, 79
        ]));
        assert!(!oeis_material(&[1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23]));
        assert!(!oeis_material(&[1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1]));
        assert_eq!(
            seed_key("code1.3.1 x3 island mask:copy birth:evens survive:evens x:256 cap:151"),
            "code1.3.1 mask:copy birth:evens survive:evens x:256"
        );
    }
}
