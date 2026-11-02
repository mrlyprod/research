use crate::census::{self, Census, CHAMPIONS};
use crate::rows::{Row, Sheet, Stop, CAP, CEILING, DEEP, TIERS};

pub const BANDS: [i128; 5] = [10, 100, 1_000, 10_000, CEILING];
pub const TAIL: i128 = 30_000;
pub const MODULI: i128 = 64;
pub const SMALL: i128 = 1_000;

fn fold<'a>(sets: impl Iterator<Item = &'a Vec<i128>>) -> Vec<u32> {
    let mut counts = vec![0u32; CEILING as usize + 1];
    for set in sets {
        for &term in set {
            counts[term as usize] += 1;
        }
    }
    counts
}

fn seen(counts: &[u32]) -> usize {
    counts[1..].iter().filter(|&&count| count > 0).count()
}

fn gap(counts: &[u32]) -> usize {
    counts[1..]
        .iter()
        .position(|&count| count == 0)
        .map_or(0, |index| index + 1)
}

fn longest(counts: &[u32]) -> usize {
    let mut best = 0;
    let mut here = 0;
    for &count in &counts[1..] {
        here = if count > 0 { here + 1 } else { 0 };
        best = best.max(here);
    }
    best
}

fn inside(window: &[i128]) -> Vec<i128> {
    let mut out: Vec<i128> = window
        .iter()
        .copied()
        .filter(|term| (1..=CEILING).contains(term))
        .collect();
    out.sort_unstable();
    out.dedup();
    out
}

fn edge(terms: &[i128]) -> Option<usize> {
    let mut previous: Option<i128> = None;
    for (index, &term) in terms.iter().enumerate() {
        if previous.is_some_and(|last| term <= last) {
            return None;
        }
        if term > CEILING {
            return Some(index);
        }
        previous = Some(term);
    }
    None
}

pub fn degree(head: &[i128]) -> Option<usize> {
    let mut layer = head.to_vec();
    for order in 0..head.len() {
        if layer.len() < 2 {
            return None;
        }
        if layer.iter().all(|&value| value == layer[0]) {
            return Some(order);
        }
        layer = layer.windows(2).map(|pair| pair[1] - pair[0]).collect();
    }
    None
}

pub fn extend(head: &[i128], length: usize) -> Option<Vec<i128>> {
    let order = degree(head)?;
    let mut table = vec![head.to_vec()];
    for _ in 0..order {
        let last: Vec<i128> = table
            .last()?
            .windows(2)
            .map(|pair| pair[1] - pair[0])
            .collect();
        table.push(last);
    }
    let mut ends: Vec<i128> = table
        .iter()
        .map(|layer| *layer.last().expect("a difference layer is nonempty"))
        .collect();
    let mut out = head.to_vec();
    while out.len() < length {
        for index in (0..order).rev() {
            ends[index] = ends[index].checked_add(ends[index + 1])?;
        }
        out.push(ends[0]);
    }
    out.truncate(length);
    Some(out)
}

pub fn rebuild(row: &Row) -> Option<Vec<i128>> {
    let terms = extend(&row.head, CAP)?;
    let window = match row.stop {
        Stop::Cap => terms,
        Stop::Ceiling => terms[..=edge(&terms)?].to_vec(),
        Stop::Budget => return None,
    };
    Some(inside(&window))
}

fn factors() -> (Vec<u32>, Vec<u32>, Vec<i128>) {
    let size = CEILING as usize + 1;
    let mut great = vec![0u32; size];
    let mut divisors = vec![0u32; size];
    for step in 1..size {
        for slot in (step..size).step_by(step) {
            divisors[slot] += 1;
        }
    }
    for value in 2..size {
        if great[value] == 0 {
            for slot in (value..size).step_by(value) {
                great[slot] = value as u32;
            }
        }
    }
    let primes = (2..size)
        .filter(|&value| great[value] == value as u32)
        .map(|value| value as i128)
        .collect();
    (great, divisors, primes)
}

fn perfect() -> Vec<i128> {
    let mut out = Vec::new();
    let mut base = 2i128;
    while base * base <= CEILING {
        let mut value = base * base;
        while value <= CEILING {
            out.push(value);
            value *= base;
        }
        base += 1;
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn mean(counts: &[u32], values: &[i128]) -> f64 {
    values
        .iter()
        .map(|&value| counts[value as usize] as f64)
        .sum::<f64>()
        / values.len() as f64
}

fn family(power: u32) -> Vec<i128> {
    let mut out = Vec::new();
    let mut root = 1i128;
    while root.pow(power) <= CEILING {
        out.push(root.pow(power));
        root += 1;
    }
    out
}

fn depth(sheet: &Sheet, book: &Census) {
    println!("DEPTH");
    let heads: Vec<Vec<i128>> = sheet.rows.iter().map(|row| inside(&row.head)).collect();
    let short = fold(heads.iter());
    let shallow = fold(sheet.rows.iter().map(|row| &row.shallow));
    let dropped = fold(sheet.rows.iter().map(|row| &row.tail));
    for (name, counts) in [("8", &short), ("32", &shallow), ("48", &book.counts)] {
        println!(
            "depth window {name} written {} missed {} first miss {}",
            seen(counts),
            CEILING as usize - seen(counts),
            gap(counts)
        );
    }
    let full = seen(&book.counts);
    println!(
        "depth reached only past term 8 {} only past term 32 {}",
        full - seen(&short),
        full - seen(&shallow)
    );
    let carried: u64 = dropped.iter().map(|&count| count as u64).sum();
    println!(
        "depth without each row's first term incidences {carried} of {} lost {} written {} never {}",
        book.incidences,
        book.incidences - carried,
        seen(&dropped),
        census::WINDOWS
            .iter()
            .map(|&window| dropped[1..=window].iter().filter(|&&count| count == 0).count().to_string())
            .collect::<Vec<_>>()
            .join(" ")
    );
    let book = Census {
        counts: dropped,
        incidences: 0,
        repeats: 0,
        low: 0,
    };
    println!(
        "depth without each row's first term leaders {}",
        census::champions(&book, 4)
            .iter()
            .map(|(value, count)| format!("{value} {count}"))
            .collect::<Vec<_>>()
            .join(" ")
    );
}

fn model(sheet: &Sheet) {
    println!("MODEL");
    let mut orders = [0usize; 7];
    for row in &sheet.rows {
        if let Some(order) = degree(&row.head) {
            orders[order] += 1;
        }
    }
    println!(
        "model rows with a head degree at most 6 {} by degree {}",
        orders.iter().sum::<usize>(),
        orders
            .iter()
            .map(|count| count.to_string())
            .collect::<Vec<_>>()
            .join(" ")
    );
    for stop in [Stop::Cap, Stop::Ceiling, Stop::Budget] {
        let batch: Vec<&Row> = sheet
            .rows
            .iter()
            .filter(|row| row.stop == stop && degree(&row.head).is_some())
            .collect();
        let built: Vec<Option<Vec<i128>>> = batch.iter().map(|row| rebuild(row)).collect();
        let tested = built.iter().filter(|slot| slot.is_some()).count();
        let pass = batch
            .iter()
            .zip(&built)
            .filter(|(row, slot)| slot.as_ref().is_some_and(|terms| *terms == row.written))
            .count();
        let mut failed: Vec<usize> = batch
            .iter()
            .zip(&built)
            .filter(|(row, slot)| slot.as_ref().is_some_and(|terms| *terms != row.written))
            .filter_map(|(row, _)| degree(&row.head))
            .collect();
        failed.sort_unstable();
        failed.dedup();
        println!(
            "model rebuild stop {} rows {} tested {tested} pass {pass} fail {} failing degrees {}",
            stop.slug(),
            batch.len(),
            tested - pass,
            failed
                .iter()
                .map(|order| order.to_string())
                .collect::<Vec<_>>()
                .join(" ")
        );
    }
}

fn deeper(sheet: &Sheet, book: &Census) {
    println!("DEEPER");
    let mut counts = book.counts.clone();
    let mut checked = 0;
    let mut used = 0;
    for row in &sheet.rows {
        if row.stop != Stop::Cap {
            continue;
        }
        let Some(built) = rebuild(row) else {
            continue;
        };
        checked += 1;
        if built != row.written {
            continue;
        }
        let Some(terms) = extend(&row.head, DEEP) else {
            continue;
        };
        used += 1;
        for term in inside(&terms[CAP..]) {
            counts[term as usize] += 1;
        }
    }
    println!(
        "deeper cap rows rebuilt {checked} extended {used} of {}",
        sheet
            .rows
            .iter()
            .filter(|row| row.stop == Stop::Cap)
            .count()
    );
    println!(
        "deeper window {DEEP} written at least {} first miss {} longest written run at least {}",
        seen(&counts),
        gap(&counts),
        longest(&counts)
    );
    let squares = family(2);
    println!(
        "deeper squares at least {} of {} first missed square {}",
        squares
            .iter()
            .filter(|&&value| counts[value as usize] > 0)
            .count(),
        squares.len(),
        squares
            .iter()
            .find(|&&value| counts[value as usize] == 0)
            .copied()
            .unwrap_or(0)
    );
}

fn arithmetic(sheet: &Sheet, book: &Census) {
    println!("ARITHMETIC");
    let counts = &book.counts;
    for power in 2u32..=6 {
        let batch = family(power);
        let written: Vec<i128> = batch
            .iter()
            .copied()
            .filter(|&value| counts[value as usize] > 0)
            .collect();
        println!(
            "arith power {power} written {} of {} first missed {} largest written {}",
            written.len(),
            batch.len(),
            batch
                .iter()
                .find(|&&value| counts[value as usize] == 0)
                .copied()
                .unwrap_or(0),
            written.last().copied().unwrap_or(0)
        );
    }
    let powers = perfect();
    let carried: u64 = powers
        .iter()
        .map(|&value| counts[value as usize] as u64)
        .sum();
    let share = carried as f64 / book.incidences as f64;
    let density = powers.len() as f64 / CEILING as f64;
    println!(
        "arith perfect powers {} incidences {carried} share {share:.4} against a density {density:.6} ratio {:.2}",
        powers.len(),
        share / density
    );
    let squares = family(2);
    let top = squares
        .iter()
        .copied()
        .rfind(|&value| counts[value as usize] > 0)
        .unwrap_or(0);
    let writers = census::writers(sheet, top);
    println!(
        "arith largest written square {top} rows {} {}",
        writers.len(),
        writers
            .iter()
            .map(|row| row.name.as_str())
            .collect::<Vec<_>>()
            .join(" ")
    );
    println!(
        "arith square frontier {}",
        (96i128..=100)
            .map(|root| format!("{root} rows {}", counts[(root * root) as usize]))
            .collect::<Vec<_>>()
            .join(" ")
    );
    let (great, divisors, primes) = factors();
    let written: Vec<i128> = primes
        .iter()
        .copied()
        .filter(|&value| counts[value as usize] > 0)
        .collect();
    println!(
        "arith primes written {} of {} first missed {} above 10000 {} of {}",
        written.len(),
        primes.len(),
        primes
            .iter()
            .find(|&&value| counts[value as usize] == 0)
            .copied()
            .unwrap_or(0),
        written.iter().filter(|&&value| value > 10_000).count(),
        primes.iter().filter(|&&value| value > 10_000).count()
    );
    let mut classes = 0;
    let mut empty = 0;
    for modulus in 2..=MODULI {
        for residue in 0..modulus {
            classes += 1;
            let first = 10_000 + (residue - 10_000).rem_euclid(modulus);
            if !(first..=CEILING)
                .step_by(modulus as usize)
                .any(|value| counts[value as usize] > 0)
            {
                empty += 1;
            }
        }
    }
    println!("arith residue classes mod 2..{MODULI} on 10000..{CEILING} {classes} with no written integer {empty}");
    for modulus in [6i128, 12] {
        let mut tally = vec![0usize; modulus as usize];
        for value in 10_000..=CEILING {
            if counts[value as usize] > 0 {
                tally[(value % modulus) as usize] += 1;
            }
        }
        let high = tally
            .iter()
            .max()
            .copied()
            .expect("a residue tally is nonempty");
        let low = tally
            .iter()
            .min()
            .copied()
            .expect("a residue tally is nonempty");
        println!(
            "arith mod {modulus} written by residue {} high {high} low {low} ratio {:.2}",
            tally
                .iter()
                .map(|count| count.to_string())
                .collect::<Vec<_>>()
                .join(" "),
            high as f64 / low as f64
        );
    }
    let mut floor = 1i128;
    for &roof in &BANDS {
        let band: Vec<i128> = (10_000..=CEILING)
            .filter(|&value| {
                i128::from(great[value as usize]) > floor
                    && i128::from(great[value as usize]) <= roof
            })
            .collect();
        println!(
            "arith greatest prime factor in {floor}..{roof} size {} written share {:.4}",
            band.len(),
            band.iter()
                .filter(|&&value| counts[value as usize] > 0)
                .count() as f64
                / band.len() as f64
        );
        floor = roof;
    }
    let all: Vec<i128> = (1..=SMALL).collect();
    let square: Vec<i128> = squares
        .iter()
        .copied()
        .filter(|&value| value <= SMALL)
        .collect();
    let rich: Vec<i128> = all
        .iter()
        .copied()
        .filter(|&value| divisors[value as usize] >= 8)
        .collect();
    let power: Vec<i128> = powers
        .iter()
        .copied()
        .filter(|&value| value <= SMALL)
        .collect();
    println!(
        "arith mean rows on 1..{SMALL} all {:.2} squares {:.2} at least eight divisors {:.2} perfect powers {:.2}",
        mean(counts, &all),
        mean(counts, &square),
        mean(counts, &rich),
        mean(counts, &power)
    );
}

fn spectrum(sheet: &Sheet, book: &Census) {
    println!("SPECTRUM");
    let counts = &book.counts;
    let mut heights: Vec<u32> = counts[1..]
        .iter()
        .copied()
        .filter(|&count| count > 0)
        .collect();
    heights.sort_unstable();
    heights.dedup();
    println!(
        "spectrum distinct multiplicities {} max {}",
        heights.len(),
        heights.last().copied().unwrap_or(0)
    );
    let champions = census::champions(book, CHAMPIONS);
    let carried: u64 = champions.iter().map(|&(_, count)| count as u64).sum();
    let mut set: Vec<usize> = champions.iter().map(|&(value, _)| value).collect();
    set.sort_unstable();
    println!(
        "spectrum top {CHAMPIONS} incidences {carried} of {} share {:.4} largest champion {}",
        book.incidences,
        carried as f64 / book.incidences as f64,
        set.last().copied().unwrap_or(0)
    );
    println!(
        "spectrum champion set ascending {}",
        set.iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(" ")
    );
    let above = |level: u32| counts[1..].iter().filter(|&&count| count >= level).count();
    let ratio = above(2) as f64 / above(1) as f64;
    println!(
        "spectrum S(1) {} S(2) {} ratio {ratio:.4} geometric S(64) {:.3e} observed {}",
        above(1),
        above(2),
        above(1) as f64 * ratio.powi(63),
        above(64)
    );
    let whole: Vec<i128> = (1..=CEILING)
        .filter(|&value| counts[value as usize] > 0)
        .collect();
    let high: Vec<i128> = whole
        .iter()
        .copied()
        .filter(|&value| value >= TAIL)
        .collect();
    println!(
        "spectrum written {} above {TAIL} {}",
        whole.len(),
        high.len()
    );
    for tier in TIERS {
        let mine = fold(
            sheet
                .rows
                .iter()
                .filter(|row| row.tier == tier)
                .map(|row| &row.written),
        );
        let rest = fold(
            sheet
                .rows
                .iter()
                .filter(|row| row.tier != tier)
                .map(|row| &row.written),
        );
        println!(
            "spectrum tier {} covers {} exclusive {} above {TAIL} {}",
            tier.slug(),
            seen(&mine),
            whole
                .iter()
                .filter(|&&value| mine[value as usize] > 0 && rest[value as usize] == 0)
                .count(),
            high.iter()
                .filter(|&&value| mine[value as usize] > 0)
                .count()
        );
    }
    let mut families: Vec<Vec<i128>> = sheet
        .rows
        .iter()
        .map(|row| {
            row.written
                .iter()
                .copied()
                .filter(|&value| value >= TAIL)
                .collect::<Vec<i128>>()
        })
        .filter(|set| !set.is_empty())
        .collect();
    families.sort();
    families.dedup();
    let mut owners = vec![0u32; CEILING as usize + 1];
    for set in &families {
        for &value in set {
            owners[value as usize] += 1;
        }
    }
    println!(
        "spectrum written sets above {TAIL} {} owning a tail integer alone {} covering {} of {}",
        families.len(),
        families
            .iter()
            .filter(|set| set.iter().any(|&value| owners[value as usize] == 1))
            .count(),
        high.iter()
            .filter(|&&value| owners[value as usize] == 1)
            .count(),
        high.len()
    );
    for (column, value) in [
        ("euler.side", 1i128),
        ("peak.side", 12),
        ("heights.side", 9),
        ("heights.side", 33),
    ] {
        let batch: Vec<&Row> = sheet
            .rows
            .iter()
            .filter(|row| row.name.ends_with(column))
            .collect();
        println!(
            "spectrum column {column} writes {value} in {} of {} rows",
            batch
                .iter()
                .filter(|row| row.written.binary_search(&value).is_ok())
                .count(),
            batch.len()
        );
    }
    let batch: Vec<&Row> = sheet
        .rows
        .iter()
        .filter(|row| row.name.ends_with("heights.side"))
        .collect();
    let mut tally: Vec<(i128, usize)> = (1..=100i128)
        .map(|value| {
            (
                value,
                batch
                    .iter()
                    .filter(|row| row.written.binary_search(&value).is_ok())
                    .count(),
            )
        })
        .collect();
    tally.sort_by(|left, right| right.1.cmp(&left.1).then(left.0.cmp(&right.0)));
    tally.truncate(8);
    let mut leaders: Vec<i128> = tally.iter().map(|&(value, _)| value).collect();
    leaders.sort_unstable();
    println!(
        "spectrum heights.side leaders to 100 {} every one a multiple of eight plus one {}",
        leaders
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(" "),
        leaders.iter().all(|value| (value - 1) % 8 == 0)
    );
}

pub fn report(sheet: &Sheet, book: &Census) {
    depth(sheet, book);
    println!();
    model(sheet);
    println!();
    deeper(sheet, book);
    println!();
    arithmetic(sheet, book);
    println!();
    spectrum(sheet, book);
}
