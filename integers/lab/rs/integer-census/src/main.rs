mod census;
mod hunt;
mod rows;
mod tables;

use census::{Census, CHAMPIONS, MISSES, WINDOWS};
use ledger::{designs, Axis, Closed, Cost, Measure, Tier, SPACES, TERMS};
use rows::{counted, Row, Sheet, Stop, CAP, CEILING, CELLS};
use std::env;
use std::path::{Path, PathBuf};
use std::time::Instant;

fn definition() -> Vec<String> {
    vec![
        "A registry row is one (design, measure, axis) key of `ledger::keys` over the closed, convolved, side and level tiers and the first twelve measures, the registry the census was run on.".to_string(),
        format!("A row's rendered window is its first `min({CAP}, B)` terms, `B` the leading terms whose footprint fits {CELLS} cells, under the ledger's own budget of {CELLS} cells a term."),
        "A term's footprint is its render cost: 1 cell for a measure closed on both axes, `number^dimension + level * span` for a convolved measure, `number^(dimension * level)` for a grid measure and for the face counts, which close on the level axis but are read here to the same window they rendered.".to_string(),
        format!("A row whose rendered terms are strictly increasing stops at the first term above {CEILING}; the count of rows truncated this way is printed, never assumed to lose nothing."),
        format!("Row `R` writes `n` iff `n` is a term of `R` inside `R`'s rendered window and `1 <= n <= {CEILING}`."),
        "Multiplicity counts rows, not `(row, index)` pairs: a row writing `n` at several indices counts once.".to_string(),
        "An integer `n` appears iff some row writes it, and is missed iff no row writes it.".to_string(),
    ]
}

fn expected(tier: Tier) -> usize {
    let cost = match tier {
        Tier::Closed => Cost::Closed,
        Tier::Convolved => Cost::Convolved,
        _ => Cost::Grid,
    };
    let axes: &[Axis] = match tier {
        Tier::SideGrid => &[Axis::Side],
        Tier::LevelGrid => &[Axis::Level],
        _ => &[Axis::Level, Axis::Side],
    };
    SPACES
        .iter()
        .map(|&(dimension, base)| {
            let designs = designs(dimension, base)
                .expect("the ledger spaces are walkable")
                .len();
            let rows: usize = Measure::ALL
                .iter()
                .filter(|measure| counted(**measure) && measure.applies(dimension, base))
                .map(|measure| {
                    axes.iter()
                        .filter(|&&axis| measure.cost(axis) == cost)
                        .count()
                })
                .sum();
            designs * rows
        })
        .sum()
}

fn verdict(agree: bool) -> &'static str {
    if agree {
        "PASS"
    } else {
        "FAIL"
    }
}

fn list(values: &[i128]) -> String {
    values
        .iter()
        .map(|value| value.to_string())
        .collect::<Vec<_>>()
        .join(" ")
}

fn stops(batch: &[&Row], stop: Stop) -> usize {
    batch.iter().filter(|row| row.stop == stop).count()
}

fn registry(sheet: &Sheet) -> Vec<String> {
    let mut out = Vec::new();
    for &(tier, count) in &sheet.tiers {
        let batch: Vec<&Row> = sheet.rows.iter().filter(|row| row.tier == tier).collect();
        out.push(format!(
            "registry tier {} rows {} derived {} {} ceiling {} cap {} budget {} silent {}",
            tier.slug(),
            count,
            expected(tier),
            verdict(count == expected(tier)),
            stops(&batch, Stop::Ceiling),
            stops(&batch, Stop::Cap),
            stops(&batch, Stop::Budget),
            batch.iter().filter(|row| row.written.is_empty()).count()
        ));
    }
    let total: usize = sheet.tiers.iter().map(|&(_, count)| count).sum();
    let all: Vec<&Row> = sheet.rows.iter().collect();
    out.push(format!(
        "registry rows {total} rendered {} unread {}",
        sheet.rows.len(),
        sheet.unread
    ));
    out.push(format!(
        "registry stop ceiling {} cap {} budget {} silent {}",
        stops(&all, Stop::Ceiling),
        stops(&all, Stop::Cap),
        stops(&all, Stop::Budget),
        all.iter().filter(|row| row.written.is_empty()).count()
    ));
    out
}

fn checks(sheet: &Sheet, book: &Census) {
    println!("CHECKS");
    let total: u64 = book.counts.iter().map(|&count| count as u64).sum();
    let walked: u64 = sheet.rows.iter().map(|row| row.written.len() as u64).sum();
    println!(
        "checks incidences histogram {total} row walk {walked} {}",
        verdict(total == walked)
    );
    let mut tried = 0u64;
    let mut wrong = 0u64;
    for row in &sheet.rows {
        let Some(form) = &row.form else {
            continue;
        };
        match form {
            Closed::Recurrence(coefficients) => {
                for index in coefficients.len()..row.head.len() {
                    if let Some(value) = rows::replay(coefficients, &row.head, index) {
                        tried += 1;
                        wrong += u64::from(value != row.head[index]);
                    }
                }
            }
            _ => {
                for (index, &term) in row.head.iter().enumerate() {
                    if let Some(value) = rows::predict(form, index) {
                        tried += 1;
                        wrong += u64::from(value != term);
                    }
                }
            }
        }
    }
    println!(
        "checks closed forms against rendered terms {tried} mismatches {wrong} {}",
        verdict(wrong == 0)
    );
    for window in WINDOWS {
        let (never, once, many) = census::split(book, window);
        println!(
            "checks window {window} never {never} once {once} multiple {many} sum {} {}",
            never + once + many,
            verdict(never + once + many == window)
        );
    }
    let classics: [(&str, &[i128]); 5] = [
        (
            "sequence_dim=2_code=7_measure=fills_axis=side",
            &[8, 21, 40, 65],
        ),
        (
            "sequence_dim=2_code=7_measure=fills_axis=level",
            &[8, 64, 512],
        ),
        (
            "sequence_dim=2_code=7_measure=voids_axis=level",
            &[1, 17, 217],
        ),
        (
            "sequence_dim=3_code=23_measure=fills_axis=side",
            &[20, 81, 208, 425],
        ),
        (
            "sequence_dim=3_code=23_measure=surface_axis=level",
            &[72, 1056, 18048],
        ),
    ];
    for (name, head) in classics {
        let row = sheet
            .rows
            .iter()
            .find(|row| row.name == name)
            .expect("the classic is a registry row");
        let seen = row.head.len() >= head.len() && row.head[..head.len()] == *head;
        let written = head
            .iter()
            .all(|&term| term > CEILING || book.counts[term as usize] > 0);
        println!(
            "checks classic {name} head {} written {} {}",
            verdict(seen),
            verdict(written),
            list(&row.head)
        );
    }
    let surface = sheet
        .rows
        .iter()
        .find(|row| row.name == "sequence_dim=3_code=23_measure=surface_axis=level")
        .expect("the sponge surface is a registry row");
    let obeys = (2..surface.head.len()).all(|index| {
        surface.head[index] == 28 * surface.head[index - 1] - 160 * surface.head[index - 2]
    });
    println!(
        "checks recurrence a(L) = 28 a(L-1) - 160 a(L-2) on {} terms {}",
        surface.head.len(),
        verdict(obeys)
    );
    let (top, count) = census::champions(book, 1)[0];
    let walk = census::writers(sheet, top as i128).len();
    println!(
        "checks champion {top} histogram {count} row walk {walk} {}",
        verdict(count as usize == walk)
    );
    println!(
        "checks incidences rows {} pairs {} terms at or below zero {}",
        book.incidences,
        book.incidences + book.repeats,
        book.low
    );
}

fn tables(book: &Census) {
    println!("WINDOWS");
    for window in WINDOWS {
        let (never, once, many) = census::split(book, window);
        println!(
            "window {window} never {never} once {once} multiple {many} written {} share {:.4}",
            once + many,
            (once + many) as f64 / window as f64
        );
    }
    println!();
    println!("CHAMPIONS");
    for (rank, (value, count)) in census::champions(book, CHAMPIONS).into_iter().enumerate() {
        println!("champion {} {value} rows {count}", rank + 1);
    }
    println!();
    println!("MISSES");
    let first = census::misses(book, MISSES);
    println!(
        "misses first {MISSES} {}",
        first
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(" ")
    );
    println!();
    println!("DECADES");
    for band in census::bands(book) {
        println!(
            "decade {}..{} width {} missed {} density {:.6}",
            band.first,
            band.last,
            band.width(),
            band.missed,
            band.density()
        );
    }
}

fn emit(out: &Path, sheet: &Sheet, book: &Census, lines: &[String]) {
    std::fs::create_dir_all(out).expect("the output directory is writable");
    let table: Vec<Vec<String>> = sheet
        .rows
        .iter()
        .map(|row| {
            vec![
                row.name.clone(),
                row.tier.slug().to_string(),
                row.form
                    .as_ref()
                    .map_or_else(|| "none".to_string(), |form| form.text()),
                list(&row.head),
                list(&row.written),
                row.stop.slug().to_string(),
            ]
        })
        .collect();
    tables::write_csv(
        &out.join("rows.csv"),
        &["key", "tier", "closed", "head", "written", "stop"],
        &table,
    );
    let multiset: Vec<Vec<String>> = (1..=CEILING as usize)
        .map(|value| vec![value.to_string(), book.counts[value].to_string()])
        .collect();
    tables::write_csv(&out.join("multiset.csv"), &["integer", "rows"], &multiset);
    let mut page = vec![
        "# Integer Census Manifest".to_string(),
        String::new(),
        "## DEFINITION".to_string(),
        String::new(),
    ];
    page.extend(definition().iter().map(|line| format!("- {line}")));
    page.push(String::new());
    page.push("## REGISTRY".to_string());
    page.push(String::new());
    page.extend(lines.iter().map(|line| format!("- {line}")));
    page.push(String::new());
    page.push("## CENSUS".to_string());
    page.push(String::new());
    for window in WINDOWS {
        let (never, once, many) = census::split(book, window);
        page.push(format!(
            "- window `1..={window}`: never {never}, once {once}, multiple {many}."
        ));
    }
    page.push(format!(
        "- incidences: {} `(row, integer)` pairs, {} `(row, index, integer)` pairs, {} terms at or below zero.",
        book.incidences,
        book.incidences + book.repeats,
        book.low
    ));
    page.push(String::new());
    page.push("## FILES".to_string());
    page.push(String::new());
    page.push(format!(
        "- `rows.csv`: one line a registry row, {} lines and a header.",
        sheet.rows.len()
    ));
    page.push(format!("- `rows.csv` fields: `key`, `tier`, `closed` (the closed form or `none`), `head` (the first {TERMS} rendered terms, space separated), `written` (the distinct terms of the rendered window inside `1..={CEILING}`, ascending, space separated, empty when the row writes none), `stop` (`ceiling`, `cap` or `budget`)."));
    page.push(format!("- `multiset.csv`: `integer,rows` for every `n` in `1..={CEILING}`, {} lines and a header, `rows` the row multiplicity.", CEILING));
    page.push("- `multiset.csv` is the fold of the `written` column of `rows.csv`, so `rows.csv` rebuilds it; the fold is not invertible the other way.".to_string());
    tables::write_lines(&out.join("MANIFEST.md"), &page);
}

fn main() {
    let out = env::args()
        .nth(1)
        .map(PathBuf::from)
        .unwrap_or_else(|| env::temp_dir().join("integer-census"));
    println!("DEFINITION");
    for line in definition() {
        println!("{line}");
    }
    println!();
    let clock = Instant::now();
    let sheet = rows::read();
    let walk = clock.elapsed().as_secs_f64();
    let book = census::build(&sheet);
    println!("REGISTRY");
    let lines = registry(&sheet);
    for line in &lines {
        println!("{line}");
    }
    println!("registry walk {walk:.1}s on {} threads", rows::THREADS);
    println!();
    tables(&book);
    println!();
    hunt::report(&sheet, &book);
    println!();
    checks(&sheet, &book);
    println!();
    emit(&out, &sheet, &book, &lines);
    println!("FILES");
    println!(
        "files rows.csv multiset.csv MANIFEST.md in {}",
        out.display()
    );
    println!("files run {:.1}s", clock.elapsed().as_secs_f64());
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_census_reads_its_registry() {
        let mut total = 0;
        for tier in rows::TIERS {
            let count = ledger::keys(tier)
                .into_iter()
                .filter(|key| counted(key.measure))
                .count();
            assert_eq!(count, expected(tier), "{}", tier.slug());
            total += count;
        }
        assert_eq!(total, 18066);
    }
}
