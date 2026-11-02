use super::records::record_by_id;
use super::{
    keys, sequence, Axis, Key, Measure, Record, Row, Tier, BUDGET, RECORDS, SPACES, TERMS,
};
use mrly_core::name::Named;
use mrly_math::name::Sequence as SequenceName;
use std::collections::{BTreeMap, BTreeSet, HashSet};
use std::fs;
use std::path::Path;

const PROSE: &str = include_str!("prose.md");

fn label(key: &Key) -> &'static str {
    match (key.dimension, key.code) {
        (2, 1) | (3, 1) => "dust",
        (2, 3) => "tree",
        (2, 6) | (3, 22) => "star",
        (2, 7) => "carpet",
        (2, 8) | (3, 128) => "point",
        (2, 9) | (3, 129) => "void",
        (2, 10) | (2, 12) | (3, 136) | (3, 160) | (3, 192) => "line",
        (2, 14) | (3, 232) => "net",
        (2, 11) => "corner and centre",
        (2, 15) | (3, 255) => "solid",
        (3, 23) => "sponge",
        _ => "design",
    }
}

fn link(id: &str) -> String {
    format!("[{id}](https://oeis.org/{id})")
}

fn spelled(row: &Row) -> String {
    let terms: Vec<String> = row
        .sequence
        .terms
        .iter()
        .map(|term| term.to_string())
        .collect();
    let mut text = format!("`{}`", terms.join(", "));
    if row.sequence.capped {
        text.push_str(" to the cell budget");
    }
    text
}

fn table(axis: Axis) -> String {
    let index = match axis {
        Axis::Level => "level 1",
        Axis::Side => "`k = 2`",
        Axis::Window => "`k = 1`",
    };
    let mut out = format!("| id | design | key | closed form | terms from {index} | record | shift | status |\n|---|---|---|---|---|---|---|---|\n");
    let mut keyed: Vec<(Key, &Record)> = RECORDS
        .iter()
        .filter_map(|record| {
            record
                .key
                .filter(|key| key.axis == axis)
                .map(|key| (key, record))
        })
        .collect();
    keyed.sort_by_key(|(key, _)| (key.dimension, key.code, key.measure));
    for (key, record) in keyed {
        let row = sequence(&key, TERMS, BUDGET).expect("every keyed record reads");
        let (matched, shift) = row.record.expect("every keyed record matches");
        assert_eq!(matched.id, record.id, "the row reads its own record");
        let closed = row
            .sequence
            .closed
            .as_ref()
            .map_or_else(|| "none".to_string(), |form| format!("`{}`", form.text()));
        out.push_str(&format!(
            "| {} | {} | `{}` | {} | {} | {} | {} | **{}** |\n",
            key.id(),
            label(&key),
            key.name(),
            closed,
            spelled(&row),
            link(record.id),
            shift,
            record.status.text()
        ));
    }
    out.trim_end().to_string()
}

fn tokens(text: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut rest = text;
    while let Some(at) = rest.find("sequence_") {
        let tail = &rest[at..];
        let end = tail
            .find(|c: char| !(c.is_ascii_lowercase() || c.is_ascii_digit() || c == '_' || c == '='))
            .unwrap_or(tail.len());
        out.push(tail[..end].to_string());
        rest = &tail[end..];
    }
    out
}

fn citations() -> BTreeMap<String, BTreeSet<String>> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let mut topics: Vec<_> = fs::read_dir(&root)
        .map(|entries| entries.flatten().map(|entry| entry.path()).collect())
        .unwrap_or_default();
    topics.retain(|path| path.is_dir() && path.file_name().is_some_and(|name| name != "shared"));
    topics.sort();
    let mut names: BTreeMap<String, BTreeSet<String>> = BTreeMap::new();
    for folder in ["claims", "notes"] {
        let mut paths: Vec<_> = Vec::new();
        for topic in &topics {
            let Ok(entries) = fs::read_dir(topic.join(folder)) else {
                continue;
            };
            paths.extend(
                entries
                    .flatten()
                    .map(|entry| entry.path())
                    .filter(|path| path.extension().is_some_and(|ext| ext == "md"))
                    .filter(|path| path.file_name().is_some_and(|name| name != "README.md")),
            );
        }
        paths.sort();
        for path in paths {
            let Ok(text) = fs::read_to_string(&path) else {
                continue;
            };
            let stem = path
                .file_stem()
                .and_then(|s| s.to_str())
                .unwrap_or_default();
            for token in tokens(&text) {
                names
                    .entry(token)
                    .or_default()
                    .insert(format!("{folder}/{stem}"));
            }
        }
    }
    names
}

fn cited() -> String {
    let keyed: HashSet<String> = RECORDS
        .iter()
        .filter_map(|record| record.key.map(|key| key.name()))
        .collect();
    let mut out =
        String::from("| id | key | closed form | terms | cited by |\n|---|---|---|---|---|\n");
    for (name, by) in citations() {
        if keyed.contains(&name) {
            continue;
        }
        let named = SequenceName::from_file(&name).expect("every cited sequence name reads");
        let key = Key::new(
            named.code,
            named.dim,
            named.base,
            Measure::parse(&named.measure).expect("a cited measure is a measure"),
            Axis::parse(&named.axis).expect("a cited axis is an axis"),
        );
        let row = sequence(&key, TERMS, BUDGET).expect("every cited sequence reads");
        let closed = row
            .sequence
            .closed
            .as_ref()
            .map_or_else(|| "none".to_string(), |form| format!("`{}`", form.text()));
        let files: Vec<String> = by.into_iter().map(|file| format!("`{file}`")).collect();
        out.push_str(&format!(
            "| `{}` | `{}` | {} | {} | {} |\n",
            key.id(),
            key.name(),
            closed,
            spelled(&row),
            files.join(", ")
        ));
    }
    out.trim_end().to_string()
}

fn prose(text: &str) -> String {
    text.chars()
        .map(|c| match c {
            '*' | '_' | '`' | '[' | ']' | '|' | '\\' => format!("\\{c}"),
            _ => c.to_string(),
        })
        .collect()
}

fn records() -> String {
    let mut out = String::from(
        "| id | record | name | offset | first terms | key | shift | status |\n|---|---|---|---|---|---|---|---|\n",
    );
    for record in RECORDS {
        let (id, key, shift) = match record.key {
            Some(key) => (
                key.id(),
                format!("`{}`", key.name()),
                record.shift.to_string(),
            ),
            None => (String::new(), String::new(), String::new()),
        };
        out.push_str(&format!(
            "| {} | {} | {} | {} | `{}` | {} | {} | **{}** |\n",
            id,
            link(record.id),
            prose(record.name),
            record.offset,
            record.terms,
            key,
            shift,
            record.status.text()
        ));
    }
    out.trim_end().to_string()
}

fn tally() -> String {
    let designs: usize = SPACES
        .iter()
        .map(|&(dimension, base)| {
            super::designs(dimension, base)
                .expect("the ledger spaces are walkable")
                .len()
        })
        .sum();
    format!(
        "The registry walks {designs} designs across {} dimension and base pairs and holds {} closed rows and {} convolved rows of {TERMS} terms each; the grid tiers render on demand within a budget of {BUDGET} cells a term.",
        SPACES.len(),
        keys(Tier::Closed).len(),
        keys(Tier::Convolved).len()
    )
}

fn fill(tag: &str) -> String {
    let record = |id: &str| record_by_id(id).expect("the placeholder names a record");
    match tag.split_once(' ') {
        Some(("terms", id)) => format!("`{}`", record(id).terms),
        Some(("formula", id)) => format!("`{}`, offset {}", record(id).formula, record(id).offset),
        _ => match tag {
            "sides" => table(Axis::Side),
            "levels" => table(Axis::Level),
            "cited" => cited(),
            "records" => records(),
            "tally" => tally(),
            other => panic!("unknown placeholder {other}"),
        },
    }
}

/// Renders the sequences page: the prose with every generated table and number filled in.
pub fn markdown() -> String {
    let mut out = String::new();
    let mut rest = PROSE;
    while let Some(open) = rest.find("{{") {
        out.push_str(&rest[..open]);
        let close = open + rest[open..].find("}}").expect("every placeholder closes");
        out.push_str(&fill(rest[open + 2..close].trim()));
        rest = &rest[close + 2..];
    }
    out.push_str(rest);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::Path;

    #[test]
    fn the_page_on_disk_is_the_rendered_page() {
        let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../sequences.md");
        let disk = fs::read_to_string(&path).expect("the page is on disk");
        let page = markdown();
        assert!(!page.contains("{{"));
        assert!(page.contains(
            "| 8a9e4ce8 | carpet | `sequence_dim=2_code=7_measure=fills_axis=side` | `3k^2 - 2k` | `8, 21, 40, 65, 96, 133, 176, 225`"
        ));
        assert_eq!(
            page, disk,
            "run `cargo run --manifest-path research/shared/ledger/Cargo.toml`"
        );
    }
}
