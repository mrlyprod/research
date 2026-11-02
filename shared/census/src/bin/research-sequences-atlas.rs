use census::save;
use mrly_core::error::{value_error, Result};
use mrly_core::{json, Json};
use mrly_core::Error;
use std::collections::{BTreeMap, BTreeSet};
use std::path::PathBuf;

const NAME: &str = "research-sequences-atlas";

fn atlas() -> Result<Json> {
    let file = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../data/research/shared/sequences/atlas.json");
    let text = std::fs::read_to_string(&file)
        .map_err(|e| Error::Value(format!("cannot read {file:?}: {e}")))?;
    text.parse::<Json>()
        .map_err(|e| Error::Value(format!("cannot parse {file:?}: {e}")))
}

fn growth(terms: &[Json]) -> f64 {
    let top = terms
        .iter()
        .filter_map(Json::as_f64)
        .map(f64::abs)
        .fold(0.0, f64::max);
    (1.0 + top).ln() / terms.len().max(1) as f64
}

fn state(class: &Json) -> u8 {
    let names: BTreeSet<&str> = class["names"]
        .as_array()
        .map(|names| names.iter().filter_map(Json::as_str).collect())
        .unwrap_or_default();
    if class["searched"] != Json::Bool(true) {
        3
    } else if names.len() > 1 {
        2
    } else if class["absent"] == Json::Bool(true) {
        0
    } else {
        1
    }
}

fn slice(code: u32, axis: usize, parity: u32) -> i128 {
    let cells: Vec<u32> = (0..8u32)
        .filter(|&i| (code >> i) & 1 == 1 && (i >> axis) & 1 == parity)
        .collect();
    let rest = |i: u32| -> Vec<u32> {
        (0..3)
            .filter(|&a| a != axis)
            .map(|a| (i >> a) & 1)
            .collect()
    };
    match cells.len() {
        4 => 1,
        2 if rest(cells[0])
            .iter()
            .zip(rest(cells[1]))
            .all(|(a, b)| *a != b) =>
        {
            -1
        }
        _ => 0,
    }
}

fn predicted(code: u32, side: i128) -> i128 {
    (0..3)
        .map(|axis| (side + 1) / 2 * slice(code, axis, 0) + (side - 1) / 2 * slice(code, axis, 1))
        .sum::<i128>()
        * (side - 1)
        * (side - 1)
}

fn edges_against_surface(rows: &[Json]) -> Result<()> {
    let mut found: BTreeMap<(u32, &str), Vec<i128>> = BTreeMap::new();
    for row in rows {
        if row["generator"] == Json::from("the ledger page") {
            continue;
        }
        let Some(name) = row["name"].as_str() else {
            continue;
        };
        let Some(rest) = name.strip_prefix("sequence_dim=3_code=") else {
            continue;
        };
        let Some((code, measure)) = rest.split_once("_measure=") else {
            continue;
        };
        let measure = match measure {
            "edges_axis=side" => "edges",
            "surface_axis=side" => "surface",
            _ => continue,
        };
        let Ok(code) = code.parse::<u32>() else {
            continue;
        };
        let terms = row["terms"]
            .as_array()
            .map(|t| t.iter().filter_map(Json::as_i64).map(i128::from).collect())
            .unwrap_or_default();
        found.insert((code, measure), terms);
    }
    let (mut designs, mut twice, mut rule, mut terms) = (0, 0, 0, 0);
    for (&(code, measure), edges) in &found {
        if measure != "edges" {
            continue;
        }
        let Some(surface) = found.get(&(code, "surface")) else {
            continue;
        };
        designs += 1;
        let count = edges.len().min(surface.len());
        let side = |k: usize| 2 * k as i128 + 3;
        for k in 0..count {
            if edges[k] - 2 * surface[k] != predicted(code, side(k)) {
                return value_error(format!(
                    "code {code} at side {}: E - 2S = {}, the slice rule says {}.",
                    side(k),
                    edges[k] - 2 * surface[k],
                    predicted(code, side(k))
                ));
            }
        }
        terms += count;
        if (0..count).all(|k| edges[k] == 2 * surface[k]) {
            twice += 1;
        }
        if (0..count).all(|k| predicted(code, side(k)) == 0) {
            rule += 1;
        }
    }
    println!(
        "edges against surface on the side axis: {designs} designs of dim 3 base 2, {terms} terms, E - 2S as the slice rule says on every term; E = 2S for {twice}, the rule predicts {rule}"
    );
    Ok(())
}

fn measures(rows: &[Json]) -> BTreeMap<(String, String, String), Vec<i128>> {
    let mut out = BTreeMap::new();
    for row in rows {
        if row["generator"] == Json::from("the ledger page") {
            continue;
        }
        let Some(name) = row["name"].as_str() else {
            continue;
        };
        let Some((design, rest)) = name.split_once("_measure=") else {
            continue;
        };
        let Some((measure, axis)) = rest.split_once("_axis=") else {
            continue;
        };
        if !design.starts_with("sequence_dim=") {
            continue;
        }
        let terms = row["terms"]
            .as_array()
            .map(|t| t.iter().filter_map(Json::as_i64).map(i128::from).collect())
            .unwrap_or_default();
        out.insert(
            (design.to_string(), measure.to_string(), axis.to_string()),
            terms,
        );
    }
    out
}

fn double_counts(rows: &[Json]) -> Result<()> {
    let found = measures(rows);
    let laws: [(&str, &str, fn(i128, i128) -> Option<i128>); 2] = [
        ("sequence_dim=2_", "edges", |f, s| {
            (s % 2 == 0).then_some(2 * f + s / 2)
        }),
        ("sequence_dim=3_", "faces", |f, s| {
            ((6 * f + s) % 2 == 0).then_some((6 * f + s) / 2)
        }),
    ];
    for (prefix, measure, law) in laws {
        for axis in ["level", "side"] {
            let (mut designs, mut terms) = (0, 0);
            for ((design, m, a), counted) in &found {
                if !design.starts_with(prefix) || m != measure || a != axis {
                    continue;
                }
                let key = |m: &str| (design.clone(), m.to_string(), a.clone());
                let (Some(fills), Some(surface)) =
                    (found.get(&key("fills")), found.get(&key("surface")))
                else {
                    continue;
                };
                let count = counted.len().min(fills.len()).min(surface.len());
                for k in 0..count {
                    if law(fills[k], surface[k]) != Some(counted[k]) {
                        return value_error(format!(
                            "{design} {measure}/{axis} breaks its double count at term {k}."
                        ));
                    }
                }
                designs += 1;
                terms += count;
            }
            println!("{measure}/{axis} by the double count: {designs} designs, {terms} terms, every one agrees");
        }
    }
    Ok(())
}

fn main() -> Result<()> {
    let atlas = atlas()?;
    let Some(classes) = atlas["classes"].as_array() else {
        return value_error("the atlas holds no classes.");
    };
    let Some(rows) = atlas["rows"].as_array() else {
        return value_error("the atlas holds no rows.");
    };
    edges_against_surface(rows)?;
    double_counts(rows)?;
    let mut order: Vec<(f64, usize, u8)> = classes
        .iter()
        .enumerate()
        .map(|(i, class)| {
            let terms = class["terms"].as_array().map(Vec::as_slice).unwrap_or(&[]);
            (growth(terms), i, state(class))
        })
        .collect();
    order.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    let side = (order.len() as f64).sqrt().ceil() as usize;
    let cells: String = order
        .iter()
        .map(|&(_, _, s)| char::from(b'0' + s))
        .collect();
    let tally: Vec<usize> = (0..4u8)
        .map(|s| order.iter().filter(|o| o.2 == s).count())
        .collect();
    println!(
        "classes {} side {side} absent {} found {} named twice {} unsearched {}",
        order.len(),
        tally[0],
        tally[1],
        tally[2],
        tally[3]
    );
    save(NAME, &json!({ "side": side, "cells": cells }))
}
