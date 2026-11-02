use crate::cache::data_dir;
use crate::laws::class;
use crate::rows::{each_row, Curve, Row};
use mrly_core::{json, Json};
use mrly_core::Result;
use std::collections::BTreeMap;
use std::fs;

pub const PER_CELL: usize = 8;
pub const MOST: usize = 3000;

fn curve_json(c: &Curve) -> Json {
    json!({ "first": c.first, "every": c.every, "last": c.last, "max": [c.max, c.max_gen] })
}

fn row_json(row: &Row, fate: &str, class: &str, family: &str) -> Json {
    let mut out = json!({
        "key": row.name,
        "fate": fate,
        "class": class,
        "family": family,
        "gen": row.gen,
        "growth": row.growth,
        "entropy": row.entropy,
        "churn": row.churn,
        "symmetry": row.symmetry,
        "heat": [row.visited, row.peak],
        "pop": curve_json(&row.pop),
        "box": curve_json(&row.bx),
    });
    let map = out.as_object_mut().unwrap();
    if let Some(p) = row.period {
        map.insert("period".into(), json!(p));
    }
    if let Some((dr, dc)) = row.d {
        map.insert("d".into(), json!([dr, dc]));
    }
    if let Some(w) = row.wall_by {
        map.insert("wall_by".into(), json!(w));
    }
    if let Some(s) = row.settled {
        map.insert("settled".into(), json!(s));
    }
    if !row.cuts.is_empty() {
        map.insert("cuts".into(), json!(row.cuts));
    }
    if !row.tilings.is_empty() {
        map.insert("tilings".into(), json!(row.tilings));
    }
    if row.flips {
        map.insert("flips".into(), json!(true));
        map.insert("background".into(), json!(row.background));
    }
    out
}

pub fn write() -> Result<()> {
    let mut shelves: BTreeMap<(String, String, String), Vec<Row>> = BTreeMap::new();
    each_row(|row| {
        let Ok(recipe) = row.recipe() else {
            return;
        };
        let cell = (row.stop.clone(), class(&recipe), recipe.seed.token());
        let shelf = shelves.entry(cell).or_default();
        if shelf.len() < PER_CELL {
            shelf.push(row);
            shelf.sort_by_key(|r| r.index);
        } else if shelf.last().is_some_and(|last| row.index < last.index) {
            shelf.pop();
            shelf.push(row);
            shelf.sort_by_key(|r| r.index);
        }
    })?;
    let mut rows: Vec<Json> = Vec::new();
    for ((fate, class, family), shelf) in &shelves {
        for row in shelf {
            rows.push(row_json(row, fate, class, family));
        }
    }
    rows.truncate(MOST);
    let value = json!({ "shelves": shelves.len(), "per_cell": PER_CELL, "rows": rows });
    let path = data_dir().join("library.json");
    fs::create_dir_all(data_dir()).map_err(|e| mrly_core::Error::Value(e.to_string()))?;
    fs::write(&path, format!("{value}\n")).map_err(|e| mrly_core::Error::Value(e.to_string()))?;
    println!(
        "library: {} shelf cells (fate, class, family), {} rows at most {PER_CELL} a cell, {} bytes, {}",
        shelves.len(),
        rows.len(),
        value.to_string().len() + 1,
        path.display()
    );
    Ok(())
}
