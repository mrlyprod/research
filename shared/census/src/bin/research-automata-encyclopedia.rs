use census::save;
use mrly_core::error::{value_error, Result};
use mrly_core::{json, Json};
use mrly_core::Error;
use std::path::PathBuf;

const NAME: &str = "research-automata-encyclopedia";
const STOPS: [&str; 6] = ["dead", "still", "loop", "mover", "wall", "limit"];

fn tally() -> Result<Json> {
    let file = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../../data/research/automata/lab/rs/encyclopedia/tally.json");
    let text = std::fs::read_to_string(&file)
        .map_err(|e| Error::Value(format!("cannot read {file:?}: {e}")))?;
    text.parse::<Json>()
        .map_err(|e| Error::Value(format!("cannot parse {file:?}: {e}")))
}

fn names(value: &Json) -> Vec<String> {
    value
        .as_array()
        .map(|a| {
            a.iter()
                .filter_map(|v| v.as_str().map(str::to_string))
                .collect()
        })
        .unwrap_or_default()
}

fn shares(cell: &Json) -> Result<[u64; 3]> {
    let counts: Vec<u64> = cell
        .as_array()
        .map(|a| a.iter().filter_map(Json::as_u64).collect())
        .unwrap_or_default();
    if counts.len() != STOPS.len() {
        return value_error(format!(
            "a tally cell holds {} stops, want {}.",
            counts.len(),
            STOPS.len()
        ));
    }
    Ok([
        counts[0],
        counts[1] + counts[2] + counts[3],
        counts[4] + counts[5],
    ])
}

fn main() -> Result<()> {
    let tally = tally()?;
    if names(&tally["stops"]) != STOPS {
        return value_error("the tally stops are not the six the figure reads.");
    }
    let masks = names(&tally["masks"]);
    let rules = names(&tally["rules"]);
    let mut grid: Vec<Vec<[u64; 3]>> = Vec::new();
    for row in tally["grid"]
        .as_array()
        .map(Vec::as_slice)
        .unwrap_or_default()
    {
        let cells = row.as_array().map(Vec::as_slice).unwrap_or_default();
        if cells.len() != masks.len() {
            return value_error("a tally row does not span the masks.");
        }
        grid.push(cells.iter().map(shares).collect::<Result<Vec<_>>>()?);
    }
    if grid.len() != rules.len() || grid.is_empty() {
        return value_error("the tally grid does not span the rules.");
    }
    let total: u64 = grid.iter().flatten().flatten().sum();
    save(
        NAME,
        &json!({
            "rows": tally["rows"],
            "total": tally["total"],
            "cut": tally["cut"],
            "shares": ["dead", "periodic", "open"],
            "masks": masks,
            "rules": rules,
            "grid": grid,
            "counted": total,
        }),
    )
}
