use crate::cache::{chunks_dir, index_of, key_of};
use crate::recipe::Recipe;
use mrly_core::Json;
use mrly_core::Result;
use std::fs;
use std::io::{BufRead, BufReader};

#[derive(Clone, Debug)]
pub struct Curve {
    pub first: Vec<u32>,
    pub every: Vec<u32>,
    pub last: u32,
    pub max: u32,
    pub max_gen: usize,
}

#[derive(Clone, Debug)]
pub struct Row {
    pub index: usize,
    pub name: String,
    pub stop: String,
    pub gen: usize,
    pub period: Option<usize>,
    pub d: Option<(i64, i64)>,
    pub wall_by: Option<usize>,
    pub settled: Option<usize>,
    pub growth: f64,
    pub entropy: i64,
    pub churn: f64,
    pub symmetry: usize,
    pub cuts: Vec<usize>,
    pub tilings: Vec<usize>,
    pub visited: usize,
    pub peak: u32,
    pub flips: bool,
    pub background: u8,
    pub pop: Curve,
    pub bx: Curve,
}

fn list(json: &Json) -> Vec<u32> {
    json.as_array()
        .map(|a| {
            a.iter()
                .filter_map(|v| v.as_u64().map(|v| v as u32))
                .collect()
        })
        .unwrap_or_default()
}

fn usizes(json: &Json) -> Vec<usize> {
    json.as_array()
        .map(|a| {
            a.iter()
                .filter_map(|v| v.as_u64().map(|v| v as usize))
                .collect()
        })
        .unwrap_or_default()
}

fn curve(json: &Json) -> Option<Curve> {
    Some(Curve {
        first: list(&json[0]),
        every: list(&json[1]),
        last: json[2].as_u64()? as u32,
        max: json[3][0].as_u64()? as u32,
        max_gen: json[3][1].as_u64()? as usize,
    })
}

impl Row {
    pub fn of(line: &str) -> Option<Row> {
        let name = key_of(line)?.to_string();
        let index = index_of(line)?;
        let readings = line.split_once("\"readings\":")?.1;
        let readings = readings.strip_suffix('}').unwrap_or(readings);
        let json: Json = readings.parse().ok()?;
        let heat = usizes(&json["heat"]);
        Some(Row {
            index,
            name,
            stop: json["stop"].as_str()?.to_string(),
            gen: json["gen"].as_u64()? as usize,
            period: json["period"].as_u64().map(|v| v as usize),
            d: json["d"]
                .as_array()
                .and_then(|a| Some((a.first()?.as_i64()?, a.get(1)?.as_i64()?))),
            wall_by: json["wall_by"].as_u64().map(|v| v as usize),
            settled: json["settled"].as_u64().map(|v| v as usize),
            growth: json["growth"].as_f64().unwrap_or(0.0),
            entropy: json["entropy"].as_i64().unwrap_or(0),
            churn: json["churn"].as_f64().unwrap_or(0.0),
            symmetry: json["symmetry"].as_u64().unwrap_or(0) as usize,
            cuts: usizes(&json["cuts"]),
            tilings: usizes(&json["tilings"]),
            visited: heat.first().copied().unwrap_or(0),
            peak: heat.get(1).copied().unwrap_or(0) as u32,
            flips: json["flips"].as_bool().unwrap_or(false),
            background: json["background"].as_u64().unwrap_or(0) as u8,
            pop: curve(&json["pop"])?,
            bx: curve(&json["box"])?,
        })
    }

    pub fn recipe(&self) -> Result<Recipe> {
        Recipe::parse(&self.name)
    }
}

pub fn each_line(mut f: impl FnMut(&str)) -> Result<usize> {
    let mut files: Vec<_> = fs::read_dir(chunks_dir())
        .map_err(|e| mrly_core::Error::Value(format!("no chunks: {e}")))?
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.extension().is_some_and(|e| e == "jsonl"))
        .collect();
    files.sort();
    let mut n = 0;
    for path in files {
        let file = fs::File::open(&path).map_err(|e| mrly_core::Error::Value(e.to_string()))?;
        for line in BufReader::new(file)
            .lines()
            .map_while(std::result::Result::ok)
        {
            f(&line);
            n += 1;
        }
    }
    Ok(n)
}

pub fn each_row(mut f: impl FnMut(Row)) -> Result<usize> {
    each_line(|line| {
        if let Some(row) = Row::of(line) {
            f(row);
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache::line;
    use crate::{run_named, short};

    #[test]
    fn a_short_line_reads_back_as_a_row() {
        let name = "dust.3.1 x1 wallpaper mask:carpet.3.1 birth:evens survive:evens x:256";
        let (recipe, run) = run_named(name).unwrap();
        let text = line(&recipe.name(), 7, &short(&run).unwrap());
        let row = Row::of(&text).unwrap();
        assert_eq!(row.index, 7);
        assert_eq!(row.name, name);
        assert_eq!(row.stop, "still");
        assert_eq!(row.pop.first, vec![4, 5, 5]);
        assert_eq!(row.pop.max, 5);
        assert!(!row.flips);
        assert_eq!(row.recipe().unwrap(), recipe);
    }
}
