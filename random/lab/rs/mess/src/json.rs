use crate::census::Entry;
use crate::frontier::Front;
use std::fs;
use std::path::PathBuf;

pub fn data_dir() -> PathBuf {
    let here = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let repo = here
        .ancestors()
        .nth(5)
        .expect("the study sits five levels under the repo");
    repo.join("data/research/random/lab/rs/mess")
}

fn list(values: impl Iterator<Item = String>) -> String {
    format!("[{}]", values.collect::<Vec<_>>().join(","))
}

pub fn write(entries: &[Entry]) {
    let body: Vec<String> = entries
        .iter()
        .map(|e| {
            let rows = list(e.rows.iter().map(|p| {
                list(
                    [
                        p.k as i128,
                        p.degree as i128,
                        p.ones as i128,
                        p.nonlinearity as i128,
                        p.walsh_max as i128,
                        p.autocorrelation_max as i128,
                        p.differential as i128,
                        p.immunity as i128,
                    ]
                    .iter()
                    .map(|x| x.to_string()),
                )
            }));
            format!(
                "{{\"name\":\"{}\",\"oeis\":{},\"g\":{},\"c\":{},\"rows\":{}}}",
                e.name,
                e.oeis
                    .as_ref()
                    .map_or("null".to_string(), |id| format!("\"{id}\"")),
                list(e.g.iter().map(|x| x.to_string())),
                list(e.c.iter().map(|x| x.to_string())),
                rows
            )
        })
        .collect();
    let dir = data_dir();
    fs::create_dir_all(&dir).expect("the data folder is made");
    fs::write(dir.join("census.json"), format!("[{}]", body.join(",\n")))
        .expect("the census is written");
    println!("census.json written, {} sources, columns k degree ones nl walsh_max autocorrelation_max differential immunity", entries.len());
}

pub fn write_frontier(fronts: &[Front]) {
    let body: Vec<String> = fronts
        .iter()
        .map(|front| {
            let rows = list(front.rows.iter().map(|r| {
                list(
                    [
                        r.count,
                        r.best_walsh,
                        r.best_count,
                        r.balanced_count,
                        r.balanced_best_walsh,
                        r.balanced_best_count,
                    ]
                    .iter()
                    .map(|x| x.to_string()),
                )
            }));
            format!("{{\"k\":{},\"rows\":{}}}", front.k, rows)
        })
        .collect();
    let dir = data_dir();
    fs::create_dir_all(&dir).expect("the data folder is made");
    fs::write(dir.join("frontier.json"), format!("[{}]", body.join(",\n")))
        .expect("the frontier is written");
    println!("frontier.json written, k = 1..{}, columns count best_walsh best_count balanced_count balanced_best_walsh balanced_best_count per degree 0..k", fronts.len());
}
