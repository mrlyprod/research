use std::collections::HashSet;
use std::fs;
use std::io::Write;
use std::path::PathBuf;

pub fn data_dir() -> PathBuf {
    let here = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let repo = here.ancestors().nth(4).unwrap().to_path_buf();
    repo.join("data").join(here.strip_prefix(&repo).unwrap())
}

pub fn chunks_dir() -> PathBuf {
    data_dir().join("chunks")
}

pub fn chunk_file(from: usize, to: usize) -> PathBuf {
    chunks_dir().join(format!("chunk-{from}-{to}.jsonl"))
}

pub fn key_of(line: &str) -> Option<&str> {
    let rest = line.strip_prefix("{\"key\":\"")?;
    let end = rest.find('"')?;
    Some(&rest[..end])
}

pub fn done_keys() -> HashSet<String> {
    let mut keys = HashSet::new();
    let Ok(entries) = fs::read_dir(chunks_dir()) else {
        return keys;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.extension().is_some_and(|e| e == "jsonl") {
            if let Ok(text) = fs::read_to_string(&path) {
                keys.extend(text.lines().filter_map(key_of).map(str::to_string));
            }
        }
    }
    keys
}

pub fn free_percent() -> f64 {
    let out = std::process::Command::new("df")
        .arg("-k")
        .arg(data_dir().ancestors().nth(4).unwrap())
        .output();
    let Ok(out) = out else {
        return 0.0;
    };
    let text = String::from_utf8_lossy(&out.stdout);
    let Some(row) = text.lines().nth(1) else {
        return 0.0;
    };
    let words: Vec<&str> = row.split_whitespace().collect();
    let (Some(used), Some(free)) = (
        words.get(2).and_then(|w| w.parse::<f64>().ok()),
        words.get(3).and_then(|w| w.parse::<f64>().ok()),
    ) else {
        return 0.0;
    };
    100.0 * free / (used + free)
}

pub fn line(key: &str, index: usize, readings: &mrly_core::Json) -> String {
    format!("{{\"key\":\"{key}\",\"index\":{index},\"readings\":{readings}}}")
}

pub fn index_of(line: &str) -> Option<usize> {
    let rest = line.split_once("\",\"index\":")?.1;
    let end = rest.find(',')?;
    rest[..end].parse().ok()
}

pub fn append(path: &PathBuf, lines: &[String]) -> std::io::Result<()> {
    fs::create_dir_all(path.parent().unwrap())?;
    let mut file = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(path)?;
    for line in lines {
        writeln!(file, "{line}")?;
    }
    file.flush()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_cache_line_starts_with_its_key() {
        let line = line("carpet.3.1 x1 island", 7, &mrly_core::json!({ "a": 1 }));
        assert_eq!(index_of(&line), Some(7));
        assert_eq!(key_of(&line), Some("carpet.3.1 x1 island"));
        assert!(data_dir().ends_with("data/automata/lab/rs/encyclopedia"));
    }
}
