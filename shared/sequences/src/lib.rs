pub mod designs;
pub mod extra;
pub mod merge;
pub mod oeis;
pub mod page;
pub mod row;
pub mod write;

use mrly_core::error::Result;
use mrly_math::bang::baseq;
use row::Row;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

pub const THREADS: usize = 4;

pub fn root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../..")
}

pub fn io<T>(result: std::io::Result<T>) -> Result<T> {
    result.or_else(|e| mrly_core::error::value_error(e.to_string()))
}

pub fn design_rows() -> Result<Vec<Row>> {
    let jobs = designs::jobs()?;
    let groups: Vec<((usize, usize), Vec<Vec<usize>>)> = designs::SPACES
        .iter()
        .map(|&(d, b)| ((d, b), baseq::group(b, d)))
        .collect();
    let next = AtomicUsize::new(0);
    let done: Mutex<Vec<(usize, Result<Vec<Row>>)>> = Mutex::new(Vec::new());
    std::thread::scope(|scope| {
        for _ in 0..THREADS {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                let Some(&(d, b, code)) = jobs.get(i) else {
                    break;
                };
                let group = &groups.iter().find(|g| g.0 == (d, b)).unwrap().1;
                let rows = designs::orbit_rows(d, b, code, group);
                done.lock().unwrap().push((i, rows));
            });
        }
    });
    let mut done = done.into_inner().unwrap();
    done.sort_by_key(|d| d.0);
    let mut out = Vec::new();
    for (_, rows) in done {
        out.extend(rows?);
    }
    Ok(out)
}

pub fn page_rows() -> Result<Vec<Row>> {
    let text = io(std::fs::read_to_string(
        root().join("research/shared/sequences.md"),
    ))?;
    Ok(page::rows(&text))
}

pub fn walk() -> Result<Vec<Row>> {
    let mut rows = design_rows()?;
    rows.extend(extra::all());
    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_walk_pins_its_rows() {
        let rows = walk().unwrap();
        assert_eq!(rows.len(), 99152);
    }

    #[test]
    fn orbit_mates_split_at_base_two_and_merge_under_shift() {
        let rows = designs::orbit_rows(2, 2, 7, &baseq::group(2, 2)).unwrap();
        let fills: Vec<&Row> = rows.iter().filter(|r| r.measure == "fills/level").collect();
        let firsts: Vec<i128> = fills.iter().map(|r| r.terms[0]).collect();
        assert_eq!(firsts, [8, 7, 5]);
        let shifted = Row::new("tail", "test", "test", fills[0].terms[2..].to_vec());
        let all = [fills[0].clone(), shifted];
        let classes = merge::classes(&all);
        assert_eq!(classes.len(), 1);
        assert_eq!(classes[0].members, [(0, 0), (1, 2)]);
    }
}
