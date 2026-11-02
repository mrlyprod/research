use crate::ledger::Book;
use crate::registry::Verdict;
use std::collections::HashSet;

pub struct Run {
    pub key: &'static str,
    pub outcome: Option<Result<(), String>>,
}

#[derive(Default)]
pub struct Sheet {
    rows: usize,
    parts: usize,
    untagged: usize,
    witnessless: usize,
    tags: [usize; 4],
    counts: [usize; 5],
    failures: Vec<(String, String)>,
    orphans: Vec<String>,
    missing: Vec<String>,
}

impl Sheet {
    pub fn count(&self, verdict: Verdict) -> usize {
        self.counts[verdict as usize]
    }

    fn record(&mut self, verdict: Verdict, key: &str, why: &str) {
        self.counts[verdict as usize] += 1;
        match verdict {
            Verdict::Failed => self.failures.push((key.to_string(), why.to_string())),
            Verdict::Orphan => self.orphans.push(key.to_string()),
            Verdict::Unchecked => self.missing.push(key.to_string()),
            _ => {}
        }
    }
}

// TALLY

pub fn tally(book: &Book, runs: &[Run]) -> Sheet {
    let mut sheet = Sheet {
        rows: book.rows.len(),
        parts: book.parts,
        untagged: book.untagged,
        witnessless: book
            .rows
            .iter()
            .filter(|row| row.witnesses.is_empty())
            .count(),
        ..Sheet::default()
    };
    for row in &book.rows {
        sheet.tags[row.tag as usize] += 1;
    }
    let claimed: HashSet<&str> = book.rows.iter().map(|row| row.key.as_str()).collect();
    let run: HashSet<&str> = runs.iter().map(|entry| entry.key).collect();
    for entry in runs {
        let (verdict, why) = match (claimed.contains(entry.key), &entry.outcome) {
            (false, _) => (Verdict::Orphan, String::new()),
            (true, None) => (Verdict::Skipped, String::new()),
            (true, Some(Ok(()))) => (Verdict::Checked, String::new()),
            (true, Some(Err(why))) => (Verdict::Failed, why.clone()),
        };
        sheet.record(verdict, entry.key, &why);
    }
    for row in &book.rows {
        if !run.contains(row.key.as_str()) {
            sheet.record(Verdict::Unchecked, &row.key, "");
        }
    }
    sheet
}

// DASHBOARD

pub fn print(sheet: &Sheet) {
    println!(
        "LEDGER {} rows in {} parts, {} untagged, {} witnessless",
        sheet.rows, sheet.parts, sheet.untagged, sheet.witnessless
    );
    println!(
        "TAGS {} proved, {} verified, {} conjecture, {} refuted",
        sheet.tags[0], sheet.tags[1], sheet.tags[2], sheet.tags[3]
    );
    println!(
        "CLAIMS {} checked, {} failed, {} unchecked, {} skipped, {} orphan",
        sheet.count(Verdict::Checked),
        sheet.count(Verdict::Failed),
        sheet.count(Verdict::Unchecked),
        sheet.count(Verdict::Skipped),
        sheet.count(Verdict::Orphan)
    );
    for (key, why) in &sheet.failures {
        println!("FAILED {key}: {why}");
    }
    for key in &sheet.orphans {
        println!("ORPHAN {key}");
    }
    for key in sheet.missing.iter().take(10) {
        println!("UNCHECKED {key}");
    }
}
