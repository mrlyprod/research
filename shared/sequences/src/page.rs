use crate::row::Row;

fn backticked(cell: &str) -> Option<&str> {
    let cell = cell.trim();
    cell.strip_prefix('`')?.strip_suffix('`')
}

fn numbers(text: &str) -> Option<Vec<i128>> {
    let terms: Option<Vec<i128>> = text.split(',').map(|t| t.trim().parse().ok()).collect();
    terms.filter(|t| t.len() >= 3)
}

fn record(cell: &str) -> Option<String> {
    let start = cell.find("[A")? + 1;
    let id = cell.get(start..start + 7)?;
    (id[1..].bytes().all(|b| b.is_ascii_digit())).then(|| id.to_string())
}

fn hex(cell: &str) -> bool {
    cell.len() == 8 && cell.bytes().all(|b| b.is_ascii_hexdigit())
}

pub fn rows(text: &str) -> Vec<Row> {
    let mut out = Vec::new();
    for line in text.lines().filter(|l| l.starts_with('|')) {
        let cells: Vec<&str> = line.split('|').collect();
        let Some(terms) = cells.iter().find_map(|c| backticked(c).and_then(numbers)) else {
            continue;
        };
        let key = cells
            .iter()
            .filter_map(|c| backticked(c))
            .find(|c| c.starts_with("sequence_"));
        let id = cells.iter().find_map(|c| record(c));
        let closed = key.and_then(|_| {
            cells
                .iter()
                .filter_map(|c| backticked(c))
                .find(|c| !c.starts_with("sequence_") && numbers(c).is_none() && !hex(c))
                .map(str::to_string)
        });
        let name = match (key, &id) {
            (Some(key), _) => key.to_string(),
            (None, Some(id)) => format!("record_{id}"),
            (None, None) => continue,
        };
        let generator = if key.is_some() {
            "the ledger page"
        } else {
            "the ledger page, an OEIS record"
        };
        out.push(
            Row::new(name, generator, "ledger", terms)
                .closed(closed)
                .stated(id),
        );
    }
    out
}
