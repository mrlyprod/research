use crate::merge::{Class, Relation};
use crate::oeis::{Hit, Search};
use crate::row::Row;

pub fn quote(text: &str) -> String {
    let mut out = String::from("\"");
    for c in text.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn optional(text: Option<&str>) -> String {
    text.map(quote).unwrap_or_else(|| "null".to_string())
}

fn numbers(terms: &[i128]) -> String {
    let inner: Vec<String> = terms.iter().map(i128::to_string).collect();
    format!("[{}]", inner.join(","))
}

fn strings<'a>(items: impl Iterator<Item = &'a str>) -> String {
    let inner: Vec<String> = items.map(quote).collect();
    format!("[{}]", inner.join(","))
}

fn hits(list: &[Hit], search: &Search) -> String {
    let inner: Vec<String> = list
        .iter()
        .map(|h| {
            format!(
                "{{\"id\":{},\"name\":{},\"transform\":{},\"offset\":{},\"window\":{},\"overlap\":{}}}",
                quote(&h.id),
                optional(search.names.get(&h.id).map(String::as_str)),
                quote(h.transform),
                h.offset,
                quote(&h.window),
                h.overlap
            )
        })
        .collect();
    format!("[{}]", inner.join(","))
}

pub fn class_of(classes: &[Class], rows: usize) -> Vec<(usize, i32)> {
    let mut out = vec![(0, 0); rows];
    for (c, class) in classes.iter().enumerate() {
        for &(i, s) in &class.members {
            out[i] = (c, s);
        }
    }
    out
}

pub fn json(rows: &[Row], classes: &[Class], relations: &[Relation], search: &Search) -> String {
    let owner = class_of(classes, rows.len());
    let empty = Vec::new();
    let found = |c: usize| search.hits.get(c).unwrap_or(&empty);
    let row_lines: Vec<String> = rows
        .iter()
        .enumerate()
        .map(|(i, r)| {
            let (c, s) = owner[i];
            format!(
                "{{\"name\":{},\"generator\":{},\"measure\":{},\"terms\":{},\"closed\":{},\"stated\":{},\"class\":{c},\"shift\":{s},\"oeis\":{}}}",
                quote(&r.name),
                quote(&r.generator),
                quote(&r.measure),
                numbers(&r.terms),
                optional(r.closed.as_deref()),
                optional(r.stated.as_deref()),
                hits(found(c), search)
            )
        })
        .collect();
    let class_lines: Vec<String> = classes
        .iter()
        .enumerate()
        .map(|(c, class)| {
            format!(
                "{{\"class\":{c},\"names\":{},\"terms\":{},\"labels\":{},\"explanation\":{},\"searched\":{},\"absent\":{},\"total\":{},\"hits\":{}}}",
                strings(class.members.iter().map(|m| rows[m.0].name.as_str())),
                numbers(&rows[class.rep].terms),
                strings(class.labels.iter().map(String::as_str)),
                optional(class.explanation.as_deref()),
                search.searched.get(c).copied().unwrap_or(false),
                search.searched.get(c).is_some() && search.absent(c),
                search.totals.get(c).copied().unwrap_or(0),
                hits(found(c), search)
            )
        })
        .collect();
    let relation_lines: Vec<String> = relations
        .iter()
        .map(|r| {
            format!(
                "{{\"from\":{},\"transform\":{},\"offset\":{},\"to\":{},\"overlap\":{}}}",
                r.from,
                quote(r.transform),
                r.offset,
                r.to,
                r.overlap
            )
        })
        .collect();
    let failed_lines: Vec<String> = search
        .failed
        .iter()
        .map(|(i, id)| {
            format!(
                "{{\"name\":{},\"id\":{}}}",
                quote(&rows[*i].name),
                quote(id)
            )
        })
        .collect();
    format!(
        "{{\"rows\":[\n{}\n],\"classes\":[\n{}\n],\"relations\":[\n{}\n],\"failed\":[\n{}\n]}}\n",
        row_lines.join(",\n"),
        class_lines.join(",\n"),
        relation_lines.join(",\n"),
        failed_lines.join(",\n")
    )
}

pub fn markdown(rows: &[Row], classes: &[Class], search: &Search) -> String {
    let owner = class_of(classes, rows.len());
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &b| {
        (&rows[a].generator, &rows[a].name).cmp(&(&rows[b].generator, &rows[b].name))
    });
    let mut out = String::from(
        "# Sequence atlas\n\n| generator | name | terms | closed form | OEIS | class |\n| --- | --- | --- | --- | --- | --- |\n",
    );
    for i in order {
        let r = &rows[i];
        let c = owner[i].0;
        let oeis = match search.hits.get(c).and_then(|h| h.first()) {
            Some(h) => format!("{} ({})", h.id, h.transform),
            None if search.searched.get(c).copied().unwrap_or(false) => "none".to_string(),
            None => "short".to_string(),
        };
        let terms: Vec<String> = r.terms.iter().take(8).map(i128::to_string).collect();
        out += &format!(
            "| `{}` | `{}` | {} | {} | {} | {} |\n",
            r.generator,
            r.name,
            terms.join(", "),
            r.closed
                .as_deref()
                .map(|c| format!("`{c}`"))
                .unwrap_or_else(|| "none".to_string()),
            oeis,
            c
        );
    }
    out
}
