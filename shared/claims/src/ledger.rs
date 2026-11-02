use std::collections::HashMap;

#[derive(Clone, Copy)]
pub enum Tag {
    Proved,
    Verified,
    Conjecture,
    Refuted,
}

pub struct Row {
    pub key: String,
    pub tag: Tag,
    pub witnesses: Vec<String>,
}

pub struct Book {
    pub rows: Vec<Row>,
    pub parts: usize,
    pub untagged: usize,
}

const TAGS: [(&str, Tag); 4] = [
    ("[Proved] ", Tag::Proved),
    ("[Verified] ", Tag::Verified),
    ("[Conjecture] ", Tag::Conjecture),
    ("[Refuted] ", Tag::Refuted),
];

const SEED: usize = 6;

struct Draft {
    section: String,
    words: Vec<String>,
    tag: Tag,
    witnesses: Vec<String>,
    line: usize,
    text: String,
}

// PARSE

pub fn parse(text: &str) -> Result<Book, String> {
    let mut parts = 0;
    let mut untagged = 0;
    let mut section = String::new();
    let mut drafts = Vec::new();
    for (index, line) in text.lines().enumerate() {
        if let Some(heading) = line.strip_prefix("# ") {
            parts += 1;
            section = slug(heading);
        } else if let Some(body) = line.strip_prefix("- ") {
            if section.is_empty() {
                continue;
            }
            match tagged(body) {
                Some((tag, rest)) => drafts.push(draft(&section, tag, rest, index + 1, line)),
                None => untagged += 1,
            }
        }
    }
    Ok(Book {
        rows: key(drafts)?,
        parts,
        untagged,
    })
}

fn tagged(body: &str) -> Option<(Tag, &str)> {
    let body = dated(body)?;
    TAGS.iter()
        .find_map(|(mark, tag)| body.strip_prefix(mark).map(|rest| (*tag, rest)))
}

fn dated(body: &str) -> Option<&str> {
    let (date, rest) = body.split_once(' ')?;
    let digits = date.bytes().filter(u8::is_ascii_digit).count();
    (date.len() == 10 && digits == 8 && date.matches('-').count() == 2).then_some(rest)
}

fn draft(section: &str, tag: Tag, rest: &str, line: usize, text: &str) -> Draft {
    let (claim, witnesses) = match rest.split_once(" Witness: ") {
        Some((claim, field)) => (claim, witnesses(field)),
        None => (rest, Vec::new()),
    };
    Draft {
        section: section.to_string(),
        words: claim.split_whitespace().map(str::to_string).collect(),
        tag,
        witnesses,
        line,
        text: text.to_string(),
    }
}

fn witnesses(field: &str) -> Vec<String> {
    field
        .trim_end_matches('.')
        .split([',', ';'])
        .map(|token| token.trim().trim_matches('`').to_string())
        .filter(|token| !token.is_empty())
        .collect()
}

// KEY

fn key(drafts: Vec<Draft>) -> Result<Vec<Row>, String> {
    let mut widths: Vec<usize> = drafts
        .iter()
        .map(|draft| SEED.min(draft.words.len()))
        .collect();
    loop {
        let mut seen: HashMap<String, Vec<usize>> = HashMap::new();
        for (index, draft) in drafts.iter().enumerate() {
            seen.entry(name(draft, widths[index]))
                .or_default()
                .push(index);
        }
        let mut clash = false;
        for (name, group) in seen {
            if group.len() < 2 {
                continue;
            }
            clash = true;
            let mut room = false;
            for &index in &group {
                if widths[index] < drafts[index].words.len() {
                    widths[index] += 1;
                    room = true;
                }
            }
            if !room {
                let lines: Vec<String> = group
                    .iter()
                    .map(|&index| format!("{}: {}", drafts[index].line, drafts[index].text))
                    .collect();
                return Err(format!("duplicate key {name}\n{}", lines.join("\n")));
            }
        }
        if !clash {
            break;
        }
    }
    Ok(drafts
        .into_iter()
        .zip(widths)
        .map(|(draft, width)| Row {
            key: name(&draft, width),
            tag: draft.tag,
            witnesses: draft.witnesses,
        })
        .collect())
}

fn name(draft: &Draft, width: usize) -> String {
    format!(
        "{}/{}",
        draft.section,
        slug(&draft.words[..width].join(" "))
    )
}

fn slug(text: &str) -> String {
    let mut out = String::new();
    for letter in text.chars() {
        if letter.is_ascii_alphanumeric() {
            out.push(letter.to_ascii_lowercase());
        } else if !out.ends_with('-') {
            out.push('-');
        }
    }
    out.trim_matches('-').to_string()
}
