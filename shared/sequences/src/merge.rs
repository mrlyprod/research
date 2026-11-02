use crate::row::{constant, transform, Row, RELATIONS};
use std::collections::{BTreeSet, HashMap};

pub const LONG: usize = 8;
pub const SHORT: usize = 4;
pub const OVERLAP: usize = 7;

#[derive(Clone, Debug)]
pub struct Class {
    pub rep: usize,
    pub members: Vec<(usize, i32)>,
    pub labels: BTreeSet<String>,
    pub explanation: Option<String>,
    pub kind: &'static str,
}

impl Class {
    pub fn cross(&self) -> bool {
        self.labels.len() > 1
    }
}

pub fn simple(terms: &[i128]) -> bool {
    if terms.len() < 3 {
        return true;
    }
    let linear = terms
        .windows(3)
        .all(|w| w[2].checked_sub(w[1]) == w[1].checked_sub(w[0]));
    let geometric = terms[0] != 0
        && terms[1] % terms[0] == 0
        && terms
            .windows(2)
            .all(|w| w[0].checked_mul(terms[1] / terms[0]) == Some(w[1]));
    constant(terms) || linear || geometric
}

fn design(name: &str) -> Option<&str> {
    name.strip_prefix("sequence_")?.split("_measure=").next()
}

fn kind(
    rows: &[Row],
    members: &[usize],
    labels: &BTreeSet<String>,
    rep: usize,
    explained: bool,
) -> &'static str {
    if members.len() == 1 {
        return "single";
    }
    if labels.len() == 1 || explained {
        return "explained";
    }
    if rows[rep].terms.len() < 3 || simple(&rows[rep].terms) {
        return "trivial";
    }
    let mut seen: HashMap<&str, &str> = HashMap::new();
    for &i in members {
        if let Some(d) = design(&rows[i].name) {
            match seen.get(d) {
                Some(&m) if m != rows[i].measure => return "design",
                _ => {
                    seen.insert(d, &rows[i].measure);
                }
            }
        }
    }
    if labels.iter().any(|l| !l.contains('/') && l != "ledger") {
        return "topic";
    }
    "codes"
}

#[derive(Clone, Debug)]
pub struct Relation {
    pub near: bool,
    pub from: usize,
    pub transform: &'static str,
    pub offset: usize,
    pub to: usize,
    pub overlap: usize,
}

fn find(parent: &mut [usize], mut x: usize) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

fn union(parent: &mut [usize], a: usize, b: usize) {
    let (a, b) = (find(parent, a), find(parent, b));
    if a != b {
        parent[a.max(b)] = a.min(b);
    }
}

fn agree(a: &[i128], b: &[i128]) -> usize {
    let overlap = a.len().min(b.len());
    if a[..overlap] == b[..overlap] {
        overlap
    } else {
        0
    }
}

fn shift(rep: &[i128], member: &[i128]) -> i32 {
    for s in 0..=6usize {
        if s < rep.len() && agree(&rep[s..], member) >= SHORT.min(member.len()) {
            return s as i32;
        }
        if s < member.len() && agree(rep, &member[s..]) >= SHORT.min(member.len() - s) {
            return -(s as i32);
        }
    }
    0
}

pub fn explain(labels: &BTreeSet<String>) -> Option<String> {
    let has = |l: &str| labels.contains(l);
    let only = |set: &[&str]| labels.iter().all(|l| set.contains(&l.as_str()));
    let text = if labels.len() == 1 {
        "one measure: codes or orbit mates share the count"
    } else if labels.iter().filter(|l| *l != "ledger").count() <= 1 {
        "one measure, and the ledger page cites its record"
    } else if only(&["source", "source-counts"]) {
        "counts drops the zeros and ones of the same source"
    } else if only(&["arcs-law", "arcs-loops"]) {
        "the proved law gives the drawn count"
    } else if only(&["surface/level", "perimeter/level"])
        || only(&["surface/side", "perimeter/side"])
    {
        "in the plane the surface is the perimeter"
    } else if has("source") || has("source-counts") {
        if labels
            .iter()
            .all(|l| l.starts_with("source") || l.ends_with("/side") || l == "ledger")
        {
            "a design Source is the side axis read from side 1"
        } else {
            return None;
        }
    } else {
        return None;
    };
    Some(text.to_string())
}

pub fn classes(rows: &[Row]) -> Vec<Class> {
    let n = rows.len();
    let mut parent: Vec<usize> = (0..n).collect();
    let mut starts: HashMap<&[i128], Vec<usize>> = HashMap::new();
    let mut windows: HashMap<&[i128], Vec<(usize, usize)>> = HashMap::new();
    for (i, row) in rows.iter().enumerate() {
        if row.terms.len() >= LONG {
            starts.entry(&row.terms[..SHORT]).or_default().push(i);
            for t in 0..4 {
                windows
                    .entry(&row.terms[t..t + SHORT])
                    .or_default()
                    .push((i, t));
            }
        }
    }
    for (a, row) in rows.iter().enumerate() {
        if row.terms.len() < LONG {
            continue;
        }
        for s in 0..4 {
            for &b in starts.get(&row.terms[s..s + SHORT]).into_iter().flatten() {
                if b != a && agree(&row.terms[s..], &rows[b].terms) >= LONG - 3 {
                    union(&mut parent, a, b);
                }
            }
        }
    }
    let mut exact: HashMap<&[i128], usize> = HashMap::new();
    for (a, row) in rows.iter().enumerate() {
        let len = row.terms.len();
        if !(SHORT..LONG).contains(&len) {
            continue;
        }
        let mut roots = BTreeSet::new();
        for &(b, t) in windows.get(&row.terms[..SHORT]).into_iter().flatten() {
            let long = &rows[b].terms;
            if long.len() >= t + len && long[t..t + len] == row.terms[..] {
                roots.insert(find(&mut parent, b));
            }
        }
        if roots.len() == 1 {
            union(&mut parent, a, *roots.iter().next().unwrap());
        } else if roots.is_empty() {
            match exact.get(&row.terms[..]) {
                Some(&b) => union(&mut parent, a, b),
                None => {
                    exact.insert(&row.terms, a);
                }
            }
        }
    }
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..n {
        groups.entry(find(&mut parent, i)).or_default().push(i);
    }
    let mut out: Vec<Class> = groups
        .into_values()
        .map(|members| {
            let rep = *members
                .iter()
                .max_by_key(|&&i| (rows[i].terms.len(), std::cmp::Reverse(i)))
                .unwrap();
            let designed: BTreeSet<&str> = members
                .iter()
                .filter(|&&i| rows[i].measure != "ledger")
                .map(|&i| rows[i].name.as_str())
                .collect();
            let labels: BTreeSet<String> = members
                .iter()
                .filter(|&&i| {
                    rows[i].measure != "ledger" || !designed.contains(rows[i].name.as_str())
                })
                .map(|&i| rows[i].measure.clone())
                .collect();
            let explanation = explain(&labels);
            let kind = kind(rows, &members, &labels, rep, explanation.is_some());
            let members = members
                .iter()
                .map(|&i| (i, shift(&rows[rep].terms, &rows[i].terms)))
                .collect();
            Class {
                rep,
                members,
                labels,
                explanation,
                kind,
            }
        })
        .collect();
    out.sort_by_key(|c| c.members.iter().map(|m| m.0).min().unwrap());
    out
}

pub fn relations(rows: &[Row], classes: &[Class]) -> Vec<Relation> {
    let mut index: HashMap<&[i128], Vec<usize>> = HashMap::new();
    for (c, class) in classes.iter().enumerate() {
        let terms = &rows[class.rep].terms;
        if terms.len() >= LONG && !simple(terms) {
            index.entry(&terms[..5]).or_default().push(c);
        }
    }
    let mut out = Vec::new();
    for (c, class) in classes.iter().enumerate() {
        let terms = &rows[class.rep].terms;
        if terms.len() < LONG || simple(terms) {
            continue;
        }
        for name in &RELATIONS {
            let Some(image) = transform(name, terms) else {
                continue;
            };
            for offset in 0..4 {
                if image.len() < offset + OVERLAP {
                    break;
                }
                for &d in index.get(&image[offset..offset + 5]).into_iter().flatten() {
                    if d == c || classes[d].labels == class.labels {
                        continue;
                    }
                    let overlap = agree(&image[offset..], &rows[classes[d].rep].terms);
                    if overlap >= OVERLAP {
                        let near = [c, d].iter().any(|&k| {
                            classes[k]
                                .labels
                                .iter()
                                .any(|l| !l.contains('/') && l != "ledger")
                        }) || classes[c].members.iter().any(|&(i, _)| {
                            let a = design(&rows[i].name);
                            a.is_some()
                                && classes[d]
                                    .members
                                    .iter()
                                    .any(|&(j, _)| design(&rows[j].name) == a)
                        });
                        out.push(Relation {
                            near,
                            from: c,
                            transform: name,
                            offset,
                            to: d,
                            overlap,
                        });
                    }
                }
            }
        }
    }
    out
}
