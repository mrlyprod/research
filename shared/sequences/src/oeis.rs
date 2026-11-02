use crate::merge::Class;
use crate::row::{transform, Row, TRANSFORMS};
use mrly_core::error::{value_error, Result};
use std::collections::{HashMap, HashSet};

pub const WINDOW: usize = 5;
pub const OFFSETS: usize = 4;
pub const KEEP: usize = 8;

#[derive(Clone, Debug)]
pub struct Hit {
    pub id: String,
    pub transform: &'static str,
    pub offset: usize,
    pub window: String,
    pub overlap: usize,
}

#[derive(Clone, Debug, Default)]
pub struct Search {
    pub hits: Vec<Vec<Hit>>,
    pub totals: Vec<usize>,
    pub searched: Vec<bool>,
    pub failed: Vec<(usize, String)>,
    pub names: HashMap<String, String>,
}

impl Search {
    pub fn absent(&self, class: usize) -> bool {
        self.searched[class] && self.totals[class] == 0
    }
}

struct Request {
    class: usize,
    transform: usize,
    offset: usize,
}

fn bounds(terms: &str) -> Vec<usize> {
    terms
        .bytes()
        .enumerate()
        .filter(|&(_, b)| b == b',')
        .map(|(i, _)| i)
        .collect()
}

fn run(image: &[String], line: &str, commas: &[usize], start: usize) -> usize {
    let mut overlap = 0;
    for (j, term) in image.iter().enumerate() {
        let i = start + j;
        if i + 1 >= commas.len() {
            break;
        }
        if &line[commas[i] + 1..commas[i + 1]] != term.as_str() {
            return 0;
        }
        overlap += 1;
    }
    overlap
}

fn contains(terms: &[i128], line: &str) -> bool {
    let commas = bounds(line);
    let image: Vec<String> = terms.iter().map(i128::to_string).collect();
    let need = image.len().min(WINDOW);
    (0..3.min(image.len())).any(|s| {
        (0..commas.len().saturating_sub(1))
            .any(|p| run(&image[s..], line, &commas, p) >= need.min(image.len() - s))
    })
}

pub fn search(rows: &[Row], classes: &[Class], stripped: &str, names: &str) -> Result<Search> {
    let mut images: Vec<Vec<Option<Vec<String>>>> = Vec::with_capacity(classes.len());
    let mut requests: HashMap<String, Vec<Request>> = HashMap::new();
    let mut searched = vec![false; classes.len()];
    for (c, class) in classes.iter().enumerate() {
        let terms = &rows[class.rep].terms;
        let mut per = Vec::new();
        for (t, name) in TRANSFORMS.iter().enumerate() {
            let image = transform(name, terms).filter(|i| i.len() >= WINDOW);
            if let Some(image) = &image {
                searched[c] = true;
                for offset in 0..OFFSETS {
                    if image.len() < offset + WINDOW {
                        break;
                    }
                    let window: String = image[offset..offset + WINDOW]
                        .iter()
                        .map(|v| format!(",{v}"))
                        .collect::<String>()
                        + ",";
                    requests.entry(window).or_default().push(Request {
                        class: c,
                        transform: t,
                        offset,
                    });
                }
            }
            per.push(image.map(|i| i.iter().map(i128::to_string).collect()));
        }
        images.push(per);
    }
    let stated: HashSet<&str> = rows.iter().filter_map(|r| r.stated.as_deref()).collect();
    let mut lines: HashMap<String, String> = HashMap::new();
    let mut hits: Vec<Vec<Hit>> = vec![Vec::new(); classes.len()];
    let mut totals = vec![0usize; classes.len()];
    let mut kept = vec![[0usize; TRANSFORMS.len()]; classes.len()];
    let mut scanned = 0usize;
    for line in stripped.lines() {
        if line.starts_with('#') || line.len() < 9 {
            continue;
        }
        scanned += 1;
        let (id, terms) = line.split_at(7);
        let terms = terms.trim_start();
        if stated.contains(id) {
            lines.insert(id.to_string(), terms.to_string());
        }
        let commas = bounds(terms);
        if commas.len() <= WINDOW {
            continue;
        }
        let mut seen: HashSet<(usize, usize)> = HashSet::new();
        for i in 0..commas.len() - WINDOW {
            let Some(found) = requests.get(&terms[commas[i]..=commas[i + WINDOW]]) else {
                continue;
            };
            for request in found {
                let image = images[request.class][request.transform].as_ref().unwrap();
                let overlap = run(&image[request.offset..], terms, &commas, i);
                if overlap < WINDOW || !seen.insert((request.class, request.transform)) {
                    continue;
                }
                totals[request.class] += 1;
                let count = &mut kept[request.class][request.transform];
                if *count < KEEP {
                    *count += 1;
                    hits[request.class].push(Hit {
                        id: id.to_string(),
                        transform: TRANSFORMS[request.transform],
                        offset: request.offset,
                        window: terms[commas[i]..=commas[i + WINDOW]].to_string(),
                        overlap,
                    });
                }
            }
        }
    }
    if scanned == 0 {
        return value_error("the stripped dump holds no records.");
    }
    for list in &mut hits {
        list.sort_by_key(|h| {
            (
                TRANSFORMS.iter().position(|t| *t == h.transform),
                std::cmp::Reverse(h.overlap),
                h.id.clone(),
            )
        });
        list.truncate(KEEP);
    }
    let failed: Vec<(usize, String)> = rows
        .iter()
        .enumerate()
        .filter_map(|(i, r)| {
            let id = r.stated.as_ref()?;
            let ok =
                r.terms.len() >= 3 && lines.get(id).is_some_and(|line| contains(&r.terms, line));
            (!ok).then(|| (i, id.clone()))
        })
        .collect();
    let mut wanted: HashSet<&str> = stated;
    for list in &hits {
        wanted.extend(list.iter().map(|h| h.id.as_str()));
    }
    let names = names
        .lines()
        .filter(|l| l.len() > 8 && wanted.contains(&l[..7]))
        .map(|l| (l[..7].to_string(), l[8..].to_string()))
        .collect();
    Ok(Search {
        hits,
        totals,
        searched,
        failed,
        names,
    })
}
