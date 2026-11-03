pub mod markdown;
pub mod records;
pub mod terms;
pub mod windows;

pub use markdown::markdown;
pub use mrly_math::counts::{fill_polynomial, Closed};
pub use mrly_sequences::{Axis, Cost, Key, Measure, Sequence, Tag};
pub use records::{identify, Record, RECORDS};
pub use terms::{closed, terms};

use mrly_core::error::{value_error, Result};
use mrly_math::bang::baseq;
use std::collections::BTreeMap;
use std::sync::{Mutex, OnceLock};

type Walked = BTreeMap<(usize, usize), &'static [u128]>;

pub const SPACES: [(usize, usize); 9] = [
    (1, 2),
    (2, 2),
    (3, 2),
    (4, 2),
    (1, 3),
    (2, 3),
    (1, 4),
    (2, 4),
    (1, 5),
];

pub const BUDGET: u128 = 500_000;

pub const TERMS: usize = 8;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Tier {
    Closed,
    Convolved,
    SideGrid,
    LevelGrid,
    WindowGrid,
}

impl Tier {
    pub const ALL: [Tier; 5] = [
        Tier::Closed,
        Tier::Convolved,
        Tier::SideGrid,
        Tier::LevelGrid,
        Tier::WindowGrid,
    ];
    pub fn slug(self) -> &'static str {
        match self {
            Tier::Closed => "closed",
            Tier::Convolved => "convolved",
            Tier::SideGrid => "side",
            Tier::LevelGrid => "level",
            Tier::WindowGrid => "window",
        }
    }
    pub fn parse(slug: &str) -> Result<Tier> {
        Tier::ALL
            .into_iter()
            .find(|tier| tier.slug() == slug)
            .map_or_else(|| value_error(format!("unknown tier {slug:?}.")), Ok)
    }
    fn cost(self) -> Cost {
        match self {
            Tier::Closed => Cost::Closed,
            Tier::Convolved => Cost::Convolved,
            _ => Cost::Grid,
        }
    }
    pub fn axes(self) -> &'static [Axis] {
        match self {
            Tier::SideGrid => &[Axis::Side],
            Tier::LevelGrid => &[Axis::Level],
            Tier::WindowGrid => &[Axis::Window],
            _ => &[Axis::Level, Axis::Side],
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Row {
    pub sequence: Sequence,
    pub record: Option<(&'static Record, i32)>,
}

impl Row {
    fn mentions(&self, needle: &str) -> bool {
        self.sequence.key.name().contains(needle)
            || self.record.is_some_and(|(record, _)| {
                record.id.to_lowercase().contains(needle)
                    || record.name.to_lowercase().contains(needle)
            })
    }
}

pub fn designs(dimension: usize, base: usize) -> Result<&'static [u128]> {
    static CACHE: OnceLock<Mutex<Walked>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(BTreeMap::new()));
    let mut guard = cache.lock().expect("the design cache is not poisoned");
    if let Some(codes) = guard.get(&(dimension, base)) {
        return Ok(codes);
    }
    let codes: Vec<u128> = baseq::representatives(base, dimension)?
        .into_iter()
        .map(|(code, _)| code.get())
        .collect();
    let leaked: &'static [u128] = Box::leak(codes.into_boxed_slice());
    guard.insert((dimension, base), leaked);
    Ok(leaked)
}

pub fn keys(tier: Tier) -> Vec<Key> {
    let mut out = Vec::new();
    for (dimension, base) in SPACES {
        let codes = designs(dimension, base).expect("the ledger spaces are walkable");
        for &code in codes {
            for measure in Measure::ALL {
                if !measure.applies(dimension, base) {
                    continue;
                }
                for &axis in tier.axes().iter().filter(|&&axis| {
                    measure.axes().contains(&axis) && measure.cost(axis) == tier.cost()
                }) {
                    out.push(Key::new(code, dimension, base, measure, axis));
                }
            }
        }
    }
    out
}

fn attach(key: &Key, terms: &[i128]) -> Option<(&'static Record, i32)> {
    if terms.len() < 4 || terms.iter().all(|&term| term == terms[0]) {
        return None;
    }
    let found = identify(terms);
    found
        .iter()
        .find(|(record, _)| record.key == Some(*key))
        .or(found.first())
        .map(|&(record, shift)| (record, shift - key.axis.start()))
}

pub fn sequence(key: &Key, count: usize, cells: u128) -> Result<Row> {
    let (terms, capped) = terms(key, count, cells)?;
    let closed = closed(key)?;
    let record = attach(key, &terms);
    let tag = record.map(|(record, _)| {
        if record.key == Some(*key) {
            record.status
        } else {
            Tag::Conjecture
        }
    });
    Ok(Row {
        sequence: Sequence {
            key: *key,
            terms,
            capped,
            closed,
            tag,
        },
        record,
    })
}

pub fn catalog(tier: Tier, count: usize) -> Vec<Row> {
    keys(tier)
        .iter()
        .filter_map(|key| sequence(key, count, BUDGET).ok())
        .collect()
}

pub fn numbers(text: &str) -> Option<Vec<i128>> {
    let tokens: Vec<&str> = text
        .split(|c: char| c == ',' || c.is_whitespace())
        .filter(|token| !token.is_empty())
        .collect();
    if tokens.is_empty() {
        return None;
    }
    tokens.iter().map(|token| token.parse().ok()).collect()
}

pub fn search(catalog: &[Row], query: &str) -> Vec<usize> {
    let query = query.trim();
    if query.is_empty() {
        return (0..catalog.len()).collect();
    }
    if let Some(window) = numbers(query) {
        return catalog
            .iter()
            .enumerate()
            .filter(|(_, row)| {
                row.sequence
                    .terms
                    .windows(window.len())
                    .any(|w| w == window.as_slice())
            })
            .map(|(index, _)| index)
            .collect();
    }
    let needle = query.to_lowercase();
    catalog
        .iter()
        .enumerate()
        .filter(|(_, row)| row.mentions(&needle))
        .map(|(index, _)| index)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_classics_read_their_records() {
        let carpet = sequence(&Key::new(7, 2, 2, Measure::Fills, Axis::Side), 4, BUDGET).unwrap();
        assert_eq!(carpet.sequence.terms, [8, 21, 40, 65]);
        assert_eq!(carpet.record.map(|(r, s)| (r.id, s)), Some(("A000567", 0)));
        assert_eq!(carpet.sequence.tag, Some(Tag::Proved));
        assert_eq!(carpet.sequence.closed.unwrap().text(), "3k^2 - 2k");
        let tree = sequence(&Key::new(3, 2, 2, Measure::Fills, Axis::Side), 4, BUDGET).unwrap();
        assert_eq!(tree.sequence.terms[..3], [6, 15, 28]);
        assert_eq!(tree.record.map(|(r, s)| (r.id, s)), Some(("A000384", 0)));
        let sponge = sequence(
            &Key::new(23, 3, 2, Measure::Surface, Axis::Level),
            3,
            BUDGET,
        )
        .unwrap();
        assert_eq!(sponge.sequence.terms, [72, 1056, 18048]);
        assert_eq!(
            sponge.sequence.closed.unwrap().text(),
            "a(level) = 28 a(level-1) - 160 a(level-2)"
        );
        let slice = sequence(
            &Key::new(23, 3, 2, Measure::Triangles, Axis::Level),
            8,
            BUDGET,
        )
        .unwrap();
        assert_eq!(slice.sequence.terms, [42, 306, 2250, 16578]);
        assert!(slice.sequence.capped);
        assert_eq!(slice.record.map(|(r, s)| (r.id, s)), Some(("A299916", 1)));
        let void = sequence(&Key::new(9, 2, 2, Measure::Voids, Axis::Side), 3, BUDGET).unwrap();
        assert_eq!(void.sequence.closed.unwrap().text(), "2k^2 - 2k");
        assert_eq!(void.sequence.terms, [4, 12, 24]);
    }

    #[test]
    fn the_keys_of_the_closed_tier_cover_every_space() {
        let closed = keys(Tier::Closed);
        assert_eq!(closed.len(), 1282 * 3 * 2 + (6 + 26 + 805) * 3 + 22 * 4);
        assert_eq!(keys(Tier::Convolved).len(), (1282 - 3 - 4 - 6 - 8) * 2 * 2);
        assert_eq!(keys(Tier::WindowGrid).len(), 6 + 26 + 805);
        assert!(keys(Tier::LevelGrid)
            .iter()
            .all(|key| key.measure != Measure::Windows));
        assert!(closed
            .iter()
            .all(|key| key.measure.cost(key.axis) == Cost::Closed));
        assert!(designs(3, 3).is_err());
        assert_eq!(Tier::parse("side").unwrap(), Tier::SideGrid);
    }

    #[test]
    fn the_search_finds_terms_and_names() {
        let rows: Vec<Row> = [
            Key::new(7, 2, 2, Measure::Fills, Axis::Level),
            Key::new(7, 2, 2, Measure::Surface, Axis::Level),
            Key::new(23, 3, 2, Measure::Fills, Axis::Level),
        ]
        .iter()
        .map(|key| sequence(key, 5, BUDGET).unwrap())
        .collect();
        assert_eq!(search(&rows, "64, 512"), [0]);
        assert_eq!(search(&rows, "80 496"), [1]);
        assert_eq!(search(&rows, "dim=3_code=23"), [2]);
        assert_eq!(search(&rows, "A381517"), [1]);
        assert_eq!(search(&rows, "surface"), [1]);
        assert_eq!(search(&rows, ""), [0, 1, 2]);
        assert!(search(&rows, "5, 6, 7").is_empty());
        assert_eq!(numbers("1, 2 3"), Some(vec![1, 2, 3]));
        assert_eq!(numbers("1, x"), None);
    }
}
