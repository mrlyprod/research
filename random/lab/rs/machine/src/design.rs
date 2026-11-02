use std::fs;
use std::path::PathBuf;

// CATALOG

pub const HOODS: [(&str, u32); 3] = [("carpet", 495), ("net", 186), ("void", 341)];
pub const SPREADS: [&str; 6] = ["id", "shift", "drop", "slide", "turn", "flip"];
pub const LETTERS: std::ops::RangeInclusive<u32> = 1..=14;

pub fn offsets(code: u32) -> Vec<(i64, i64)> {
    (0..9)
        .filter(|&b| code >> b & 1 == 1 && b != 4)
        .map(|b| (b as i64 / 3 - 1, b as i64 % 3 - 1))
        .collect()
}

pub fn spread_at(name: &str, n: usize, i: usize, j: usize) -> (usize, usize) {
    match name {
        "id" => (i, j),
        "shift" => (i, (j + 1) % n),
        "drop" => ((i + 1) % n, j),
        "slide" => ((i + 1) % n, (j + 1) % n),
        "turn" => (j, n - 1 - i),
        "flip" => (j, i),
        _ => unreachable!(),
    }
}

#[derive(Clone, PartialEq)]
pub enum Spread {
    Map(&'static str),
    Sum([(i64, i64); 3]),
}

impl Spread {
    pub fn name(&self) -> String {
        match self {
            Spread::Map(s) => s.to_string(),
            Spread::Sum(v) => format!(
                "sum{}",
                v.iter()
                    .map(|(a, b)| format!("({a},{b})"))
                    .collect::<String>()
            ),
        }
    }
    pub fn sources(&self, n: usize) -> (Vec<[u32; 3]>, usize) {
        let m = n as i64;
        match self {
            Spread::Map(name) => {
                let mut src = vec![[0u32; 3]; n * n];
                for p in 0..n * n {
                    let (i, j) = spread_at(name, n, p / n, p % n);
                    src[i * n + j][0] = p as u32;
                }
                (src, 1)
            }
            Spread::Sum(v) => {
                let src = (0..n * n)
                    .map(|q| {
                        let (i, j) = ((q / n) as i64, (q % n) as i64);
                        v.map(|(a, b)| ((i - a).rem_euclid(m) * m + (j - b).rem_euclid(m)) as u32)
                    })
                    .collect();
                (src, 3)
            }
        }
    }
}

pub fn triples(radius: i64) -> Vec<[(i64, i64); 3]> {
    let box_: Vec<(i64, i64)> = (-radius..=radius)
        .flat_map(|a| (-radius..=radius).map(move |b| (a, b)))
        .collect();
    let mut out = Vec::new();
    for x in 0..box_.len() {
        for y in x + 1..box_.len() {
            for z in y + 1..box_.len() {
                out.push([box_[x], box_[y], box_[z]]);
            }
        }
    }
    out
}

pub fn power(letter: u32, level: u32) -> Vec<bool> {
    let n = 1usize << level;
    (0..n * n)
        .map(|p| {
            let (i, j) = (p / n, p % n);
            (0..level).all(|l| letter >> (2 * (i >> l & 1) + (j >> l & 1)) & 1 == 1)
        })
        .collect()
}

#[derive(Clone)]
pub struct Split {
    pub name: String,
    pub letter: u32,
    pub anti: bool,
}

impl Split {
    pub fn cells(&self, level: u32) -> Vec<bool> {
        power(self.letter, level)
            .into_iter()
            .map(|a| a != self.anti)
            .collect()
    }
}

pub fn splits() -> Vec<Split> {
    [false, true]
        .into_iter()
        .flat_map(|anti| {
            LETTERS.map(move |letter| Split {
                name: format!("{}{letter}^L", if anti { "anti " } else { "" }),
                letter,
                anti,
            })
        })
        .collect()
}

// SOURCES

const WEIGHTS: [(&str, &[u32]); 8] = [
    ("carpet", &[0, 1, 1]),
    ("net", &[1, 1, 2]),
    ("tree", &[0, 1]),
    ("void", &[0, 2]),
    ("point", &[2]),
    ("dust", &[0]),
    ("line", &[1, 2]),
    ("star", &[1, 1]),
];

fn fill(corners: &[u32], side: u64) -> u64 {
    let (even, odd) = ((side + 1) / 2, side / 2);
    corners.iter().map(|&w| even.pow(2 - w) * odd.pow(w)).sum()
}

fn design_set(corners: &[u32], void: bool, limit: u64) -> Vec<u64> {
    let mut out = Vec::new();
    let mut side = 1;
    while side <= limit + 3 {
        let f = fill(corners, side);
        let value = if void { side * side - f } else { f };
        if value > limit {
            break;
        }
        out.push(value);
        side += 2;
    }
    out
}

fn is_prime(v: u64) -> bool {
    v >= 2 && (2..v).take_while(|d| d * d <= v).all(|d| v % d != 0)
}

pub fn named_sources(limit: u64) -> Vec<(String, Vec<u64>)> {
    let mut fib = vec![0u64, 1];
    while fib[fib.len() - 1] + fib[fib.len() - 2] <= limit {
        fib.push(fib[fib.len() - 1] + fib[fib.len() - 2]);
    }
    let mut out: Vec<(String, Vec<u64>)> = vec![
        ("evens".into(), (0..=limit).step_by(2).collect()),
        ("odds".into(), (1..=limit).step_by(2).collect()),
        (
            "primes".into(),
            (0..=limit).filter(|&v| is_prime(v)).collect(),
        ),
        (
            "binary".into(),
            (0..64)
                .map(|i| 1u64 << i)
                .take_while(|&v| v <= limit)
                .collect(),
        ),
        ("fibonacci".into(), fib),
        (
            "grid_squares".into(),
            (1..=limit + 1)
                .step_by(2)
                .map(|s| s * s)
                .take_while(|&v| v <= limit)
                .collect(),
        ),
    ];
    for (name, corners) in WEIGHTS {
        out.push((format!("{name}_fills"), design_set(corners, false, limit)));
        out.push((format!("{name}_voids"), design_set(corners, true, limit)));
    }
    out
}

pub fn indicator(terms: &[u64], k: usize) -> Vec<u8> {
    let mut g = vec![0u8; k + 1];
    for &t in terms {
        if (t as usize) <= k {
            g[t as usize] = 1;
        }
    }
    g
}

pub fn census_check(k: usize) -> Option<(usize, usize)> {
    let here = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let repo = here.ancestors().nth(5)?;
    let text =
        fs::read_to_string(repo.join("data/research/random/lab/rs/mess/census.json")).ok()?;
    let mine = named_sources(64);
    let (mut read, mut agree) = (0, 0);
    for chunk in text.split("{\"name\":\"").skip(1) {
        let name = &chunk[..chunk.find('"')?];
        let at = chunk.find("\"g\":[")? + 5;
        let g: Vec<u8> = chunk[at..at + chunk[at..].find(']')?]
            .split(',')
            .map(|v| v.trim().parse().unwrap_or(9))
            .take(k + 1)
            .collect();
        read += 1;
        if mine
            .iter()
            .any(|(n, terms)| n == name && indicator(terms, k) == g)
        {
            agree += 1;
        }
    }
    Some((read, agree))
}

// RULES

#[derive(Clone)]
pub enum Mess {
    Feistel(Vec<u8>),
    Life { birth: Vec<u8>, survive: Vec<u8> },
}

#[derive(Clone)]
pub struct Rule {
    pub name: String,
    pub mess: Mess,
}

pub fn bits(g: &[u8]) -> String {
    g.iter().map(|b| char::from(b'0' + b)).collect()
}

pub fn parity_distance(g: &[u8]) -> usize {
    g.iter()
        .enumerate()
        .filter(|&(j, &v)| v as usize != j % 2)
        .count()
}

pub fn sources_at(k: usize) -> Vec<(String, Vec<u8>)> {
    let mut out: Vec<(String, Vec<u8>)> = Vec::new();
    for (name, terms) in named_sources(64) {
        let g = indicator(&terms, k);
        match out.iter_mut().find(|(_, h)| *h == g) {
            Some((names, _)) => {
                names.push('/');
                names.push_str(&name);
            }
            None => out.push((name, g)),
        }
    }
    out
}

pub fn dense_at(k: usize) -> Vec<Vec<u8>> {
    let (lo, hi) = ((k + 1) / 2, (k + 2) / 2);
    (0u32..1 << (k + 1))
        .map(|m| (0..=k).map(|j| (m >> j & 1) as u8).collect::<Vec<u8>>())
        .filter(|g| (lo..=hi).contains(&parity_distance(g)))
        .collect()
}

pub fn rules_at(k: usize) -> (Vec<Rule>, usize, usize) {
    let named = sources_at(k);
    let dense = dense_at(k);
    let mut out: Vec<Rule> = named
        .iter()
        .map(|(name, g)| Rule {
            name: name.clone(),
            mess: Mess::Feistel(g.clone()),
        })
        .collect();
    for g in &dense {
        if !named.iter().any(|(_, h)| h == g) {
            out.push(Rule {
                name: format!("g={}", bits(g)),
                mess: Mess::Feistel(g.clone()),
            });
        }
    }
    (out, named.len(), dense.len())
}

pub fn obstruction(k: usize) -> Rule {
    let named = named_sources(64);
    let get = |want: &str| indicator(&named.iter().find(|(n, _)| n == want).unwrap().1, k);
    Rule {
        name: "life B=primes S=odds".into(),
        mess: Mess::Life {
            birth: get("primes"),
            survive: get("odds"),
        },
    }
}

pub fn parity_agreement(birth: &[u8], survive: &[u8]) -> f64 {
    let k = birth.len() - 1;
    let total = (1u64 << k) as f64;
    let mut agree = 0.0;
    let mut binom = 1u64;
    for j in 0..=k {
        let par = (j % 2) as u8;
        let hits = (birth[j] == par) as u64 + (survive[j] == par) as u64;
        agree += hits as f64 * binom as f64 / (2.0 * total);
        binom = binom * (k - j) as u64 / (j + 1) as u64;
    }
    agree
}

pub fn anf(g: &[u8]) -> [u64; 9] {
    let mut sel = [0u64; 9];
    for (w, s) in sel.iter_mut().enumerate() {
        if w < g.len() {
            let c = (0..=w).filter(|&j| j & w == j).fold(0, |c, j| c ^ g[j]);
            *s = if c == 1 { !0 } else { 0 };
        }
    }
    sel
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primes_birth_odds_survive_agrees_with_parity_on_476_of_512() {
        if let Mess::Life { birth, survive } = obstruction(8).mess {
            assert!((parity_agreement(&birth, &survive) - 476.0 / 512.0).abs() < 1e-12);
        }
    }
}
