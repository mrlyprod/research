use std::time::{Duration, Instant};

pub const TERMS: usize = 16;
pub const SECONDS: Duration = Duration::from_secs(2);

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Row {
    pub name: String,
    pub generator: String,
    pub measure: String,
    pub terms: Vec<i128>,
    pub closed: Option<String>,
    pub stated: Option<String>,
}

impl Row {
    pub fn new(name: impl Into<String>, generator: &str, measure: &str, terms: Vec<i128>) -> Row {
        Row {
            name: name.into(),
            generator: generator.to_string(),
            measure: measure.to_string(),
            terms,
            closed: None,
            stated: None,
        }
    }

    pub fn closed(mut self, closed: Option<String>) -> Row {
        self.closed = closed;
        self
    }

    pub fn stated(mut self, stated: Option<String>) -> Row {
        self.stated = stated;
        self
    }
}

pub fn timed<F: FnMut(usize) -> Option<i128>>(count: usize, mut term: F) -> Vec<i128> {
    let start = Instant::now();
    let mut out = Vec::new();
    for index in 0..count {
        match term(index) {
            Some(value) => out.push(value),
            None => break,
        }
        if start.elapsed() > SECONDS {
            break;
        }
    }
    out
}

pub const RELATIONS: [&str; 8] = ["-a", "diff", "psum", "a-1", "a+1", "a/2", "even", "odd"];
pub const TRANSFORMS: [&str; 8] = ["a", "-a", "diff", "psum", "a-1", "a+1", "a/2", "a/gcd"];

pub fn transform(name: &str, terms: &[i128]) -> Option<Vec<i128>> {
    match name {
        "a" => Some(terms.to_vec()),
        "-a" => terms.iter().map(|t| t.checked_neg()).collect(),
        "diff" => terms
            .windows(2)
            .map(|w| w[1].checked_sub(w[0]))
            .collect::<Option<Vec<i128>>>(),
        "psum" => {
            let mut total = 0i128;
            terms
                .iter()
                .map(|&t| {
                    total = total.checked_add(t)?;
                    Some(total)
                })
                .collect()
        }
        "a-1" => terms.iter().map(|t| t.checked_sub(1)).collect(),
        "a+1" => terms.iter().map(|t| t.checked_add(1)).collect(),
        "even" => Some(terms.iter().step_by(2).copied().collect()),
        "odd" => Some(terms.iter().skip(1).step_by(2).copied().collect()),
        "a/2" => terms
            .iter()
            .map(|t| (t % 2 == 0).then_some(t / 2))
            .collect(),
        "a/gcd" => {
            let g = terms.iter().fold(0u128, |g, t| gcd(g, t.unsigned_abs()));
            (g > 1).then(|| terms.iter().map(|t| t / g as i128).collect())
        }
        _ => None,
    }
}

fn gcd(a: u128, b: u128) -> u128 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

pub fn constant(terms: &[i128]) -> bool {
    terms.windows(2).all(|w| w[0] == w[1])
}
