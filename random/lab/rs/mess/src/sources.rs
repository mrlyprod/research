use mrly_life::source::{counts, Source};

pub const KMAX: usize = 64;

pub struct Named {
    pub name: String,
    pub oeis: Option<String>,
    pub source: Source,
    pub zeros: bool,
    pub ones: bool,
    pub g: Vec<u8>,
}

fn indicator(source: Source, limit: usize, zeros: bool, ones: bool) -> Vec<u8> {
    let terms = counts(source, limit, zeros, ones).expect("a fixed source builds");
    let mut g = vec![0u8; limit + 1];
    for term in &terms {
        g[*term] = 1;
    }
    g
}

pub fn fixed() -> Vec<Named> {
    variant(true, true)
}

pub fn variant(zeros: bool, ones: bool) -> Vec<Named> {
    Source::all()
        .into_iter()
        .filter(|source| !source.is_random())
        .map(|source| Named {
            name: source.name(),
            oeis: source.oeis(),
            source,
            zeros,
            ones,
            g: indicator(source, KMAX, zeros, ones),
        })
        .collect()
}

pub fn prefix(named: &Named, k: usize) -> Vec<u8> {
    let g = indicator(named.source, k, named.zeros, named.ones);
    assert_eq!(
        g[..],
        named.g[..=k],
        "{} is not prefix stable at {k}",
        named.name
    );
    g
}

pub fn designs() -> Vec<usize> {
    let mut fills: Vec<usize> = (1..=9).collect();
    fills.extend((1..=9).map(|f| f * f));
    fills.sort_unstable();
    fills.dedup();
    fills.retain(|&k| k <= KMAX);
    fills
}
