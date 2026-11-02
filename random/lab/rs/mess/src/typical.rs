use crate::mess::{balanced, ceiling, profile, Tables};
use mrly_life::source::{sequence, Source};

pub const SEED: u64 = 0x6D65_7373;
pub const DRAWS: usize = 2000;
pub const KS: [usize; 6] = [8, 16, 24, 32, 48, 64];

struct Split(u64);

impl Split {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

pub struct Typical {
    pub k: usize,
    pub median_ratio: f64,
    pub min_ratio: f64,
    pub max_ratio: f64,
    pub full_degree: f64,
    pub balanced: f64,
    pub median_degree: i32,
}

pub fn run(k: usize, tree: bool) -> Typical {
    let tables = Tables::new(k);
    let mut rng = Split(SEED ^ k as u64);
    let top = ceiling(k) as f64;
    let mut ratios = Vec::new();
    let mut degrees = Vec::new();
    let mut full = 0usize;
    let mut even = 0usize;
    for draw in 0..DRAWS {
        let g: Vec<u8> = if tree {
            let mut g = vec![0u8; k + 1];
            for term in sequence(Source::Random(draw as u64), k).expect("a random source builds") {
                g[term] = 1;
            }
            g
        } else {
            (0..=k).map(|_| (rng.next() >> 63) as u8).collect()
        };
        let p = profile(&tables, &g);
        ratios.push(p.nonlinearity as f64 / top);
        degrees.push(p.degree);
        full += (p.degree == k as i32) as usize;
        even += balanced(&p) as usize;
    }
    ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
    degrees.sort_unstable();
    Typical {
        k,
        median_ratio: ratios[DRAWS / 2],
        min_ratio: ratios[0],
        max_ratio: ratios[DRAWS - 1],
        full_degree: full as f64 / DRAWS as f64,
        balanced: even as f64 / DRAWS as f64,
        median_degree: degrees[DRAWS / 2],
    }
}
