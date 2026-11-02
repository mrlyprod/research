use mrly_math::bang::factory::{code_to_corners, corners_to_code, residue_corners};
use mrly_math::bang::Code;
use mrly_num::factor::{divisors, factorize, gcd, mobius, radical};
use mrly_num::series::zeta;
use std::collections::BTreeSet;

const BUDGET: usize = 200_000;
const TOLERANCE: f64 = 0.06;
const EXACT_LEVEL: u32 = 4;
const ZETA_TERMS: usize = 10_000;
const HEAD: usize = 6;

pub struct Line {
    pub code: Code,
    pub k: usize,
    pub index: usize,
    pub spanning: bool,
    pub bracket: f64,
    pub predicted: f64,
    pub measured: f64,
    pub level: u32,
    pub exact: bool,
    pub within: bool,
    pub monotone: bool,
    pub head: Vec<u64>,
}

pub struct Family {
    pub label: String,
    pub base: usize,
    pub lines: Vec<Line>,
}

impl Family {
    fn spanning(&self) -> usize {
        self.lines.iter().filter(|line| line.spanning).count()
    }

    fn distinct(&self, spanning_only: bool) -> usize {
        let heads: BTreeSet<&Vec<u64>> = self
            .lines
            .iter()
            .filter(|line| line.spanning || !spanning_only)
            .map(|line| &line.head)
            .collect();
        heads.len()
    }

    fn flagged(&self) -> usize {
        self.lines
            .iter()
            .filter(|line| line.spanning && (!line.within || !line.exact))
            .count()
    }

    fn exact_failures(&self) -> usize {
        self.lines.iter().filter(|line| !line.exact).count()
    }
}

fn squarefree_divisors(base: usize) -> Vec<usize> {
    divisors(radical(base) as u64)
        .into_iter()
        .map(|d| d as usize)
        .collect()
}

fn hits(corners: &[Vec<u8>], divisor: usize) -> usize {
    corners
        .iter()
        .filter(|corner| corner.iter().all(|&r| (r as usize).is_multiple_of(divisor)))
        .count()
}

pub fn bracket(corners: &[Vec<u8>], base: usize) -> f64 {
    let sum: f64 = squarefree_divisors(base)
        .iter()
        .map(|&e| f64::from(mobius(e)) * hits(corners, e) as f64)
        .sum();
    sum / corners.len() as f64
}

pub fn foreign(base: usize, dimension: usize) -> f64 {
    let mut out = 1.0 / zeta(dimension as f64, ZETA_TERMS).expect("the dimension is above one");
    for (prime, _) in factorize(base) {
        out /= 1.0 - (prime as f64).powi(-(dimension as i32));
    }
    out
}

fn determinant(rows: &[Vec<i64>]) -> i64 {
    if rows.len() == 1 {
        return rows[0][0];
    }
    let mut total = 0;
    for (column, head) in rows[0].iter().enumerate() {
        let minor: Vec<Vec<i64>> = rows[1..]
            .iter()
            .map(|row| {
                row.iter()
                    .enumerate()
                    .filter(|(i, _)| *i != column)
                    .map(|(_, v)| *v)
                    .collect()
            })
            .collect();
        let sign = if column % 2 == 0 { 1 } else { -1 };
        total += sign * head * determinant(&minor);
    }
    total
}

fn choose(start: usize, count: usize, size: usize) -> Vec<Vec<usize>> {
    if size == 0 {
        return vec![Vec::new()];
    }
    let mut out = Vec::new();
    for first in start..count {
        for mut rest in choose(first + 1, count, size - 1) {
            rest.insert(0, first);
            out.push(rest);
        }
    }
    out
}

pub fn index(corners: &[Vec<u8>], dimension: usize) -> usize {
    let first = &corners[0];
    let diffs: Vec<Vec<i64>> = corners[1..]
        .iter()
        .map(|corner| {
            corner
                .iter()
                .zip(first)
                .map(|(a, b)| i64::from(*a) - i64::from(*b))
                .collect()
        })
        .collect();
    let mut out = 0usize;
    for rows in choose(0, diffs.len(), dimension) {
        let minor: Vec<Vec<i64>> = rows.iter().map(|&r| diffs[r].clone()).collect();
        out = gcd(out as u128, determinant(&minor).unsigned_abs() as u128) as usize;
        if out == 1 {
            return out;
        }
    }
    out
}

fn expand(columns: &mut [Vec<usize>], corners: &[Vec<u8>], base: usize) {
    for (axis, column) in columns.iter_mut().enumerate() {
        let mut next = Vec::with_capacity(column.len() * corners.len());
        for corner in corners {
            next.extend(column.iter().map(|&x| base * x + corner[axis] as usize));
        }
        *column = next;
    }
}

fn common(columns: &[Vec<usize>]) -> Vec<usize> {
    (0..columns[0].len())
        .map(|i| {
            columns
                .iter()
                .fold(0, |g, column| gcd(g as u128, column[i] as u128) as usize)
        })
        .collect()
}

fn measure(corners: &[Vec<u8>], base: usize, dimension: usize) -> (Vec<u64>, bool) {
    let k = corners.len();
    let mut columns = vec![vec![0usize]; dimension];
    let mut terms = Vec::new();
    let mut exact = None;
    let mut level = 0u32;
    while k.pow(level + 1) <= BUDGET {
        level += 1;
        expand(&mut columns, corners, base);
        let shared = common(&columns);
        terms.push(shared.iter().filter(|&&g| g == 1).count() as u64);
        if level == EXACT_LEVEL {
            exact = Some(
                squarefree_divisors(base)
                    .iter()
                    .filter(|&&e| e > 1)
                    .all(|&e| {
                        let seen = shared.iter().filter(|&&g| g % e == 0).count();
                        seen == hits(corners, e) * k.pow(level - 1)
                    }),
            );
        }
    }
    (
        terms,
        exact.expect("the point budget reaches the exact level"),
    )
}

pub fn survey(
    base: usize,
    dimension: usize,
    codes: impl Iterator<Item = Code>,
    label: &str,
) -> Family {
    let mut lines = Vec::new();
    for code in codes {
        let corners = code_to_corners(code, dimension, base).expect("the code fits its cells");
        let k = corners.len();
        if k < 2 {
            continue;
        }
        let (terms, exact) = measure(&corners, base, dimension);
        let level = terms.len() as u32;
        let measured = terms[terms.len() - 1] as f64 / k.pow(level) as f64;
        let bracket = bracket(&corners, base);
        let predicted = bracket * foreign(base, dimension);
        let gaps: Vec<f64> = terms
            .iter()
            .enumerate()
            .map(|(step, &count)| (count as f64 / k.pow(step as u32 + 1) as f64 - predicted).abs())
            .collect();
        let monotone = gaps.windows(2).all(|pair| pair[1] <= pair[0]);
        let index = index(&corners, dimension);
        lines.push(Line {
            code,
            k,
            index,
            spanning: index == 1,
            bracket,
            predicted,
            measured,
            level,
            exact,
            within: (measured - predicted).abs() < TOLERANCE,
            monotone,
            head: terms.iter().take(HEAD).copied().collect(),
        });
    }
    Family {
        label: label.to_string(),
        base,
        lines,
    }
}

fn sponge() -> Code {
    let filled: Vec<Vec<u8>> = residue_corners(3, 3)
        .into_iter()
        .filter(|corner| corner.iter().filter(|&&r| r == 1).count() <= 1)
        .collect();
    corners_to_code(&filled, 3, 3)
}

fn base6() -> Vec<Code> {
    let samples: [[[u8; 2]; 8]; 2] = [
        [
            [0, 0],
            [1, 1],
            [2, 3],
            [3, 2],
            [4, 5],
            [5, 4],
            [1, 0],
            [0, 1],
        ],
        [
            [0, 0],
            [2, 0],
            [4, 0],
            [0, 3],
            [1, 1],
            [5, 5],
            [1, 2],
            [2, 1],
        ],
    ];
    samples
        .iter()
        .map(|sample| {
            let filled: Vec<Vec<u8>> = sample.iter().map(|corner| corner.to_vec()).collect();
            corners_to_code(&filled, 2, 6)
        })
        .collect()
}

pub fn report() {
    let families = [
        survey(2, 2, (1..16u128).map(Code::from), "base 2 D 2"),
        survey(2, 3, (1..256u128).map(Code::from), "base 2 D 3"),
        survey(3, 2, (1..512u128).map(Code::from), "base 3 D 2"),
        survey(3, 3, [sponge()].into_iter(), "menger sponge"),
        survey(6, 2, base6().into_iter(), "base 6 samples"),
    ];
    println!("coprime census: designs with k >= 2, point budget {BUDGET}, tolerance {TOLERANCE}, exact identity at n = {EXACT_LEVEL}");
    println!(
        "{:<16}{:>9}{:>10}{:>10}{:>19}{:>9}{:>13}",
        "family", "designs", "spanning", "distinct", "distinct spanning", "flagged", "exact fails"
    );
    for family in &families {
        println!(
            "{:<16}{:>9}{:>10}{:>10}{:>19}{:>9}{:>13}",
            family.label,
            family.lines.len(),
            family.spanning(),
            family.distinct(false),
            family.distinct(true),
            family.flagged(),
            family.exact_failures()
        );
    }
    println!(
        "{:<16}{:>9}{:>10}{:>10}{:>19}{:>9}{:>13}",
        "total",
        families.iter().map(|f| f.lines.len()).sum::<usize>(),
        families.iter().map(Family::spanning).sum::<usize>(),
        "-",
        "-",
        families.iter().map(Family::flagged).sum::<usize>(),
        families.iter().map(Family::exact_failures).sum::<usize>()
    );
    let by = |keep: &dyn Fn(&Line, usize) -> bool| {
        families
            .iter()
            .flat_map(|f| f.lines.iter().map(move |line| (line, f.base)))
            .filter(|(line, base)| keep(line, *base))
            .count()
    };
    println!(
        "spanning by dimension log_q(k): above {}, exactly one {}, below {}; index 3 at k = q {}",
        by(&|line, base| line.spanning && line.k > base),
        by(&|line, base| line.spanning && line.k == base),
        by(&|line, base| line.spanning && line.k < base),
        by(&|line, base| line.index == 3 && line.k == base)
    );
    println!(
        "spanning lines whose gap to the predicted density widens at some level: {}",
        by(&|line, _| line.spanning && !line.monotone)
    );
    for family in &families[3..] {
        for line in &family.lines {
            let corners =
                code_to_corners(line.code, 2 + usize::from(family.base == 3), family.base)
                    .expect("the named code fits its cells");
            let counts: Vec<String> = squarefree_divisors(family.base)
                .iter()
                .filter(|&&e| e > 1)
                .map(|&e| format!("k_{e} {}", hits(&corners, e)))
                .collect();
            println!(
                "{}: code {} k {} {} bracket {:.4} pred {:.6} meas {:.6} n {} terms {:?}",
                family.label,
                line.code,
                line.k,
                counts.join(" "),
                line.bracket,
                line.predicted,
                line.measured,
                line.level,
                line.head
            );
        }
    }
    for k in [2usize, 3] {
        let level = families
            .iter()
            .flat_map(|f| f.lines.iter())
            .find(|line| line.k == k)
            .map(|line| line.level)
            .expect("a design with that k exists");
        println!("deepest level at k = {k}: n = {level}");
    }
}
