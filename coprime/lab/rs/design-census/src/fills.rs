use crate::quasi::{degree, fit, fraction, leading, text};
use crate::tables::write_csv;
use mrly_math::bang::baseq::{
    canonical, distinct_designs, group, group_order, orbit, representatives, WALK_LIMIT,
};
use mrly_math::bang::factory::{corners_to_code, levels_code, residue_corners};
use mrly_math::bang::universe;
use mrly_math::bang::Code;
use mrly_math::counts::{fill, void};
use mrly_math::rules::render;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

const SIDES: usize = 12;
const CASES: [(usize, usize); 4] = [(2, 2), (2, 3), (3, 2), (3, 3)];
const SHOWN: usize = 6;

pub fn cell_index(cell: &[u8], base: usize) -> usize {
    cell.iter().fold(0, |acc, &r| acc * base + r as usize)
}

fn named(base: usize, dimension: usize) -> Vec<(&'static str, Code)> {
    let cells = residue_corners(dimension, base);
    let select = |keep: &dyn Fn(&[u8]) -> bool| {
        let filled: Vec<Vec<u8>> = cells.iter().filter(|cell| keep(cell)).cloned().collect();
        corners_to_code(&filled, dimension, base)
    };
    let free = if dimension == 2 { 0 } else { dimension - 1 };
    vec![
        ("carpet", levels_code(dimension, base, &[0, 1])),
        (
            "net",
            select(&|cell| cell.iter().map(|&r| r as usize).sum::<usize>() + 1 >= dimension),
        ),
        ("void", select(&|cell| cell.iter().all(|&r| r == cell[0]))),
        (
            "tree",
            select(&|cell| (0..dimension).filter(|&a| a != free).all(|a| cell[a] == 0)),
        ),
    ]
}

fn burnside(base: usize, dimension: usize) -> u128 {
    distinct_designs(base, dimension).expect("the Burnside average is an integer")
}

pub fn orbits_report() {
    println!("base 2 cube group: order 2^D D!, designs up to symmetry two ways");
    for dimension in 1..=4 {
        let walk = representatives(2, dimension)
            .expect("the walk stays under the code limit")
            .len();
        let order =
            group_order(2, dimension).expect("dimension four is under the factorial ceiling");
        let burnside = burnside(2, dimension);
        println!("D {dimension}: order {order} orbit walk {walk} Burnside {burnside}");
    }
}

pub struct Row {
    pub base: usize,
    pub dimension: usize,
    pub code: Code,
    pub label: String,
    pub orbit: usize,
    pub popcount: u32,
    pub gf2: Option<i32>,
    pub gfq: i64,
    pub mobius: i64,
    pub poly: String,
    pub degree: String,
    pub lead: String,
    pub fills: Vec<u128>,
    pub voids: Vec<u128>,
}

fn rendered(code: Code, number: usize, dimension: usize, base: usize) -> u128 {
    let filled: Vec<Vec<u8>> = residue_corners(dimension, base)
        .into_iter()
        .filter(|r| code.get() >> cell_index(r, base) & 1 == 1)
        .collect();
    let tile = render(&filled, number, dimension, base).expect("the tile renders");
    u128::from(tile.sum())
}

fn live_degree(cells: &[Vec<u8>], live: impl Fn(usize) -> bool) -> i64 {
    cells
        .iter()
        .enumerate()
        .filter(|(index, _)| live(*index))
        .map(|(_, cell)| cell.iter().map(|&r| i64::from(r)).sum::<i64>())
        .max()
        .unwrap_or(-1)
}

fn mobius_degree(code: Code, cells: &[Vec<u8>], dimension: usize, base: usize) -> i64 {
    let mut coeff: Vec<i64> = (0..cells.len())
        .map(|i| (code.get() >> i & 1) as i64)
        .collect();
    for axis in 0..dimension {
        for value in (1..base as u8).rev() {
            for (index, cell) in cells.iter().enumerate() {
                if cell[axis] == value {
                    let mut lower = cell.clone();
                    lower[axis] -= 1;
                    coeff[index] -= coeff[cell_index(&lower, base)];
                }
            }
        }
    }
    live_degree(cells, |index| coeff[index] != 0)
}

fn inverse_modular(value: i64, modulus: i64) -> i64 {
    (1..modulus)
        .find(|candidate| candidate * value.rem_euclid(modulus) % modulus == 1)
        .expect("the pivot is invertible in a prime field")
}

fn inverse_vandermonde(q: usize) -> Vec<Vec<i64>> {
    let modulus = q as i64;
    let mut rows: Vec<Vec<i64>> = (0..q)
        .map(|i| {
            let mut row: Vec<i64> = (0..q)
                .map(|j| (i as i64).pow(j as u32).rem_euclid(modulus))
                .collect();
            row.extend((0..q).map(|c| i64::from(c == i)));
            row
        })
        .collect();
    for column in 0..q {
        let pivot = (column..q)
            .find(|&r| rows[r][column] != 0)
            .expect("the Vandermonde matrix is nonsingular over a prime field");
        rows.swap(column, pivot);
        let scale = inverse_modular(rows[column][column], modulus);
        for value in rows[column].iter_mut() {
            *value = *value * scale % modulus;
        }
        for r in 0..q {
            if r != column && rows[r][column] != 0 {
                let factor = rows[r][column];
                let pivot_row = rows[column].clone();
                for (value, above) in rows[r].iter_mut().zip(&pivot_row) {
                    *value = (*value - factor * above).rem_euclid(modulus);
                }
            }
        }
    }
    rows.into_iter().map(|row| row[q..].to_vec()).collect()
}

fn gfq_degree(code: Code, cells: &[Vec<u8>], dimension: usize, q: usize) -> i64 {
    let inverse = inverse_vandermonde(q);
    let mut coeff: Vec<i64> = (0..cells.len())
        .map(|i| (code.get() >> i & 1) as i64)
        .collect();
    for axis in 0..dimension {
        let mut next = vec![0i64; cells.len()];
        for cell in cells.iter().filter(|cell| cell[axis] == 0) {
            let line: Vec<i64> = (0..q)
                .map(|t| {
                    let mut key = cell.clone();
                    key[axis] = t as u8;
                    coeff[cell_index(&key, q)]
                })
                .collect();
            for (exponent, row) in inverse.iter().enumerate() {
                let acc: i64 = row.iter().zip(&line).map(|(a, b)| a * b).sum();
                let mut key = cell.clone();
                key[axis] = exponent as u8;
                next[cell_index(&key, q)] = acc.rem_euclid(q as i64);
            }
        }
        coeff = next;
    }
    live_degree(cells, |index| coeff[index] != 0)
}

fn row(base: usize, dimension: usize, code: Code, orbit: usize, label: String) -> Row {
    let cells = residue_corners(dimension, base);
    let count =
        |n: usize, level: u32| fill(code, n, dimension, level, base).expect("the fill counts");
    let fills: Vec<u128> = (1..=SIDES).map(|n| count(n, 1)).collect();
    let voids: Vec<u128> = (1..=SIDES)
        .map(|n| void(code, n, dimension, 1, base).expect("the void counts"))
        .collect();
    for n in 1..=SIDES {
        assert!(
            rendered(code, n, dimension, base) == fills[n - 1],
            "the closed form matches the rendered tile"
        );
    }
    for n in [base, base + 1, 2 * base] {
        for level in [2u32, 3] {
            assert!(
                count(n, level) == fills[n - 1].pow(level),
                "the level law holds"
            );
        }
    }
    let (polys, collapses) = fit(code, dimension, base);
    let top = polys.iter().filter_map(|poly| degree(poly)).max();
    let poly = if collapses {
        text(&polys[0])
    } else {
        (0..base)
            .map(|class| format!("n\u{2261}{class}(mod {base}): {}", text(&polys[class])))
            .collect::<Vec<String>>()
            .join(" ; ")
    };
    Row {
        base,
        dimension,
        code,
        label,
        orbit,
        popcount: code.get().count_ones(),
        gf2: (base == 2).then(|| universe::degree(code, dimension)),
        gfq: gfq_degree(code, &cells, dimension, base),
        mobius: mobius_degree(code, &cells, dimension, base),
        poly,
        degree: top.map_or(String::from("-oo"), |d| d.to_string()),
        lead: fraction(&leading(&polys[0])),
        fills,
        voids,
    }
}

fn case(base: usize, dimension: usize) -> (Vec<Row>, bool) {
    let group = group(base, dimension);
    let full = base.pow(dimension as u32) <= WALK_LIMIT;
    if full {
        let labels: BTreeMap<Code, &str> = named(base, dimension)
            .into_iter()
            .map(|(name, code)| {
                (
                    canonical(&group, code).expect("the group is not empty"),
                    name,
                )
            })
            .collect();
        let rows = representatives(base, dimension)
            .expect("the walk stays under the code limit")
            .into_iter()
            .map(|(code, size)| {
                let label = labels.get(&code).copied().unwrap_or_default().to_string();
                row(base, dimension, code, size, label)
            })
            .collect();
        (rows, true)
    } else {
        let rows = named(base, dimension)
            .into_iter()
            .map(|(name, code)| {
                let representative = canonical(&group, code).expect("the group is not empty");
                let size = orbit(&group, representative).len();
                row(base, dimension, representative, size, name.to_string())
            })
            .collect();
        (rows, false)
    }
}

fn header() -> Vec<String> {
    let mut out: Vec<String> = [
        "base",
        "dim",
        "rep_code",
        "label",
        "orbit_size",
        "popcount",
        "gf2_degree",
        "gfq_degree",
        "intmobius_degree",
        "fill_poly",
        "fill_deg",
        "fill_lead",
    ]
    .iter()
    .map(|name| (*name).to_string())
    .collect();
    out.extend((1..=SIDES).map(|n| format!("fill_n{n}")));
    out.extend((1..=SIDES).map(|n| format!("void_n{n}")));
    out
}

fn record(row: &Row) -> Vec<String> {
    let mut out = vec![
        row.base.to_string(),
        row.dimension.to_string(),
        row.code.to_string(),
        row.label.clone(),
        row.orbit.to_string(),
        row.popcount.to_string(),
        row.gf2.map(|d| d.to_string()).unwrap_or_default(),
        row.gfq.to_string(),
        row.mobius.to_string(),
        row.poly.clone(),
        row.degree.clone(),
        row.lead.clone(),
    ];
    out.extend(row.fills.iter().map(u128::to_string));
    out.extend(row.voids.iter().map(u128::to_string));
    out
}

pub fn report(path: &Path) {
    println!("fill census: one row per cube-group orbit, sides 1..{SIDES}, every row rendered cell by cell");
    let mut rows: Vec<Row> = Vec::new();
    for (base, dimension) in CASES {
        let (case_rows, full) = case(base, dimension);
        let expected = burnside(base, dimension);
        let scope = if full { "full walk" } else { "named only" };
        let polys: BTreeSet<&str> = case_rows.iter().map(|row| row.poly.as_str()).collect();
        println!(
            "base {base} D {dimension}: {} designs, Burnside {expected}, {} distinct fill polynomials, {scope}",
            case_rows.len(),
            polys.len()
        );
        for row in case_rows.iter().filter(|row| !row.label.is_empty()) {
            println!(
                "  {:<6} code {:>9} orbit {:>3} popcount {:>2} degree {} lead {:<5} fill n1..n{SHOWN} {:?}",
                row.label,
                row.code,
                row.orbit,
                row.popcount,
                row.degree,
                row.lead,
                &row.fills[..SHOWN]
            );
        }
        rows.extend(case_rows);
    }
    let labels = rows.iter().filter(|row| !row.label.is_empty()).count();
    println!(
        "{} designs, {} columns, {labels} labelled rows, written sequences.csv",
        rows.len(),
        header().len()
    );
    let records: Vec<Vec<String>> = rows.iter().map(record).collect();
    write_csv(path, &header(), &records);
}
