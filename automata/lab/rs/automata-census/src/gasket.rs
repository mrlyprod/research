use crate::rules::{output, single_seed, Diagram};
use mrly_core::tensor::Tensor;

const CANDIDATES: [usize; 4] = [7, 11, 13, 14];
const DEPTH: usize = 8;

fn tile(code: usize, level: usize) -> Vec<Vec<u8>> {
    let mut side = 1usize;
    let mut cells = vec![vec![1u8]];
    for _ in 0..level {
        let next = side * 2;
        let mut out = vec![vec![0u8; next]; next];
        for row in 0..side {
            for column in 0..side {
                if cells[row][column] == 1 {
                    for i in 0..4usize {
                        if (code >> i) & 1 == 1 {
                            out[2 * row + i / 2][2 * column + i % 2] = 1;
                        }
                    }
                }
            }
        }
        side = next;
        cells = out;
    }
    cells
}

fn crate_tile(code: usize, level: usize) -> Vec<Vec<u8>> {
    let seed = Tensor::of(
        (0..4).map(|i| ((code >> i) & 1) as u8).collect(),
        vec![2, 2],
    )
    .expect("four bits fill a two by two");
    let grid = seed.fractal(level);
    let side = 1usize << level;
    (0..side)
        .map(|row| {
            (0..side)
                .map(|column| grid.at(grid.index(&[row, column])) as u8)
                .collect()
        })
        .collect()
}

fn matches(diagram: &Diagram, code: usize, reading: &str) -> bool {
    for level in 1..=DEPTH {
        let grid = tile(code, level);
        let side = 1usize << level;
        for (t, row) in grid.iter().enumerate() {
            for (j, &cell) in row.iter().enumerate() {
                let offset = match reading {
                    "right" => j as i64,
                    "left" => j as i64 - (side as i64 - 1),
                    _ => 2 * j as i64 - t as i64,
                };
                if diagram.signed(t, offset) != cell {
                    return false;
                }
            }
        }
    }
    true
}

fn hits(rule: usize, reading: &str) -> Vec<usize> {
    let diagram = single_seed(rule, 1 << DEPTH);
    CANDIDATES
        .into_iter()
        .filter(|&code| matches(&diagram, code, reading))
        .collect()
}

fn polynomials(rows: usize) -> Vec<Vec<u8>> {
    let width = 2 * rows + 1;
    let mut out = Vec::with_capacity(rows);
    let mut row = vec![0u8; width];
    row[0] = 1;
    for _ in 0..rows {
        out.push(row.clone());
        let mut next = vec![0u8; width];
        for i in 0..width {
            next[i] = row[i];
            if i >= 1 {
                next[i] ^= row[i - 1];
            }
            if i >= 2 {
                next[i] ^= row[i - 2];
            }
        }
        row = next;
    }
    out
}

fn population_totals(depth: usize) -> (Vec<u64>, Vec<u64>) {
    let rows = 1usize << depth;
    let width = 2 * rows + 1;
    let mut row = vec![0u8; width];
    row[0] = 1;
    let mut totals = Vec::with_capacity(depth + 1);
    let mut pairs = Vec::with_capacity(depth + 1);
    let mut running = 0u64;
    let mut adjacent = 0u64;
    let mut next_mark = 1usize;
    for t in 0..rows {
        running += row.iter().filter(|&&c| c == 1).count() as u64;
        adjacent += row.windows(2).filter(|w| w[0] == 1 && w[1] == 1).count() as u64;
        if t + 1 == next_mark {
            totals.push(running);
            pairs.push(adjacent);
            next_mark *= 2;
        }
        let mut next = vec![0u8; width];
        for i in 0..width {
            next[i] = row[i];
            if i >= 1 {
                next[i] ^= row[i - 1];
            }
            if i >= 2 {
                next[i] ^= row[i - 2];
            }
        }
        row = next;
    }
    (totals, pairs)
}

pub fn report() {
    println!("THE GASKET IDENTITY");
    for level in 1..=DEPTH {
        for code in CANDIDATES {
            assert_eq!(
                tile(code, level),
                crate_tile(code, level),
                "the study renderer and Tensor::fractal disagree at code {code} level {level}"
            );
        }
    }
    println!("both renderers agree on codes {CANDIDATES:?} at levels 1..{DEPTH}");
    for (rule, reading, want) in [
        (60usize, "right", 13usize),
        (102, "left", 14),
        (90, "shear", 13),
    ] {
        let found = hits(rule, reading);
        println!("rule {rule} read {reading}: matching codes {found:?}");
        assert_eq!(
            found,
            vec![want],
            "rule {rule} does not match exactly bang dim 2, code {want}"
        );
    }
    println!("rule 60 is bang dim 2, code 13; rule 102 is bang dim 2, code 14; rule 90 sheared by j = (t+i)/2 is bang dim 2, code 13, cell for cell to level {DEPTH}");
    println!("RULE 150");
    let rows = polynomials(129);
    let diagram = single_seed(150, 128);
    for (t, row) in rows.iter().enumerate() {
        for (j, &cell) in row.iter().enumerate().take(t + 1) {
            assert_eq!(
                diagram.signed(t, j as i64 - t as i64),
                cell,
                "the GF(2) row polynomial and the evolved diagram disagree at ({t},{j})"
            );
        }
    }
    println!("the row polynomial (1+x+x^2)^t matches the evolved diagram on rows 0..128");
    for t in 0..64usize {
        for j in 0..rows[2 * t].len() {
            let want = if j % 2 == 0 { rows[t][j / 2] } else { 0 };
            assert_eq!(
                rows[2 * t][j],
                want,
                "row {} is not row {t} spread by two",
                2 * t
            );
        }
        for j in 0..rows[2 * t + 1].len() {
            let mut want = rows[2 * t][j];
            if j >= 1 {
                want ^= rows[2 * t][j - 1];
            }
            if j >= 2 {
                want ^= rows[2 * t][j - 2];
            }
            assert_eq!(
                rows[2 * t + 1][j],
                want,
                "row {} is not row {} xor its two unit shifts",
                2 * t + 1,
                2 * t
            );
        }
    }
    println!("(1+x+x^2)^(2n) = (1+x^2+x^4)^n on rows 0..128: row 2t is row t spread by two, row 2t+1 is row 2t xor its two unit shifts");
    let populations: Vec<usize> = (0..=64)
        .map(|t| rows[t].iter().filter(|&&c| c == 1).count())
        .collect();
    println!("row populations t = 0..64");
    println!(
        "{}",
        populations
            .iter()
            .map(|p| p.to_string())
            .collect::<Vec<_>>()
            .join(",")
    );
    println!("this is OEIS A071053, the number of ON cells at generation t of rule 150 from a single ON cell");
    let (totals, pairs) = population_totals(12);
    let root = 5f64.sqrt();
    let constant = (5.0 + 3.0 * root) / 10.0;
    println!("k P(k) log2(P)/k P/(1+sqrt5)^k");
    for (k, total) in totals.iter().enumerate().skip(1) {
        println!(
            "{} {} {:.9} {:.9}",
            k,
            total,
            (*total as f64).log2() / k as f64,
            *total as f64 / (1.0 + root).powi(k as i32)
        );
    }
    println!("k P(k) B(k)");
    for k in 0..totals.len() {
        println!("{k} {} {}", totals[k], pairs[k]);
    }
    for k in 0..totals.len() - 1 {
        assert_eq!(
            totals[k + 1],
            4 * totals[k] - 2 * pairs[k],
            "P(k+1) = 4 P(k) - 2 B(k) breaks at k = {k}"
        );
        assert_eq!(
            pairs[k + 1],
            2 * totals[k] - 2 * pairs[k],
            "B(k+1) = 2 P(k) - 2 B(k) breaks at k = {k}"
        );
    }
    println!("the two-term system P(k+1) = 4 P(k) - 2 B(k), B(k+1) = 2 P(k) - 2 B(k) holds at k = 0..{}, matrix [[4,-2],[2,-2]] of trace 2 and determinant -4", totals.len() - 2);
    for k in 2..totals.len() {
        assert_eq!(
            totals[k],
            2 * totals[k - 1] + 4 * totals[k - 2],
            "the cumulative population breaks P(k+1) = 2 P(k) + 4 P(k-1) at k = {k}"
        );
    }
    let mut exact = vec![1u64, 4];
    while exact.len() <= totals.len() {
        let n = exact.len();
        exact.push(2 * exact[n - 1] + 4 * exact[n - 2]);
    }
    assert_eq!(
        totals[..],
        exact[..totals.len()],
        "the cumulative population leaves the closed form"
    );
    for (k, total) in totals.iter().enumerate() {
        let e = k as i32;
        let closed = constant * (1.0 + root).powi(e) + (1.0 - constant) * (1.0 - root).powi(e);
        assert!(
            (closed - *total as f64).abs() < 1e-3,
            "the closed form misses P({e}) = {total}"
        );
    }
    let mut fib = vec![0u64, 1];
    while fib.len() < totals.len() + 3 {
        let n = fib.len();
        fib.push(fib[n - 1] + fib[n - 2]);
    }
    for (k, total) in totals.iter().enumerate() {
        assert_eq!(*total, (1u64 << k) * fib[k + 2], "P({k}) is not 2^k F(k+2)");
    }
    println!("P(k) = 2 P(k-1) + 4 P(k-2) with P(0) = 1, P(1) = 4, so P(k) = c (1+sqrt5)^k + (1-c) (1-sqrt5)^k with c = (5+3 sqrt5)/10 = {constant:.12}");
    println!("equivalently P(k) = 2^k F(k+2), asserted at k = 0..12");
    println!(
        "the exponent is exactly log2(1+sqrt5) = {:.7}",
        (1.0 + root).log2()
    );
    println!("P(k) over rows 0..2^k - 1 is OEIS A087206, whose %N carries the recurrence and whose %F carries 2^n Fibonacci(n+2) and the (1 +- sqrt5) form; A071053 %F states Sum_{{k = 0..2^n - 1}} a(k) = A087206(n). The study re-derives them.");
    assert_eq!(output(150, 1, 1, 1), 1, "rule 150 is not the xor rule");
}
