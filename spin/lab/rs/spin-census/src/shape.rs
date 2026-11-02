use crate::design::{plane, BASE};
use crate::orbit;
use std::collections::HashMap;

pub struct Classes {
    pub index: Vec<u32>,
    pub count: usize,
    pub cells: usize,
}

pub fn d4(n: usize) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = Vec::new();
    for flip in 0..2 {
        for turn in 0..4 {
            let map: Vec<usize> = (0..n * n)
                .map(|f| {
                    let (mut r, mut c) = (f / n, f % n);
                    if flip == 1 {
                        std::mem::swap(&mut r, &mut c);
                    }
                    for _ in 0..turn {
                        let next = (c, n - 1 - r);
                        r = next.0;
                        c = next.1;
                    }
                    r * n + c
                })
                .collect();
            if !out.contains(&map) {
                out.push(map);
            }
        }
    }
    out
}

pub fn oct(n: usize) -> Vec<Vec<usize>> {
    let axes = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let mut out: Vec<Vec<usize>> = Vec::new();
    for axis in axes {
        for signs in 0..8u32 {
            let map: Vec<usize> = (0..n * n * n)
                .map(|f| {
                    let p = [f / (n * n), (f / n) % n, f % n];
                    let mut q = [p[axis[0]], p[axis[1]], p[axis[2]]];
                    for (k, v) in q.iter_mut().enumerate() {
                        if signs >> k & 1 == 1 {
                            *v = n - 1 - *v;
                        }
                    }
                    (q[0] * n + q[1]) * n + q[2]
                })
                .collect();
            if !out.contains(&map) {
                out.push(map);
            }
        }
    }
    out
}

pub fn classes(cells: usize, group: &[Vec<usize>]) -> Classes {
    let mut index = vec![u32::MAX; cells * cells];
    let mut count = 0u32;
    for j in 0..cells {
        for l in j..cells {
            if index[j * cells + l] != u32::MAX {
                continue;
            }
            for map in group {
                let (a, b) = (map[j], map[l]);
                index[a * cells + b] = count;
                index[b * cells + a] = count;
            }
            count += 1;
        }
    }
    Classes {
        index,
        count: count as usize,
        cells,
    }
}

pub fn burnside(cells: usize, group: &[Vec<usize>]) -> u128 {
    let mut total = 0u128;
    for map in group {
        let mut seen = vec![false; cells];
        let mut cycles = 0u32;
        for j in 0..cells {
            if seen[j] {
                continue;
            }
            cycles += 1;
            let mut x = j;
            while !seen[x] {
                seen[x] = true;
                x = map[x];
            }
        }
        total += 1u128 << cycles;
    }
    total / group.len() as u128
}

pub fn census(points: &[usize], table: &Classes) -> Vec<u32> {
    let mut out = vec![0u32; table.count];
    for a in 0..points.len() {
        for b in a..points.len() {
            out[table.index[points[a] * table.cells + points[b]] as usize] += 1;
        }
    }
    out
}

pub fn canon(bits: u128, cells: usize, group: &[Vec<usize>]) -> u128 {
    let mut best = u128::MAX;
    for map in group {
        let mut x = 0u128;
        for (j, &target) in map.iter().enumerate().take(cells) {
            if bits >> j & 1 == 1 {
                x |= 1u128 << target;
            }
        }
        best = best.min(x);
    }
    best
}

fn cells_of(bits: u128, cells: usize) -> Vec<usize> {
    (0..cells).filter(|j| bits >> j & 1 == 1).collect()
}

fn kron_plane(bits: u128, q: usize) -> Vec<usize> {
    let base = cells_of(bits, q * q);
    let side = q * q;
    let mut out = Vec::with_capacity(base.len() * base.len());
    for &p in &base {
        for &s in &base {
            let r = (p / q) * q + s / q;
            let c = (p % q) * q + s % q;
            out.push(r * side + c);
        }
    }
    out.sort_unstable();
    out
}

fn picture(bits: u128) -> String {
    (0..3)
        .map(|r| {
            (0..3)
                .map(|c| {
                    if bits >> (r * 3 + c) & 1 == 1 {
                        '#'
                    } else {
                        '.'
                    }
                })
                .collect::<String>()
        })
        .collect::<Vec<_>>()
        .join("/")
}

fn solve(matrix: &mut [Vec<f64>], rhs: &mut [f64]) -> Option<Vec<f64>> {
    let n = rhs.len();
    for col in 0..n {
        let mut pivot = col;
        for row in col..n {
            if matrix[row][col].abs() > matrix[pivot][col].abs() {
                pivot = row;
            }
        }
        if matrix[pivot][col].abs() < 1e-9 {
            return None;
        }
        matrix.swap(col, pivot);
        rhs.swap(col, pivot);
        let lead = matrix[col].clone();
        let lead_rhs = rhs[col];
        for (row, (line, value)) in matrix.iter_mut().zip(rhs.iter_mut()).enumerate() {
            if row == col {
                continue;
            }
            let f = line[col] / lead[col];
            for (cell, &p) in line[col..n].iter_mut().zip(&lead[col..n]) {
                *cell -= f * p;
            }
            *value -= f * lead_rhs;
        }
    }
    Some((0..n).map(|i| rhs[i] / matrix[i][i]).collect())
}

fn rank(rows: &[Vec<f64>], tolerance: f64) -> usize {
    let mut work: Vec<Vec<f64>> = rows.to_vec();
    let width = work[0].len();
    let mut got = 0usize;
    for col in 0..width {
        let Some(pivot) = (got..work.len()).find(|&row| work[row][col].abs() > tolerance) else {
            continue;
        };
        work.swap(got, pivot);
        let lead = work[got].clone();
        for (row, line) in work.iter_mut().enumerate() {
            if row == got {
                continue;
            }
            let f = line[col] / lead[col];
            for (cell, &p) in line[col..width].iter_mut().zip(&lead[col..width]) {
                *cell -= f * p;
            }
        }
        got += 1;
        if got == work.len() {
            break;
        }
    }
    got
}

pub fn control() {
    println!();
    println!("SHAPE READING, THE 511 AS CONTROL");
    let g1 = d4(BASE);
    let g2 = d4(BASE * BASE);
    let t1 = classes(9, &g1);
    let t2 = classes(81, &g2);
    println!(
        "  pair classes level 1 {} level 2 {}; Burnside orbits {} of which {} nonempty",
        t1.count,
        t2.count,
        burnside(9, &g1),
        burnside(9, &g1) - 1
    );
    let mut mismatch = 0usize;
    for code in 1..512u128 {
        let grid = plane(code, BASE, 2);
        let live: Vec<usize> = (0..grid.size()).filter(|&f| grid.at(f) != 0).collect();
        if live != kron_plane(code, BASE) {
            mismatch += 1;
        }
    }
    println!("  level-2 renders disagreeing with the Kronecker square: {mismatch}");
    let mut reps: Vec<u128> = (1..512u128).map(|c| canon(c, 9, &g1)).collect();
    reps.sort_unstable();
    reps.dedup();
    println!("  nonempty orbits by canonical form: {}", reps.len());
    let mut first: HashMap<Vec<u32>, Vec<u128>> = HashMap::new();
    let mut second: HashMap<Vec<u32>, Vec<u128>> = HashMap::new();
    for &code in &reps {
        first
            .entry(census(&cells_of(code, 9), &t1))
            .or_default()
            .push(code);
        second
            .entry(census(&kron_plane(code, BASE), &t2))
            .or_default()
            .push(code);
    }
    println!(
        "  distinct pair censuses: level 1 {} level 2 {}",
        first.len(),
        second.len()
    );
    let mut clashes: Vec<Vec<u128>> = first.values().filter(|v| v.len() > 1).cloned().collect();
    clashes.sort();
    let spin1 = orbit::census(1, 1024, 12);
    let spin2 = orbit::census(2, 768, 12);
    for group in &clashes {
        let (a, b) = (group[0], group[1]);
        println!(
            "    level-1 clash {a} {} against {b} {} fill {} spectrum gap level 1 {:.2e} level 2 {:.2e}",
            picture(a),
            picture(b),
            a.count_ones(),
            orbit::gap(&spin1[a as usize], &spin1[b as usize]),
            orbit::gap(&spin2[a as usize], &spin2[b as usize])
        );
    }
    let mut buckets: Vec<usize> = Vec::new();
    for code in 1..512usize {
        if !buckets
            .iter()
            .any(|&head| orbit::agree(&spin1[head], &spin1[code], 1e-9))
        {
            buckets.push(code);
        }
    }
    println!(
        "  spectra from level 1 alone: {} against {} pair censuses",
        buckets.len(),
        first.len()
    );
    let rows: Vec<Vec<f64>> = (1..512usize)
        .map(|code| {
            census(&cells_of(code as u128, 9), &t1)
                .iter()
                .map(|&v| v as f64)
                .collect()
        })
        .collect();
    let mut chosen: Vec<usize> = Vec::new();
    for i in 0..rows.len() {
        let mut trial: Vec<Vec<f64>> = chosen.iter().map(|&j| rows[j].clone()).collect();
        trial.push(rows[i].clone());
        if rank(&trial, 1e-9) == trial.len() {
            chosen.push(i);
        }
        if chosen.len() == t1.count {
            break;
        }
    }
    let mut weights: Vec<Vec<f64>> = Vec::new();
    let mut worst = 0.0f64;
    let columns: Vec<Vec<f64>> = (0..13)
        .map(|m| spin1.iter().map(|row| row[m]).collect())
        .collect();
    for (m, column) in columns.iter().enumerate() {
        let mut matrix: Vec<Vec<f64>> = chosen.iter().map(|&j| rows[j].clone()).collect();
        let mut rhs: Vec<f64> = chosen.iter().map(|&j| column[j + 1]).collect();
        let Some(w) = solve(&mut matrix, &mut rhs) else {
            println!("  the level-1 system is singular at order {m}");
            return;
        };
        let scale = column[1..512]
            .iter()
            .map(|v| v.abs())
            .fold(0.0f64, f64::max)
            .max(1e-300);
        for (i, row) in rows.iter().enumerate() {
            let predicted: f64 = row.iter().zip(&w).map(|(a, b)| a * b).sum();
            worst = worst.max((predicted - column[i + 1]).abs() / scale);
        }
        let top = w.iter().fold(0.0f64, |a, b| a.max(b.abs())).max(1e-300);
        weights.push(w.iter().map(|v| v / top).collect());
    }
    let odd: Vec<Vec<f64>> = weights
        .iter()
        .enumerate()
        .filter(|(m, _)| m % 2 == 1)
        .map(|(_, w)| w.clone())
        .collect();
    let even: Vec<Vec<f64>> = weights
        .iter()
        .enumerate()
        .filter(|(m, _)| m % 2 == 0)
        .map(|(_, w)| w.clone())
        .collect();
    println!(
        "  the level-1 spectrum is a linear functional of the pair census: independent censuses solved {}, worst relative residual over 511 codes by 13 orders {:.2e}, coefficient rank {} of {}",
        chosen.len(),
        worst,
        rank(&weights, 1e-8),
        t1.count
    );
    println!(
        "  half-turning one member of a pair caps the {} odd orders at the 3-dimensional antisymmetric part and they reach {}; the {} even orders reach {}",
        odd.len(),
        rank(&odd, 1e-8),
        even.len(),
        rank(&even, 1e-8)
    );
}

fn chunk_tables(cells: usize, group: &[Vec<usize>]) -> Vec<Vec<Vec<u32>>> {
    let chunks = cells.div_ceil(8);
    group
        .iter()
        .map(|map| {
            (0..chunks)
                .map(|c| {
                    (0..256u32)
                        .map(|v| {
                            let mut x = 0u32;
                            for b in 0..8 {
                                let j = c * 8 + b;
                                if j < cells && v >> b & 1 == 1 {
                                    x |= 1u32 << map[j];
                                }
                            }
                            x
                        })
                        .collect()
                })
                .collect()
        })
        .collect()
}

fn image(code: u32, tabs: &[Vec<u32>]) -> u32 {
    let mut x = 0u32;
    for (c, tab) in tabs.iter().enumerate() {
        x |= tab[(code >> (8 * c) & 255) as usize];
    }
    x
}

fn pack(counts: &[u32], width: u32) -> (u128, u128) {
    let per = 128 / width as usize;
    let (mut lo, mut hi) = (0u128, 0u128);
    for (c, &v) in counts.iter().enumerate() {
        assert!(v < 1 << width, "a class count overflows the key");
        if c < per {
            lo |= (v as u128) << (c as u32 * width);
        } else {
            hi |= (v as u128) << ((c - per) as u32 * width);
        }
    }
    (lo, hi)
}

fn pair_list(table: &Classes) -> Vec<(u32, u16)> {
    let mut out = Vec::new();
    for j in 0..table.cells {
        for l in j + 1..table.cells {
            out.push((
                (1u32 << j) | (1u32 << l),
                table.index[j * table.cells + l] as u16,
            ));
        }
    }
    out
}

fn level_one_sweep(
    cells: usize,
    group: &[Vec<usize>],
    table: &Classes,
    width: u32,
) -> Vec<(u128, u128, u32)> {
    let tabs = chunk_tables(cells, group);
    let pairs = pair_list(table);
    let diag: Vec<u16> = (0..cells)
        .map(|j| table.index[j * cells + j] as u16)
        .collect();
    let mut out: Vec<(u128, u128, u32)> = Vec::new();
    let mut counts = vec![0u32; table.count];
    for code in 1..(1u32 << cells) {
        if tabs.iter().any(|t| image(code, t) < code) {
            continue;
        }
        counts.iter_mut().for_each(|v| *v = 0);
        for j in 0..cells {
            if code >> j & 1 == 1 {
                counts[diag[j] as usize] += 1;
            }
        }
        for &(mask, class) in &pairs {
            if code & mask == mask {
                counts[class as usize] += 1;
            }
        }
        let (lo, hi) = pack(&counts, width);
        out.push((lo, hi, code));
    }
    out.sort_unstable();
    out
}

fn groups(sweep: &[(u128, u128, u32)]) -> Vec<Vec<u32>> {
    let mut out: Vec<Vec<u32>> = Vec::new();
    let mut k = 0usize;
    while k < sweep.len() {
        let mut j = k + 1;
        while j < sweep.len() && sweep[j].0 == sweep[k].0 && sweep[j].1 == sweep[k].1 {
            j += 1;
        }
        if j - k > 1 {
            out.push(sweep[k..j].iter().map(|e| e.2).collect());
        }
        k = j;
    }
    out
}

fn mix(x: u64) -> u64 {
    let mut z = x.wrapping_add(0x9E3779B97F4A7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

fn digest(points: &[usize], table: &Classes, scratch: &mut [u32], touched: &mut Vec<u32>) -> u128 {
    touched.clear();
    for a in 0..points.len() {
        for b in a..points.len() {
            let c = table.index[points[a] * table.cells + points[b]] as usize;
            if scratch[c] == 0 {
                touched.push(c as u32);
            }
            scratch[c] += 1;
        }
    }
    let (mut one, mut two) = (0u64, 0u64);
    for &c in touched.iter() {
        let n = scratch[c as usize];
        scratch[c as usize] = 0;
        let k = mix(((c as u64) << 34) ^ n as u64);
        one = one.wrapping_add(k);
        two ^= mix(k ^ 0x51_7C_C1_B7_27_22_0A_95);
    }
    (one as u128) << 64 | two as u128
}

fn exact(points: &[usize], table: &Classes) -> Vec<(u32, u32)> {
    let mut seen: Vec<(u32, u32)> = Vec::new();
    let mut tally = vec![0u32; table.count];
    for a in 0..points.len() {
        for b in a..points.len() {
            tally[table.index[points[a] * table.cells + points[b]] as usize] += 1;
        }
    }
    for (c, &n) in tally.iter().enumerate() {
        if n > 0 {
            seen.push((c as u32, n));
        }
    }
    seen
}

fn kron_cube(bits: u32, q: usize) -> Vec<usize> {
    let base: Vec<usize> = (0..q * q * q).filter(|j| bits >> j & 1 == 1).collect();
    let side = q * q;
    let mut out = Vec::with_capacity(base.len() * base.len());
    for &p in &base {
        for &s in &base {
            let a = (p / (q * q)) * q + s / (q * q);
            let b = ((p / q) % q) * q + (s / q) % q;
            let c = (p % q) * q + s % q;
            out.push((a * side + b) * side + c);
        }
    }
    out.sort_unstable();
    out
}

fn level_two_window(
    clashes: &[Vec<u32>],
    table: &Classes,
    render: impl Fn(u32) -> Vec<usize>,
    budget: u128,
) -> (u32, usize, usize, Vec<(u32, u32)>) {
    let mut strata: Vec<Vec<&Vec<u32>>> = vec![Vec::new(); 64];
    for group in clashes {
        let weight = group.iter().map(|c| c.count_ones()).max().unwrap_or(0);
        strata[weight as usize].push(group);
    }
    let (mut top, mut done, mut spent) = (0u32, 0usize, 0u128);
    let mut survivors: Vec<(u32, u32)> = Vec::new();
    let mut scratch = vec![0u32; table.count];
    let mut touched: Vec<u32> = Vec::new();
    for (weight, stratum) in strata.iter().enumerate() {
        if stratum.is_empty() {
            continue;
        }
        let bill: u128 = stratum
            .iter()
            .map(|group| {
                group
                    .iter()
                    .map(|c| {
                        let n = c.count_ones() as u128 * c.count_ones() as u128;
                        n * (n + 1) / 2
                    })
                    .sum::<u128>()
            })
            .sum();
        if spent + bill > budget {
            break;
        }
        spent += bill;
        top = weight as u32;
        done += stratum.len();
        for group in stratum {
            let seen: Vec<u128> = group
                .iter()
                .map(|&c| digest(&render(c), table, &mut scratch, &mut touched))
                .collect();
            for a in 0..group.len() {
                for b in a + 1..group.len() {
                    if seen[a] == seen[b]
                        && exact(&render(group[a]), table) == exact(&render(group[b]), table)
                    {
                        survivors.push((group[a], group[b]));
                    }
                }
            }
        }
    }
    (top, done, clashes.len(), survivors)
}

pub fn base_five() {
    println!();
    println!("SHAPE READING, BASE FIVE PLANE");
    let g1 = d4(5);
    let g2 = d4(25);
    let t1 = classes(25, &g1);
    let t2 = classes(625, &g2);
    let orbits = burnside(25, &g1);
    println!(
        "  pair classes level 1 {} level 2 {}; Burnside orbits {} of which {} nonempty",
        t1.count,
        t2.count,
        orbits,
        orbits - 1
    );
    let sweep = level_one_sweep(25, &g1, &t1, 4);
    println!(
        "  canonical codes swept over all 2^25: {}; distinct level-1 pair censuses: {}",
        sweep.len(),
        {
            let mut k = 0usize;
            for j in 0..sweep.len() {
                if j == 0 || (sweep[j].0, sweep[j].1) != (sweep[j - 1].0, sweep[j - 1].1) {
                    k += 1;
                }
            }
            k
        }
    );
    let clashes = groups(&sweep);
    let orbits_in: usize = clashes.iter().map(|g| g.len()).sum();
    println!(
        "  level-1 collisions: {} censuses shared by {} orbits, the heaviest group holding {}",
        clashes.len(),
        orbits_in,
        clashes.iter().map(|g| g.len()).max().unwrap_or(0)
    );
    let (top, done, all, alive) =
        level_two_window(&clashes, &t2, |c| kron_plane(c as u128, 5), 24_000_000_000);
    println!(
        "  level-2 window: every colliding census of weight at most {top}, {done} of {all}; pairs surviving the level-2 pair census: {}",
        alive.len()
    );
    for (a, b) in alive.iter().take(4) {
        println!("    surviving pair {a} against {b}");
    }
}

pub fn cube_three() {
    println!();
    println!("SHAPE READING, BASE THREE AT D = 3");
    let g1 = oct(3);
    let g2 = oct(9);
    let t1 = classes(27, &g1);
    let t2 = classes(729, &g2);
    let orbits = burnside(27, &g1);
    println!(
        "  cube group order {}; pair classes level 1 {} level 2 {}; Burnside orbits {} of which {} nonempty",
        g1.len(),
        t1.count,
        t2.count,
        orbits,
        orbits - 1
    );
    let sweep = level_one_sweep(27, &g1, &t1, 6);
    println!(
        "  canonical codes swept over all 2^27: {}; distinct level-1 pair censuses: {}",
        sweep.len(),
        {
            let mut k = 0usize;
            for j in 0..sweep.len() {
                if j == 0 || (sweep[j].0, sweep[j].1) != (sweep[j - 1].0, sweep[j - 1].1) {
                    k += 1;
                }
            }
            k
        }
    );
    let clashes = groups(&sweep);
    let orbits_in: usize = clashes.iter().map(|g| g.len()).sum();
    println!(
        "  level-1 collisions: {} censuses shared by {} orbits, the heaviest group holding {}",
        clashes.len(),
        orbits_in,
        clashes.iter().map(|g| g.len()).max().unwrap_or(0)
    );
    let (top, done, all, alive) =
        level_two_window(&clashes, &t2, |c| kron_cube(c, 3), 200_000_000_000);
    println!(
        "  level-2 window: every colliding census of weight at most {top}, {done} of {all}; pairs surviving the level-2 pair census: {}",
        alive.len()
    );
    for (a, b) in alive.iter().take(4) {
        println!("    surviving pair {a} against {b}");
    }
}
