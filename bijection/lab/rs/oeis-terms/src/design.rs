use num_bigint::BigUint;
use std::collections::HashSet;

const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

type Matrix = [[i64; 3]; 3];

const IDENTITY: Matrix = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];

#[derive(Clone, Copy)]
pub struct Element {
    perm: [usize; 3],
    sign: [i64; 3],
    shift: [i64; 3],
}

pub fn group(n: usize) -> Vec<Element> {
    let mut out = Vec::with_capacity(48 * n * n * n);
    for perm in PERMS {
        for bits in 0..8i64 {
            let sign = [
                1 - 2 * (bits >> 2 & 1),
                1 - 2 * (bits >> 1 & 1),
                1 - 2 * (bits & 1),
            ];
            for s0 in 0..n as i64 {
                for s1 in 0..n as i64 {
                    for s2 in 0..n as i64 {
                        out.push(Element {
                            perm,
                            sign,
                            shift: [s0, s1, s2],
                        });
                    }
                }
            }
        }
    }
    out
}

fn image(n: usize, g: &Element, x: [i64; 3]) -> usize {
    let m = n as i64;
    let mut index = 0usize;
    for t in 0..3 {
        let y = (g.sign[t] * x[g.perm[t]] + g.shift[t]).rem_euclid(m);
        index = index * n + y as usize;
    }
    index
}

fn images(n: usize, g: &Element, img: &mut [u32]) {
    let mut index = 0usize;
    for x0 in 0..n as i64 {
        for x1 in 0..n as i64 {
            for x2 in 0..n as i64 {
                img[index] = image(n, g, [x0, x1, x2]) as u32;
                index += 1;
            }
        }
    }
}

fn walk(img: &[u32], seen: &mut [bool]) -> usize {
    seen.fill(false);
    let mut count = 0;
    for start in 0..img.len() {
        if seen[start] {
            continue;
        }
        count += 1;
        let mut j = start;
        while !seen[j] {
            seen[j] = true;
            j = img[j] as usize;
        }
    }
    count
}

fn burnside(n: usize, histogram: &[u64]) -> BigUint {
    let mut total = BigUint::from(0u64);
    for (cycles, tally) in histogram.iter().enumerate() {
        if *tally != 0 {
            total += BigUint::from(*tally) << cycles;
        }
    }
    let order = BigUint::from(48 * (n as u64).pow(3));
    assert!((&total % &order) == BigUint::from(0u64));
    total / order
}

pub fn by_cycles(n: usize) -> BigUint {
    let cells = n * n * n;
    let mut histogram = vec![0u64; cells + 1];
    let mut img = vec![0u32; cells];
    let mut seen = vec![false; cells];
    for g in group(n) {
        images(n, &g, &mut img);
        histogram[walk(&img, &mut seen)] += 1;
    }
    burnside(n, &histogram)
}

fn matrix(g: &Element) -> Matrix {
    let mut rows: Matrix = [[0; 3]; 3];
    for t in 0..3 {
        rows[t][g.perm[t]] = g.sign[t];
    }
    rows
}

fn matmul(a: &Matrix, b: &Matrix, m: i64) -> Matrix {
    let mut out: Matrix = [[0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum::<i64>().rem_euclid(m);
        }
    }
    out
}

fn affine(a: &Matrix, v: [i64; 3], t: [i64; 3], m: i64) -> [i64; 3] {
    let mut out = [0i64; 3];
    for i in 0..3 {
        out[i] = ((0..3).map(|k| a[i][k] * v[k]).sum::<i64>() + t[i]).rem_euclid(m);
    }
    out
}

fn fixed_points(a: &Matrix, t: [i64; 3], m: i64) -> u64 {
    let mut count = 0u64;
    for x0 in 0..m {
        for x1 in 0..m {
            for x2 in 0..m {
                let x = [x0, x1, x2];
                if affine(a, x, t, m) == x {
                    count += 1;
                }
            }
        }
    }
    count
}

fn cycles_by_powers(g: &Element, m: i64) -> usize {
    let base = matrix(g);
    let identity = matmul(&IDENTITY, &IDENTITY, m);
    let mut power = matmul(&base, &IDENTITY, m);
    let mut offset = g.shift;
    let mut total = 0u64;
    let mut k = 0u64;
    loop {
        k += 1;
        total += fixed_points(&power, offset, m);
        if power == identity && offset == [0, 0, 0] {
            break;
        }
        power = matmul(&base, &power, m);
        offset = affine(&base, offset, g.shift, m);
    }
    assert!(total.is_multiple_of(k));
    (total / k) as usize
}

pub fn by_affine(n: usize) -> BigUint {
    let cells = n * n * n;
    let mut histogram = vec![0u64; cells + 1];
    for g in group(n) {
        histogram[cycles_by_powers(&g, n as i64)] += 1;
    }
    burnside(n, &histogram)
}

pub fn by_orbits(n: usize) -> u64 {
    let cells = n * n * n;
    let mut img = vec![0u32; cells];
    let mut maps: HashSet<Vec<u32>> = HashSet::new();
    for g in group(n) {
        images(n, &g, &mut img);
        maps.insert(img.clone());
    }
    let mut canon: HashSet<u64> = HashSet::new();
    for mask in 0u64..1u64 << cells {
        let mut best = mask;
        for map in &maps {
            let mut moved = 0u64;
            for (i, target) in map.iter().enumerate() {
                if mask >> i & 1 == 1 {
                    moved |= 1u64 << target;
                }
            }
            best = best.min(moved);
        }
        canon.insert(best);
    }
    canon.len() as u64
}
