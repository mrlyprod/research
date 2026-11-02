use crate::word::{observe, render, Obs, CODES};
use std::collections::HashMap;

pub struct Row {
    pub total: usize,
    pub fill: usize,
    pub diagonal: usize,
    pub boundary: usize,
    pub perimeter: usize,
    pub components: usize,
    pub euler: usize,
    pub holes: usize,
    pub profile: usize,
    pub peak: usize,
    pub support: usize,
}

pub fn words(codes: &[u8], length: usize) -> Vec<Vec<u8>> {
    let mut out: Vec<Vec<u8>> = vec![Vec::new()];
    for _ in 0..length {
        let mut next: Vec<Vec<u8>> = Vec::new();
        for word in out.iter() {
            for &code in codes.iter() {
                let mut child = word.clone();
                child.push(code);
                next.push(child);
            }
        }
        out = next;
    }
    out
}

fn peak(obs: &Obs) -> u64 {
    *obs.profile.iter().max().expect("a profile has a row")
}

fn support(obs: &Obs) -> usize {
    obs.profile.iter().filter(|value| **value > 0).count()
}

pub fn sensitivity(codes: &[u8], length: usize) -> Row {
    let mut groups: HashMap<Vec<u8>, Vec<Vec<u8>>> = HashMap::new();
    for word in words(codes, length) {
        let mut key = word.clone();
        key.sort_unstable();
        groups.entry(key).or_default().push(word);
    }
    let mut row = Row {
        total: 0,
        fill: 0,
        diagonal: 0,
        boundary: 0,
        perimeter: 0,
        components: 0,
        euler: 0,
        holes: 0,
        profile: 0,
        peak: 0,
        support: 0,
    };
    for (_, family) in groups.iter() {
        if family.len() < 2 {
            continue;
        }
        row.total += 1;
        let seen: Vec<Obs> = family.iter().map(|word| observe(word)).collect();
        let head = &seen[0];
        let split = |test: &dyn Fn(&Obs, &Obs) -> bool| seen.iter().any(|obs| test(head, obs));
        if split(&|a, b| a.fill != b.fill) {
            row.fill += 1;
        }
        if split(&|a, b| a.diagonal != b.diagonal) {
            row.diagonal += 1;
        }
        if split(&|a, b| a.boundary != b.boundary) {
            row.boundary += 1;
        }
        if split(&|a, b| a.perimeter != b.perimeter) {
            row.perimeter += 1;
        }
        if split(&|a, b| a.components != b.components) {
            row.components += 1;
        }
        if split(&|a, b| a.euler != b.euler) {
            row.euler += 1;
        }
        if split(&|a, b| a.holes != b.holes) {
            row.holes += 1;
        }
        if split(&|a, b| a.profile != b.profile) {
            row.profile += 1;
        }
        if split(&|a, b| peak(a) != peak(b)) {
            row.peak += 1;
        }
        if split(&|a, b| support(a) != support(b)) {
            row.support += 1;
        }
    }
    row
}

fn index(word: &[u8]) -> usize {
    word.iter()
        .fold(0usize, |at, code| at * 15 + (*code as usize - 1))
}

pub fn library_search() -> (usize, Vec<Vec<u8>>) {
    let mut table = vec![[0i64; 4]; 15 * 15 * 15];
    for word in words(&CODES, 3) {
        let obs = observe(&word);
        table[index(&word)] = [
            obs.boundary as i64,
            obs.components as i64,
            obs.euler,
            obs.holes as i64,
        ];
    }
    let target = [36usize, 188, 188, 100];
    let mut hits: Vec<Vec<u8>> = Vec::new();
    let mut scanned = 0usize;
    let mut pick: Vec<u8> = Vec::new();
    choose(&CODES, 10, 0, &mut pick, &mut |subset: &[u8]| {
        scanned += 1;
        let mut counts = [0usize; 4];
        let mut total = 0usize;
        for a in 0..10 {
            for b in a..10 {
                for c in b..10 {
                    let letters = [subset[a], subset[b], subset[c]];
                    let mut orders: Vec<[u8; 3]> = Vec::new();
                    for order in [
                        [0usize, 1, 2],
                        [0, 2, 1],
                        [1, 0, 2],
                        [1, 2, 0],
                        [2, 0, 1],
                        [2, 1, 0],
                    ] {
                        let word = [letters[order[0]], letters[order[1]], letters[order[2]]];
                        if !orders.contains(&word) {
                            orders.push(word);
                        }
                    }
                    if orders.len() < 2 {
                        continue;
                    }
                    total += 1;
                    for slot in 0..4 {
                        let head = table[index(&orders[0])][slot];
                        if orders.iter().any(|word| table[index(word)][slot] != head) {
                            counts[slot] += 1;
                        }
                    }
                }
            }
        }
        if total == 210 && counts == target {
            hits.push(subset.to_vec());
        }
    });
    (scanned, hits)
}

fn choose(pool: &[u8], size: usize, at: usize, pick: &mut Vec<u8>, visit: &mut dyn FnMut(&[u8])) {
    if pick.len() == size {
        visit(pick);
        return;
    }
    if at >= pool.len() || pool.len() - at < size - pick.len() {
        return;
    }
    pick.push(pool[at]);
    choose(pool, size, at + 1, pick, visit);
    pick.pop();
    choose(pool, size, at + 1, pick, visit);
}

pub fn diagonal_factors() -> usize {
    let mut bad = 0usize;
    for word in words(&CODES, 3) {
        let product: u64 = word
            .iter()
            .map(|code| render(&[*code]).diagonal())
            .product();
        if render(&word).diagonal() != product {
            bad += 1;
        }
    }
    bad
}

pub fn contacts_multiply() -> usize {
    let mut bad = 0usize;
    for word in words(&CODES, 3) {
        let mut rows = 1u64;
        let mut cols = 1u64;
        for code in word.iter() {
            let (h, v) = render(&[*code]).contacts();
            rows *= h;
            cols *= v;
        }
        let (h, v) = render(&word).contacts();
        if h != rows || v != cols {
            bad += 1;
        }
    }
    bad
}

pub fn boundary_pairs() -> (usize, usize) {
    let mut bad = 0usize;
    let mut pairs = 0usize;
    for a in 0..16u8 {
        for b in 0..16u8 {
            pairs += 1;
            let left = pair_grid(a, b);
            let right = pair_grid(b, a);
            if left != right {
                bad += 1;
            }
        }
    }
    (pairs, bad)
}

fn pair_grid(a: u8, b: u8) -> (u64, u64) {
    let side = 4usize;
    let mut cells = vec![false; side * side];
    for (ra, ca) in crate::word::corners(a) {
        for (rb, cb) in crate::word::corners(b) {
            cells[(2 * ra + rb) * side + 2 * ca + cb] = true;
        }
    }
    let grid = crate::word::Grid { side, cells };
    (grid.boundary(), grid.interior())
}

pub fn pair_component_table(pool: &[u8]) -> Vec<(u8, u8, u64, u64)> {
    let mut out = Vec::new();
    for (i, &a) in pool.iter().enumerate() {
        for &b in pool.iter().skip(i + 1) {
            let left = render(&[a, b]).components();
            let right = render(&[b, a]).components();
            out.push((a, b, left, right));
        }
    }
    out
}

pub const PERIODS: [(&[u8], usize); 6] = [
    (&[3, 6], 2),
    (&[3, 6], 3),
    (&[7, 9], 3),
    (&[6, 9], 3),
    (&[3, 5, 6], 2),
    (&[7, 11, 13], 2),
];

pub fn block_reduction() -> (usize, usize) {
    use mrly_math::bang::{magic, MagicLayer};
    use mrly_math::name::Bang;
    let mut cases = 0usize;
    let mut bad = 0usize;
    for (period, repeats) in PERIODS.iter() {
        let layers: Vec<MagicLayer> = period
            .iter()
            .map(|code| MagicLayer::new(Bang::new(*code as u128, 2, 2), 2))
            .collect();
        let composite = magic(&layers).expect("the factory accepts a plane period");
        let power = composite.fractal(*repeats);
        let mut flat: Vec<u8> = Vec::new();
        for _ in 0..*repeats {
            flat.extend_from_slice(period);
        }
        let grid = render(&flat);
        cases += 1;
        let same = power.shape == vec![grid.side, grid.side]
            && power
                .bytes()
                .expect("the factory renders bytes")
                .iter()
                .zip(grid.cells.iter())
                .all(|(byte, cell)| (*byte == 1) == *cell);
        if !same {
            bad += 1;
        }
    }
    (cases, bad)
}

pub fn crate_agreement(length: usize) -> (usize, usize) {
    use mrly_math::bang::{magic, MagicLayer};
    use mrly_math::name::Bang;
    let mut checked = 0usize;
    let mut bad = 0usize;
    for word in words(&CODES, length) {
        let layers: Vec<MagicLayer> = word
            .iter()
            .map(|code| MagicLayer::new(Bang::new(*code as u128, 2, 2), 2))
            .collect();
        let tensor = magic(&layers).expect("the factory accepts a plane word");
        let grid = render(&word);
        checked += 1;
        let same = tensor.shape == vec![grid.side, grid.side]
            && tensor
                .bytes()
                .expect("the factory renders bytes")
                .iter()
                .zip(grid.cells.iter())
                .all(|(byte, cell)| (*byte == 1) == *cell);
        if !same {
            bad += 1;
        }
    }
    (checked, bad)
}
