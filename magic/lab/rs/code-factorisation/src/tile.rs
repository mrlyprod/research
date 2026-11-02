use std::collections::BTreeSet;

#[derive(Clone, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub struct Tile {
    pub side: usize,
    pub cells: Vec<bool>,
}

impl Tile {
    pub fn new(side: usize) -> Tile {
        Tile {
            side,
            cells: vec![false; side * side],
        }
    }

    pub fn at(&self, r: usize, c: usize) -> bool {
        self.cells[r * self.side + c]
    }

    pub fn set(&mut self, r: usize, c: usize) {
        let side = self.side;
        self.cells[r * side + c] = true;
    }

    pub fn fill(&self) -> usize {
        self.cells.iter().filter(|cell| **cell).count()
    }

    pub fn empty(&self) -> bool {
        !self.cells.iter().any(|cell| *cell)
    }

    pub fn pack(&self) -> u64 {
        let mut key = 0u64;
        for (index, cell) in self.cells.iter().enumerate() {
            if *cell {
                key |= 1u64 << index;
            }
        }
        key
    }

    pub fn key(&self) -> [u64; 3] {
        let mut key = [0u64; 3];
        for (index, cell) in self.cells.iter().enumerate() {
            if *cell {
                key[index / 64] |= 1u64 << (index % 64);
            }
        }
        key
    }

    pub fn support(&self) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        for r in 0..self.side {
            for c in 0..self.side {
                if self.at(r, c) {
                    out.push((r, c));
                }
            }
        }
        out
    }

    pub fn text(&self) -> String {
        let cells: Vec<String> = self
            .support()
            .iter()
            .map(|(r, c)| format!("({r},{c})"))
            .collect();
        format!("[{}]{{{}}}", self.side, cells.join(","))
    }
}

pub fn mask_tile(mask: u64, side: usize) -> Tile {
    let mut tile = Tile::new(side);
    for index in 0..side * side {
        if (mask >> index) & 1 == 1 {
            tile.set(index / side, index % side);
        }
    }
    tile
}

pub fn unpack(key: u64, side: usize) -> Tile {
    mask_tile(key, side)
}

pub fn kron(outer: &Tile, inner: &Tile) -> Tile {
    let side = outer.side * inner.side;
    let mut out = Tile::new(side);
    let n = inner.side;
    for i in 0..outer.side {
        for j in 0..outer.side {
            if !outer.at(i, j) {
                continue;
            }
            for p in 0..n {
                for q in 0..n {
                    if inner.at(p, q) {
                        out.set(i * n + p, j * n + q);
                    }
                }
            }
        }
    }
    out
}

pub fn divisors(n: usize) -> Vec<usize> {
    (1..=n).filter(|d| n.is_multiple_of(*d)).collect()
}

pub fn split(tile: &Tile, d: usize) -> Option<(Tile, Tile)> {
    if !tile.side.is_multiple_of(d) {
        return None;
    }
    let n = tile.side / d;
    let mut outer = Tile::new(d);
    let mut inner: Option<Tile> = None;
    for i in 0..d {
        for j in 0..d {
            let mut block = Tile::new(n);
            let mut live = false;
            for p in 0..n {
                for q in 0..n {
                    if tile.at(i * n + p, j * n + q) {
                        block.set(p, q);
                        live = true;
                    }
                }
            }
            if !live {
                continue;
            }
            outer.set(i, j);
            match &inner {
                None => inner = Some(block),
                Some(first) => {
                    if *first != block {
                        return None;
                    }
                }
            }
        }
    }
    inner.map(|block| (outer, block))
}

pub fn cuts(tile: &Tile) -> Vec<usize> {
    divisors(tile.side)
        .into_iter()
        .filter(|d| split(tile, *d).is_some())
        .collect()
}

pub fn irreducible(tile: &Tile) -> bool {
    if tile.side < 2 || tile.empty() {
        return false;
    }
    !divisors(tile.side)
        .into_iter()
        .any(|d| d > 1 && d < tile.side && split(tile, d).is_some())
}

pub fn separable(tile: &Tile) -> bool {
    let mut rows = BTreeSet::new();
    let mut cols = BTreeSet::new();
    for (r, c) in tile.support() {
        rows.insert(r);
        cols.insert(c);
    }
    rows.len() * cols.len() == tile.fill()
}

pub fn factorisations(tile: &Tile) -> Vec<Vec<Tile>> {
    let mut out: Vec<Vec<Tile>> = Vec::new();
    let mut found = false;
    for d in divisors(tile.side) {
        if d == 1 || d == tile.side {
            continue;
        }
        if let Some((outer, inner)) = split(tile, d) {
            found = true;
            for left in factorisations(&outer) {
                for right in factorisations(&inner) {
                    let mut whole = left.clone();
                    whole.extend(right.iter().cloned());
                    if !out.contains(&whole) {
                        out.push(whole);
                    }
                }
            }
        }
    }
    if !found {
        out.push(vec![tile.clone()]);
    }
    out
}

pub fn profile(word: &[Tile]) -> Vec<usize> {
    word.iter().map(|tile| tile.side).collect()
}

pub fn chain(word: &[Tile]) -> Vec<usize> {
    let mut out = vec![1usize];
    let mut run = 1usize;
    for tile in word {
        run *= tile.side;
        out.push(run);
    }
    out
}

pub fn totally_ordered(set: &BTreeSet<usize>) -> bool {
    let list: Vec<usize> = set.iter().copied().collect();
    for i in 0..list.len() {
        for j in i + 1..list.len() {
            if !list[j].is_multiple_of(list[i]) {
                return false;
            }
        }
    }
    true
}

pub fn incomparable(list: &[usize]) -> bool {
    for i in 0..list.len() {
        for j in i + 1..list.len() {
            if !list[j].is_multiple_of(list[i]) && !list[i].is_multiple_of(list[j]) {
                return true;
            }
        }
    }
    false
}

pub fn gcd(a: usize, b: usize) -> usize {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

pub fn lcm(a: usize, b: usize) -> usize {
    a / gcd(a, b) * b
}

pub fn line_kron(outer: u128, outer_side: usize, inner: u128, inner_side: usize) -> u128 {
    let mut out = 0u128;
    for i in 0..outer_side {
        if (outer >> i) & 1 == 0 {
            continue;
        }
        for p in 0..inner_side {
            if (inner >> p) & 1 == 1 {
                out |= 1u128 << (i * inner_side + p);
            }
        }
    }
    out
}

pub fn line_split(mask: u128, side: usize, d: usize) -> Option<(u128, u128)> {
    if !side.is_multiple_of(d) {
        return None;
    }
    let n = side / d;
    let full = if n >= 128 {
        u128::MAX
    } else {
        (1u128 << n) - 1
    };
    let mut outer = 0u128;
    let mut inner: Option<u128> = None;
    for i in 0..d {
        let block = (mask >> (i * n)) & full;
        if block == 0 {
            continue;
        }
        outer |= 1u128 << i;
        match inner {
            None => inner = Some(block),
            Some(first) => {
                if first != block {
                    return None;
                }
            }
        }
    }
    inner.map(|block| (outer, block))
}

pub fn line_cuts(mask: u128, side: usize) -> Vec<usize> {
    divisors(side)
        .into_iter()
        .filter(|d| line_split(mask, side, *d).is_some())
        .collect()
}

pub fn line_reducible(mask: u128, side: usize) -> bool {
    divisors(side)
        .into_iter()
        .any(|d| d > 1 && d < side && line_split(mask, side, d).is_some())
}

pub fn line_text(mask: u128, side: usize) -> String {
    let cells: Vec<String> = (0..side)
        .filter(|i| (mask >> i) & 1 == 1)
        .map(|i| format!("{i}"))
        .collect();
    format!("{{{}}}", cells.join(","))
}
