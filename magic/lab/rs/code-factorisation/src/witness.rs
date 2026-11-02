use crate::tile::{
    cuts, factorisations, gcd, incomparable, kron, lcm, line_cuts, line_kron, mask_tile, profile,
    unpack, Tile,
};
use std::collections::BTreeSet;

pub fn from_cells(side: usize, cells: &[(usize, usize)]) -> Tile {
    let mut tile = Tile::new(side);
    for (r, c) in cells {
        tile.set(*r, *c);
    }
    tile
}

pub fn identity(side: usize) -> Tile {
    let cells: Vec<(usize, usize)> = (0..side).map(|i| (i, i)).collect();
    from_cells(side, &cells)
}

pub fn antidiagonal(side: usize) -> Tile {
    let cells: Vec<(usize, usize)> = (0..side).map(|i| (i, side - 1 - i)).collect();
    from_cells(side, &cells)
}

pub fn gcd_closed(list: &[usize]) -> bool {
    let set: BTreeSet<usize> = list.iter().copied().collect();
    for a in &set {
        for b in &set {
            if !set.contains(&gcd(*a, *b)) {
                return false;
            }
        }
    }
    true
}

pub fn lcm_closed(list: &[usize], side: usize) -> bool {
    let set: BTreeSet<usize> = list.iter().copied().collect();
    for a in &set {
        for b in &set {
            let value = lcm(*a, *b);
            if value <= side && !set.contains(&value) {
                return false;
            }
        }
    }
    true
}

pub struct LineClosure {
    pub gcd_violations: usize,
    pub lcm_violations: usize,
    pub first_lcm: Option<(usize, u128, Vec<usize>)>,
}

pub fn line_closure(max_side: usize) -> LineClosure {
    let mut out = LineClosure {
        gcd_violations: 0,
        lcm_violations: 0,
        first_lcm: None,
    };
    for side in 1..=max_side {
        for mask in 1u128..1u128 << side {
            let list = line_cuts(mask, side);
            if !gcd_closed(&list) {
                out.gcd_violations += 1;
            }
            if !lcm_closed(&list, side) {
                out.lcm_violations += 1;
                if out.first_lcm.is_none() {
                    out.first_lcm = Some((side, mask, list.clone()));
                }
            }
        }
    }
    out
}

pub fn line_commuting(m: usize, n: usize) -> Vec<(u128, u128)> {
    let mut out = Vec::new();
    for a in 1u128..1u128 << m {
        for b in 1u128..1u128 << n {
            if line_kron(a, m, b, n) == line_kron(b, n, a, m) {
                out.push((a, b));
            }
        }
    }
    out
}

pub struct TwelveSweep {
    pub tiles: usize,
    pub gcd_violations: usize,
    pub lcm_violations: usize,
    pub mismatches: usize,
    pub multiple: usize,
    pub unequal_length: usize,
    pub max_factorisations: usize,
}

pub fn twelve_sweep(six: &[u64]) -> TwelveSweep {
    let small: Vec<Tile> = (1u32..16).map(|c| mask_tile(c as u64, 2)).collect();
    let mut out = TwelveSweep {
        tiles: 0,
        gcd_violations: 0,
        lcm_violations: 0,
        mismatches: 0,
        multiple: 0,
        unequal_length: 0,
        max_factorisations: 0,
    };
    let mut seen: std::collections::HashSet<[u64; 3]> = std::collections::HashSet::new();
    for key in six {
        let inner = unpack(*key, 6);
        for code in &small {
            for whole in [kron(code, &inner), kron(&inner, code)] {
                if !seen.insert(whole.key()) {
                    continue;
                }
                out.tiles += 1;
                let list = cuts(&whole);
                if !gcd_closed(&list) {
                    out.gcd_violations += 1;
                }
                if !lcm_closed(&list, 12) {
                    out.lcm_violations += 1;
                }
                let words = factorisations(&whole);
                out.max_factorisations = out.max_factorisations.max(words.len());
                let many = words.len() > 1;
                if many {
                    out.multiple += 1;
                }
                if many != incomparable(&list) {
                    out.mismatches += 1;
                }
                let lengths: BTreeSet<usize> = words.iter().map(|word| word.len()).collect();
                if lengths.len() > 1 {
                    out.unequal_length += 1;
                }
            }
        }
    }
    out
}

pub fn profiles_text(tile: &Tile) -> String {
    let mut rows: Vec<String> = factorisations(tile)
        .iter()
        .map(|word| {
            let sides: Vec<String> = profile(word).iter().map(|s| format!("{s}")).collect();
            format!("({})", sides.join(" x "))
        })
        .collect();
    rows.sort();
    rows.join(" ")
}
