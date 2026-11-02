use crate::tile::{kron, mask_tile, separable, split, unpack, Tile};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};

pub fn codes(base: usize) -> Vec<u32> {
    let total = 1u32 << (base * base);
    (1..total).collect()
}

pub struct SideFour {
    pub products: usize,
    pub distinct: usize,
    pub max_preimage: usize,
    pub by_product: HashSet<u64>,
    pub by_block: HashSet<u64>,
}

pub fn side_four() -> SideFour {
    let mut counts: HashMap<u64, usize> = HashMap::new();
    let mut products = 0usize;
    for a in codes(2) {
        let outer = mask_tile(a as u64, 2);
        for b in codes(2) {
            let inner = mask_tile(b as u64, 2);
            let whole = kron(&outer, &inner);
            *counts.entry(whole.pack()).or_insert(0) += 1;
            products += 1;
        }
    }
    let by_product: HashSet<u64> = counts.keys().copied().collect();
    let max_preimage = counts.values().copied().max().unwrap_or(0);
    let distinct = by_product.len();
    let mut by_block: HashSet<u64> = HashSet::new();
    for mask in 1u64..1 << 16 {
        let tile = mask_tile(mask, 4);
        if split(&tile, 2).is_some() {
            by_block.insert(mask);
        }
    }
    SideFour {
        products,
        distinct,
        max_preimage,
        by_product,
        by_block,
    }
}

pub struct SideSix {
    pub image23: HashMap<u64, (u32, u32)>,
    pub image32: HashMap<u64, (u32, u32)>,
    pub collisions23: usize,
    pub collisions32: usize,
    pub cross: Vec<u64>,
    pub reducible: Vec<u64>,
}

pub fn side_six() -> SideSix {
    let mut image23: HashMap<u64, (u32, u32)> = HashMap::new();
    let mut collisions23 = 0usize;
    for a in codes(2) {
        let outer = mask_tile(a as u64, 2);
        for b in codes(3) {
            let inner = mask_tile(b as u64, 3);
            let key = kron(&outer, &inner).pack();
            if image23.insert(key, (a, b)).is_some() {
                collisions23 += 1;
            }
        }
    }
    let mut image32: HashMap<u64, (u32, u32)> = HashMap::new();
    let mut collisions32 = 0usize;
    for a in codes(3) {
        let outer = mask_tile(a as u64, 3);
        for b in codes(2) {
            let inner = mask_tile(b as u64, 2);
            let key = kron(&outer, &inner).pack();
            if image32.insert(key, (a, b)).is_some() {
                collisions32 += 1;
            }
        }
    }
    let mut cross: Vec<u64> = image23
        .keys()
        .copied()
        .filter(|key| image32.contains_key(key))
        .collect();
    cross.sort_unstable();
    let mut union: BTreeSet<u64> = image23.keys().copied().collect();
    union.extend(image32.keys().copied());
    let reducible: Vec<u64> = union.into_iter().collect();
    SideSix {
        image23,
        image32,
        collisions23,
        collisions32,
        cross,
        reducible,
    }
}

pub struct CrossAnatomy {
    pub separable: usize,
    pub non_separable: usize,
    pub commuting: usize,
    pub rewriting: usize,
    pub non_separable_non_commuting: usize,
    pub non_separable_commuting: Vec<(u32, u32)>,
    pub fills: BTreeMap<usize, usize>,
    pub fills_non_separable: BTreeMap<usize, usize>,
    pub outer_fills: BTreeMap<(usize, usize), usize>,
    pub one_cell_any: BTreeMap<usize, usize>,
    pub one_cell_any_non_separable: BTreeMap<usize, usize>,
    pub one_cell_outer: BTreeMap<usize, usize>,
    pub one_cell_outer_non_separable: BTreeMap<usize, usize>,
    pub commuting_pairs: Vec<(u32, u32)>,
    pub separable_set: HashSet<u64>,
    pub set_non_commuting: HashSet<u64>,
    pub set_one_any: HashSet<u64>,
    pub set_one_outer: HashSet<u64>,
}

pub fn cross_anatomy(six: &SideSix) -> CrossAnatomy {
    let mut out = CrossAnatomy {
        separable: 0,
        non_separable: 0,
        commuting: 0,
        rewriting: 0,
        non_separable_non_commuting: 0,
        non_separable_commuting: Vec::new(),
        fills: BTreeMap::new(),
        fills_non_separable: BTreeMap::new(),
        outer_fills: BTreeMap::new(),
        one_cell_any: BTreeMap::new(),
        one_cell_any_non_separable: BTreeMap::new(),
        one_cell_outer: BTreeMap::new(),
        one_cell_outer_non_separable: BTreeMap::new(),
        commuting_pairs: Vec::new(),
        separable_set: HashSet::new(),
        set_non_commuting: HashSet::new(),
        set_one_any: HashSet::new(),
        set_one_outer: HashSet::new(),
    };
    for key in &six.cross {
        let tile = unpack(*key, 6);
        let (a2, b3) = six.image23[key];
        let (x3, y2) = six.image32[key];
        let flat = separable(&tile);
        let commutes = a2 == y2 && b3 == x3;
        let fill_a2 = mask_tile(a2 as u64, 2).fill();
        let fill_b3 = mask_tile(b3 as u64, 3).fill();
        let fill_x3 = mask_tile(x3 as u64, 3).fill();
        let fill_y2 = mask_tile(y2 as u64, 2).fill();
        if flat {
            out.separable += 1;
            out.separable_set.insert(*key);
        } else {
            out.non_separable += 1;
        }
        if commutes {
            out.commuting += 1;
            out.commuting_pairs.push((a2, b3));
        } else {
            out.rewriting += 1;
        }
        if !flat && !commutes {
            out.non_separable_non_commuting += 1;
        }
        if !flat && commutes {
            out.non_separable_commuting.push((a2, b3));
        }
        *out.fills.entry(tile.fill()).or_insert(0) += 1;
        if !flat {
            *out.fills_non_separable.entry(tile.fill()).or_insert(0) += 1;
            let low = fill_a2.min(fill_x3);
            let high = fill_a2.max(fill_x3);
            *out.outer_fills.entry((low, high)).or_insert(0) += 1;
        }
        let any =
            usize::from(fill_a2 == 1 || fill_b3 == 1) + usize::from(fill_x3 == 1 || fill_y2 == 1);
        let outer = usize::from(fill_a2 == 1) + usize::from(fill_x3 == 1);
        *out.one_cell_any.entry(any).or_insert(0) += 1;
        *out.one_cell_outer.entry(outer).or_insert(0) += 1;
        if !flat {
            *out.one_cell_any_non_separable.entry(any).or_insert(0) += 1;
            *out.one_cell_outer_non_separable.entry(outer).or_insert(0) += 1;
            if !commutes {
                out.set_non_commuting.insert(*key);
            }
            if any > 0 {
                out.set_one_any.insert(*key);
            }
            if outer > 0 {
                out.set_one_outer.insert(*key);
            }
        }
    }
    out.commuting_pairs.sort_unstable();
    out.non_separable_commuting.sort_unstable();
    out
}

pub fn two_radix_lines() -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = Vec::new();
    for mask in 1u128..1 << 6 {
        let up = crate::tile::line_split(mask, 6, 2).is_some();
        let down = crate::tile::line_split(mask, 6, 3).is_some();
        if up && down {
            out.push((0..6).filter(|i| (mask >> i) & 1 == 1).collect());
        }
    }
    out
}

pub fn rectangle_set() -> HashSet<u64> {
    let lines: Vec<u128> = (1u128..1 << 6)
        .filter(|mask| {
            crate::tile::line_split(*mask, 6, 2).is_some()
                && crate::tile::line_split(*mask, 6, 3).is_some()
        })
        .collect();
    let mut out = HashSet::new();
    for rows in &lines {
        for cols in &lines {
            let mut tile = Tile::new(6);
            for r in 0..6 {
                if (rows >> r) & 1 == 0 {
                    continue;
                }
                for c in 0..6 {
                    if (cols >> c) & 1 == 1 {
                        tile.set(r, c);
                    }
                }
            }
            out.insert(tile.pack());
        }
    }
    out
}

pub struct SideEight {
    pub image24: usize,
    pub image42: usize,
    pub intersection: usize,
    pub triples_match: bool,
}

pub fn side_eight() -> SideEight {
    let small: Vec<Tile> = codes(2)
        .into_iter()
        .map(|c| mask_tile(c as u64, 2))
        .collect();
    let mut image24: HashSet<u64> = HashSet::new();
    let mut image42: HashSet<u64> = HashSet::new();
    for mask in 1u64..1 << 16 {
        let big = mask_tile(mask, 4);
        for tile in &small {
            image24.insert(kron(tile, &big).pack());
            image42.insert(kron(&big, tile).pack());
        }
    }
    let intersection: HashSet<u64> = image24.intersection(&image42).copied().collect();
    let mut triples: HashSet<u64> = HashSet::new();
    for x in &small {
        for y in &small {
            for z in &small {
                triples.insert(kron(&kron(x, y), z).pack());
            }
        }
    }
    SideEight {
        image24: image24.len(),
        image42: image42.len(),
        intersection: intersection.len(),
        triples_match: triples == intersection,
    }
}

pub struct SideNine {
    pub products: usize,
    pub distinct: usize,
    pub collisions: usize,
}

pub fn side_nine() -> SideNine {
    let mut seen: HashSet<u128> = HashSet::new();
    let mut products = 0usize;
    let mut collisions = 0usize;
    let tiles: Vec<Tile> = codes(3)
        .into_iter()
        .map(|c| mask_tile(c as u64, 3))
        .collect();
    for outer in &tiles {
        for inner in &tiles {
            let whole = kron(outer, inner);
            let mut key = 0u128;
            for (index, cell) in whole.cells.iter().enumerate() {
                if *cell {
                    key |= 1u128 << index;
                }
            }
            products += 1;
            if !seen.insert(key) {
                collisions += 1;
            }
        }
    }
    SideNine {
        products,
        distinct: seen.len(),
        collisions,
    }
}

pub struct WordTwelve {
    pub per_shape: usize,
    pub distinct: [usize; 3],
    pub pairs: [usize; 3],
    pub triple: usize,
    pub union: usize,
    pub witness_shapes: [bool; 3],
}

pub fn word_twelve(witness: &Tile) -> WordTwelve {
    let two: Vec<Tile> = codes(2)
        .into_iter()
        .map(|c| mask_tile(c as u64, 2))
        .collect();
    let three: Vec<Tile> = codes(3)
        .into_iter()
        .map(|c| mask_tile(c as u64, 3))
        .collect();
    let mut sets: Vec<HashSet<[u64; 3]>> = Vec::new();
    let shapes: [[usize; 3]; 3] = [[2, 2, 3], [2, 3, 2], [3, 2, 2]];
    let mut per_shape = 0usize;
    for shape in shapes {
        let mut set: HashSet<[u64; 3]> = HashSet::new();
        let mut count = 0usize;
        let pick = |side: usize| -> &Vec<Tile> {
            if side == 2 {
                &two
            } else {
                &three
            }
        };
        for a in pick(shape[0]) {
            for b in pick(shape[1]) {
                let left = kron(a, b);
                for c in pick(shape[2]) {
                    set.insert(kron(&left, c).key());
                    count += 1;
                }
            }
        }
        per_shape = count;
        sets.push(set);
    }
    let distinct = [sets[0].len(), sets[1].len(), sets[2].len()];
    let pairs = [
        sets[0].intersection(&sets[1]).count(),
        sets[1].intersection(&sets[2]).count(),
        sets[0].intersection(&sets[2]).count(),
    ];
    let triple = sets[0]
        .iter()
        .filter(|key| sets[1].contains(*key) && sets[2].contains(*key))
        .count();
    let union = distinct[0] + distinct[1] + distinct[2] - pairs[0] - pairs[1] - pairs[2] + triple;
    let key = witness.key();
    let witness_shapes = [
        sets[0].contains(&key),
        sets[1].contains(&key),
        sets[2].contains(&key),
    ];
    WordTwelve {
        per_shape,
        distinct,
        pairs,
        triple,
        union,
        witness_shapes,
    }
}
