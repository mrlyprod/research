use crate::tile::Tile;

pub fn residue_tile(code: u32, base: usize, side: usize) -> Tile {
    let mut tile = Tile::new(side);
    for r in 0..side {
        for c in 0..side {
            let index = (r % base) * base + (c % base);
            if (code >> index) & 1 == 1 {
                tile.set(r, c);
            }
        }
    }
    tile
}

pub fn collisions(side: usize) -> Vec<(u32, u32)> {
    let two: Vec<(u32, Tile)> = (1u32..16).map(|c| (c, residue_tile(c, 2, side))).collect();
    let three: Vec<(u32, Tile)> = (1u32..512).map(|c| (c, residue_tile(c, 3, side))).collect();
    let mut out = Vec::new();
    for (a, left) in &two {
        for (b, right) in &three {
            if left == right {
                out.push((*a, *b));
            }
        }
    }
    out
}

pub fn partners(code: u32, side: usize) -> Vec<u32> {
    let truth = residue_tile(code, 2, side);
    (1u32..512)
        .filter(|c| residue_tile(*c, 3, side) == truth)
        .collect()
}

pub fn self_power(code: u32, base: usize, level: usize) -> Tile {
    let seed = residue_tile(code, base, base);
    let mut out = seed.clone();
    for _ in 1..level {
        out = crate::tile::kron(&out, &seed);
    }
    out
}
