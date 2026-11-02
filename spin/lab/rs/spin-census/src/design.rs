use mrly_core::Tensor;
use mrly_math::bang::Code;
use std::collections::BTreeSet;

pub const BASE: usize = 3;

pub fn plane(code: u128, base: usize, level: usize) -> Tensor {
    mrly_math::two::create(Code::from(code), base, level, 0, base)
        .expect("a plane design renders")
        .types()
        .clone()
}

pub fn floats(grid: &Tensor) -> Vec<f32> {
    (0..grid.size()).map(|i| grid.at(i) as f32).collect()
}

pub fn bit_cells() -> Vec<(usize, usize)> {
    (0..9)
        .map(|bit| {
            let grid = plane(1u128 << bit, BASE, 1);
            let flat = grid
                .bytes()
                .expect("a design is bytes")
                .iter()
                .position(|cell| *cell != 0)
                .expect("one filled cell");
            (flat / BASE, flat % BASE)
        })
        .collect()
}

pub fn square_group() -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = Vec::new();
    for flip in 0..2 {
        for turn in 0..4 {
            let map: Vec<usize> = (0..9)
                .map(|flat| {
                    let (mut r, mut c) = (flat / BASE, flat % BASE);
                    if flip == 1 {
                        std::mem::swap(&mut r, &mut c);
                    }
                    for _ in 0..turn {
                        let next = (c, BASE - 1 - r);
                        r = next.0;
                        c = next.1;
                    }
                    r * BASE + c
                })
                .collect();
            if !out.contains(&map) {
                out.push(map);
            }
        }
    }
    out
}

pub fn carry(map: &[usize], code: u128, table: &[(usize, usize)]) -> u128 {
    let index: Vec<usize> = table.iter().map(|(r, c)| r * BASE + c).collect();
    let mut out = 0u128;
    for bit in 0..9 {
        if code >> bit & 1 == 1 {
            let image = map[index[bit]];
            let target = index
                .iter()
                .position(|flat| *flat == image)
                .expect("a cell");
            out |= 1u128 << target;
        }
    }
    out
}

pub fn orbit(group: &[Vec<usize>], code: u128, table: &[(usize, usize)]) -> BTreeSet<u128> {
    group.iter().map(|map| carry(map, code, table)).collect()
}

pub fn canonical(group: &[Vec<usize>], code: u128, table: &[(usize, usize)]) -> u128 {
    *orbit(group, code, table).iter().next().expect("an orbit")
}
