use std::collections::BTreeMap;

pub struct Dsu {
    parent: Vec<u32>,
}

impl Dsu {
    pub fn new(size: usize) -> Dsu {
        Dsu {
            parent: (0..size as u32).collect(),
        }
    }

    pub fn find(&mut self, mut i: u32) -> u32 {
        while self.parent[i as usize] != i {
            let up = self.parent[self.parent[i as usize] as usize];
            self.parent[i as usize] = up;
            i = up;
        }
        i
    }

    pub fn union(&mut self, a: u32, b: u32) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        self.parent[ra as usize] = rb;
        true
    }
}

pub fn arcs(side: usize, on: &[bool]) -> (u64, u64) {
    let rows = side * (side + 1);
    let h = |x: usize, y: usize| (y * side + x) as u32;
    let v = |x: usize, y: usize| (rows + y * (side + 1) + x) as u32;
    let mut dsu = Dsu::new(2 * rows);
    let mut degree = vec![0u8; 2 * rows];
    for y in 0..side {
        for x in 0..side {
            let (bottom, top, left, right) = (h(x, y), h(x, y + 1), v(x, y), v(x + 1, y));
            let pairs = if on[y * side + x] {
                [(left, bottom), (right, top)]
            } else {
                [(bottom, right), (left, top)]
            };
            for (a, b) in pairs {
                degree[a as usize] += 1;
                degree[b as usize] += 1;
                dsu.union(a, b);
            }
        }
    }
    let mut closed = vec![1u8; 2 * rows];
    let mut root = vec![false; 2 * rows];
    for (i, &deg) in degree.iter().enumerate() {
        let r = dsu.find(i as u32) as usize;
        root[r] = true;
        if deg != 2 {
            closed[r] = 0;
        }
    }
    let curves = root.iter().filter(|&&r| r).count() as u64;
    let loops = (0..2 * rows).filter(|&i| root[i] && closed[i] == 1).count() as u64;
    (loops, curves - loops)
}

pub fn mirrors(side: usize, on: &[bool]) -> i64 {
    let w = side + 1;
    let p = |x: usize, y: usize| (y * w + x) as u32;
    let mut dsu = Dsu::new(w * w);
    let mut comp = (w * w) as i64;
    for y in 0..side {
        for x in 0..side {
            let (a, b) = if on[y * side + x] {
                (p(x + 1, y), p(x, y + 1))
            } else {
                (p(x, y), p(x + 1, y + 1))
            };
            if dsu.union(a, b) {
                comp -= 1;
            }
        }
    }
    comp - 2 * side as i64 - 1
}

pub struct Block {
    pub side: usize,
    pub partner: Vec<u32>,
    pub loops: u128,
    pub lengths: BTreeMap<u64, u64>,
}

pub fn unit(on: bool) -> Block {
    let partner = if on {
        vec![3, 2, 1, 0]
    } else {
        vec![1, 0, 3, 2]
    };
    Block {
        side: 1,
        partner,
        loops: 0,
        lengths: BTreeMap::new(),
    }
}

pub fn void_partner(s: usize, p: usize) -> usize {
    let t = p % s;
    match p / s {
        0 => 2 * s - 1 - t,
        1 => s - 1 - t,
        2 => 4 * s - 1 - t,
        _ => 3 * s - 1 - t,
    }
}

pub fn void_block(s: usize) -> Block {
    Block {
        side: s,
        partner: (0..4 * s).map(|p| void_partner(s, p) as u32).collect(),
        loops: 0,
        lengths: BTreeMap::new(),
    }
}

pub fn glue(n: usize, tile: &[bool], kept: &Block, keep: bool) -> Block {
    let s = kept.side;
    let big = n * s;
    let s4 = 4 * s;
    let mut seen = vec![0u64; (n * n * s4).div_ceil(64)];
    let mut partner = if keep {
        vec![u32::MAX; 4 * big]
    } else {
        Vec::new()
    };
    let inside = |b: usize, p: usize| {
        if tile[b] {
            kept.partner[p] as usize
        } else {
            void_partner(s, p)
        }
    };
    let outer = |i: usize, j: usize, p: usize| -> Option<usize> {
        let t = p % s;
        match p / s {
            0 => (i == 0).then_some(j * s + t),
            1 => (j == n - 1).then_some(big + i * s + t),
            2 => (i == n - 1).then_some(2 * big + j * s + t),
            _ => (j == 0).then_some(3 * big + i * s + t),
        }
    };
    let across = |i: usize, j: usize, p: usize| -> (usize, usize, usize) {
        let t = p % s;
        match p / s {
            0 => (i - 1, j, 2 * s + t),
            1 => (i, j + 1, 3 * s + t),
            2 => (i + 1, j, t),
            _ => (i, j - 1, s + t),
        }
    };
    let key = |i: usize, j: usize, p: usize| (i * n + j) * s4 + p;
    let mark = |seen: &mut Vec<u64>, k: usize| seen[k / 64] |= 1 << (k % 64);
    let is_seen = |seen: &Vec<u64>, k: usize| seen[k / 64] >> (k % 64) & 1 == 1;
    for i in 0..n {
        for j in 0..n {
            for p in 0..s4 {
                let Some(start) = outer(i, j, p) else {
                    continue;
                };
                if is_seen(&seen, key(i, j, p)) {
                    continue;
                }
                let (mut a, mut b, mut c) = (i, j, p);
                loop {
                    mark(&mut seen, key(a, b, c));
                    let q = inside(a * n + b, c);
                    mark(&mut seen, key(a, b, q));
                    if let Some(end) = outer(a, b, q) {
                        if keep {
                            partner[start] = end as u32;
                            partner[end] = start as u32;
                        }
                        break;
                    }
                    (a, b, c) = across(a, b, q);
                }
            }
        }
    }
    let mut cycles = 0u128;
    let mut lengths = BTreeMap::new();
    for i in 0..n {
        for j in 0..n {
            for p in 0..s4 {
                if is_seen(&seen, key(i, j, p)) {
                    continue;
                }
                cycles += 1;
                let mut length = 0u64;
                let (mut a, mut b, mut c) = (i, j, p);
                loop {
                    length += 1;
                    mark(&mut seen, key(a, b, c));
                    let q = inside(a * n + b, c);
                    mark(&mut seen, key(a, b, q));
                    (a, b, c) = across(a, b, q);
                    if (a, b, c) == (i, j, p) {
                        break;
                    }
                }
                *lengths.entry(length).or_insert(0) += 1;
            }
        }
    }
    let inner = tile.iter().filter(|&&on| on).count() as u128 * kept.loops;
    Block {
        side: big,
        partner,
        loops: inner + cycles,
        lengths,
    }
}

pub fn series(n: usize, tile: &[bool], top: usize) -> Vec<u128> {
    let mut kept = unit(true);
    let mut out = vec![0];
    for level in 0..top {
        kept = glue(n, tile, &kept, level + 1 < top);
        out.push(kept.loops);
    }
    out
}

pub fn blocks(n: usize, tile: &[bool], top: usize) -> Vec<Block> {
    let mut out = vec![unit(true)];
    for level in 0..top {
        let next = glue(n, tile, &out[level], true);
        out.push(next);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::Design;

    #[test]
    fn three_counts_agree_at_base_two() {
        for code in 0..16u128 {
            let d = Design::full(code, 2);
            let blocks = series(2, &d.tile, 5);
            for (level, &want) in blocks.iter().enumerate() {
                let (side, on) = d.cells(level);
                let (loops, strands) = arcs(side, &on);
                assert_eq!(loops as u128, want);
                assert_eq!(mirrors(side, &on), loops as i64);
                assert_eq!(strands, 2 * side as u64);
            }
        }
    }

    #[test]
    fn carpet_law() {
        let s = series(3, &Design::new(7, 3, 2).tile, 8);
        for (n, &x) in s.iter().enumerate() {
            let n = n as u32;
            assert_eq!(
                x as i128,
                (8i128.pow(n) - 1) / 7 - 3i128.pow(n) + n as i128 + 1
            );
        }
    }

    #[test]
    fn void_block_glues_to_its_formula() {
        for n in 2..6 {
            let glued = glue(n, &vec![true; n * n], &void_block(n), true);
            assert_eq!(glued.partner, void_block(n * n).partner);
            assert_eq!(glued.loops, 0);
        }
    }
}
