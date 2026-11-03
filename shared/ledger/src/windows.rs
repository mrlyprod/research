use std::collections::{BTreeSet, HashSet};

#[derive(Clone)]
struct Picture {
    side: usize,
    cells: Vec<u8>,
}

impl Picture {
    fn at(&self, row: usize, col: usize) -> u8 {
        self.cells[row * self.side + col]
    }
    fn pairs(&self) -> BTreeSet<[u8; 4]> {
        let mut out = BTreeSet::new();
        for row in 0..self.side.saturating_sub(1) {
            for col in 0..self.side - 1 {
                out.insert([
                    self.at(row, col),
                    self.at(row, col + 1),
                    self.at(row + 1, col),
                    self.at(row + 1, col + 1),
                ]);
            }
        }
        out
    }
    fn substitute(&self, tile: &Picture) -> Picture {
        let s = tile.side;
        let side = self.side * s;
        let mut cells = vec![0u8; side * side];
        for row in 0..side {
            for col in 0..side {
                cells[row * side + col] = self.at(row / s, col / s) & tile.at(row % s, col % s);
            }
        }
        Picture { side, cells }
    }
    fn blocks(&self, k: usize, seen: &mut HashSet<Vec<u64>>) {
        let wide = self.side + 1 - k;
        for row in 0..wide {
            for col in 0..wide {
                let mut key = vec![0u64; (k * k).div_ceil(64)];
                for t in 0..k * k {
                    if self.at(row + t / k, col + t % k) != 0 {
                        key[t / 64] |= 1 << (t % 64);
                    }
                }
                seen.insert(key);
            }
        }
    }
}

fn language(tile: &Picture) -> Vec<Picture> {
    let mut union = BTreeSet::new();
    let mut walked: Vec<BTreeSet<[u8; 4]>> = Vec::new();
    let mut level = tile.pairs();
    while !walked.contains(&level) {
        union.extend(level.iter().copied());
        let next = level
            .iter()
            .flat_map(|window| {
                Picture {
                    side: 2,
                    cells: window.to_vec(),
                }
                .substitute(tile)
                .pairs()
            })
            .collect();
        walked.push(level);
        level = next;
    }
    union
        .into_iter()
        .map(|window| Picture {
            side: 2,
            cells: window.to_vec(),
        })
        .collect()
}

pub fn census(cells: &[u8], side: usize, count: usize, budget: u128) -> (Vec<i128>, bool) {
    let tile = Picture {
        side,
        cells: cells.iter().map(|&cell| u8::from(cell != 0)).collect(),
    };
    let mut pictures = language(&tile);
    let mut reach = 2;
    let mut out = Vec::with_capacity(count);
    for k in 1..=count {
        while reach / 2 + 1 < k {
            pictures = pictures.iter().map(|p| p.substitute(&tile)).collect();
            reach *= side;
        }
        let area = (reach as u128).pow(2) * pictures.len() as u128;
        if area > budget {
            return (out, true);
        }
        let mut seen = HashSet::new();
        for picture in &pictures {
            picture.blocks(k, &mut seen);
        }
        out.push(seen.len() as i128);
    }
    (out, false)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn brute(render: &Picture, top: usize) -> Vec<i128> {
        (1..=top)
            .map(|k| {
                let mut seen = HashSet::new();
                render.blocks(k, &mut seen);
                seen.len() as i128
            })
            .collect()
    }

    #[test]
    fn the_census_matches_a_brute_count_on_a_deep_render() {
        for (side, code) in [(3usize, 495u32), (3, 186), (2, 7), (3, 84), (4, 0x9669)] {
            let tile = Picture {
                side,
                cells: (0..side * side).map(|c| ((code >> c) & 1) as u8).collect(),
            };
            let mut render = tile.clone();
            while render.side * side <= 300 {
                render = render.substitute(&tile);
            }
            let (got, capped) = census(&tile.cells, side, 10, 1 << 20);
            assert!(!capped);
            assert_eq!(got, brute(&render, 10), "side {side} code {code}");
        }
    }
}
