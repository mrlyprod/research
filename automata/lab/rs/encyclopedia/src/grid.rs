use mrly_core::tensor::Tensor;
use mrly_math::two::Cell2d;
use mrly_core::Result;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Grid {
    pub h: usize,
    pub w: usize,
    pub bits: Vec<u8>,
}

impl Grid {
    pub fn new(h: usize, w: usize) -> Grid {
        Grid {
            h,
            w,
            bits: vec![0; h * w],
        }
    }

    pub fn full(h: usize, w: usize) -> Grid {
        Grid {
            h,
            w,
            bits: vec![1; h * w],
        }
    }

    pub fn of(h: usize, w: usize, bits: Vec<u8>) -> Grid {
        assert_eq!(bits.len(), h * w, "grid bits must fill the shape");
        Grid { h, w, bits }
    }

    pub fn from_rows(rows: &[&str]) -> Grid {
        let h = rows.len();
        let w = rows.first().map_or(0, |r| r.len());
        let bits = rows
            .iter()
            .flat_map(|r| r.bytes().map(|b| u8::from(b != b'.' && b != b'0')))
            .collect();
        Grid::of(h, w, bits)
    }

    pub fn from_tensor(t: &Tensor) -> Result<Grid> {
        Ok(Grid::of(t.shape[0], t.shape[1], t.bytes()?.to_vec()))
    }

    pub fn from_cell(c: &Cell2d) -> Result<Grid> {
        Grid::from_tensor(c.types())
    }

    pub fn tensor(&self) -> Result<Tensor> {
        Tensor::of(self.bits.clone(), vec![self.h, self.w])
    }

    pub fn cell(&self) -> Result<Cell2d> {
        Cell2d::new(self.tensor()?)
    }

    pub fn get(&self, r: usize, c: usize) -> u8 {
        self.bits[r * self.w + c]
    }

    pub fn set(&mut self, r: usize, c: usize, v: u8) {
        self.bits[r * self.w + c] = v;
    }

    pub fn population(&self) -> usize {
        self.bits.iter().filter(|&&b| b != 0).count()
    }

    pub fn is_empty(&self) -> bool {
        self.bits.iter().all(|&b| b == 0)
    }

    pub fn bbox(&self) -> Option<(usize, usize, usize, usize)> {
        let (mut r0, mut c0, mut r1, mut c1) = (usize::MAX, usize::MAX, 0, 0);
        for r in 0..self.h {
            let row = &self.bits[r * self.w..(r + 1) * self.w];
            let Some(first) = row.iter().position(|&b| b != 0) else {
                continue;
            };
            let last = row.iter().rposition(|&b| b != 0).unwrap();
            r0 = r0.min(r);
            r1 = r;
            c0 = c0.min(first);
            c1 = c1.max(last);
        }
        (r0 != usize::MAX).then_some((r0, c0, r1, c1))
    }

    pub fn window(&self, r0: usize, c0: usize, h: usize, w: usize) -> Grid {
        let mut out = Grid::new(h, w);
        for r in 0..h {
            out.bits[r * w..(r + 1) * w]
                .copy_from_slice(&self.bits[(r0 + r) * self.w + c0..(r0 + r) * self.w + c0 + w]);
        }
        out
    }

    pub fn crop(&self) -> Grid {
        match self.bbox() {
            Some((r0, c0, r1, c1)) => self.window(r0, c0, r1 - r0 + 1, c1 - c0 + 1),
            None => Grid::new(0, 0),
        }
    }

    pub fn embed(&self, h: usize, w: usize, r0: usize, c0: usize) -> Grid {
        let mut out = Grid::new(h, w);
        for r in 0..self.h {
            out.bits[(r0 + r) * w + c0..(r0 + r) * w + c0 + self.w]
                .copy_from_slice(&self.bits[r * self.w..(r + 1) * self.w]);
        }
        out
    }

    pub fn centred(&self, side: usize) -> Grid {
        self.embed(side, side, (side - self.h) / 2, (side - self.w) / 2)
    }

    pub fn kron(&self, other: &Grid) -> Grid {
        let (h, w) = (self.h * other.h, self.w * other.w);
        let mut out = Grid::new(h, w);
        for r in 0..h {
            for c in 0..w {
                out.bits[r * w + c] =
                    self.get(r / other.h, c / other.w) & other.get(r % other.h, c % other.w);
            }
        }
        out
    }

    pub fn tile(&self, t: usize) -> Grid {
        Grid::full(t, t).kron(self)
    }

    pub fn invert(&self) -> Grid {
        Grid::of(self.h, self.w, self.bits.iter().map(|&b| 1 - b).collect())
    }

    pub fn popped(&self) -> Grid {
        let mut out = self.clone();
        out.set((self.h - 1) / 2, (self.w - 1) / 2, 0);
        out
    }

    pub fn transpose(&self) -> Grid {
        let mut out = Grid::new(self.w, self.h);
        for r in 0..self.h {
            for c in 0..self.w {
                out.set(c, r, self.get(r, c));
            }
        }
        out
    }

    pub fn flip_rows(&self) -> Grid {
        let mut out = Grid::new(self.h, self.w);
        for r in 0..self.h {
            out.bits[r * self.w..(r + 1) * self.w]
                .copy_from_slice(&self.bits[(self.h - 1 - r) * self.w..(self.h - r) * self.w]);
        }
        out
    }

    pub fn flip_cols(&self) -> Grid {
        let mut out = self.clone();
        for r in 0..self.h {
            out.bits[r * self.w..(r + 1) * self.w].reverse();
        }
        out
    }

    pub fn symmetry(&self) -> usize {
        if self.h == 0 || self.w == 0 {
            return 8;
        }
        let a = self.flip_rows();
        let b = self.flip_cols();
        let mut order = 1;
        order += usize::from(a == *self);
        order += usize::from(b == *self);
        order += usize::from(a.flip_cols() == *self);
        if self.h == self.w {
            let t = self.transpose();
            order += usize::from(t == *self);
            order += usize::from(t.flip_rows() == *self);
            order += usize::from(t.flip_cols() == *self);
            order += usize::from(t.flip_rows().flip_cols() == *self);
        }
        order
    }

    pub fn split(&self, d: usize) -> Option<(Grid, Grid)> {
        if d == 0 || !self.h.is_multiple_of(d) || !self.w.is_multiple_of(d) {
            return None;
        }
        let (qh, qw) = (self.h / d, self.w / d);
        let mut quotient = Grid::new(qh, qw);
        let mut tile: Option<Grid> = None;
        for i in 0..qh {
            for j in 0..qw {
                let block = self.window(i * d, j * d, d, d);
                if block.is_empty() {
                    continue;
                }
                match &tile {
                    None => tile = Some(block),
                    Some(t) if *t == block => {}
                    Some(_) => return None,
                }
                quotient.set(i, j, 1);
            }
        }
        tile.map(|t| (quotient, t))
    }

    pub fn cuts(&self) -> (Vec<usize>, Vec<usize>) {
        let (mut proper, mut tilings) = (Vec::new(), Vec::new());
        let top = self.h.max(self.w);
        for d in 2..top {
            if let Some((quotient, _)) = self.split(d) {
                if quotient.population() == quotient.bits.len() {
                    tilings.push(d);
                } else {
                    proper.push(d);
                }
            }
        }
        (proper, tilings)
    }

    pub fn hex_rows(&self) -> Vec<String> {
        (0..self.h)
            .map(|r| {
                let row = &self.bits[r * self.w..(r + 1) * self.w];
                row.chunks(4)
                    .map(|nib| {
                        let v = nib
                            .iter()
                            .enumerate()
                            .fold(0u8, |acc, (i, &b)| acc | (b << (3 - i)));
                        char::from_digit(u32::from(v), 16).unwrap()
                    })
                    .collect()
            })
            .collect()
    }

    pub fn hash(&self) -> u128 {
        let mut a: u64 = 0xcbf29ce484222325;
        let mut b: u64 = 0x9e3779b97f4a7c15;
        for chunk in self.bits.chunks(64) {
            let mut word = 0u64;
            for (i, &bit) in chunk.iter().enumerate() {
                word |= u64::from(bit) << i;
            }
            a = (a ^ word).wrapping_mul(0x100000001b3);
            b = (b.rotate_left(29) ^ word).wrapping_mul(0xbf58476d1ce4e5b9);
        }
        a = (a ^ (self.w as u64)).wrapping_mul(0x100000001b3);
        b = (b ^ (self.h as u64)).wrapping_mul(0x94d049bb133111eb);
        (u128::from(a) << 64) | u128::from(b)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_tiled_carpet_splits_at_its_tile_and_a_shifted_one_does_not() {
        let carpet = Grid::from_rows(&["111", "101", "111"]);
        let tiled = carpet.tile(3);
        let (proper, tilings) = tiled.cuts();
        assert_eq!(tilings, vec![3]);
        assert!(proper.is_empty());
        let mut holed = tiled.clone();
        for r in 3..6 {
            for c in 3..6 {
                holed.set(r, c, 0);
            }
        }
        assert_eq!(holed.cuts(), (vec![3], vec![]));
        assert_eq!(carpet.symmetry(), 8);
        assert_eq!(Grid::from_rows(&["110", "000", "000"]).symmetry(), 1);
        assert_eq!(carpet.hex_rows(), vec!["e", "a", "e"]);
        assert_eq!(carpet.crop().bbox(), Some((0, 0, 2, 2)));
    }
}
