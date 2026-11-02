use crate::design::Shape;
use crate::grid::Grid;
use mrly_core::error::value_error;
use mrly_core::Result;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MaskSpec {
    Design {
        shape: Shape,
        side: usize,
        level: usize,
    },
    Copy {
        inverted: bool,
    },
}

impl MaskSpec {
    pub fn token(self) -> String {
        match self {
            MaskSpec::Design { shape, side, level } => format!("{}.{side}.{level}", shape.token()),
            MaskSpec::Copy { inverted: false } => "copy".to_string(),
            MaskSpec::Copy { inverted: true } => "copyinv".to_string(),
        }
    }

    pub fn parse(token: &str) -> Result<MaskSpec> {
        match token {
            "copy" => Ok(MaskSpec::Copy { inverted: false }),
            "copyinv" => Ok(MaskSpec::Copy { inverted: true }),
            _ => {
                let (shape, side, level) = parse_design(token)?;
                Ok(MaskSpec::Design { shape, side, level })
            }
        }
    }
}

pub fn parse_design(token: &str) -> Result<(Shape, usize, usize)> {
    let parts: Vec<&str> = token.split('.').collect();
    if parts.len() != 3 {
        return value_error(format!("design token {token:?} wants shape.side.level."));
    }
    let shape = Shape::parse(parts[0])?;
    let (Ok(side), Ok(level)) = (parts[1].parse::<usize>(), parts[2].parse::<usize>()) else {
        return value_error(format!("design token {token:?} wants numbers."));
    };
    Ok((shape, side, level))
}

#[derive(Clone, Debug)]
pub struct Factor {
    pub n: usize,
    pub taps: Vec<(i64, i64)>,
}

impl Factor {
    pub fn of(grid: &Grid) -> Factor {
        assert_eq!(grid.h, grid.w, "a factor is square");
        let n = grid.h;
        let centre = (n as i64 - 1) / 2;
        let mut taps = Vec::new();
        for r in 0..n {
            for c in 0..n {
                if grid.get(r, c) != 0 {
                    taps.push((r as i64 - centre, c as i64 - centre));
                }
            }
        }
        Factor { n, taps }
    }

    pub fn ones(n: usize) -> Factor {
        Factor::of(&Grid::full(n, n))
    }
}

#[derive(Clone, Debug)]
pub struct Kernel {
    pub side: usize,
    pub reach: usize,
    pub budget: usize,
    pub taps: usize,
    terms: Vec<(i32, Vec<Factor>)>,
    centre: i32,
}

impl Kernel {
    pub fn of(terms: Vec<(i32, Vec<Factor>)>, centre_value: u8) -> Result<Kernel> {
        let side = terms
            .first()
            .map(|(_, chain)| chain.iter().map(|f| f.n).product::<usize>())
            .unwrap_or(1);
        if terms
            .iter()
            .any(|(_, chain)| chain.iter().map(|f| f.n).product::<usize>() != side)
        {
            return value_error("every term of a kernel spans the same side.");
        }
        if side.is_multiple_of(2) {
            return value_error("a kernel side is odd.");
        }
        let taps = terms
            .iter()
            .map(|(_, chain)| chain.iter().map(|f| f.taps.len()).sum::<usize>())
            .sum();
        let mut kernel = Kernel {
            side,
            reach: (side - 1) / 2,
            budget: 0,
            taps,
            terms,
            centre: i32::from(centre_value),
        };
        let dense = kernel.dense();
        if dense.get((side - 1) / 2, (side - 1) / 2) != 0 {
            return value_error("the kernel centre must be popped.");
        }
        kernel.budget = dense.population();
        Ok(kernel)
    }

    pub fn from_grid(mask: &Grid) -> Result<Kernel> {
        if mask.h != mask.w {
            return value_error("a mask is square.");
        }
        let centre = mask.get((mask.h - 1) / 2, (mask.w - 1) / 2);
        Kernel::of(vec![(1, vec![Factor::of(mask)])], centre)
    }

    pub fn design(shape: Shape, side: usize, level: usize) -> Result<Kernel> {
        let base = shape.build(side, 1)?;
        let full = shape.build(side, level)?;
        let centre = full.get((full.h - 1) / 2, (full.w - 1) / 2);
        let chain: Vec<Factor> = (0..level).map(|_| Factor::of(&base)).collect();
        Kernel::of(vec![(1, chain)], centre)
    }

    pub fn copy(
        shape: Shape,
        side: usize,
        level: usize,
        tiling: usize,
        inverted: bool,
    ) -> Result<Kernel> {
        let base = shape.build(side, 1)?;
        let board = shape.build(side, level)?.tile(tiling);
        let centre = board.get((board.h - 1) / 2, (board.w - 1) / 2);
        let mut chain = Vec::new();
        if tiling > 1 {
            chain.push(Factor::ones(tiling));
        }
        chain.extend((0..level).map(|_| Factor::of(&base)));
        if !inverted {
            return Kernel::of(vec![(1, chain)], centre);
        }
        let mut full = Vec::new();
        if tiling > 1 {
            full.push(Factor::ones(tiling));
        }
        full.extend((0..level).map(|_| Factor::ones(side)));
        Kernel::of(vec![(1, full), (-1, chain)], 1 - centre)
    }

    pub fn dense(&self) -> Grid {
        let n = self.side;
        let mut acc = vec![0i32; n * n];
        for (sign, chain) in &self.terms {
            let mut grid = Grid::full(1, 1);
            for factor in chain {
                let mut g = Grid::new(factor.n, factor.n);
                let centre = (factor.n as i64 - 1) / 2;
                for &(dr, dc) in &factor.taps {
                    g.set((dr + centre) as usize, (dc + centre) as usize, 1);
                }
                grid = grid.kron(&g);
            }
            for (a, &b) in acc.iter_mut().zip(&grid.bits) {
                *a += sign * i32::from(b);
            }
        }
        acc[(n / 2) * n + n / 2] -= self.centre;
        Grid::of(
            n,
            n,
            acc.iter()
                .map(|&v| {
                    assert!((0..=1).contains(&v), "a kernel is a 0/1 mask");
                    v as u8
                })
                .collect(),
        )
    }

    pub fn count(&self, grid: &Grid, wrap: bool) -> Vec<i32> {
        if wrap {
            return self.chain(grid, true);
        }
        let r = self.reach;
        let grown = self.chain(&grid.embed(grid.h + 2 * r, grid.w + 2 * r, r, r), false);
        let w = grid.w + 2 * r;
        let mut out = Vec::with_capacity(grid.h * grid.w);
        for row in r..r + grid.h {
            out.extend_from_slice(&grown[row * w + r..row * w + r + grid.w]);
        }
        out
    }

    pub fn count_grown(&self, grid: &Grid) -> (Grid, Vec<i32>) {
        let r = self.reach;
        let grown = grid.embed(grid.h + 2 * r, grid.w + 2 * r, r, r);
        let counts = self.chain(&grown, false);
        (grown, counts)
    }

    fn chain(&self, grid: &Grid, wrap: bool) -> Vec<i32> {
        let (h, w) = (grid.h, grid.w);
        let base: Vec<i32> = grid.bits.iter().map(|&b| i32::from(b)).collect();
        let mut out = vec![0i32; h * w];
        for (sign, chain) in &self.terms {
            let mut acc = base.clone();
            let mut dilation = 1i64;
            for factor in chain.iter().rev() {
                acc = convolve(&acc, h, w, factor, dilation, wrap);
                dilation *= factor.n as i64;
            }
            for (o, a) in out.iter_mut().zip(&acc) {
                *o += sign * a;
            }
        }
        if self.centre != 0 {
            for (o, b) in out.iter_mut().zip(&base) {
                *o -= self.centre * b;
            }
        }
        out
    }
}

fn convolve(src: &[i32], h: usize, w: usize, factor: &Factor, d: i64, wrap: bool) -> Vec<i32> {
    let mut out = vec![0i32; h * w];
    for &(tr, tc) in &factor.taps {
        let (dr, dc) = (tr * d, tc * d);
        if wrap {
            let (hh, ww) = (h as i64, w as i64);
            for r in 0..h {
                let sr = (r as i64 + dr).rem_euclid(hh) as usize;
                for c in 0..w {
                    let sc = (c as i64 + dc).rem_euclid(ww) as usize;
                    out[r * w + c] += src[sr * w + sc];
                }
            }
        } else {
            shift_add(&mut out, src, h, w, dr, dc);
        }
    }
    out
}

fn shift_add(out: &mut [i32], src: &[i32], h: usize, w: usize, dr: i64, dc: i64) {
    let (hh, ww) = (h as i64, w as i64);
    let (r_lo, r_hi) = ((-dr).max(0), (hh - dr).min(hh));
    let (c_lo, c_hi) = ((-dc).max(0), (ww - dc).min(ww));
    if r_lo >= r_hi || c_lo >= c_hi {
        return;
    }
    let (c_lo, c_hi) = (c_lo as usize, c_hi as usize);
    for r in r_lo as usize..r_hi as usize {
        let sr = (r as i64 + dr) as usize;
        let sc = (c_lo as i64 + dc) as usize;
        let o = &mut out[r * w + c_lo..r * w + c_hi];
        let s = &src[sr * w + sc..sr * w + sc + (c_hi - c_lo)];
        for (a, b) in o.iter_mut().zip(s) {
            *a += *b;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::standard_shapes;
    use mrly_core::Rng;
    use mrly_life::{design_mask, moore, next_grid, Boundary};
    use mrly_math::bang::Code;

    fn noise(h: usize, w: usize, seed: u64) -> Grid {
        let mut rng = Rng::new(seed);
        Grid::of(
            h,
            w,
            (0..h * w).map(|_| u8::from(rng.chance(0.4))).collect(),
        )
    }

    fn mrly_counts(grid: &Grid, mask: &Grid, wrap: bool) -> Vec<i32> {
        let boundary = if wrap {
            Boundary::Wrap
        } else {
            Boundary::Constant
        };
        let mut out = vec![0i32; grid.h * grid.w];
        for n in 0..=mask.population() {
            let fired = next_grid(
                &grid.cell().unwrap(),
                &[n],
                &[n],
                &mask.tensor().unwrap(),
                boundary,
            )
            .unwrap();
            let fired = Grid::from_cell(&fired).unwrap();
            for (o, &b) in out.iter_mut().zip(&fired.bits) {
                if b != 0 {
                    *o = n as i32;
                }
            }
        }
        out
    }

    #[test]
    fn every_grammar_kernel_is_the_popped_design_and_counts_as_the_engine_does() {
        let shapes = standard_shapes();
        let mut checked = 0;
        for &shape in &shapes {
            for side in [3, 5, 7, 9] {
                for level in [1, 2] {
                    let kernel = Kernel::design(shape, side, level).unwrap();
                    assert_eq!(kernel.dense(), shape.build(side, level).unwrap().popped());
                    if let Shape::Code(code) = shape {
                        let reference = design_mask(2, Code::from(code), side, level).unwrap();
                        assert_eq!(kernel.dense(), Grid::from_tensor(&reference).unwrap());
                    }
                    checked += 1;
                }
            }
            for tiling in [1, 3, 5] {
                for inverted in [false, true] {
                    let kernel = Kernel::copy(shape, 3, 2, tiling, inverted).unwrap();
                    let board = shape.build(3, 2).unwrap().tile(tiling);
                    let want = if inverted { board.invert() } else { board };
                    assert_eq!(kernel.dense(), want.popped());
                }
            }
        }
        assert_eq!(checked, 15 * 8);
    }

    #[test]
    fn kernel_counts_match_the_engine_on_both_boundaries() {
        let grid = noise(31, 31, 7);
        let small = noise(13, 13, 8);
        let cases = [
            Kernel::design(Shape::parse("vtree").unwrap(), 3, 2).unwrap(),
            Kernel::design(Shape::parse("code3").unwrap(), 5, 1).unwrap(),
            Kernel::copy(Shape::parse("net").unwrap(), 3, 1, 3, true).unwrap(),
            Kernel::copy(Shape::parse("htree").unwrap(), 5, 1, 3, false).unwrap(),
            Kernel::from_grid(&Grid::from_cell(&moore().unwrap()).unwrap()).unwrap(),
        ];
        for kernel in &cases {
            let mask = kernel.dense();
            for wrap in [false, true] {
                assert_eq!(kernel.count(&grid, wrap), mrly_counts(&grid, &mask, wrap));
                assert_eq!(kernel.count(&small, wrap), mrly_counts(&small, &mask, wrap));
            }
        }
    }
}
