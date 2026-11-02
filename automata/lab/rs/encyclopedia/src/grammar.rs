use crate::design::{distinct, standard_shapes, Shape};
use crate::grid::Grid;
use crate::mask::MaskSpec;
use crate::recipe::{Mode, Recipe, Side};
use mrly_life::Source;
use mrly_core::Result;

#[derive(Clone, Debug)]
pub struct Grammar {
    pub shapes: Vec<Shape>,
    pub sides: Vec<usize>,
    pub levels: Vec<usize>,
    pub tilings: Vec<usize>,
    pub wallpaper: bool,
    pub mask_sides: Vec<usize>,
    pub mask_levels: Vec<usize>,
    pub copies: bool,
    pub sources: Vec<Source>,
    pub rolls: Vec<(bool, bool)>,
    pub life: bool,
    pub pairs: Vec<(Source, Source)>,
    pub x: usize,
    pub cap: usize,
}

pub const ROLLS: [(bool, bool); 4] = [(false, false), (true, false), (false, true), (true, true)];

pub const PAIRS: [(Source, Source); 8] = [
    (Source::CarpetFills, Source::CarpetVoids),
    (Source::NetFills, Source::NetVoids),
    (Source::TreeFills, Source::TreeVoids),
    (Source::VoidFills, Source::VoidVoids),
    (Source::PointFills, Source::PointVoids),
    (Source::DustFills, Source::DustVoids),
    (Source::LineFills, Source::LineVoids),
    (Source::StarFills, Source::StarVoids),
];

impl Grammar {
    pub fn standard(x: usize, cap: usize) -> Grammar {
        Grammar {
            shapes: standard_shapes(),
            sides: vec![3, 5, 7],
            levels: vec![1, 2],
            tilings: vec![1, 3, 5],
            wallpaper: true,
            mask_sides: vec![3, 5, 7, 9],
            mask_levels: vec![1, 2],
            copies: true,
            sources: Source::all().to_vec(),
            rolls: ROLLS.to_vec(),
            life: true,
            pairs: PAIRS.to_vec(),
            x,
            cap,
        }
    }

    pub fn seeds(&self) -> Result<Vec<(Shape, usize, usize, Grid)>> {
        let mut out = Vec::new();
        for &side in &self.sides {
            for &level in &self.levels {
                for (shape, grid) in distinct(&self.shapes, side, level)? {
                    out.push((shape, side, level, grid));
                }
            }
        }
        Ok(out)
    }

    pub fn masks(&self) -> Result<Vec<(MaskSpec, Grid)>> {
        let mut out: Vec<(MaskSpec, Grid)> = Vec::new();
        for &side in &self.mask_sides {
            for &level in &self.mask_levels {
                for (shape, grid) in distinct(&self.shapes, side, level)? {
                    let popped = grid.popped();
                    if popped.is_empty() || out.iter().any(|(_, g)| *g == popped) {
                        continue;
                    }
                    out.push((MaskSpec::Design { shape, side, level }, popped));
                }
            }
        }
        Ok(out)
    }

    pub fn rules(&self) -> Vec<(Side, Side)> {
        let mut out = Vec::new();
        for &seq in &self.sources {
            for &(zeros, ones) in &self.rolls {
                out.push((Side::drawn(seq, zeros, ones), Side::drawn(seq, zeros, ones)));
            }
        }
        if self.life {
            out.push((Side::list(&[3]), Side::list(&[2, 3])));
        }
        for &(fills, voids) in &self.pairs {
            for &(zeros, ones) in &self.rolls {
                out.push((
                    Side::drawn(fills, zeros, ones),
                    Side::drawn(voids, zeros, ones),
                ));
            }
        }
        out
    }

    fn copies_for(&self, board: &Grid, masks: &[(MaskSpec, Grid)]) -> Vec<MaskSpec> {
        if !self.copies {
            return Vec::new();
        }
        let mut out = Vec::new();
        for inverted in [false, true] {
            let grid = if inverted {
                board.invert().popped()
            } else {
                board.popped()
            };
            if !grid.is_empty() && !masks.iter().any(|(_, g)| *g == grid) {
                out.push(MaskSpec::Copy { inverted });
            }
        }
        out
    }

    pub fn walk(&self, from: usize, to: usize, mut f: impl FnMut(usize, Recipe)) -> Result<usize> {
        let seeds = self.seeds()?;
        let masks = self.masks()?;
        let rules = self.rules();
        let mut index = 0;
        let mut slots: Vec<(usize, Mode, Vec<MaskSpec>)> = Vec::new();
        for (shape, side, level, grid) in &seeds {
            slots.clear();
            for &tiling in &self.tilings {
                let mut list: Vec<MaskSpec> = masks.iter().map(|(m, _)| *m).collect();
                list.extend(self.copies_for(&grid.tile(tiling), &masks));
                slots.push((tiling, Mode::Island, list));
            }
            if self.wallpaper {
                let mut list: Vec<MaskSpec> = masks.iter().map(|(m, _)| *m).collect();
                list.extend(self.copies_for(grid, &masks));
                slots.push((1, Mode::Wallpaper, list));
            }
            for (tiling, mode, list) in &slots {
                let block = list.len() * rules.len();
                if index + block <= from || index >= to {
                    index += block;
                    continue;
                }
                for mask in list {
                    for (birth, survive) in &rules {
                        if index >= from && index < to {
                            f(
                                index,
                                Recipe {
                                    seed: *shape,
                                    side: *side,
                                    level: *level,
                                    tiling: *tiling,
                                    mode: *mode,
                                    mask: *mask,
                                    birth: birth.clone(),
                                    survive: survive.clone(),
                                    x: self.x,
                                    cap: if *mode == Mode::Island { self.cap } else { 0 },
                                },
                            );
                        }
                        index += 1;
                    }
                }
            }
        }
        Ok(index)
    }

    pub fn count(&self) -> Result<usize> {
        self.walk(0, 0, |_, _| {})
    }

    pub fn recipes(&self, from: usize, to: usize) -> Result<Vec<Recipe>> {
        let mut out = Vec::new();
        self.walk(from, to, |_, recipe| out.push(recipe))?;
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_grammar_counts_its_recipes_and_every_name_regenerates_its_recipe() {
        let mut small = Grammar::standard(64, 101);
        small.sides = vec![3];
        small.levels = vec![1];
        small.mask_sides = vec![3];
        small.mask_levels = vec![1];
        small.sources = vec![Source::Binary];
        small.rolls = vec![(false, false)];
        small.pairs = vec![];
        let seeds = small.seeds().unwrap();
        let masks = small.masks().unwrap();
        let rules = small.rules();
        assert_eq!(rules.len(), 2);
        let all = small.recipes(0, usize::MAX).unwrap();
        assert_eq!(all.len(), small.count().unwrap());
        assert!(all.len() >= seeds.len() * (3 + 1) * masks.len() * rules.len());
        let names: std::collections::HashSet<String> = all.iter().map(Recipe::name).collect();
        assert_eq!(names.len(), all.len());
        for recipe in &all {
            assert_eq!(Recipe::parse(&recipe.name()).unwrap(), *recipe);
        }
        let slice = small.recipes(5, 9).unwrap();
        assert_eq!(slice, all[5..9].to_vec());
        let copies = all
            .iter()
            .filter(|r| matches!(r.mask, MaskSpec::Copy { .. }))
            .count();
        assert!(copies > 0);
        assert!(all
            .iter()
            .all(|r| r.tiling > 1 || r.mask != MaskSpec::Copy { inverted: false }));
    }
}
