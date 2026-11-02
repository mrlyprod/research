use crate::design::{classes, Shape};
use crate::engine::{torus, Fate, Rule};
use crate::grid::Grid;
use crate::mask::Kernel;
use crate::recipe::Side;
use mrly_life::Source;
use mrly_core::Result;

pub type Table = Vec<(String, [usize; 4])>;

pub const CANVAS: usize = 27;
pub const GENERATIONS: usize = 64;
pub const CLAIMED_NUMBERS: [usize; 4] = [106, 24, 56, 132];
pub const ATLAS_TODAY: [(&str, [usize; 4]); 24] = [
    ("evens", [20, 4, 7, 22]),
    ("odds", [44, 1, 0, 8]),
    ("random_0", [7, 3, 15, 28]),
    ("primes", [4, 3, 18, 28]),
    ("binary", [20, 4, 8, 21]),
    ("fibonacci", [4, 3, 17, 29]),
    ("grid_squares", [53, 0, 0, 0]),
    ("carpet_fills", [53, 0, 0, 0]),
    ("carpet_voids", [46, 0, 0, 7]),
    ("net_fills", [52, 1, 0, 0]),
    ("net_voids", [46, 0, 0, 7]),
    ("tree_fills", [53, 0, 0, 0]),
    ("tree_voids", [45, 1, 0, 7]),
    ("void_fills", [53, 0, 0, 0]),
    ("void_voids", [47, 0, 0, 6]),
    ("point_fills", [46, 0, 0, 7]),
    ("point_voids", [53, 0, 0, 0]),
    ("dust_fills", [46, 0, 0, 7]),
    ("dust_voids", [52, 1, 0, 0]),
    ("line_fills", [45, 1, 0, 7]),
    ("line_voids", [53, 0, 0, 0]),
    ("star_fills", [47, 0, 0, 6]),
    ("star_voids", [53, 0, 0, 0]),
    ("b3s23", [31, 17, 4, 1]),
];
pub const CLAIMED_DESIGNS: [usize; 4] = [448, 2, 0, 27];
pub const CLAIMED_LIFE: [usize; 4] = [31, 17, 4, 1];

pub const OLD_DESIGNS: [Source; 9] = [
    Source::GridSquares,
    Source::CarpetFills,
    Source::CarpetVoids,
    Source::NetFills,
    Source::NetVoids,
    Source::TreeFills,
    Source::TreeVoids,
    Source::VoidFills,
    Source::VoidVoids,
];

pub fn rules() -> Vec<(String, Side, Side)> {
    let mut out: Vec<(String, Side, Side)> = Source::all()
        .into_iter()
        .map(|seq| {
            let side = Side::drawn(seq, false, false);
            (seq.name(), side.clone(), side)
        })
        .collect();
    out.push(("b3s23".to_string(), Side::list(&[3]), Side::list(&[2, 3])));
    out
}

pub fn masks(board: &Grid) -> Result<Vec<Grid>> {
    let mut out: Vec<Grid> = Vec::new();
    let mut push = |grid: Grid| {
        if !out.contains(&grid) {
            out.push(grid);
        }
    };
    push(Shape::Code(7).build(3, 1)?.popped());
    for code in classes() {
        if code != 7 {
            push(Shape::Code(code).build(3, 1)?.popped());
        }
    }
    push(board.popped());
    push(board.invert().popped());
    Ok(out)
}

pub fn mini() -> Result<(usize, Table)> {
    let rules = rules();
    let mut table: Table = rules
        .iter()
        .map(|(name, _, _)| (name.clone(), [0; 4]))
        .collect();
    let mut runs = 0;
    for code in classes() {
        let seed = Shape::Code(code).build(3, 1)?;
        for tessellation in [1, 3] {
            let board = seed.tile(tessellation);
            for mask in masks(&board)? {
                let kernel = Kernel::from_grid(&mask)?;
                for (i, (_, birth, survive)) in rules.iter().enumerate() {
                    let rule = Rule::new(
                        &birth.values(kernel.budget)?,
                        &survive.values(kernel.budget)?,
                        kernel.budget,
                    );
                    let (fate, _) = torus(&board, &kernel, &rule, CANVAS, GENERATIONS);
                    let slot = match fate {
                        Fate::Dead => 0,
                        Fate::Alive => 1,
                        Fate::Loop => 2,
                        Fate::Timeout => 3,
                    };
                    table[i].1[slot] += 1;
                    runs += 1;
                }
            }
        }
    }
    Ok((runs, table))
}

pub fn sum(table: &[(String, [usize; 4])], names: &[String]) -> [usize; 4] {
    let mut out = [0; 4];
    for (name, counts) in table {
        if names.contains(name) {
            for k in 0..4 {
                out[k] += counts[k];
            }
        }
    }
    out
}

pub struct Summary {
    pub runs: usize,
    pub numbers: [usize; 4],
    pub old_designs: [usize; 4],
    pub all_designs: [usize; 4],
    pub life: [usize; 4],
    pub table: Table,
}

pub fn summary() -> Result<Summary> {
    let (runs, table) = mini()?;
    let names = |list: &[Source]| list.iter().map(|s| s.name()).collect::<Vec<_>>();
    Ok(Summary {
        runs,
        numbers: sum(&table, &names(&Source::numbers())),
        old_designs: sum(&table, &names(&OLD_DESIGNS)),
        all_designs: sum(&table, &names(&Source::designs())),
        life: sum(&table, &["b3s23".to_string()]),
        table,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_mini_control_replays_the_old_fate_table() {
        let s = summary().unwrap();
        println!("mini control: {} runs of the 848", s.runs);
        println!("numbers {:?} claimed {:?}", s.numbers, CLAIMED_NUMBERS);
        println!(
            "old designs {:?} claimed {:?}",
            s.old_designs, CLAIMED_DESIGNS
        );
        println!("all designs {:?}", s.all_designs);
        println!("life {:?} claimed {:?}", s.life, CLAIMED_LIFE);
        assert_eq!(s.runs, 1272);
        for (name, counts) in ATLAS_TODAY {
            let row = s.table.iter().find(|(n, _)| n == name).unwrap();
            assert_eq!(row.1, counts, "{name}");
        }
        assert_eq!(s.life, CLAIMED_LIFE);
        assert_eq!(s.old_designs, CLAIMED_DESIGNS);
        assert_ne!(s.numbers, CLAIMED_NUMBERS);
    }
}
