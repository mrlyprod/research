use census::save;
use mrly_core::error::Result;
use mrly_core::json;

const NAME: &str = "paper-spin-harmonics";
const SIDE: usize = 3;
const CODE: u32 = 45;

type Cell = (usize, usize);
type Chord = (Cell, Cell, i64);

fn cells() -> Vec<Cell> {
    (0..SIDE * SIDE)
        .filter(|j| CODE >> j & 1 == 1)
        .map(|j| (j / SIDE, j % SIDE))
        .collect()
}

fn chords(cells: &[Cell]) -> Vec<Chord> {
    let mut out = Vec::new();
    for (a, p) in cells.iter().enumerate() {
        for q in &cells[a + 1..] {
            let dr = p.0 as i64 - q.0 as i64;
            let dc = p.1 as i64 - q.1 as i64;
            out.push((*p, *q, dr * dr + dc * dc));
        }
    }
    out
}

fn main() -> Result<()> {
    let cells = cells();
    let chords = chords(&cells);
    let long = chords.iter().filter(|c| c.2 == 4).count();
    assert_eq!(cells.len(), 4, "code 45 fills four cells");
    assert_eq!(chords.len(), 6, "four centres carry six chords");
    assert_eq!(long, 2, "two chords have squared length 4/9");
    save(
        NAME,
        &json!({
            "side": SIDE,
            "cells": cells.iter().map(|c| [c.0, c.1]).collect::<Vec<_>>(),
            "chords": chords
                .iter()
                .map(|c| [c.0 .0 as i64, c.0 .1 as i64, c.1 .0 as i64, c.1 .1 as i64, c.2])
                .collect::<Vec<_>>(),
        }),
    )
}
