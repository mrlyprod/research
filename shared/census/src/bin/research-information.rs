use census::save;
use mrly_core::error::Result;
use mrly_core::json;
use mrly_math::bang::Code;
use mrly_math::two::designs;

const NAME: &str = "research-information";
const CODE: u128 = 495;
const SPLIT: usize = 9;
const SIDE: usize = SPLIT * SPLIT;

// MASKS

fn render(level: usize) -> Result<Vec<Vec<bool>>> {
    let cells = designs::create(Code::from(CODE), 3, level, 0, 3)?;
    let types = cells.types();
    let side = cells.width();
    Ok((0..side)
        .map(|row| {
            (0..side)
                .map(|col| types.get(&[row, col]).is_ok_and(|v| v != 0))
                .collect()
        })
        .collect())
}

fn rearrange(cells: &[Vec<bool>]) -> Vec<Vec<bool>> {
    let mut out = vec![vec![false; SIDE]; SIDE];
    for block_row in 0..SPLIT {
        for block_col in 0..SPLIT {
            for row in 0..SPLIT {
                for col in 0..SPLIT {
                    out[block_row * SPLIT + block_col][row * SPLIT + col] =
                        cells[block_row * SPLIT + row][block_col * SPLIT + col];
                }
            }
        }
    }
    out
}

fn ones(mask: &[Vec<bool>]) -> usize {
    mask.iter().flatten().filter(|on| **on).count()
}

// WRITE

fn main() -> Result<()> {
    let whole = render(4)?;
    let factor = render(2)?;
    assert_eq!(whole.len(), SIDE);
    assert_eq!(factor.len(), SPLIT);
    assert_eq!(ones(&whole), 4096);
    assert_eq!(ones(&factor), 64);

    let rearranged = rearrange(&whole);
    assert_eq!(rearranged.len(), SIDE);
    assert_eq!(rearranged[0].len(), SIDE);
    assert_eq!(ones(&rearranged), 4096);

    let vector: Vec<bool> = factor.iter().flatten().copied().collect();
    assert_eq!(vector.len(), SIDE);
    for (row, line) in rearranged.iter().enumerate() {
        for (col, on) in line.iter().enumerate() {
            assert_eq!(*on, vector[row] && vector[col]);
        }
    }

    let rows: Vec<String> = rearranged
        .iter()
        .map(|line| line.iter().map(|on| if *on { '1' } else { '0' }).collect())
        .collect();
    save(NAME, &json!({ "rearranged": rows }))
}
