use census::save;
use mrly_core::error::Result;
use mrly_core::json;

const NAME: &str = "paper-slice-sign-even-half";
const THREE: [usize; 20] = [
    0, 1, 2, 4, 6, 8, 11, 14, 18, 22, 27, 32, 38, 44, 50, 57, 64, 72, 81, 89,
];
const OTHERS: [[usize; 3]; 4] = [[5, 4, 16], [7, 4, 26], [9, 4, 0], [11, 4, 0]];

fn depth(first: usize, second: usize, dim: usize) -> usize {
    if second > 0 && dim >= second {
        2
    } else if dim >= first {
        1
    } else {
        0
    }
}

fn main() -> Result<()> {
    let peak = *THREE.last().unwrap();
    assert_eq!(THREE.len(), 20);
    assert_eq!(peak, 89);
    let others: Vec<Vec<usize>> = OTHERS
        .iter()
        .map(|row| {
            (0..THREE.len())
                .map(|index| depth(row[1], row[2], 2 + 2 * index))
                .collect()
        })
        .collect();
    save(NAME, &json!({"three": THREE, "others": others}))
}
