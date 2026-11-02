use census::save;
use mrly_core::error::Result;
use mrly_core::{json, Rng};
use mrly_math::bang::Code;
use mrly_math::two::designs;
use mrly_math::two::Cell2d;
use std::collections::VecDeque;

const NAME: &str = "research-walks";
const LEVEL: usize = 4;
const STEPS: usize = 20000;
const SEED: u64 = 7;

fn sites(cell: &Cell2d) -> (usize, Vec<bool>) {
    let side = cell.width();
    let types = cell.types();
    let mut mask = vec![false; side * side];
    for row in 0..side {
        for col in 0..side {
            mask[row * side + col] = types.get(&[row, col]).is_ok_and(|v| v != 0);
        }
    }
    (side, mask)
}

fn around(side: usize, at: usize) -> Vec<usize> {
    let (row, col) = (at / side, at % side);
    let mut out = Vec::with_capacity(4);
    if row > 0 {
        out.push(at - side);
    }
    if row + 1 < side {
        out.push(at + side);
    }
    if col > 0 {
        out.push(at - 1);
    }
    if col + 1 < side {
        out.push(at + 1);
    }
    out
}

fn giant(side: usize, mask: &[bool]) -> Vec<usize> {
    let mut seen = vec![false; mask.len()];
    let mut best = Vec::new();
    for start in 0..mask.len() {
        if !mask[start] || seen[start] {
            continue;
        }
        let mut part = Vec::new();
        let mut queue = VecDeque::from([start]);
        seen[start] = true;
        while let Some(at) = queue.pop_front() {
            part.push(at);
            for next in around(side, at) {
                if mask[next] && !seen[next] {
                    seen[next] = true;
                    queue.push_back(next);
                }
            }
        }
        if part.len() > best.len() {
            best = part;
        }
    }
    best
}

fn middle(side: usize, part: &[usize]) -> usize {
    let mid = (side as f64 - 1.0) / 2.0;
    *part
        .iter()
        .min_by(|a, b| {
            let reach = |at: &usize| {
                let (row, col) = ((at / side) as f64, (at % side) as f64);
                (row - mid).hypot(col - mid)
            };
            reach(a).partial_cmp(&reach(b)).unwrap()
        })
        .unwrap()
}

fn walk(side: usize, mask: &[bool], start: usize, seed: u64) -> Vec<usize> {
    let mut rng = Rng::new(seed);
    let mut at = start;
    let mut path = vec![at];
    for _ in 0..STEPS {
        let (row, col) = (at / side, at % side);
        let next = match rng.below(4) {
            0 if row > 0 => at - side,
            1 if row + 1 < side => at + side,
            2 if col > 0 => at - 1,
            3 if col + 1 < side => at + 1,
            _ => at,
        };
        if next != at && mask[next] {
            at = next;
            path.push(at);
        }
    }
    path
}

fn visited(code: u128) -> Result<Vec<usize>> {
    let cell = designs::create(Code::from(code), 3, LEVEL, 0, 3)?;
    let (side, mask) = sites(&cell);
    let part = giant(side, &mask);
    let trail = walk(side, &mask, middle(side, &part), SEED);
    let mut seen = vec![false; mask.len()];
    let mut order = Vec::new();
    for at in trail {
        if !seen[at] {
            seen[at] = true;
            order.push(at);
        }
    }
    assert_eq!(
        mask.iter().filter(|f| **f).count(),
        7usize.pow(LEVEL as u32)
    );
    Ok(order)
}

fn main() -> Result<()> {
    let a = visited(127)?;
    let b = visited(239)?;
    save(NAME, &json!({"127": a, "239": b}))
}
