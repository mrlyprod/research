use census::save;
use mrly_core::error::Result;
use mrly_core::{json, Rng};
use mrly_math::bang::factory;
use mrly_math::bang::Code;

const NAME: &str = "research-dimers";
const LEVEL: usize = 3;
const SIDE: usize = 27;
const CELLS: usize = 512;
const SEED: u64 = 3;
const FLIPS: usize = 400_000;
const NONE: usize = usize::MAX;

fn neighbours(filled: &[bool], cell: usize) -> Vec<usize> {
    let (r, c) = (cell / SIDE, cell % SIDE);
    let mut out = Vec::with_capacity(4);
    if r > 0 {
        out.push(cell - SIDE);
    }
    if c + 1 < SIDE {
        out.push(cell + 1);
    }
    if r + 1 < SIDE {
        out.push(cell + SIDE);
    }
    if c > 0 {
        out.push(cell - 1);
    }
    out.retain(|&n| filled[n]);
    out
}

fn augment(links: &[Vec<usize>], mate: &mut [usize], seen: &mut [bool], black: usize) -> bool {
    for &white in &links[black] {
        if std::mem::replace(&mut seen[white], true) {
            continue;
        }
        if mate[white] == NONE || augment(links, mate, seen, mate[white]) {
            mate[white] = black;
            mate[black] = white;
            return true;
        }
    }
    false
}

fn flip(filled: &[bool], mate: &mut [usize], r: usize, c: usize) {
    let a = r * SIDE + c;
    let (b, d, e) = (a + 1, a + SIDE, a + SIDE + 1);
    if ![a, b, d, e].iter().all(|&x| filled[x]) {
        return;
    }
    if mate[a] == b && mate[d] == e {
        (mate[a], mate[d], mate[b], mate[e]) = (d, a, e, b);
    } else if mate[a] == d && mate[b] == e {
        (mate[a], mate[b], mate[d], mate[e]) = (b, a, e, d);
    }
}

fn main() -> Result<()> {
    let design = factory::create(Code::from(495u128), 3, 2, 3, LEVEL)?;
    assert_eq!(design.shape, vec![SIDE, SIDE]);
    let filled: Vec<bool> = (0..SIDE * SIDE).map(|f| design.at(f) == 1).collect();
    assert_eq!(filled.iter().filter(|&&on| on).count(), CELLS);
    let links: Vec<Vec<usize>> = (0..SIDE * SIDE).map(|i| neighbours(&filled, i)).collect();
    let blacks: Vec<usize> = (0..SIDE * SIDE)
        .filter(|&i| filled[i] && (i / SIDE + i % SIDE).is_multiple_of(2))
        .collect();
    assert_eq!(2 * blacks.len(), CELLS);

    let mut mate = vec![NONE; SIDE * SIDE];
    for &black in &blacks {
        let mut seen = vec![false; SIDE * SIDE];
        assert!(augment(&links, &mut mate, &mut seen, black));
    }
    let mut rng = Rng::new(SEED);
    for _ in 0..FLIPS {
        flip(&filled, &mut mate, rng.below(SIDE - 1), rng.below(SIDE - 1));
    }
    for cell in 0..SIDE * SIDE {
        if filled[cell] {
            assert!(links[cell].contains(&mate[cell]));
            assert_eq!(mate[mate[cell]], cell);
        } else {
            assert_eq!(mate[cell], NONE);
        }
    }

    let dominoes: Vec<_> = blacks
        .iter()
        .map(|&black| {
            let white = mate[black];
            let (lo, hi) = (black.min(white), black.max(white));
            json!([lo / SIDE, lo % SIDE, usize::from(hi == lo + 1)])
        })
        .collect();
    assert_eq!(2 * dominoes.len(), CELLS);
    save(NAME, &json!({"side": SIDE, "dominoes": dominoes}))
}
