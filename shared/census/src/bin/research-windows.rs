use census::save;
use mrly_core::error::Result;
use mrly_core::json;
use mrly_math::bang::factory;
use mrly_math::bang::Code;
use std::collections::BTreeSet;

const NAME: &str = "research-windows";
const LEVEL: usize = 3;
const SIDE: usize = 27;
const SCAN: usize = 5;
const COUNTS: [usize; 4] = [2, 10, 40, 74];
const RING: usize = 11;

type Window = Vec<u8>;

fn render(level: usize) -> Result<(usize, Vec<u8>)> {
    let design = factory::create(Code::from(495u128), 3, 2, 3, level)?;
    let side = 3usize.pow(level as u32);
    assert_eq!(design.shape, vec![side, side]);
    assert_eq!(design.sum(), 8u64.pow(level as u32));
    let parity = factory::create(Code::from(7u128), 3, 2, 2, level)?;
    assert_eq!(design.bytes()?, parity.bytes()?);
    Ok((side, design.bytes()?.to_vec()))
}

fn windows(side: usize, cells: &[u8], k: usize) -> BTreeSet<Window> {
    let mut out = BTreeSet::new();
    for r in 0..=side - k {
        for c in 0..=side - k {
            out.insert(
                (0..k * k)
                    .map(|i| cells[(r + i / k) * side + c + i % k])
                    .collect(),
            );
        }
    }
    out
}

fn turn(w: &[u8]) -> Window {
    (0..9).map(|i| w[(2 - i % 3) * 3 + i / 3]).collect()
}

fn orbit(w: &[u8]) -> Vec<Window> {
    let mut out = vec![w.to_vec()];
    loop {
        let next = turn(&out[out.len() - 1]);
        if next == out[0] {
            return out;
        }
        out.push(next);
    }
}

fn ones(w: &[u8]) -> usize {
    w.iter().filter(|&&v| v == 1).count()
}

fn low(w: &[u8]) -> usize {
    (0..9).filter(|&i| w[i] == 1).map(|i| i / 3).sum()
}

fn slots(three: &BTreeSet<Window>) -> Vec<Window> {
    let mut seen = BTreeSet::new();
    let mut reps = Vec::new();
    let mut fixed = Vec::new();
    for w in three {
        if seen.contains(w) {
            continue;
        }
        let orb = orbit(w);
        seen.extend(orb.iter().cloned());
        if orb.len() == 4 {
            reps.push(orb.into_iter().max_by_key(|o| (low(o), o.clone())).unwrap());
        } else {
            fixed.extend(orb);
        }
    }
    reps.sort_by_key(|w| (ones(w), w.clone()));
    assert_eq!(reps.len(), RING - 2);
    assert_eq!(fixed.len(), 4);
    let pick = |n: usize, head: u8| -> Window {
        let found: Vec<&Window> = fixed
            .iter()
            .filter(|w| ones(w) == n && w[0] == head)
            .collect();
        assert_eq!(found.len(), 1);
        found[0].clone()
    };
    let mut slots = vec![pick(7, 1), pick(7, 0), pick(8, 1), pick(0, 0)];
    assert_eq!(turn(&slots[0]), slots[1]);
    for rep in &reps {
        let a = turn(rep);
        let b = turn(&a);
        let c = turn(&b);
        slots.extend([rep.clone(), a, b, c]);
    }
    let drawn: BTreeSet<Window> = slots.iter().cloned().collect();
    assert_eq!(&drawn, three);
    assert_eq!(slots.len(), 4 * RING - 4);
    slots
}

fn mask(w: &[u8]) -> usize {
    assert!(w.iter().all(|&v| v < 2));
    w.iter()
        .enumerate()
        .fold(0, |acc, (i, &v)| acc | (v as usize) << i)
}

fn main() -> Result<()> {
    let (side, scan) = render(SCAN)?;
    let (before, coarse) = render(SCAN - 1)?;
    for (k, count) in COUNTS.iter().enumerate() {
        let found = windows(side, &scan, k + 1);
        assert_eq!(found.len(), *count);
        assert_eq!(found, windows(before, &coarse, k + 1));
    }
    let three = windows(side, &scan, 3);
    let ring = slots(&three);
    let (ground, _) = render(LEVEL)?;
    assert_eq!(ground, SIDE);
    let masks: Vec<usize> = ring.iter().map(|w| mask(w)).collect();
    save(NAME, &json!({ "windows": masks }))
}
