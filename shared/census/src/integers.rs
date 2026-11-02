use ledger::{keys, terms, Axis, Cost, Key, Measure, Tier};

const CAP: usize = 48;
const CELLS: u128 = 100_000;
const BLOCK: usize = 8;
const ROWS: usize = 18066;
const TIERS: [Tier; 4] = [
    Tier::Closed,
    Tier::Convolved,
    Tier::SideGrid,
    Tier::LevelGrid,
];
const MEASURES: usize = 12;

fn registry() -> Vec<Key> {
    TIERS
        .into_iter()
        .flat_map(keys)
        .filter(|key| Measure::ALL[..MEASURES].contains(&key.measure))
        .collect()
}

fn footprint(key: &Key, index: usize) -> Option<u128> {
    let (number, level) = key.axis.place(index, key.number());
    let number = number as u128;
    let dimension = key.dimension as u32;
    match key.measure.cost(Axis::Side) {
        Cost::Closed => Some(1),
        Cost::Convolved => {
            let tile = number.checked_pow(dimension)?;
            let side = number.checked_pow(level)?;
            let span = key.dimension as u128 * (side - 1) + 1;
            tile.checked_add(span.checked_mul(level as u128)?)
        }
        Cost::Grid => number.checked_pow(dimension.checked_mul(level)?),
    }
}

fn allowance(key: &Key) -> usize {
    (0..CAP)
        .take_while(|&index| footprint(key, index).is_some_and(|cells| cells <= CELLS))
        .count()
}

fn ceiling_stop(read: &[i128], ceiling: i128) -> Option<usize> {
    let mut previous: Option<i128> = None;
    for (index, &term) in read.iter().enumerate() {
        if previous.is_some_and(|last| term <= last) {
            return None;
        }
        if term > ceiling {
            return Some(index);
        }
        previous = Some(term);
    }
    None
}

fn rendered(key: &Key, ceiling: i128) -> Option<Vec<i128>> {
    let allowed = allowance(key);
    let mut count = BLOCK.min(allowed);
    loop {
        let (read, capped) = terms(key, count, CELLS).ok()?;
        if let Some(edge) = ceiling_stop(&read, ceiling) {
            return Some(read[..=edge].to_vec());
        }
        if capped || read.len() < count || count >= allowed {
            return Some(read);
        }
        count = (count * 2).min(allowed);
    }
}

fn tally(
    rows: &[Key],
    lane: usize,
    lanes: usize,
    ceiling: i128,
    window: usize,
    depths: &[usize],
) -> Vec<Vec<u32>> {
    let mut counts = vec![vec![0u32; window + 1]; depths.len()];
    for key in rows.iter().skip(lane).step_by(lanes) {
        let Some(read) = rendered(key, ceiling) else {
            continue;
        };
        for (slot, &depth) in depths.iter().enumerate() {
            let head = &read[..read.len().min(depth)];
            let mut written: Vec<usize> = head
                .iter()
                .filter(|&&term| (1..=window as i128).contains(&term))
                .map(|&term| term as usize)
                .collect();
            written.sort_unstable();
            written.dedup();
            for value in written {
                counts[slot][value] += 1;
            }
        }
    }
    counts
}

/// Counts, for each depth, how many of the 18066 ledger sequences write each value from 1 to `window` in their first `depth` terms, reading each sequence until a term passes `ceiling`; returns one tally of `window + 1` slots per depth, slot 0 always empty.
pub fn census(ceiling: i128, window: usize, depths: &[usize]) -> Vec<Vec<u32>> {
    let rows = registry();
    assert_eq!(rows.len(), ROWS);
    let lanes = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let mut counts = vec![vec![0u32; window + 1]; depths.len()];
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..lanes)
            .map(|lane| {
                let rows = &rows;
                scope.spawn(move || tally(rows, lane, lanes, ceiling, window, depths))
            })
            .collect();
        for handle in handles {
            let part = handle.join().expect("a lane of the census");
            for (total, tier) in counts.iter_mut().zip(part) {
                for (sum, count) in total.iter_mut().zip(tier) {
                    *sum += count;
                }
            }
        }
    });
    counts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_census_reads_its_registry() {
        assert_eq!(registry().len(), ROWS);
    }
}
