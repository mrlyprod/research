use ledger::{closed, keys, terms, Axis, Closed, Cost, Key, Measure, Tier, TERMS};
use std::thread;

pub const CEILING: i128 = 100_000;
pub const CAP: usize = 48;
pub const DEEP: usize = 96;
pub const SHALLOW: usize = 32;
pub const CELLS: u128 = 100_000;
pub const BLOCK: usize = 8;
pub const THREADS: usize = 4;
pub const TIERS: [Tier; 4] = [
    Tier::Closed,
    Tier::Convolved,
    Tier::SideGrid,
    Tier::LevelGrid,
];
pub const MEASURES: usize = 12;

pub fn counted(measure: Measure) -> bool {
    Measure::ALL[..MEASURES].contains(&measure)
}

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Stop {
    Ceiling,
    Cap,
    Budget,
}

impl Stop {
    pub fn slug(self) -> &'static str {
        match self {
            Stop::Ceiling => "ceiling",
            Stop::Cap => "cap",
            Stop::Budget => "budget",
        }
    }
}

pub struct Row {
    pub name: String,
    pub tier: Tier,
    pub form: Option<Closed>,
    pub head: Vec<i128>,
    pub written: Vec<i128>,
    pub shallow: Vec<i128>,
    pub tail: Vec<i128>,
    pub repeats: usize,
    pub low: usize,
    pub stop: Stop,
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

fn ceiling_stop(read: &[i128]) -> Option<usize> {
    let mut previous: Option<i128> = None;
    for (index, &term) in read.iter().enumerate() {
        if previous.is_some_and(|last| term <= last) {
            return None;
        }
        if term > CEILING {
            return Some(index);
        }
        previous = Some(term);
    }
    None
}

fn gather(window: &[i128]) -> (Vec<i128>, usize, usize) {
    let low = window.iter().filter(|&&term| term < 1).count();
    let mut inside: Vec<i128> = window
        .iter()
        .copied()
        .filter(|&term| (1..=CEILING).contains(&term))
        .collect();
    let all = inside.len();
    inside.sort_unstable();
    inside.dedup();
    let repeats = all - inside.len();
    (inside, repeats, low)
}

fn render(key: &Key, tier: Tier) -> Option<Row> {
    let allowed = allowance(key);
    let mut count = BLOCK.min(allowed);
    let mut head = Vec::new();
    let window;
    let stop;
    loop {
        let (read, capped) = terms(key, count, CELLS).ok()?;
        let short = capped || read.len() < count;
        if head.len() < TERMS.min(read.len()) {
            head = read[..TERMS.min(read.len())].to_vec();
        }
        if let Some(edge) = ceiling_stop(&read) {
            window = read[..=edge].to_vec();
            stop = Stop::Ceiling;
            break;
        }
        if short {
            window = read;
            stop = Stop::Budget;
            break;
        }
        if count >= allowed {
            window = read;
            stop = if allowed == CAP {
                Stop::Cap
            } else {
                Stop::Budget
            };
            break;
        }
        count = (count * 2).min(allowed);
    }
    let (written, repeats, low) = gather(&window);
    let (shallow, _, _) = gather(&window[..SHALLOW.min(window.len())]);
    let (tail, _, _) = gather(&window[1.min(window.len())..]);
    Some(Row {
        name: key.name(),
        tier,
        form: closed(key).ok().flatten(),
        head,
        written,
        shallow,
        tail,
        repeats,
        low,
        stop,
    })
}

pub fn predict(form: &Closed, index: usize) -> Option<i128> {
    match form {
        Closed::Power(fill) => i128::try_from(*fill).ok()?.checked_pow(index as u32 + 1),
        Closed::Difference(all, fill) => {
            let level = index as u32 + 1;
            let whole = i128::try_from(*all).ok()?.checked_pow(level)?;
            whole.checked_sub(i128::try_from(*fill).ok()?.checked_pow(level)?)
        }
        Closed::Polynomial(coefficients) => {
            let side = index as i128 + 2;
            coefficients
                .iter()
                .enumerate()
                .try_fold(0i128, |sum, (power, &c)| {
                    sum.checked_add(c.checked_mul(side.checked_pow(power as u32)?)?)
                })
        }
        Closed::Recurrence(_) => None,
    }
}

pub fn replay(coefficients: &[i128], head: &[i128], index: usize) -> Option<i128> {
    coefficients
        .iter()
        .enumerate()
        .try_fold(0i128, |sum, (back, &c)| {
            sum.checked_add(c.checked_mul(*head.get(index.checked_sub(back + 1)?)?)?)
        })
}

pub struct Sheet {
    pub rows: Vec<Row>,
    pub tiers: Vec<(Tier, usize)>,
    pub unread: usize,
}

pub fn read() -> Sheet {
    let mut tiers = Vec::new();
    let mut listed: Vec<(Key, Tier)> = Vec::new();
    for tier in TIERS {
        let batch: Vec<Key> = keys(tier)
            .into_iter()
            .filter(|key| counted(key.measure))
            .collect();
        tiers.push((tier, batch.len()));
        listed.extend(batch.into_iter().map(|key| (key, tier)));
    }
    let chunk = listed.len().div_ceil(THREADS);
    let parts: Vec<Vec<Option<Row>>> = thread::scope(|scope| {
        let handles: Vec<_> = listed
            .chunks(chunk)
            .map(|slice| {
                scope.spawn(move || slice.iter().map(|(key, tier)| render(key, *tier)).collect())
            })
            .collect();
        handles
            .into_iter()
            .map(|handle| handle.join().expect("the row walk lands"))
            .collect()
    });
    let mut rows = Vec::with_capacity(listed.len());
    let mut unread = 0;
    for part in parts {
        for slot in part {
            match slot {
                Some(row) => rows.push(row),
                None => unread += 1,
            }
        }
    }
    Sheet {
        rows,
        tiers,
        unread,
    }
}
