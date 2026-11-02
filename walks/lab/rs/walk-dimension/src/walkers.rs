use crate::design::{coords, strides};
use crate::gasket::Gasket;
use mrly_core::{Rng, Tensor};

pub const SEED: u64 = 7;
pub const WALKERS: usize = 20000;
pub const FIRST_FIT: f64 = 32.0;
pub const LAST_STEP: f64 = 1.2e5;
pub const RUNGS: usize = 110;

pub fn ladder() -> Vec<i64> {
    let step = LAST_STEP.log10() / (RUNGS - 1) as f64;
    let mut out: Vec<i64> = (0..RUNGS)
        .map(|index| 10f64.powf(index as f64 * step) as i64)
        .collect();
    out.sort_unstable();
    out.dedup();
    out
}

pub struct Trace {
    pub times: Vec<f64>,
    pub readings: Vec<f64>,
}

pub fn grid_walk(grid: &Tensor, rng: &mut Rng, walkers: usize, side_cap: f64) -> Trace {
    let shape = grid.shape.clone();
    let dims = shape.len();
    let stride = strides(&shape);
    let mut cells: Vec<usize> = (0..grid.size())
        .filter(|&flat| grid.at(flat) != 0)
        .collect();
    let mut low = vec![usize::MAX; dims];
    let mut high = vec![0usize; dims];
    for flat in &cells {
        for (axis, at) in coords(*flat, &shape).iter().enumerate() {
            low[axis] = low[axis].min(*at);
            high[axis] = high[axis].max(*at);
        }
    }
    let bulk: Vec<usize> = cells
        .iter()
        .copied()
        .filter(|flat| {
            coords(*flat, &shape).iter().enumerate().all(|(axis, at)| {
                let span = high[axis] - low[axis] + 1;
                span <= 16 || (*at >= low[axis] + span / 4 && *at <= high[axis] - span / 4)
            })
        })
        .collect();
    if bulk.len() >= 200 {
        cells = bulk;
    }
    let mut start: Vec<Vec<i64>> = Vec::with_capacity(walkers);
    let mut place: Vec<Vec<i64>> = Vec::with_capacity(walkers);
    let mut seat: Vec<usize> = Vec::with_capacity(walkers);
    for _ in 0..walkers {
        let flat = cells[rng.below(cells.len())];
        let at: Vec<i64> = coords(flat, &shape).iter().map(|v| *v as i64).collect();
        start.push(at.clone());
        place.push(at);
        seat.push(flat);
    }
    let narrow = shape.iter().min().copied().unwrap_or(0) as f64;
    let cap = (narrow / side_cap).powi(2);
    let ladder = ladder();
    let mut trace = Trace {
        times: Vec::new(),
        readings: Vec::new(),
    };
    let mut rung = 0;
    let mut time = 0i64;
    while rung < ladder.len() {
        time += 1;
        for walker in 0..walkers {
            let drawn = rng.below(2 * dims);
            let axis = drawn / 2;
            let step: i64 = if drawn.is_multiple_of(2) { 1 } else { -1 };
            let moved = place[walker][axis] + step;
            if moved < 0 || moved >= shape[axis] as i64 {
                continue;
            }
            let next = (seat[walker] as i64 + step * stride[axis] as i64) as usize;
            if grid.at(next) != 0 {
                place[walker][axis] = moved;
                seat[walker] = next;
            }
        }
        if time == ladder[rung] {
            let total: i64 = (0..walkers)
                .map(|walker| {
                    (0..dims)
                        .map(|axis| {
                            let gap = place[walker][axis] - start[walker][axis];
                            gap * gap
                        })
                        .sum::<i64>()
                })
                .sum();
            let reading = total as f64 / walkers as f64;
            trace.times.push(time as f64);
            trace.readings.push(reading);
            rung += 1;
            if reading > cap {
                break;
            }
        }
    }
    trace
}

pub fn gasket_walk(gasket: &Gasket, rng: &mut Rng, walkers: usize) -> Trace {
    let nodes = gasket.points.len();
    let table: Vec<[i64; 4]> = gasket
        .graph
        .adjacency
        .iter()
        .map(|row| {
            let mut slots = [-1i64; 4];
            for (slot, other) in row.iter().enumerate().take(4) {
                slots[slot] = *other as i64;
            }
            slots
        })
        .collect();
    let start: Vec<usize> = (0..walkers).map(|_| rng.below(nodes)).collect();
    let mut place = start.clone();
    let reach = gasket
        .points
        .iter()
        .flat_map(|point| [point.0, point.1])
        .max()
        .unwrap_or(0) as f64;
    let cap = (reach / 4.0).powi(2);
    let ladder = ladder();
    let mut trace = Trace {
        times: Vec::new(),
        readings: Vec::new(),
    };
    let mut rung = 0;
    let mut time = 0i64;
    while rung < ladder.len() {
        time += 1;
        for walker in 0..walkers {
            let next = table[place[walker]][rng.below(4)];
            if next >= 0 {
                place[walker] = next as usize;
            }
        }
        if time == ladder[rung] {
            let total: i64 = (0..walkers)
                .map(|walker| {
                    let here = gasket.points[place[walker]];
                    let there = gasket.points[start[walker]];
                    (here.0 - there.0).pow(2) + (here.1 - there.1).pow(2)
                })
                .sum();
            let reading = total as f64 / walkers as f64;
            trace.times.push(time as f64);
            trace.readings.push(reading);
            rung += 1;
            if reading > cap {
                break;
            }
        }
    }
    trace
}

pub fn slope(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let cov: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let var: f64 = x.iter().map(|a| (a - mx) * (a - mx)).sum();
    cov / var
}

pub fn fit(trace: &Trace, ceiling: f64) -> (f64, f64) {
    let mut x: Vec<f64> = Vec::new();
    let mut y: Vec<f64> = Vec::new();
    for (time, reading) in trace.times.iter().zip(&trace.readings) {
        if *time >= FIRST_FIT && *reading < ceiling {
            x.push(time.ln());
            y.push(reading.ln());
        }
    }
    if x.len() < 6 {
        return (f64::NAN, f64::NAN);
    }
    let half = x.len() / 2;
    let first = 2.0 / slope(&x[..half], &y[..half]);
    let second = 2.0 / slope(&x[half..], &y[half..]);
    (2.0 / slope(&x, &y), (first - second).abs())
}
