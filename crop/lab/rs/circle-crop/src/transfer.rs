use crate::{list, ones_table, root_floor, trim};

// SHELL LADDER

pub(crate) struct Rung {
    pub(crate) high: Vec<u64>,
    start: Vec<usize>,
    size: usize,
}

pub(crate) fn rung(radius: u64, level: u32) -> Rung {
    let scale = 3u64.pow(level);
    let square = radius * radius;
    let top = radius / scale;
    let mut high = Vec::with_capacity(top as usize + 2);
    for i in 0..=top + 1 {
        let reach = scale * i;
        if reach * reach > square {
            high.push(0);
        } else {
            high.push(root_floor(square - reach * reach) / scale);
        }
    }
    let mut start = Vec::with_capacity(top as usize + 1);
    let mut size = 0usize;
    for i in 0..=top as usize {
        start.push(size);
        size += (high[i] - high[i + 1] + 1) as usize;
    }
    assert_eq!(size as u64, 2 * top + 1);
    Rung { high, start, size }
}

impl Rung {
    fn top(&self) -> usize {
        self.high.len() - 2
    }

    fn place(&self, i: usize, y: u64) -> usize {
        self.start[i] + (y - self.high[i + 1]) as usize
    }
}

// LOCAL PATTERN

fn pattern(here: &Rung, below: &Rung) -> Vec<u16> {
    let mut out = vec![0u16; here.size];
    let mut children = 0u64;
    for i in 0..=here.top() {
        let lo = here.high[i + 1];
        let hi = here.high[i];
        for y in lo..=hi {
            let mut mask = 0u16;
            for k in 0..3usize {
                let column = 3 * i + k;
                if column > below.top() {
                    continue;
                }
                let inner_lo = below.high[column + 1];
                let inner_hi = below.high[column];
                for t in 0..3u64 {
                    let row = 3 * y + t;
                    if row >= inner_lo && row <= inner_hi {
                        mask |= 1 << (3 * k as u16 + t as u16);
                    }
                }
            }
            assert!(mask != 0);
            children += u64::from(mask.count_ones());
            out[here.place(i, y)] = mask;
        }
    }
    assert_eq!(children, below.size as u64);
    out
}

// EXACT RATIOS

fn below(left: (u128, u128), right: (u128, u128)) -> bool {
    left.0 * right.1 < right.0 * left.1
}

fn ratio(pair: (u128, u128)) -> f64 {
    pair.0 as f64 / pair.1 as f64
}

// PERRON BRACKET

fn weights(matrix: &[Vec<u64>], rows: &[u64]) -> Vec<u128> {
    let size = matrix.len();
    let mut vector = vec![1.0f64; size];
    for _ in 0..400 {
        let mut next = vec![0.0f64; size];
        for source in 0..size {
            if rows[source] == 0 {
                continue;
            }
            let mut total = 0.0f64;
            for target in 0..size {
                total += matrix[source][target] as f64 * vector[target];
            }
            next[source] = total / rows[source] as f64;
        }
        let peak = next.iter().cloned().fold(0.0f64, f64::max);
        if peak <= 0.0 {
            break;
        }
        for slot in next.iter_mut() {
            *slot /= peak;
        }
        vector = next;
    }
    vector
        .iter()
        .map(|value| ((value * 1.0e9).round() as u128).max(1))
        .collect()
}

fn bracket(matrix: &[Vec<u64>], rows: &[u64]) -> ((u128, u128), (u128, u128)) {
    let size = matrix.len();
    let vector = weights(matrix, rows);
    let mut low: Option<(u128, u128)> = None;
    let mut high: Option<(u128, u128)> = None;
    for source in 0..size {
        let mut total = 0u128;
        for target in 0..size {
            total += u128::from(matrix[source][target]) * vector[target];
        }
        let cell = (total, u128::from(rows[source]) * vector[source]);
        if low.is_none() || below(cell, low.unwrap()) {
            low = Some(cell);
        }
        if high.is_none() || below(high.unwrap(), cell) {
            high = Some(cell);
        }
    }
    (low.unwrap(), high.unwrap())
}

fn living(matrix: &[Vec<u64>]) -> Vec<usize> {
    let size = matrix.len();
    let mut live: Vec<bool> = (0..size)
        .map(|source| matrix[source].iter().any(|&count| count > 0))
        .collect();
    loop {
        let mut moved = false;
        for source in 0..size {
            if !live[source] {
                continue;
            }
            if !(0..size).any(|target| live[target] && matrix[source][target] > 0) {
                live[source] = false;
                moved = true;
            }
        }
        if !moved {
            break;
        }
    }
    (0..size).filter(|&source| live[source]).collect()
}

fn block(matrix: &[Vec<u64>], keep: &[usize]) -> Vec<Vec<u64>> {
    keep.iter()
        .map(|&source| keep.iter().map(|&target| matrix[source][target]).collect())
        .collect()
}

// DOEBLIN

fn doeblin(matrix: &[Vec<u64>], steps: usize) -> Vec<f64> {
    const UNIT: u128 = 1 << 40;
    let size = matrix.len();
    let rows: Vec<u128> = matrix
        .iter()
        .map(|row| row.iter().map(|&count| u128::from(count)).sum())
        .collect();
    let base: Vec<Vec<u128>> = (0..size)
        .map(|source| {
            (0..size)
                .map(|target| u128::from(matrix[source][target]) * UNIT / rows[source])
                .collect()
        })
        .collect();
    let mut power = base.clone();
    let mut out = Vec::new();
    for step in 1..=steps {
        if step > 1 {
            let mut next = vec![vec![0u128; size]; size];
            for source in 0..size {
                for middle in 0..size {
                    if power[source][middle] == 0 {
                        continue;
                    }
                    for target in 0..size {
                        next[source][target] += power[source][middle] * base[middle][target] / UNIT;
                    }
                }
            }
            power = next;
        }
        let mass: u128 = (0..size)
            .map(|target| (0..size).map(|source| power[source][target]).min().unwrap())
            .sum();
        out.push(mass as f64 / UNIT as f64);
    }
    out
}

// DOBRUSHIN

fn dobrushin(matrix: &[Vec<u64>]) -> (u128, u128) {
    let size = matrix.len();
    let rows: Vec<u64> = matrix.iter().map(|row| row.iter().sum()).collect();
    let mut worst = (0u128, 1u128);
    for first in 0..size {
        if rows[first] == 0 {
            continue;
        }
        for second in first + 1..size {
            if rows[second] == 0 {
                continue;
            }
            let (one, two) = (u128::from(rows[first]), u128::from(rows[second]));
            let mut gap = 0u128;
            for (&a, &b) in matrix[first].iter().zip(&matrix[second]).take(size) {
                let left = u128::from(a) * two;
                let right = u128::from(b) * one;
                gap += left.max(right) - left.min(right);
            }
            let cell = (gap, 2 * one * two);
            if below(worst, cell) {
                worst = cell;
            }
        }
    }
    worst
}

// OPERATOR

struct Operator {
    order: Vec<u16>,
    slot: Vec<usize>,
    all: Vec<Vec<u64>>,
    keep: Vec<Vec<u64>>,
    rows: Vec<u64>,
    wide_all: Vec<Vec<u64>>,
    wide_keep: Vec<Vec<u64>>,
    wide_rows: Vec<u64>,
    wide_live: usize,
    centres: u64,
    parents: u64,
}

fn operator(radius: u64, deep: u32, shallow: u32) -> Operator {
    let mut rungs: Vec<Rung> = Vec::new();
    for level in 0..=shallow + 1 {
        rungs.push(rung(radius, level));
    }
    let mut masks: Vec<Vec<u16>> = vec![Vec::new()];
    for level in 1..=shallow + 1 {
        masks.push(pattern(&rungs[level as usize], &rungs[level as usize - 1]));
    }
    let mut seen = [false; 512];
    for level in deep..=shallow + 1 {
        for &mask in masks[level as usize].iter() {
            seen[mask as usize] = true;
        }
    }
    let order: Vec<u16> = (0..512u16).filter(|&mask| seen[mask as usize]).collect();
    let mut slot = vec![usize::MAX; 512];
    for (index, &mask) in order.iter().enumerate() {
        slot[mask as usize] = index;
    }
    let size = order.len();
    let mut all = vec![vec![0u64; size]; size];
    let mut keep = vec![vec![0u64; size]; size];
    let mut rows = vec![0u64; size];
    let broad = size * 9;
    let mut wide_all = vec![vec![0u64; broad]; broad];
    let mut wide_keep = vec![vec![0u64; broad]; broad];
    let mut wide_rows = vec![0u64; broad];
    let mut centres = 0u64;
    let mut parents = 0u64;
    for level in deep..=shallow {
        let here = &rungs[level as usize + 1];
        let inner = &rungs[level as usize];
        let outer_mask = &masks[level as usize + 1];
        let inner_mask = &masks[level as usize];
        for i in 0..=here.top() {
            let lo = here.high[i + 1];
            let hi = here.high[i];
            for y in lo..=hi {
                let mask = outer_mask[here.place(i, y)];
                let source = slot[mask as usize];
                let seat = 3 * (i % 3) + (y % 3) as usize;
                let wide_source = source * 9 + seat;
                rows[source] += 1;
                wide_rows[wide_source] += 1;
                parents += 1;
                if mask & (1 << 4) != 0 {
                    centres += 1;
                }
                for bit in 0..9u16 {
                    if mask & (1 << bit) == 0 {
                        continue;
                    }
                    let column = 3 * i + (bit / 3) as usize;
                    let row = 3 * y + u64::from(bit % 3);
                    let target = slot[inner_mask[inner.place(column, row)] as usize];
                    let wide_target = target * 9 + 3 * (bit / 3) as usize + (bit % 3) as usize;
                    all[source][target] += 1;
                    wide_all[wide_source][wide_target] += 1;
                    if bit != 4 {
                        keep[source][target] += 1;
                        wide_keep[wide_source][wide_target] += 1;
                    }
                }
            }
        }
    }
    let seen: Vec<usize> = (0..broad).filter(|&index| wide_rows[index] > 0).collect();
    let wide_live = seen.len();
    let wide_all = block(&wide_all, &seen);
    let wide_keep = block(&wide_keep, &seen);
    let wide_rows: Vec<u64> = seen.iter().map(|&index| wide_rows[index]).collect();
    Operator {
        order,
        slot,
        all,
        keep,
        rows,
        wide_all,
        wide_keep,
        wide_rows,
        wide_live,
        centres,
        parents,
    }
}

// SURVIVAL LADDER

struct Ladder {
    filled: u64,
    total: u64,
    margin: Vec<f64>,
    condition: Vec<f64>,
    gain: Vec<f64>,
}

fn ladder(radius: u64, depth: usize, ones: &[u32]) -> Ladder {
    assert!(3u64.pow(depth as u32) > radius);
    let square = radius * radius;
    let mut filled = 0u64;
    let mut bad = vec![0u64; depth];
    let mut crest = vec![0u64; depth];
    let mut total = 0u64;
    for i in 0..=radius {
        let hi = root_floor(square - i * i);
        let step = (i + 1) * (i + 1);
        let lo = if step >= square {
            0
        } else {
            root_floor(square - step)
        };
        let first = ones[i as usize];
        for y in lo..=hi {
            let mask = first & ones[y as usize];
            total += 1;
            if mask == 0 {
                filled += 1;
                continue;
            }
            let mut rest = mask;
            while rest != 0 {
                bad[rest.trailing_zeros() as usize] += 1;
                rest &= rest - 1;
            }
            crest[(31 - mask.leading_zeros()) as usize] += 1;
        }
    }
    assert_eq!(total, 2 * radius + 1);
    let mut nested = vec![0u64; depth + 1];
    nested[0] = total;
    for (cut, slot) in nested.iter_mut().enumerate().skip(1) {
        *slot = filled + crest.iter().take(depth - cut).sum::<u64>();
    }
    assert_eq!(nested[depth], filled);
    let mass = total as f64;
    let mut margin = Vec::with_capacity(depth);
    let mut condition = Vec::with_capacity(depth);
    let mut gain = Vec::with_capacity(depth);
    let mut product = 1.0f64;
    for (position, &missed) in bad.iter().enumerate() {
        let cut = depth - position;
        let free = (total - missed) as f64 / mass;
        let tight = nested[cut] as f64 / nested[cut - 1] as f64;
        margin.push(free);
        condition.push(tight);
        gain.push(tight / free);
        product *= tight / free;
    }
    let direct = (filled as f64 / mass) / margin.iter().product::<f64>();
    assert!((product - direct).abs() <= 1.0e-12 * direct);
    assert_eq!(gain[depth - 1], 1.0);
    Ladder {
        filled,
        total,
        margin,
        condition,
        gain,
    }
}

fn survival(label: &str, radius: u64, depth: usize) {
    let ones = ones_table(depth);
    let read = ladder(radius, depth, &ones);
    let (filled, total) = (read.filled, read.total);
    let product: f64 = read.gain.iter().product();
    let logs: Vec<f64> = read.gain.iter().map(|value| value.ln()).collect();
    let sum: f64 = logs.iter().map(|value| value.abs()).sum();
    let tail: f64 = logs[..depth.saturating_sub(3)]
        .iter()
        .map(|value| value.abs())
        .sum();
    println!(
        "circle-crop ladder {label} r={radius} L={depth} C={filled} Cfull={total} psi={:.6} logsum={:.6} logtail={:.6} share={:.6}",
        product,
        trim(sum, true),
        trim(tail, true),
        trim(tail / sum, true)
    );
    println!(
        "circle-crop ladder {label} r={radius} marginal={} conditional={} gain={}",
        list(&read.margin),
        list(&read.condition),
        list(&read.gain)
    );
    let scaled: Vec<f64> = logs.iter().map(|value| value * 1.0e3).collect();
    println!(
        "circle-crop ladder {label} r={radius} loggain_x1000={}",
        list(&scaled)
    );
    let top = read.margin[depth - 1];
    let next = read.margin[depth - 2];
    println!(
        "circle-crop ladder {label} r={radius} g0={:.12} g1={:.12} top_marginal={:.12} next_marginal={:.12} both_centres_crossed={}",
        read.gain[depth - 1],
        read.gain[depth - 2],
        top,
        next,
        top < 1.0 && next < 1.0
    );
}

// PROFILE

struct Gather {
    depth: usize,
    rows: u64,
    logs: Vec<f64>,
    sums: Vec<f64>,
    sizes: Vec<f64>,
    counts: Vec<u64>,
    flat: u64,
}

fn gather(window: u32, step: u64) -> Gather {
    let depth = window as usize + 1;
    let ones = ones_table(depth);
    let low = 3u64.pow(window);
    let high = 3u64.pow(window + 1) - 1;
    let mut sums = vec![0.0f64; depth];
    let mut sizes = vec![0.0f64; depth];
    let mut counts = vec![0u64; depth];
    let mut logs = Vec::new();
    let mut rows = 0u64;
    let mut flat = 0u64;
    let mut radius = low;
    while radius <= high {
        let read = ladder(radius, depth, &ones);
        let mut total = 0.0f64;
        for position in 0..depth {
            let value = read.gain[position].ln();
            let rank = depth - 1 - position;
            sums[rank] += value;
            sizes[rank] += value.abs();
            counts[rank] += 1;
            total += value;
        }
        if read.gain[depth - 2] == 1.0 {
            flat += 1;
        }
        logs.push(total);
        rows += 1;
        radius += step;
    }
    Gather {
        depth,
        rows,
        logs,
        sums,
        sizes,
        counts,
        flat,
    }
}

fn bulk(read: &Gather) -> f64 {
    (0..read.depth)
        .map(|rank| read.sizes[rank] / read.counts[rank] as f64 * 1.0e3)
        .sum()
}

fn profile(label: &str, window: u32, samples: u64) -> (f64, f64, f64) {
    let low = 3u64.pow(window);
    let high = 3u64.pow(window + 1) - 1;
    let step = ((high - low + 1) / samples).max(1);
    let read = gather(window, step);
    let depth = read.depth;
    let rows = read.rows;
    let mean: f64 = read.logs.iter().sum::<f64>() / rows as f64;
    let spread = (read
        .logs
        .iter()
        .map(|value| (value - mean).powi(2))
        .sum::<f64>()
        / rows as f64)
        .sqrt();
    let worst = read.logs.iter().cloned().fold(f64::MIN, f64::max);
    let best = read.logs.iter().cloned().fold(f64::MAX, f64::min);
    let signed: Vec<f64> = (0..depth)
        .map(|rank| read.sums[rank] / read.counts[rank] as f64 * 1.0e3)
        .collect();
    let sized: Vec<f64> = (0..depth)
        .map(|rank| read.sizes[rank] / read.counts[rank] as f64 * 1.0e3)
        .collect();
    let deep: f64 = if depth > 6 {
        read.sums[6..].iter().sum::<f64>() / read.counts[6..].iter().sum::<u64>() as f64 * 1.0e3
    } else {
        0.0
    };
    let deep_size: f64 = if depth > 6 {
        read.sizes[6..].iter().sum::<f64>() / read.counts[6..].iter().sum::<u64>() as f64 * 1.0e3
    } else {
        0.0
    };
    println!(
        "circle-crop profile {label} window={low}..{high} L={depth} radii={rows} logpsi mean={:.6} sd={:.6} min={:.6} max={:.6} psi=[{:.6},{:.6}]",
        mean,
        trim(spread, true),
        trim(best, false),
        trim(worst, true),
        trim(best.exp(), false),
        trim(worst.exp(), true)
    );
    println!(
        "circle-crop profile {label} window={low}..{high} loggain_x1000={} absgain_x1000={}",
        list(&signed),
        list(&sized)
    );
    println!(
        "circle-crop profile {label} window={low}..{high} deep_signed_x1000={:.6} deep_abs_x1000={:.6} deep_ranks={}",
        deep,
        deep_size,
        depth.saturating_sub(6)
    );
    let total = bulk(&read);
    let steps: Vec<f64> = (2..depth)
        .map(|rank| sized[rank] / sized[rank - 1])
        .collect();
    println!(
        "circle-crop profile {label} window={low}..{high} decay={} tail_decay={:.6} rank1_flat={} of {rows}",
        list(&steps),
        trim(steps[1..].iter().cloned().fold(0.0f64, f64::max), true),
        read.flat
    );
    (total, trim(best.exp(), false), trim(worst.exp(), true))
}

fn tail(totals: &[f64]) -> (Vec<f64>, Vec<f64>, f64, f64) {
    let steps: Vec<f64> = (1..totals.len())
        .map(|index| totals[index] - totals[index - 1])
        .collect();
    let decay: Vec<f64> = (1..steps.len())
        .map(|index| steps[index] / steps[index - 1])
        .collect();
    let worst = decay.iter().cloned().fold(0.0f64, f64::max);
    let last = *steps.last().unwrap();
    let cap = totals.last().unwrap() + last * worst / (1.0 - worst);
    (steps, decay, worst, cap)
}

fn census(label: &str, totals: &[f64], floor: f64, roof: f64) {
    let (steps, decay, _, cap) = tail(totals);
    let band = (
        trim((-cap / 1000.0).exp(), false),
        trim((cap / 1000.0).exp(), true),
    );
    println!(
        "circle-crop profile {label} totals_x1000={} steps_x1000={} decay={} cap_x1000={:.6} psi=[{:.6},{:.6}] extremes=[{:.6},{:.6}] band_holds={}",
        list(totals),
        list(&steps),
        list(&decay),
        trim(cap, true),
        band.0,
        band.1,
        floor,
        roof,
        band.0 <= floor && roof <= band.1
    );
}

fn resample(label: &str, windows: &[u32], strides: &[u64], head: &[f64]) {
    for &step in strides {
        let mut totals: Vec<f64> = head.to_vec();
        let mut radii = Vec::new();
        for &window in windows {
            let read = gather(window, step);
            radii.push(read.rows as f64);
            totals.push(bulk(&read));
        }
        let (steps, _, worst, cap) = tail(&totals);
        println!(
            "circle-crop profile {label} stride={step} radii={} totals_x1000={} last_step_x1000={:.6} worst_decay={:.6} cap_x1000={:.6} converges={}",
            list(&radii),
            list(&totals),
            *steps.last().unwrap(),
            trim(worst, true),
            trim(cap, true),
            worst < 1.0
        );
    }
}

// REPORT

fn spectrum(label: &str, tag: &str, matrix: &[Vec<u64>], rows: &[u64], target: f64) -> (f64, f64) {
    let (low, high) = bracket(matrix, rows);
    let (bottom, top) = (trim(ratio(low), false), trim(ratio(high), true));
    let holds = bottom <= target && target <= top;
    println!(
        "circle-crop operator {label} {tag} states={} perron=[{bottom:.6},{top:.6}] target={target:.6} covers={holds}",
        matrix.len()
    );
    (bottom, top)
}

fn report(label: &str, radius: u64, deep: u32, shallow: u32) -> Vec<u16> {
    let built = operator(radius, deep, shallow);
    let size = built.order.len();
    let branches: Vec<u64> = built
        .order
        .iter()
        .map(|&mask| u64::from(mask.count_ones()))
        .collect();
    let mean: f64 = built
        .rows
        .iter()
        .zip(branches.iter())
        .map(|(&count, &arms)| (count * arms) as f64)
        .sum::<f64>()
        / built.parents as f64;
    let centre = built.centres as f64 / built.parents as f64;
    println!(
        "circle-crop operator {label} r={radius} levels={deep}..{shallow} states={size} parents={} meanb={:.6} centre={:.6} centre_x3={:.6}",
        built.parents,
        mean,
        centre,
        centre * 3.0
    );
    let widest = built
        .order
        .iter()
        .max_by_key(|&&mask| built.rows[built.slot[mask as usize]])
        .unwrap();
    println!(
        "circle-crop operator {label} r={radius} masks={} shares={}",
        built
            .order
            .iter()
            .map(|mask| format!("{mask:03o}"))
            .collect::<Vec<String>>()
            .join(","),
        list(
            &built
                .rows
                .iter()
                .map(|&count| count as f64 / built.parents as f64)
                .collect::<Vec<f64>>()
        )
    );
    println!(
        "circle-crop operator {label} r={radius} commonest={widest:03o} share={:.6}",
        built.rows[built.slot[*widest as usize]] as f64 / built.parents as f64
    );
    let whole = spectrum(label, "all", &built.all, &built.rows, 3.0);
    let alive = living(&built.keep);
    let trimmed = block(&built.keep, &alive);
    let kept_rows: Vec<u64> = alive.iter().map(|&index| built.rows[index]).collect();
    let pruned = spectrum(label, "pruned", &trimmed, &kept_rows, 8.0 / 3.0);
    let share = (
        trim(pruned.0 / whole.1, false),
        trim(pruned.1 / whole.0, true),
    );
    let drift = (
        trim((share.0 * 9.0 / 8.0).ln() / 3.0f64.ln(), false),
        trim((share.1 * 9.0 / 8.0).ln() / 3.0f64.ln(), true),
    );
    println!(
        "circle-crop operator {label} r={radius} share=[{:.6},{:.6}] target={:.6} drift=[{:.6},{:.6}]",
        share.0,
        share.1,
        8.0 / 9.0,
        drift.0,
        drift.1
    );
    let wide_whole = spectrum(label, "wide-all", &built.wide_all, &built.wide_rows, 3.0);
    let wide_alive = living(&built.wide_keep);
    let wide_trimmed = block(&built.wide_keep, &wide_alive);
    let wide_kept: Vec<u64> = wide_alive
        .iter()
        .map(|&index| built.wide_rows[index])
        .collect();
    let wide_pruned = spectrum(label, "wide-pruned", &wide_trimmed, &wide_kept, 8.0 / 3.0);
    let wide_share = (
        trim(wide_pruned.0 / wide_whole.1, false),
        trim(wide_pruned.1 / wide_whole.0, true),
    );
    println!(
        "circle-crop operator {label} r={radius} wide_states={} wide_share=[{:.6},{:.6}]",
        built.wide_live, wide_share.0, wide_share.1
    );
    let mixing = dobrushin(&built.all);
    let floors = doeblin(&built.all, 6);
    let least = floors
        .iter()
        .enumerate()
        .map(|(step, &mass)| (1.0 - mass).powf(1.0 / (step + 1) as f64))
        .fold(f64::MAX, f64::min);
    println!(
        "circle-crop operator {label} r={radius} dobrushin={:.6} doeblin={} rate={:.6}",
        trim(ratio(mixing), true),
        list(
            &floors
                .iter()
                .map(|mass| trim(*mass, false))
                .collect::<Vec<f64>>()
        ),
        trim(least, true)
    );
    built.order
}

// ALPHABET

pub(crate) fn levels(radius: u64) -> u32 {
    let mut level = 0u32;
    while 3u64.pow(level) <= radius {
        level += 1;
    }
    level
}

pub(crate) fn masks_seen(radius: u64, deep: u32, shallow: u32) -> Vec<u16> {
    let mut seen = [false; 512];
    for level in deep..=shallow {
        let here = rung(radius, level);
        let under = rung(radius, level - 1);
        for &mask in pattern(&here, &under).iter() {
            seen[mask as usize] = true;
        }
    }
    (0..512u16).filter(|&mask| seen[mask as usize]).collect()
}

pub(crate) fn show(masks: &[u16]) -> String {
    masks
        .iter()
        .map(|mask| format!("{mask:03o}"))
        .collect::<Vec<String>>()
        .join(",")
}

fn extra(seen: &[u16], base: &[u16]) -> Vec<u16> {
    seen.iter()
        .filter(|mask| !base.contains(mask))
        .cloned()
        .collect()
}

fn whole(label: &str, radius: u64, base: &[u16]) {
    let depth = levels(radius);
    let seen = masks_seen(radius, 1, depth);
    let cut = masks_seen(radius, 1, depth - 3);
    assert_eq!(cut, base);
    let more = extra(&seen, base);
    let root = pattern(&rung(radius, depth), &rung(radius, depth - 1))[0];
    println!(
        "circle-crop operator {label} r={radius} L={depth} every_level_states={} truncated_states={} extra={} root_state={root:03o}",
        seen.len(),
        cut.len(),
        show(&more)
    );
}

fn pinned(label: &str, radius: u64, base: &[u16], want: &[f64]) {
    let depth = levels(radius);
    let mut widths = Vec::new();
    let mut seats = [false; 512];
    for back in 1..=4u32 {
        let seen = masks_seen(radius, 1, depth - back);
        widths.push(seen.len() as f64);
        for &mask in extra(&seen, base).iter() {
            seats[mask as usize] = true;
        }
    }
    assert_eq!(widths, want);
    let found: Vec<u16> = (0..512u16).filter(|&mask| seats[mask as usize]).collect();
    let cut = masks_seen(radius, 1, depth - 3);
    let more = extra(&cut, base);
    let mut first = 0u32;
    if !more.is_empty() {
        for level in 1..=depth - 3 {
            if pattern(&rung(radius, level), &rung(radius, level - 1)).contains(&more[0]) {
                first = level;
                break;
            }
        }
    }
    println!(
        "circle-crop operator {label} r={radius} L={depth} states_at_L_less_1_to_4={} truncated_states={} extra={} first_level={first} all_extras={}",
        list(&widths),
        cut.len(),
        show(&more),
        show(&found)
    );
}

fn scan(label: &str, low: u64, high: u64, back: u32, base: &[u16]) {
    let mut rows = 0u64;
    let mut over = 0u64;
    let mut under = 0u64;
    let mut seats = [false; 512];
    let mut witness: Vec<u64> = Vec::new();
    for radius in low..=high {
        let depth = levels(radius);
        if depth < back + 2 {
            continue;
        }
        let seen = masks_seen(radius, 1, depth - back);
        rows += 1;
        let more = extra(&seen, base);
        if !more.is_empty() {
            over += 1;
            for &mask in more.iter() {
                seats[mask as usize] = true;
            }
            if witness.len() < 6 {
                witness.push(radius);
            }
        }
        if seen.len() - more.len() < base.len() {
            under += 1;
        }
    }
    let found: Vec<u16> = (0..512u16).filter(|&mask| seats[mask as usize]).collect();
    println!(
        "circle-crop operator {label} scan={low}..{high} levels=1..L-{back} radii={rows} over={over} under={under} extras={} first={}",
        show(&found),
        witness
            .iter()
            .map(|value| value.to_string())
            .collect::<Vec<String>>()
            .join(",")
    );
}

// MEMORY ONE

fn memory(label: &str, radius: u64, deep: u32, shallow: u32, start: u32) {
    let built = operator(radius, deep, shallow);
    let size = built.order.len();
    for source in 0..size {
        let arms = u64::from(built.order[source].count_ones());
        let row: u64 = built.all[source].iter().sum();
        assert_eq!(row, built.rows[source] * arms);
    }
    let masks = pattern(&rung(radius, start), &rung(radius, start - 1));
    let mut count = vec![0.0f64; size];
    for &mask in masks.iter() {
        count[built.slot[mask as usize]] += 1.0;
    }
    let mut model = Vec::new();
    let mut truth = Vec::new();
    let mut branch = Vec::new();
    let mut level = start;
    while level >= 1 {
        let mut next = vec![0.0f64; size];
        for (source, &weight) in count.iter().enumerate() {
            if weight == 0.0 || built.rows[source] == 0 {
                continue;
            }
            let share = weight / built.rows[source] as f64;
            for (slot, &cell) in next.iter_mut().zip(&built.all[source]) {
                *slot += share * cell as f64;
            }
        }
        model.push(next.iter().sum::<f64>());
        let big = (2 * (radius / 3u64.pow(level - 1)) + 1) as f64;
        let small = (2 * (radius / 3u64.pow(level)) + 1) as f64;
        truth.push(big);
        branch.push(big / small);
        count = next;
        level -= 1;
    }
    let error = model
        .iter()
        .zip(truth.iter())
        .map(|(read, want)| (read - want).abs())
        .fold(0.0f64, f64::max);
    let mut step = Vec::new();
    let mut slip = 0.0f64;
    for back in (1..=start).rev() {
        let seen = pattern(&rung(radius, back), &rung(radius, back - 1));
        let mut heads = vec![0.0f64; size];
        for &mask in seen.iter() {
            heads[built.slot[mask as usize]] += 1.0;
        }
        let mass: f64 = (0..size)
            .map(|source| heads[source] * f64::from(built.order[source].count_ones()))
            .sum();
        step.push(mass);
        slip = slip.max((mass - (2 * (radius / 3u64.pow(back - 1)) + 1) as f64).abs());
    }
    println!(
        "circle-crop operator {label} r={radius} rowsum_is_popcount=true start={start} onestep={} truth={} onestep_err={:.9} iterated={} iterated_maxerr={:.9} branching={}",
        list(&step),
        list(&truth),
        slip,
        list(&model),
        error,
        list(&branch)
    );
}

fn deepen(label: &str, radius: u64, top: u32) {
    let mut gaps = Vec::new();
    for back in (1..=top).rev() {
        let built = operator(radius, back, back);
        let (low, high) = bracket(&built.all, &built.rows);
        let (bottom, roof) = (trim(ratio(low), false), trim(ratio(high), true));
        let gap = if bottom <= 3.0 && 3.0 <= roof {
            0.0
        } else {
            (bottom - 3.0).abs().min((roof - 3.0).abs())
        };
        gaps.push(trim(gap, false));
        println!(
            "circle-crop operator {label} r={radius} pair={}..{back} parents={} states={} perron=[{bottom:.6},{roof:.6}] gap={:.6}",
            back + 1,
            built.parents,
            built.all.len(),
            trim(gap, false)
        );
    }
    println!(
        "circle-crop operator {label} r={radius} gap_by_pair={} monotone={}",
        list(&gaps),
        gaps.windows(2).all(|pair| pair[1] <= pair[0])
    );
}

fn signs(label: &str, radii: &[u64]) {
    let mut above = 0u32;
    let mut misses = 0u32;
    let mut least = f64::MAX;
    let mut most = 0.0f64;
    let mut rows = Vec::new();
    for &radius in radii {
        let depth = levels(radius);
        let built = operator(radius, 1, depth - 4);
        let (low, high) = bracket(&built.all, &built.rows);
        let (bottom, roof) = (trim(ratio(low), false), trim(ratio(high), true));
        if !(bottom <= 3.0 && 3.0 <= roof) {
            misses += 1;
        }
        let middle = (bottom + roof) / 2.0 - 3.0;
        if middle > 0.0 {
            above += 1;
        }
        least = least.min(middle.abs());
        most = most.max(middle.abs());
        rows.push(middle * 1.0e3);
    }
    println!(
        "circle-crop operator {label} sweep_radii={} perron_less_3_x1000={} above={above} below={} excludes_3={misses} size_x1000=[{:.6},{:.6}]",
        radii.len(),
        list(&rows),
        radii.len() as u32 - above,
        trim(least * 1.0e3, false),
        trim(most * 1.0e3, true)
    );
}

// INDEX BOUND

fn floor_shift(square: i64, back: i64, step: i64) -> i64 {
    let root = if square <= 0 {
        0
    } else {
        root_floor(square as u64) as i64
    };
    (root - back).div_euclid(step)
}

fn index_cap(reach: f64) -> f64 {
    13.60 * reach.powf(-1.0 / 3.0) + 305.08 / reach.sqrt() + 9.84 / reach
}

struct Index {
    rate: Vec<f64>,
    cap: Vec<f64>,
    slack: f64,
    drift: f64,
    logind: f64,
    live: usize,
}

fn index(radius: u64) -> Index {
    let square = radius * radius;
    let mut high = Vec::with_capacity(radius as usize + 2);
    for x in 0..=radius + 1 {
        high.push(if x * x > square {
            0
        } else {
            root_floor(square - x * x)
        });
    }
    let depth = levels(radius);
    let shell = 2 * radius + 1;
    let mut rate = Vec::new();
    let mut cap = Vec::new();
    let mut slack = 0.0f64;
    let mut drift = 0.0f64;
    let mut logind = 0.0f64;
    let mut live = 0usize;
    for level in 0..depth {
        let scale = 3u64.pow(level);
        let step = 3 * scale;
        let last = radius / scale;
        let mut cells = 0u64;
        let mut wide = 0u64;
        let mut tall = 0u64;
        let mut seats = 0u64;
        let mut seat_cells = 0u64;
        let mut seat_wide = 0u64;
        let mut seat_tall = 0u64;
        for c in 0..=last {
            let head = high[(scale * c) as usize];
            let foot = if scale * (c + 1) > radius {
                0
            } else {
                high[(scale * (c + 1)) as usize]
            };
            let low = foot / scale;
            let top = head / scale;
            let span = (top - low + 1) as usize;
            let mut runs = vec![0u64; span];
            let mut load = vec![0u64; span];
            for x in scale * c..=(scale * (c + 1) - 1).min(radius) {
                let up = high[x as usize] / scale;
                let down = high[(x + 1) as usize] / scale;
                for t in down..=up {
                    let seat = (t - low) as usize;
                    runs[seat] += 1;
                    let roof = high[x as usize].min(scale * (t + 1) - 1);
                    let base = high[(x + 1) as usize].max(scale * t);
                    load[seat] += roof - base + 1;
                }
            }
            for t in low..=top {
                let seat = (t - low) as usize;
                let rows = head.min(scale * (t + 1) - 1) - foot.max(scale * t) + 1;
                assert_eq!(load[seat], runs[seat] + rows - 1);
                cells += load[seat];
                wide += runs[seat];
                tall += rows;
                if c % 3 == 1 && t % 3 == 1 {
                    seats += 1;
                    seat_cells += load[seat];
                    seat_wide += runs[seat];
                    seat_tall += rows;
                }
            }
        }
        assert_eq!(cells, shell);
        assert_eq!(wide, last + radius + 1);
        assert_eq!(tall, last + radius + 1);
        assert_eq!(seat_wide, seat_tall);
        assert_eq!(seat_cells, 2 * seat_wide - seats);
        let mut form_wide = 0i64;
        for x in 0..=radius {
            if (x / scale) % 3 != 1 {
                continue;
            }
            let near = square as i64 - (x * x) as i64;
            let far = square as i64 - ((x + 1) * (x + 1)) as i64;
            form_wide += floor_shift(near, scale as i64, step as i64)
                - floor_shift(far, 2 * scale as i64, step as i64);
        }
        assert_eq!(form_wide, seat_wide as i64);
        let mut form_seats = 0i64;
        for c in 0..=last {
            if c % 3 != 1 {
                continue;
            }
            let near = square as i64 - (scale * c).pow(2) as i64;
            let far = square as i64 - (scale * (c + 1)).pow(2) as i64;
            form_seats += floor_shift(near, scale as i64, step as i64)
                - floor_shift(far, 2 * scale as i64, step as i64);
        }
        assert_eq!(form_seats, seats as i64);
        assert!(shell - seat_cells >= scale);
        let share = seat_cells as f64 / shell as f64;
        let reach = radius as f64 / scale as f64;
        let roof = index_cap(reach);
        let gap = (share - 1.0 / 9.0).abs();
        assert!(gap <= roof);
        if roof < 8.0 / 9.0 {
            live += 1;
        }
        slack = slack.max(gap / roof);
        drift += gap;
        logind += (1.0 - share).ln() - (8.0f64 / 9.0).ln();
        rate.push(share);
        cap.push(roof);
    }
    assert!(drift <= 781.0);
    assert!(logind.abs() <= 1191.0);
    Index {
        rate,
        cap,
        slack,
        drift,
        logind,
        live,
    }
}

fn live(label: &str, radius: u64) {
    let square = radius * radius;
    let mut seats = 0i64;
    let mut x = 1u64;
    while x <= radius {
        let near = square as i64 - (x * x) as i64;
        let far = square as i64 - ((x + 1) * (x + 1)) as i64;
        seats += floor_shift(near, 1, 3) - floor_shift(far, 2, 3);
        x += 3;
    }
    assert!(seats >= 0);
    let shell = 2 * radius + 1;
    let share = seats as f64 / shell as f64;
    let roof = index_cap(radius as f64);
    let gap = (share - 1.0 / 9.0).abs();
    assert!(roof < 8.0 / 9.0);
    assert!(gap <= roof);
    println!(
        "circle-crop index {label} live r={radius} j=0 seats={seats} p0={share:.9} gap={} cap={:.9} ratio={} headroom={:.9}",
        trim(gap, true),
        (roof * 1e9).ceil() / 1e9,
        trim(gap / roof, true),
        ((8.0 / 9.0 - roof) * 1e9).floor() / 1e9
    );
}

fn indexed(label: &str, radius: u64) -> usize {
    let read = index(radius);
    println!(
        "circle-crop index {label} r={radius} L={} live_levels={} of={} drift={} cap=781.000000 logind={} cap=1191.000000 slack={} rate={} bound={}",
        read.rate.len(),
        read.live,
        read.rate.len(),
        trim(read.drift, true),
        trim(read.logind.abs(), true),
        trim(read.slack, true),
        list(&read.rate),
        list(&read.cap)
    );
    read.live
}

// ENTRY

pub fn transfer() {
    println!("circle-crop transfer convention: every line below is the carpet, the crossing shell at level j is the whole grid's shell at the real radius r/3^j, a box's state is the 9-bit pattern of its crossed children printed in octal with bit 4 the centre, every matrix entry is an exact count of parent-child pairs over the printed levels, a gain is indexed by the depth from the top k = L - 1 - j so g_k = (T_(k+1)/T_k)/(1 - p_(L-1-k)), and a profile line samples its window at an even stride");
    for &radius in &[80u64, 242, 1000, 2186, 6560, 19682, 12345] {
        survival("carpet", radius, levels(radius) as usize);
    }
    let mut totals = Vec::new();
    let mut floor = f64::MAX;
    let mut roof = 0.0f64;
    for window in 3..=8u32 {
        let (total, best, worst) = profile("carpet", window, 96);
        totals.push(total);
        floor = floor.min(best);
        roof = roof.max(worst);
    }
    census("carpet", &totals, floor, roof);
    resample("carpet", &[6, 7, 8], &[24, 48, 72], &totals[..3]);
    let mut alphabet: Option<Vec<u16>> = None;
    for &(radius, deep, shallow) in &[(6560u64, 1u32, 4u32), (19682, 1, 5), (12345, 1, 5)] {
        let order = report("carpet", radius, deep, shallow);
        match alphabet {
            None => alphabet = Some(order),
            Some(ref first) => assert_eq!(first, &order),
        }
    }
    let base = alphabet.unwrap();
    println!(
        "circle-crop operator carpet alphabet={} states={} levels=1..L-3",
        show(&base),
        base.len()
    );
    memory("carpet", 6560, 1, 4, 5);
    deepen("carpet", 6560, 4);
    signs(
        "carpet",
        &[
            3001, 3901, 4801, 5701, 6601, 7501, 8401, 9301, 10201, 11101, 12001, 12901, 13801,
            14701, 15601, 16501, 17401, 18301,
        ],
    );
    for &radius in &[6560u64, 19682, 12345] {
        whole("carpet", radius, &base);
    }
    for &(radius, want) in &[
        (1395u64, [31.0f64, 30.0, 30.0, 30.0]),
        (1739, [31.0, 30.0, 30.0, 30.0]),
        (6570, [31.0, 30.0, 30.0, 30.0]),
        (15122, [31.0, 31.0, 31.0, 30.0]),
        (3182, [31.0, 31.0, 31.0, 31.0]),
    ] {
        pinned("carpet", radius, &base, &want);
    }
    scan("carpet", 3000, 19682, 3, &base);
    scan("carpet", 3000, 19682, 1, &base);
    let mut live_levels = 0usize;
    for &radius in &[80u64, 242, 1000, 2186, 6560, 12345, 19682] {
        live_levels += indexed("carpet", radius);
    }
    println!(
        "circle-crop index carpet totals radii=7 live_levels={live_levels} note=every level of these seven radii has cap >= 8/9, so the assert cannot fail there and the numbers are a check on the identities and not on the bound; the bound bites only at R >= 212957 and beats 1/9 only at R >= 23157375"
    );
    for &radius in &[212957u64, 531441, 2000000] {
        live("carpet", radius);
    }
}
